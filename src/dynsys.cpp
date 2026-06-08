#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <cglm/call.h>
#include <cglm/cglm.h>

#include <algorithm>
#include <array>
#include <functional>
#include <chrono>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <unistd.h>
#include <deque>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <atomic>
#include <memory>

extern "C" {
#include "arena.h"
#include "ast.h"
#include "pratt.h"
}

#include "expr_ir.h"

#include "analysis.h"
#include "expr_ir_ad.h"
#include "cas_bridge.h"

#define PNG_WRITER_IMPLEMENTATION
#include "png_writer.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

namespace {

struct Point {
  float x;
  float y;
  float z;
};

struct State {
  std::vector<double> v;
  double t = 0.0;

  State() : v(3, 0.0), t(0.0) {}
  State(double x, double y, double z, double time = 0.0) : v{x, y, z}, t(time) {}
  explicit State(size_t dim, double value = 0.0, double time = 0.0) : v(dim, value), t(time) {}
};

double state_at(const State &s, size_t index) {
  return index < s.v.size() ? s.v[index] : 0.0;
}

void set_state_at(State &s, size_t index, double value) {
  if (index >= s.v.size()) {
    s.v.resize(index + 1, 0.0);
  }
  s.v[index] = value;
}

void resize_state(State &s, size_t dim) {
  if (dim == 0) dim = 1;
  s.v.resize(dim, 0.0);
}

State make_state_like(size_t dim, double time = 0.0) {
  return State(dim, 0.0, time);
}

[[maybe_unused]] std::string join_names(const std::vector<std::string> &names) {
  std::string out;
  for (size_t i = 0; i < names.size(); ++i) {
    if (i) out += ", ";
    out += names[i];
  }
  return out;
}

struct Point2 {
  double x = 0.0;
  double y = 0.0;
};

struct KeyState {
  bool pressed = false;
  double last_pressed_time = 0.0;
};

enum class SystemMode {
  ODE,
  Map,
  IFS,  /* iterated function system: affine maps + chaos game */
};

enum class Integrator {
  Euler,
  RK2,      /* midpoint */
  Heun,     /* explicit trapezoid (RK2, alternative) */
  RK4,      /* classic 4th order */
  RK38,     /* 3/8 rule, 4th order */
  RKF45,    /* Runge-Kutta-Fehlberg, adaptive step (embedded 4/5) */
  DOPRI45,  /* Dormand-Prince, adaptive step (embedded 5/4) */
};

enum class SectionDirection {
  Any,
  Positive,
  Negative,
};

/* Infer a reasonable [lo,hi] slider range for a value that has no declared
 * range, based on its current magnitude. Symmetric around 0 for small
 * values, widened to a round number so dragging feels natural. */
void auto_slider_range(double value, double *lo, double *hi) {
  const double a = std::fabs(value);
  double span;
  if (a < 1e-9) span = 1.0;               /* zero-ish: a unit window */
  else {
    /* round the magnitude up to 1/2/5 x 10^k, then use a few of them */
    const double mag = std::pow(10.0, std::floor(std::log10(a)));
    const double n = a / mag;
    const double step = n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10;
    span = 4.0 * step * mag;
  }
  if (value >= 0) { *lo = 0.0; *hi = span; if (value > span) *hi = value * 1.5; }
  else { *hi = 0.0; *lo = -span; if (value < -span) *lo = value * 1.5; }
  /* if the value sits at an edge, pad both sides so it isn't pinned */
  if (*hi - *lo < 1e-9) { *lo = value - 1.0; *hi = value + 1.0; }
}

const char *mode_name(SystemMode mode) {
  switch (mode) {
    case SystemMode::Map: return "map";
    case SystemMode::IFS: return "ifs";
    default: return "ode";
  }
}

const char *integrator_name(Integrator integrator) {
  switch (integrator) {
  case Integrator::Euler: return "Euler";
  case Integrator::RK2: return "RK2 midpoint";
  case Integrator::Heun: return "Heun (RK2)";
  case Integrator::RK4: return "RK4";
  case Integrator::RK38: return "RK 3/8";
  case Integrator::RKF45: return "RKF45 (adaptive)";
  case Integrator::DOPRI45: return "Dormand-Prince (adaptive)";
  }
  return "unknown";
}

struct AppState {
  arena_t system_arena{};
  bool arena_ready = false;

  struct Definition {
    std::string name;
    std::vector<std::string> params;
    node_t *body = nullptr;
    int line = 0;
    bool derivative = false;
    int derivative_index = -1;
    bool next_equation = false;
    int next_index = -1;
  };

  struct Param {
    std::string name;
    node_t *default_body = nullptr;
    double value = 0.0;
    double default_value = 0.0;
    double min_value = -10.0;
    double max_value = 10.0;
    bool has_range = false;
    int line = 0;
    /* UI scratch (auto-slider): when the model declared no [min,max], we
     * synthesize an editable range from the value's magnitude. */
    double ui_lo = 0.0, ui_hi = 0.0;
    bool ui_range_init = false;
    bool ui_show_exact = false;
  };

  struct Observable {
    std::string name;
    node_t *body = nullptr;
    int line = 0;
  };

  struct PhaseTrajectory {
    std::string label;
    State current;
    std::deque<State> history;
    bool active = true;
    bool visible = true;
    uint32_t color = 0xffffffffu;
  };

  std::vector<Definition> definitions;
  std::vector<Param> params;
  std::vector<Observable> observables;
  std::vector<std::string> state_names = {"x", "y", "z"};
  std::vector<node_t *> equations;
  std::vector<node_t *> next_equations;
  std::array<node_t *, 3> plot3d_bodies = {nullptr, nullptr, nullptr};
  std::array<std::string, 3> plot3d_labels = {"x", "y", "z"};

  /* === Lowered IR programs (parallel to the AST arrays above) ===
   * Built in compile_system after parsing finishes. The hot loop
   * runs these instead of walking the AST.
   *
   * param_values is a flat copy of params[i].value kept in sync with
   * the GUI; the IR run path reads this without going through the
   * Param struct on every dispatch. */
  std::vector<dynsys::ir::Program> equation_programs;
  std::vector<dynsys::ir::Program> next_equation_programs;
  std::array<dynsys::ir::Program, 3> plot3d_programs;
  dynsys::ir::Program section_program;
  dynsys::ir::Program section_x_program;
  dynsys::ir::Program section_y_program;
  std::vector<dynsys::ir::Program> definition_programs;
  std::vector<dynsys::ir::Program> observable_programs;
  std::vector<double> param_values;
  dynsys::ir::Scratch eval_scratch;
  /* PHASE3: scratch for the forward-mode-AD executor used by the
   * analysis layer to build exact Jacobian columns and df/dp. Sized
   * to the def count whenever the system recompiles. */
  dynsys::ir::DualScratch ad_scratch;

  /* When true, hot-path callers (eval_rhs, step_map_state,
   * eval_plot3d, maybe_record_poincare, value_by_name) fall back to
   * the original eval_expr_at AST walker instead of the IR. Only
   * used for differential benchmarking from the --headless driver. */
  bool use_ast_fallback = false;

  /* === Integrator scratch states ===
   * Pre-allocated to dim, reused every step so RK4 doesn't allocate
   * 8 vector<double> per step. */
  State scratch_k1, scratch_k2, scratch_k3, scratch_k4, scratch_mid;

  SystemMode mode = SystemMode::ODE;
  /* PHASE C: when mode == IFS, the system is a list of affine maps run via
   * the chaos game (no state-vector ODE/map equations). The maps may have
   * coefficients that are EXPRESSIONS in the parameters (e.g. a rotation
   * angle), so we store parsed ASTs and evaluate them against the current
   * parameter values at render time. ifs_maps holds the last evaluated
   * numeric maps (for chaos_game / dimension). */
  std::vector<dynsys::analysis::AffineMap> ifs_maps;
  struct IFSMapExpr { node_t *a=nullptr,*b=nullptr,*c=nullptr,*d=nullptr,*e=nullptr,*f=nullptr,*p=nullptr; };
  std::vector<IFSMapExpr> ifs_map_exprs;   /* parsed coefficient expressions */
  /* per-coefficient editability: a coefficient is directly editable iff its
   * expression is a plain literal number. Index order: a,b,c,d,e,f,p. A
   * parameter-driven coefficient (e.g. s*cos(theta)) is not directly
   * editable — its slider is shown disabled and tracks the live value. */
  std::vector<std::array<bool,7>> ifs_coef_literal;
  arena_t ifs_arena{};                     /* owns the ifs_map_exprs ASTs */
  bool ifs_arena_init = false;
  /* true when EVERY coefficient is a plain constant: then add/remove-map is
   * offered (a fully hand-built IFS). Per-coefficient editing works
   * regardless, via ifs_coef_literal. */
  bool ifs_maps_editable = false;
  Integrator integrator = Integrator::RK4;
  /* PHASE B: adaptive-step controls (RKF45 / DOPRI45). The visible dt is
   * the target output step; adaptive methods subdivide internally to meet
   * this error tolerance, clamping each substep to [dt_min, dt_max]. */
  double adaptive_tol = 1e-6;
  double adaptive_dt_min = 1e-7;
  double adaptive_dt_max = 0.1;
  double last_adaptive_dt = 0.0; /* reported: last accepted substep size */

  char system_input[16384] = {0};
  std::string parse_error;
  std::string runtime_error;
  std::string analysis_message;

  State current{0.1, 0.1, 0.1, 0.0};
  State start{0.1, 0.1, 0.1, 0.0};
  double dt = 0.01;

  int num_points = 20000;
  int current_point = 0;
  int steps_per_frame = 1;
  bool paused = false;
  int selected_preset = 0;

  float center_x = 0.0f;
  float center_y = 0.0f;
  float center_z = 0.0f;
  float angle_x = 0.0f;
  float angle_y = 0.0f;
  float zoom = 1.0f;
  float scene_scale = 0.05f;
  float camera_distance = 50.0f;
  bool bridge_cam_init = false;   /* set a geometry-fitting camera on first bridge frame */
  bool orthographic_3d = false;
  bool show_axes = true;
  bool show_center_cross = true;

  int phase_x_index = 0;
  int phase_y_index = 1;
  bool phase_show_trajectory = true;
  bool phase_show_vector_field = true;
  bool phase_show_nullclines = true;
  bool phase_show_nullcline_arrows = true;
  bool phase_auto_bounds = true;
  float phase_x_min = -20.0f;
  float phase_x_max = 20.0f;
  float phase_y_min = -20.0f;
  float phase_y_max = 20.0f;
  /* "home" phase window from the preset's view2d (or a default), so a Reset
   * View button can always return to a sane framing. */
  bool home_view_set = false;
  float home_x_min = -10.0f, home_x_max = 10.0f, home_y_min = -10.0f, home_y_max = 10.0f;
  /* PHASE0/6: smoothed auto-bounds. The raw data extent is recomputed
   * each frame, but the displayed window eases toward it and only grows
   * (or shrinks slowly when the data has shrunk for a while), so a
   * wandering or converging orbit no longer makes the view jitter or
   * collapse frame-to-frame. */
  bool phase_bounds_valid = false;
  double phase_view_xmin = -1.0, phase_view_xmax = 1.0;
  double phase_view_ymin = -1.0, phase_view_ymax = 1.0;
  int phase_grid = 21;
  bool phase_equal_aspect = true;
  bool phase_normalize_vectors = true;
  bool phase_show_fixed_points = true;
  bool phase_show_separatrices = true;
  int phase_trace_limit = 5000;
  int phase_max_extra_trajectories = 64;
  std::vector<PhaseTrajectory> phase_trajectories;

  int fps = 0;
  /* Adaptive performance governor: when frames get expensive (heavy fractal /
   * IFS / large step counts), automatically shed load to keep the UI
   * responsive, and restore it when there's headroom. perf_frame_ms is an
   * exponentially-smoothed frame time; perf_throttle in [0,1] is how much we're
   * currently holding back (0 = full speed, 1 = maximally throttled). User can
   * disable from the Setup tab. */
  bool perf_governor_enabled = true;
  double perf_frame_ms = 16.0;   /* smoothed ms/frame */
  double perf_throttle = 0.0;    /* 0..1 current throttle level */
  bool perf_throttle_active = false; /* true while throttling (for the HUD note) */
  int window_width = 1100;
  int window_height = 820;
  std::string screenshot_msg;     /* transient "saved <path>" toast */
  int screenshot_msg_timer = 0;
  bool capture_request = false;   /* set by Save PNG; serviced post-render */
  bool capturing_clean = false;   /* true during the extra panel-free capture frame */

  GLuint vbo = 0;
  GLuint vao = 0;
  GLuint cbo = 0;
  GLuint shader_program = 0;
  mat4 projection{};

  /* PHASE6-UI: offscreen render target for the 3D scene, so it lives in
   * a real ImGui panel (as a texture) instead of being the fullscreen
   * GL backdrop. Sized to the panel each frame; recreated on resize. */
  GLuint scene_fbo = 0;
  GLuint scene_tex = 0;
  GLuint scene_depth = 0;
  int scene_tex_w = 0;
  int scene_tex_h = 0;
  bool scene_hovered = false; /* is the mouse over the 3D panel image? */

  /* PHASE6-UI: dimensionality model. effective dimension is auto-detected
   * (a state whose dynamics are trivial — dz=0 for ODEs, z_next=z for
   * maps — doesn't count), but the user can force a choice. */
  enum class DimOverride { Auto, Force2D, Force3D };
  DimOverride dim_override = DimOverride::Auto;
  int detected_dim = 3; /* recomputed at compile time */

  /* PHASE6-UI: which plot fills the window background. */
  /* PHASE-A: the active system exposes a set of full-area views; the
   * toolbar offers only the ones valid for its dimension/mode. */
  enum class ActiveView { Line1D, Phase2D, Scene3D, Bifurcation, Fractal, Scene3DBridge, Basins, ParamScan2D, Continuation, IFS, LimitCycle };
  ActiveView active_view = ActiveView::Phase2D;
  bool show_side_panel = true;
  float window_toolbar_h = 32.0f;
  /* on-screen rect of the controls panel (window pixels), for cropping the
   * panel out of a saved PNG so only the plotted system is captured. */
  float panel_x0 = 0, panel_y0 = 0, panel_x1 = 0, panel_y1 = 0;

  std::vector<Point> points;
  std::vector<float> colors;
  std::deque<State> history;
  int history_limit = 4000;

  bool poincare_enabled = false;
  node_t *section_body = nullptr;
  node_t *section_x_body = nullptr;
  node_t *section_y_body = nullptr;
  SectionDirection section_direction = SectionDirection::Positive;
  std::vector<Point2> poincare_points;
  int poincare_limit = 20000;
  bool have_last_section = false;
  State last_section_state{};
  double last_section_value = 0.0;

  bool lyapunov_enabled = false;
  bool lyapunov_ready = false;
  State lyapunov_shadow{};
  double lyapunov_epsilon = 1e-6;
  double lyapunov_sum = 0.0;
  long long lyapunov_samples = 0;
  double lyapunov_estimate = 0.0;

  /* PHASE B: full Lyapunov spectrum + Kaplan-Yorke dimension. */
  std::vector<double> lyap_spectrum;     /* sorted descending */
  double lyap_kaplan_yorke = 0.0;
  double lyap_spectrum_sum = 0.0;
  bool lyap_spectrum_ready = false;
  std::string lyap_spectrum_msg;
  int lyap_spec_transient = 2000;
  int lyap_spec_steps = 20000;

  /* PHASE B/C: box-counting fractal dimension of the on-screen set. */
  bool boxdim_ready = false;
  double boxdim_value = 0.0;
  double boxdim_r2 = 0.0;
  long boxdim_n_points = 0;
  std::string boxdim_msg;

  /* PHASE D step 2 (foundation): limit-cycle period & amplitude of the
   * current trajectory (oscillating ODEs). */
  bool lc_ready = false;
  double lc_period = 0.0;
  double lc_amplitude = 0.0;
  std::string lc_msg;

  /* PHASE D step 2: limit-cycle continuation diagram. Sweep a parameter,
   * measure the periodic orbit's period & amplitude at each value, and plot
   * both curves vs the parameter (amplitude growing from zero marks a Hopf
   * bifurcation). */
  char lcc_param[64] = {0};
  double lcc_p_min = 0.0, lcc_p_max = 3.0;
  int lcc_slices = 80;
  bool lcc_has_data = false;
  bool lcc_autorun_done = false;
  std::string lcc_msg;
  std::vector<double> lcc_pp;       /* parameter value per slice */
  std::vector<double> lcc_period;   /* measured period (NaN where no cycle) */
  std::vector<double> lcc_amp;      /* measured amplitude (0 where no cycle) */
  double lcc_amp_max = 1.0, lcc_period_max = 1.0; /* for axis scaling */

  /* PHASE C: IFS / chaos game view. A built-in gallery of iterated function
   * systems (fern, Sierpinski, dragon, tree), rendered via the chaos game
   * into a density texture (reusing the fractal/basin texture pipeline). */
  int ifs_selected = 0;          /* index into the built-in IFS gallery */
  long ifs_iterations = 400000;
  GLuint ifs_tex = 0;
  int ifs_tex_w = 0, ifs_tex_h = 0;
  bool ifs_dirty = true;
  double ifs_box_dim = 0.0;      /* measured dimension of the current attractor */
  bool ifs_box_dim_ready = false;

  char export_path[256] = "dynsys_export.csv";
  char bif_param[64] = "rho";
  char bif_observable[128] = "x";
  double bif_start = 0.0;
  double bif_end = 100.0;
  int bif_slices = 160;
  int bif_discard = 1000;
  int bif_keep = 80;
  std::vector<Point2> bifurcation_points;
  /* PHASE7: Lyapunov exponent vs parameter, computed on the same sweep as
   * the orbit diagram so the two line up. */
  std::vector<Point2> bifurcation_lyapunov;
  bool bif_compute_lyapunov = true;
  bool bif_show = false; /* show the bifurcation diagram view */
  /* PHASE B+: interactive view window for the bifurcation diagram. */
  bool bif_view_valid = false;
  bool bif_autorun_done = false; /* auto-run the sweep once when the view is first shown */
  /* per-slice detected period (1,2,4,...; 0 = chaotic/aperiodic) for maps,
   * used to annotate windows on the diagram */
  std::vector<Point2> bifurcation_period; /* (param, period) */
  bool bif_show_period = true;
  double bif_view_xmin = 0, bif_view_xmax = 1, bif_view_ymin = 0, bif_view_ymax = 1;

  /* PHASE C: escape-time fractal view. Iterates the system's 2D map over
   * the plane; each pixel is either an initial condition (Julia-type) or
   * a parameter value (Mandelbrot-type), colored by escape time. This is
   * the explicit bridge between dynamical systems and fractals:
   *   - Mandelbrot = parameter map of z->z^2+c (orbit of 0 stays bounded)
   *   - Julia      = same map, fixed c, sweeping the initial point
   *   - logistic bifurcation = the real-axis analog of the Mandelbrot set
   */
  enum class FractalMode { ParameterSpace /*Mandelbrot-like*/, StateSpace /*Julia-like*/, Buddhabrot /*trajectory density*/ };
  FractalMode fractal_mode = FractalMode::ParameterSpace;
  double fractal_xmin = -2.5, fractal_xmax = 1.0;   /* view window (re) */
  double fractal_ymin = -1.5, fractal_ymax = 1.5;   /* view window (im) */
  /* Buddhabrot accumulation (trajectory-density rendering) */
  std::vector<uint32_t> buddha_accum;
  int buddha_w = 0, buddha_h = 0;
  uint32_t buddha_max = 1u;
  long buddha_samples = 0;
  uint32_t buddha_rng = 22695477u;
  int fractal_max_iter = 200;
  double fractal_escape_r = 4.0;
  bool fractal_smooth = true;
  int fractal_param_cx_index = 0;  /* which params are the two plane axes */
  int fractal_param_cy_index = 1;  /* (ParameterSpace mode) */
  /* Julia constant for StateSpace mode = current values of those params. */
  bool fractal_dirty = true;       /* recompute the image when true */
  int fractal_settle = 0;          /* debounce: frames to wait after a view change */
  int fractal_prog_level = 0;      /* progressive refine: downscale factor (0 = done) */
  GLuint fractal_tex = 0;
  int fractal_tex_w = 0, fractal_tex_h = 0;
  bool fractal_view_init = false;

  /* PHASE C+: the "3D bridge" scene — the Mandelbrot set lying in its
   * plane with the period-doubling bifurcation diagram rising out of it
   * along the real axis, making the conjugacy between z^2+c and the
   * logistic map visible in one picture. Its own GL buffers (drawn as
   * GL_POINTS), independent of the live-orbit line strip. */
  GLuint bridge_vbo = 0, bridge_cbo = 0, bridge_vao = 0;
  int bridge_point_count = 0;
  bool bridge_built = false;
  int bridge_mandel_res = 420;   /* c-plane sampling resolution */
  int bridge_bif_slices = 900;   /* bifurcation parameter slices (denser curtain) */
  int bridge_bif_keep = 200;     /* attractor points kept per slice */
  bool bridge_color_by_period = true; /* color bulbs/cascade by cycle period */
  float bridge_height = 18.0f;   /* vertical scale of the bifurcation */
  bool bridge_show_mandelbrot = true;
  bool bridge_show_bifurcation = true;
  /* The 3D bridge has two modes. CLASSIC is the famous logistic/Mandelbrot
   * bridge (the quadratic family's cascade rising out of the Mandelbrot set) —
   * a special correspondence that exists ONLY for the quadratic family.
   * CURRENT_SYSTEM lifts the CURRENT system's own bifurcation diagram into 3D,
   * sweeping bif_param along x, the attractor value as height (y), and a SECOND
   * parameter (bridge_param2) along z, so any system's bifurcation structure
   * can be explored as a 3D surface of stacked diagrams. */
  enum class BridgeMode { Classic, CurrentSystem, ProjectionSolid };
  BridgeMode bridge_mode = BridgeMode::ProjectionSolid;
  /* Which polynomial family the projection-solid iterates: z^2+c (the classic
   * Mandelbrot / logistic pairing) or z^3+c (the cubic Mandelbrot / cubic-map
   * pairing). Both give one object whose footprint is that family's
   * connectedness set and whose real-axis silhouette is its bifurcation
   * cascade -- because the real map is the real restriction of the complex one. */
  enum class BridgeFamily { Quadratic, Cubic, Sine };
  BridgeFamily bridge_family = BridgeFamily::Quadratic;
  char bridge_param2[64] = "";     /* second swept parameter (z axis) */
  double bridge_p2_min = 0.0, bridge_p2_max = 1.0;
  int bridge_p2_slices = 24;       /* number of stacked diagrams along z */

  /* PHASE B/C: basins of attraction. Integrate a grid of initial
   * conditions over the current 2D view and color each pixel by which
   * attractor it converges to — the basin-of-attraction fractal. (A Julia
   * set is the escape-to-infinity special case, already in the fractal
   * view; this covers attractors at finite points / cycles, e.g. Newton
   * fractals and forced-oscillator basins.) */
  GLuint basin_tex = 0;
  int basin_tex_w = 0, basin_tex_h = 0;
  bool basin_dirty = true;
  int basin_settle = 0;
  int basin_prog_level = 0;       /* progressive refine downscale (0 = done), like the fractal */
  int basin_steps = 1500;        /* integration steps per cell */
  int basin_res = 1;             /* unused placeholder for future supersampling */
  double basin_cluster_tol = 1e-2;
  bool basin_shade_speed = true; /* modulate brightness by convergence speed */
  int basin_attractor_count = 0;
  long basin_n_converged = 0, basin_n_diverged = 0, basin_n_nonconvergent = 0;

  /* PHASE D: equilibrium continuation (the MatCont-style bifurcation
   * diagram) — drawing/view state. The sweep scalars (cont_param,
   * cont_p_min/max, cont_max_points, cont_branch) are declared with the
   * older continuation block below; here we add the flattened branch used
   * by the full-window diagram view. */
  int cont_y_index = 0;             /* which state component on the y-axis */
  bool cont_has_branch = false;
  bool cont_autorun_done = false;   /* auto-run the trace once when the view is first shown */
  bool cont_view_valid = false;     /* false => refit axes to the branch */
  double cont_view_xmin = 0, cont_view_xmax = 1, cont_view_ymin = 0, cont_view_ymax = 1;
  std::string cont_message;
  int cont_n_fold = 0, cont_n_hopf = 0;
  /* the traced branch, flattened for drawing: parallel arrays */
  std::vector<double> cont_pp;            /* parameter value per point */
  std::vector<double> cont_yy;            /* observable (state[cont_y_index]) per point */
  std::vector<unsigned char> cont_stable; /* 1 = stable, 0 = unstable */
  std::vector<int> cont_special;          /* 0 none, 1 fold, 2 hopf (per point) */
  std::vector<std::size_t> cont_break;    /* indices where the polyline should break (dir join) */

  /* PHASE B: 2-parameter scan ("shrimp" map). Color a grid over two
   * parameters by the largest Lyapunov exponent (chaotic regions bright,
   * periodic windows dark) — the shrimp-shaped periodic structures in
   * parameter space. Reuses the renormalized-shadow-orbit estimate. */
  GLuint scan_tex = 0;
  int scan_tex_w = 0, scan_tex_h = 0;
  bool scan_dirty = true;
  int scan_settle = 0;
  int scan_prog_level = 0;        /* progressive downscale (0 = done/idle), like basins */
  int scan_px_index = 0;     /* parameter on the horizontal axis */
  int scan_py_index = 1;     /* parameter on the vertical axis */
  double scan_xmin = 0, scan_xmax = 4;   /* param-x range */
  double scan_ymin = 0, scan_ymax = 4;   /* param-y range */
  int scan_transient = 200;
  int scan_iterations = 200;
  bool scan_view_init = false;
  double scan_lyap_min = 0, scan_lyap_max = 0; /* observed range, for the legend */

  bool fixed_ready = false;
  State fixed_point{};
  double fixed_residual = 0.0;
  std::vector<double> fixed_jacobian;  std::string fixed_classification;
  std::vector<std::pair<double, double>> fixed_eigenvalues;
  std::vector<PhaseTrajectory> separatrix_curves;
  int separatrix_steps = 2500;
  double separatrix_epsilon = 1e-4;

  /* PHASE1/2: general (N-D) equilibrium analysis + equilibrium
   * continuation. fixed_eigenvalues above is kept for the legacy 2D
   * panel; this richer classification works in any dimension. */
  dynsys::analysis::Classification fixed_general{};
  bool fixed_general_ready = false;
  /* Hopf first Lyapunov coefficient (normal-form criticality) at the located
   * equilibrium, computed on demand. l1<0 supercritical, l1>0 subcritical. */
  bool hopf_l1_ready = false;
  double hopf_l1 = 0.0, hopf_omega = 0.0;
  std::string hopf_l1_msg;
  /* Fold (limit point) normal-form coefficient at the located equilibrium. */
  bool fold_a_ready = false;
  double fold_a = 0.0, fold_lambda0 = 0.0;
  std::string fold_a_msg;
  /* Exact equilibria via the CAS (solve-poly / Groebner). */
  bool cas_equi_ready = false;
  bool cas_equi_is_poly = false;
  std::vector<std::string> cas_equi_lines;
  std::string cas_equi_msg;
  /* Two-parameter continuation of a fold or Hopf curve in a (p,q) plane. */
  dynsys::analysis::TwoParamCurve twopar_curve{};
  bool twopar_ready = false;
  char twopar_p2[64] = "";       /* the second parameter name (q axis) */
  std::string twopar_msg;
  /* Collocation-based periodic-orbit continuation (the BVP solver: finds the
   * cycle AND its period exactly, follows unstable cycles and folds of cycles,
   * unlike the simulate-and-measure lcc_* sweep). */
  dynsys::analysis::CycleBranch cyc_branch{};
  bool cyc_ready = false;
  std::string cyc_msg;
  /* Two-parameter fold-of-cycles (LPC) curve in the (cont_param, twopar_p2) plane. */
  dynsys::analysis::LPCCurve lpc_curve_data{};
  bool lpc_ready = false;
  std::string lpc_msg;
  /* Collocation periodic-orbit continuation (period/amplitude vs parameter). */
  dynsys::analysis::CycleBranch cycle_branch{};
  bool cycle_ready = false;
  std::string cycle_msg;
  /* PHASE E: cached exact eigenvalue report from the Sangaku CAS for the
   * located equilibrium (computed on demand via the button, not every frame).
   * cas_eig_requested guards recomputation; cas_eig_exact records whether the
   * Jacobian rationalized exactly (so an approximate result is labelled). */
  dynsys::cas::EigenReport cas_eig{};
  bool cas_eig_ready = false;
  bool cas_eig_exact = false;
  char cont_param[64] = "rho";
  double cont_p_min = -10.0;
  double cont_p_max = 110.0;
  double cont_h0 = 0.05;
  int cont_max_points = 600;
  int cont_direction = 1;
  bool cont_detect_fold = true;
  bool cont_detect_hopf = true;
  dynsys::analysis::Branch cont_branch{};
  /* A second branch produced by branch switching at a branch point, drawn
   * alongside the primary branch. */
  dynsys::analysis::Branch cont_switched_branch{};
  bool cont_switched_ready = false;
  bool cont_ready = false;

  /* Homoclinic-orbit solve (truncated BVP to a saddle). Stored for drawing in
   * the Continuation view and reporting the half-time / amplitude. */
  dynsys::analysis::HomoclinicResult homoclinic{};
  bool homoclinic_ready = false;
  std::string homoclinic_msg;
  /* two-parameter homoclinic locus (codim-1 curve in (cont_param, twopar_q)) */
  dynsys::analysis::HomoclinicCurve homoclinic_curve{};
  bool homoclinic_curve_ready = false;

  /* PHASE6-UI: auto-scanned fixed points for the phase plane (pplane
   * style). Recomputed when the view/params change, not every frame. */
  std::vector<dynsys::analysis::FixedPoint2D> phase_fixed_points;
  bool phase_auto_fixed_points = true;
  bool phase_show_manifolds = true;
  std::vector<float> custom_orbit_ic; /* typed initial conditions for "add orbit" */
  double phase_fp_scan_xmin = 0, phase_fp_scan_xmax = 0;
  double phase_fp_scan_ymin = 0, phase_fp_scan_ymax = 0;
  int phase_fp_scan_px = -1, phase_fp_scan_py = -1; /* axis pair when scanned */
  double phase_fp_scan_param_sig = 0; /* sum of params, cheap change check */

  KeyState keys[GLFW_KEY_LAST + 1]{};
  bool mouse_left_pressed = false;
  double last_mouse_x = 0.0;
  double last_mouse_y = 0.0;
};

void update_projection(AppState &app);

constexpr const char *kGlslVersion = "#version 330";

const char *kVertexShaderSrc =
    "#version 330 core\n"
    "layout(location = 0) in vec3 a_pos;\n"
    "layout(location = 1) in vec3 a_color;\n"
    "out vec3 our_color;\n"
    "uniform mat4 model;\n"
    "uniform mat4 view;\n"
    "uniform mat4 projection;\n"
    "void main()\n"
    "{\n"
    "    gl_Position = projection * view * model * vec4(a_pos, 1.0);\n"
    "    our_color = a_color;\n"
    "}\n";

const char *kFragmentShaderSrc =
    "#version 330 core\n"
    "in vec3 our_color;\n"
    "out vec4 frag_color;\n"
    "void main()\n"
    "{\n"
    "    frag_color = vec4(our_color, 1.0);\n"
    "}\n";

/* Mouse-wheel deltas vary wildly by device, and on a slow frame (e.g. while a
 * basin or param-scan is recomputing) ImGui can hand us a large ACCUMULATED
 * io.MouseWheel for one gesture. Feeding that into pow(0.8, wheel) yields an
 * enormous one-shot zoom -- the "scrolled a little, zoomed out a lot" bug.
 * Clamp every wheel notch to a small predictable step so zoom is gradual. */
static inline double clamped_wheel(double w) {
  if (w > 1.0) w = 1.0;
  if (w < -1.0) w = -1.0;
  return w;
}

void set_vec3(vec3 v, float x, float y, float z) {
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

bool parse_number_literal(const char *text, double *out) {
  if (text == nullptr || *text == '\0') {
    return false;
  }
  errno = 0;
  char *end = nullptr;
  const double value = std::strtod(text, &end);
  if (errno != 0 || end == text || *end != '\0') {
    return false;
  }
  *out = value;
  return true;
}

void set_error(char *err, size_t err_cap, const char *fmt, ...) {
  if (err == nullptr || err_cap == 0) {
    return;
  }
  va_list args;
  va_start(args, fmt);
  std::vsnprintf(err, err_cap, fmt, args);
  va_end(args);
}

std::string trim_copy(const std::string &s) {
  size_t first = 0;
  while (first < s.size() && std::isspace(static_cast<unsigned char>(s[first]))) {
    ++first;
  }
  size_t last = s.size();
  while (last > first && std::isspace(static_cast<unsigned char>(s[last - 1]))) {
    --last;
  }
  return s.substr(first, last - first);
}

bool starts_with(const std::string &s, const char *prefix) {
  const size_t n = std::strlen(prefix);
  return s.size() >= n && s.compare(0, n, prefix) == 0;
}

bool is_ident_start(char c) {
  return std::isalpha(static_cast<unsigned char>(c)) || c == '_';
}

bool is_ident_char(char c) {
  return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
}

bool is_identifier(const std::string &s) {
  if (s.empty() || !is_ident_start(s[0])) {
    return false;
  }
  for (size_t i = 1; i < s.size(); ++i) {
    if (!is_ident_char(s[i])) {
      return false;
    }
  }
  return true;
}

std::string strip_line_comment(const std::string &line) {
  size_t cut = std::string::npos;
  const size_t hash = line.find('#');
  if (hash != std::string::npos) {
    cut = hash;
  }
  const size_t slash = line.find("//");
  if (slash != std::string::npos) {
    cut = std::min(cut, slash);
  }
  if (cut == std::string::npos) {
    return line;
  }
  return line.substr(0, cut);
}

std::vector<std::string> split_lines(const char *text) {
  std::vector<std::string> lines;
  const char *start = text != nullptr ? text : "";
  const char *p = start;
  while (*p != '\0') {
    if (*p == '\n') {
      lines.emplace_back(start, static_cast<size_t>(p - start));
      start = p + 1;
    }
    ++p;
  }
  lines.emplace_back(start, static_cast<size_t>(p - start));
  return lines;
}

bool is_builtin_function_name(const std::string &name);

std::vector<std::string> split_commas(const std::string &text) {
  std::vector<std::string> out;
  size_t pos = 0;
  while (pos <= text.size()) {
    const size_t comma = text.find(',', pos);
    const size_t end = comma == std::string::npos ? text.size() : comma;
    out.push_back(trim_copy(text.substr(pos, end - pos)));
    if (comma == std::string::npos) {
      break;
    }
    pos = comma + 1;
  }
  return out;
}

bool parse_state_name_list(const std::string &text, std::vector<std::string> *out, std::string *error) {
  out->clear();
  for (const std::string &name : split_commas(text)) {
    if (name.empty()) {
      *error = "empty state-variable name";
      return false;
    }
    if (!is_identifier(name)) {
      *error = "invalid state-variable name: '" + name + "'";
      return false;
    }
    if (name == "t" || is_builtin_function_name(name)) {
      *error = "reserved state-variable name: '" + name + "'";
      return false;
    }
    if (std::find(out->begin(), out->end(), name) != out->end()) {
      *error = "duplicate state-variable name: '" + name + "'";
      return false;
    }
    out->push_back(name);
  }
  if (out->empty()) {
    *error = "state declaration must list at least one variable";
    return false;
  }
  return true;
}

bool parse_lhs(const std::string &lhs_raw, std::string *name,
               std::vector<std::string> *params, std::string *error) {
  const std::string lhs = trim_copy(lhs_raw);

  // Generic ODE derivative forms: dx/dt, dtheta/dt, dx, dtheta.
  if (starts_with(lhs, "d")) {
    std::string var;
    if (lhs.size() > 4 && lhs.compare(lhs.size() - 3, 3, "/dt") == 0) {
      var = lhs.substr(1, lhs.size() - 4);
    } else {
      var = lhs.substr(1);
    }
    if (!var.empty() && is_identifier(var)) {
      *name = "d" + var;
      params->clear();
      return true;
    }
  }

  // Generic map forms: x', theta', x_next, theta_next, next x.
  if (!lhs.empty() && lhs.back() == '\'') {
    std::string var = lhs.substr(0, lhs.size() - 1);
    if (is_identifier(var)) {
      *name = var + "_next";
      params->clear();
      return true;
    }
  }
  if (starts_with(lhs, "next ")) {
    std::string var = trim_copy(lhs.substr(5));
    if (is_identifier(var)) {
      *name = var + "_next";
      params->clear();
      return true;
    }
  }
  if (lhs.size() > 5 && lhs.compare(lhs.size() - 5, 5, "_next") == 0) {
    std::string var = lhs.substr(0, lhs.size() - 5);
    if (is_identifier(var)) {
      *name = lhs;
      params->clear();
      return true;
    }
  }

  const size_t open = lhs.find('(');
  if (open == std::string::npos) {
    if (!is_identifier(lhs)) {
      *error = "invalid definition name: '" + lhs + "'";
      return false;
    }
    *name = lhs;
    params->clear();
    return true;
  }

  const size_t close = lhs.rfind(')');
  if (close == std::string::npos || close != lhs.size() - 1) {
    *error = "function definition LHS must look like name(arg1, arg2)";
    return false;
  }

  const std::string function_name = trim_copy(lhs.substr(0, open));
  if (!is_identifier(function_name)) {
    *error = "invalid function name: '" + function_name + "'";
    return false;
  }

  std::vector<std::string> parsed_params;
  for (const std::string &param : split_commas(lhs.substr(open + 1, close - open - 1))) {
    if (param.empty()) {
      *error = "empty parameter in function definition";
      return false;
    }
    if (!is_identifier(param)) {
      *error = "invalid function parameter: '" + param + "'";
      return false;
    }
    if (std::find(parsed_params.begin(), parsed_params.end(), param) != parsed_params.end()) {
      *error = "duplicate function parameter: '" + param + "'";
      return false;
    }
    parsed_params.push_back(param);
  }

  *name = function_name;
  *params = std::move(parsed_params);
  return true;
}

int state_index_for_name(const std::vector<std::string> &state_names, const std::string &name) {
  for (size_t i = 0; i < state_names.size(); ++i) {
    if (state_names[i] == name) return static_cast<int>(i);
  }
  return -1;
}

[[maybe_unused]] int state_index_for_name(const AppState &app, const std::string &name) {
  return state_index_for_name(app.state_names, name);
}

int derivative_index_for_name(const std::vector<std::string> &state_names, const std::string &name) {
  if (name.size() < 2 || name[0] != 'd') return -1;
  return state_index_for_name(state_names, name.substr(1));
}

int next_index_for_name(const std::vector<std::string> &state_names, const std::string &name) {
  constexpr const char *suffix = "_next";
  constexpr size_t suffix_len = 5;
  if (name.size() <= suffix_len) return -1;
  if (name.compare(name.size() - suffix_len, suffix_len, suffix) != 0) return -1;
  return state_index_for_name(state_names, name.substr(0, name.size() - suffix_len));
}

bool is_reserved_value_name(const std::string &name) {
  /* Only `t` (time) is structurally reserved. `pi` and `e` are math constants
   * that a user symbol may shadow (the evaluator checks params/definitions
   * before falling back to the constants), so they are allowed as parameter
   * or state names — useful since `e` is a common attractor parameter. */
  return name == "t";
}

bool is_state_value_name(const std::vector<std::string> &state_names, const std::string &name) {
  return state_index_for_name(state_names, name) >= 0;
}

bool is_builtin_function_name(const std::string &name) {
  return name == "sin" || name == "cos" || name == "tan" ||
         name == "asin" || name == "acos" || name == "atan" ||
         name == "exp" || name == "log" || name == "log10" ||
         name == "sqrt" || name == "abs" || name == "floor" ||
         name == "ceil" || name == "sign" || name == "pow" ||
         name == "min" || name == "max" || name == "mod" ||
         name == "clamp" || name == "select";
}

const AppState::Definition *find_definition(const AppState &app,
                                            const std::string &name,
                                            size_t arity,
                                            size_t *index_out = nullptr) {
  for (size_t i = 0; i < app.definitions.size(); ++i) {
    const auto &def = app.definitions[i];
    if (def.name == name && def.params.size() == arity) {
      if (index_out != nullptr) {
        *index_out = i;
      }
      return &def;
    }
  }
  return nullptr;
}

AppState::Param *find_param(AppState &app, const std::string &name) {
  for (auto &param : app.params) {
    if (param.name == name) {
      return &param;
    }
  }
  return nullptr;
}

[[maybe_unused]] const AppState::Param *find_param(const AppState &app, const std::string &name) {
  for (const auto &param : app.params) {
    if (param.name == name) {
      return &param;
    }
  }
  return nullptr;
}

[[maybe_unused]] const AppState::Observable *find_observable(const AppState &app, const std::string &name) {
  for (const auto &obs : app.observables) {
    if (obs.name == name) {
      return &obs;
    }
  }
  return nullptr;
}

struct EvalContext {
  AppState &app;
  State state;
  std::vector<std::pair<std::string, double>> locals;
  std::vector<int> active;
  std::vector<bool> cached;
  std::vector<double> cache;
  int depth = 0;
};

bool eval_ast(EvalContext &ctx, const node_t *ast, double *out, char *err,
              size_t err_cap);

bool eval_builtin_call(EvalContext &ctx, const char *name, size_t argc,
                       node_t **args, double *out, char *err, size_t err_cap) {
  double values[3] = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < argc && i < 3; ++i) {
    if (!eval_ast(ctx, args[i], &values[i], err, err_cap)) {
      return false;
    }
  }

  if (argc == 1) {
    const double x = values[0];
    if (std::strcmp(name, "sin") == 0) { *out = std::sin(x); return true; }
    if (std::strcmp(name, "cos") == 0) { *out = std::cos(x); return true; }
    if (std::strcmp(name, "tan") == 0) { *out = std::tan(x); return true; }
    if (std::strcmp(name, "asin") == 0) { *out = std::asin(x); return true; }
    if (std::strcmp(name, "acos") == 0) { *out = std::acos(x); return true; }
    if (std::strcmp(name, "atan") == 0) { *out = std::atan(x); return true; }
    if (std::strcmp(name, "exp") == 0) { *out = std::exp(x); return true; }
    if (std::strcmp(name, "log") == 0) { *out = std::log(x); return true; }
    if (std::strcmp(name, "log10") == 0) { *out = std::log10(x); return true; }
    if (std::strcmp(name, "sqrt") == 0) { *out = std::sqrt(x); return true; }
    if (std::strcmp(name, "abs") == 0) { *out = std::fabs(x); return true; }
    if (std::strcmp(name, "floor") == 0) { *out = std::floor(x); return true; }
    if (std::strcmp(name, "ceil") == 0) { *out = std::ceil(x); return true; }
    if (std::strcmp(name, "sign") == 0) { *out = (x > 0.0) - (x < 0.0); return true; }
  }

  if (argc == 2) {
    const double a = values[0];
    const double b = values[1];
    if (std::strcmp(name, "pow") == 0) { *out = std::pow(a, b); return true; }
    if (std::strcmp(name, "min") == 0) { *out = std::fmin(a, b); return true; }
    if (std::strcmp(name, "max") == 0) { *out = std::fmax(a, b); return true; }
    if (std::strcmp(name, "mod") == 0) { *out = std::fmod(a, b); return true; }
  }

  if (argc == 3) {
    const double a = values[0];
    const double b = values[1];
    const double c = values[2];
    if (std::strcmp(name, "clamp") == 0) {
      *out = std::fmax(b, std::fmin(c, a));
      return true;
    }
    if (std::strcmp(name, "select") == 0) {
      *out = (a != 0.0) ? b : c;
      return true;
    }
  }

  return false;
}

bool eval_user_definition(EvalContext &ctx, size_t def_index, double *out,
                          char *err, size_t err_cap) {
  if (def_index >= ctx.app.definitions.size()) {
    set_error(err, err_cap, "internal definition index error");
    return false;
  }
  const auto &def = ctx.app.definitions[def_index];
  if (!def.params.empty()) {
    set_error(err, err_cap, "function '%s' needs %zu argument(s)",
              def.name.c_str(), def.params.size());
    return false;
  }
  if (ctx.cached[def_index]) {
    *out = ctx.cache[def_index];
    return true;
  }
  if (ctx.active[def_index] != 0) {
    set_error(err, err_cap, "cyclic definition involving '%s'", def.name.c_str());
    return false;
  }
  ctx.active[def_index] = 1;
  double value = 0.0;
  const bool ok = eval_ast(ctx, def.body, &value, err, err_cap);
  ctx.active[def_index] = 0;
  if (!ok) {
    return false;
  }
  ctx.cached[def_index] = true;
  ctx.cache[def_index] = value;
  *out = value;
  return true;
}

bool eval_user_call(EvalContext &ctx, const char *name, size_t argc,
                    node_t **args, double *out, char *err, size_t err_cap) {
  size_t def_index = 0;
  const auto *def = find_definition(ctx.app, name, argc, &def_index);
  if (def == nullptr || def->params.empty()) {
    return false;
  }
  if (ctx.depth > 64) {
    set_error(err, err_cap, "call depth exceeded while calling '%s'", name);
    return false;
  }
  if (ctx.active[def_index] != 0) {
    set_error(err, err_cap, "recursive function call involving '%s'", name);
    return false;
  }

  std::vector<double> values(argc, 0.0);
  for (size_t i = 0; i < argc; ++i) {
    if (!eval_ast(ctx, args[i], &values[i], err, err_cap)) {
      return false;
    }
  }

  const size_t local_base = ctx.locals.size();
  for (size_t i = 0; i < argc; ++i) {
    ctx.locals.emplace_back(def->params[i], values[i]);
  }

  ctx.active[def_index] = 1;
  ctx.depth += 1;
  const bool ok = eval_ast(ctx, def->body, out, err, err_cap);
  ctx.depth -= 1;
  ctx.active[def_index] = 0;
  ctx.locals.resize(local_base);
  return ok;
}

bool eval_call(EvalContext &ctx, const char *name, size_t argc, node_t **args,
               double *out, char *err, size_t err_cap) {
  if (eval_builtin_call(ctx, name, argc, args, out, err, err_cap)) {
    return true;
  }
  if (eval_user_call(ctx, name, argc, args, out, err, err_cap)) {
    return true;
  }
  set_error(err, err_cap, "unsupported function call: %s/%zu", name, argc);
  return false;
}

bool eval_ast(EvalContext &ctx, const node_t *ast, double *out, char *err,
              size_t err_cap) {
  if (ast == nullptr) {
    set_error(err, err_cap, "null AST");
    return false;
  }

  switch (ast->kind) {
  case NODE_CONST: {
    double value = 0.0;
    if (parse_number_literal(ast->cnst.name, &value)) {
      *out = value;
      return true;
    }
    set_error(err, err_cap, "unsupported constant: %s", ast->cnst.name);
    return false;
  }

  case NODE_VAR: {
    const char *name = ast->var.name;
    for (auto it = ctx.locals.rbegin(); it != ctx.locals.rend(); ++it) {
      if (it->first == name) {
        *out = it->second;
        return true;
      }
    }
    for (size_t i = 0; i < ctx.app.state_names.size(); ++i) {
      if (ctx.app.state_names[i] == name) {
        *out = state_at(ctx.state, i);
        return true;
      }
    }
    if (std::strcmp(name, "t") == 0) { *out = ctx.state.t; return true; }

    /* User-declared parameters and definitions take precedence over the math
     * constants pi and e, so a system may legitimately name a parameter `e`
     * (a common attractor parameter) without colliding with Euler's number.
     * `t` above stays reserved (it's the structural time variable). */
    const AppState::Param *param = find_param(ctx.app, name);
    if (param != nullptr) {
      *out = param->value;
      return true;
    }

    size_t def_index = 0;
    const auto *def = find_definition(ctx.app, name, 0, &def_index);
    if (def != nullptr) {
      if (def->derivative || def->next_equation) {
        set_error(err, err_cap, "equation '%s' cannot be used as a scalar dependency", name);
        return false;
      }
      return eval_user_definition(ctx, def_index, out, err, err_cap);
    }

    /* fall back to the math constants only if no user param/definition shadowed
     * them — so `pi`/`e` keep working in systems that don't redefine them. */
    if (std::strcmp(name, "pi") == 0) { *out = M_PI; return true; }
    if (std::strcmp(name, "e") == 0) { *out = M_E; return true; }

    set_error(err, err_cap, "unknown variable: %s", name);
    return false;
  }

  case NODE_APP:
    if (ast->app.head->kind == NODE_CONST && ast->app.head->cnst.op != nullptr) {
      const op_info_t *op = ast->app.head->cnst.op;
      if (ast->app.argc == 2) {
        double a = 0.0;
        double b = 0.0;
        if (!eval_ast(ctx, ast->app.args[0], &a, err, err_cap) ||
            !eval_ast(ctx, ast->app.args[1], &b, err, err_cap)) {
          return false;
        }
        if (op == &OP_ADD) { *out = a + b; return true; }
        if (op == &OP_SUB) { *out = a - b; return true; }
        if (op == &OP_MUL) { *out = a * b; return true; }
        if (op == &OP_DIV) { *out = a / b; return true; }
      }
      set_error(err, err_cap, "unsupported operator application: %s/%zu",
                op->name, ast->app.argc);
      return false;
    }

    if (ast->app.head->kind == NODE_VAR) {
      return eval_call(ctx, ast->app.head->var.name, ast->app.argc,
                       ast->app.args, out, err, err_cap);
    }

    set_error(err, err_cap, "unsupported application near byte range [%zu,%zu)",
              ast->span.start, ast->span.end);
    return false;

  case NODE_LAMBDA:
  case NODE_FORALL:
  case NODE_EXISTS:
    set_error(err, err_cap, "binders are not valid inside numeric systems");
    return false;
  }

  set_error(err, err_cap, "unsupported AST node");
  return false;
}

bool eval_expr_at(AppState &app, node_t *expr, const State &state, double *out,
                  char *err, size_t err_cap) {
  EvalContext ctx{app, State{}, {}, {}, {}, {}};
  ctx.state = state;
  ctx.active.assign(app.definitions.size(), 0);
  ctx.cached.assign(app.definitions.size(), false);
  ctx.cache.assign(app.definitions.size(), 0.0);
  return eval_ast(ctx, expr, out, err, err_cap);
}

/* === IR-backed evaluation ===
 *
 * The hot integration loop uses these instead of eval_expr_at. The
 * cold paths (param defaults, initial-value expressions, the
 * AST-print panel) keep going through eval_expr_at because they run
 * once at compile time or once per redraw — the strcmp tax is
 * invisible there. */

void sync_param_values(AppState &app) {
  app.param_values.resize(app.params.size());
  for (size_t i = 0; i < app.params.size(); ++i) {
    app.param_values[i] = app.params[i].value;
  }
}

bool eval_program_at(AppState &app, const dynsys::ir::Program &prog,
                     const State &state, double *out,
                     char *err, size_t err_cap) {
  dynsys::ir::scratch_reset_eval(&app.eval_scratch);
  dynsys::ir::RunContext rc;
  rc.state    = state.v.data();
  rc.n_state  = state.v.size();
  rc.t        = state.t;
  rc.params   = app.param_values.data();
  rc.n_params = app.param_values.size();
  rc.defs     = app.definition_programs.data();
  rc.n_defs   = app.definition_programs.size();
  return dynsys::ir::run(prog, rc, app.eval_scratch, out, err, err_cap);
}

bool eval_rhs(AppState &app, const State &state, State *deriv, char *err,
              size_t err_cap) {
  const size_t dim = app.state_names.size();
  if (app.equation_programs.size() != dim) {
    set_error(err, err_cap, "ODE equations are not compiled");
    return false;
  }
  resize_state(*deriv, dim);
  for (size_t i = 0; i < dim; ++i) {
    double value = 0.0;
    if (app.use_ast_fallback) {
      if (!eval_expr_at(app, app.equations[i], state, &value, err, err_cap)) return false;
    } else {
      if (!eval_program_at(app, app.equation_programs[i], state, &value, err, err_cap)) return false;
    }
    set_state_at(*deriv, i, value);
  }
  deriv->t = 1.0;
  return true;
}

/* In-place scaled add: out = base + h * k. `out` is sized to `dim`
 * (resized only when its capacity is insufficient — usual case is
 * already-sized scratch). No allocations on the hot path. */
static void add_scaled_into(const State &base, const State &k, double h,
                            size_t dim, State *out) {
  resize_state(*out, dim);
  out->t = base.t + h * k.t;
  for (size_t i = 0; i < dim; ++i) {
    set_state_at(*out, i, state_at(base, i) + h * state_at(k, i));
  }
}

[[maybe_unused]] State add_scaled(const State &s, const State &k, double h) {
  const size_t dim = std::max(s.v.size(), k.v.size());
  State out = make_state_like(dim, s.t + h * k.t);
  for (size_t i = 0; i < dim; ++i) {
    set_state_at(out, i, state_at(s, i) + h * state_at(k, i));
  }
  return out;
}

/* PHASE B: adaptive embedded Runge-Kutta step from t to t+dt.
 * Subdivides [t, t+dt] into internal substeps whose size is chosen to keep
 * the local error estimate (difference of the embedded 4th/5th order
 * solutions) under adaptive_tol, clamped to [dt_min, dt_max]. Supports
 * RKF45 (Fehlberg) and DOPRI45 (Dormand-Prince) Butcher tableaux. */
bool step_ode_adaptive(AppState &app, const State &in, State *out, char *err, size_t err_cap) {
  const size_t dim = app.state_names.size();
  const bool dopri = (app.integrator == Integrator::DOPRI45);
  const double total = app.dt;               /* signed target output step */
  const double dir = total >= 0 ? 1.0 : -1.0;
  double remaining = std::fabs(total);
  const double tol = std::max(1e-12, app.adaptive_tol);
  const double hmin = std::max(1e-12, app.adaptive_dt_min);
  const double hmax = std::max(hmin, app.adaptive_dt_max);

  std::vector<double> y(dim), ytmp(dim), y4(dim), y5(dim);
  for (size_t i = 0; i < dim; ++i) y[i] = state_at(in, i);
  double t = in.t;

  /* Butcher c-nodes and stage combinations. We evaluate via eval_rhs on a
   * State wrapper. k has up to 7 stages (DOPRI). */
  std::vector<std::vector<double>> k(7, std::vector<double>(dim));
  State stage = in; resize_state(stage, dim);
  auto rhs_at = [&](const std::vector<double> &yv, double tt, std::vector<double> &kout) -> bool {
    for (size_t i = 0; i < dim; ++i) set_state_at(stage, i, yv[i]);
    stage.t = tt;
    State d{};
    if (!eval_rhs(app, stage, &d, err, err_cap)) return false;
    resize_state(d, dim);
    for (size_t i = 0; i < dim; ++i) kout[i] = state_at(d, i);
    return true;
  };

  double h = std::min(hmax, remaining > 0 ? remaining : hmax);
  if (h < hmin) h = hmin;
  int guard = 0;
  const int max_substeps = 100000;

  while (remaining > 1e-15 && guard++ < max_substeps) {
    if (h > remaining) h = remaining;
    const double hs = dir * h;
    bool ok = true;

    if (!dopri) {
      /* Runge-Kutta-Fehlberg 4(5) */
      ok = rhs_at(y, t, k[0]);
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(0.25*k[0][i]); ok=rhs_at(ytmp,t+0.25*hs,k[1]); }
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(3.0/32*k[0][i]+9.0/32*k[1][i]); ok=rhs_at(ytmp,t+3.0/8*hs,k[2]); }
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(1932.0/2197*k[0][i]-7200.0/2197*k[1][i]+7296.0/2197*k[2][i]); ok=rhs_at(ytmp,t+12.0/13*hs,k[3]); }
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(439.0/216*k[0][i]-8.0*k[1][i]+3680.0/513*k[2][i]-845.0/4104*k[3][i]); ok=rhs_at(ytmp,t+hs,k[4]); }
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(-8.0/27*k[0][i]+2.0*k[1][i]-3544.0/2565*k[2][i]+1859.0/4104*k[3][i]-11.0/40*k[4][i]); ok=rhs_at(ytmp,t+0.5*hs,k[5]); }
      if (ok) for (size_t i=0;i<dim;++i) {
        y4[i]=y[i]+hs*(25.0/216*k[0][i]+1408.0/2565*k[2][i]+2197.0/4104*k[3][i]-1.0/5*k[4][i]);
        y5[i]=y[i]+hs*(16.0/135*k[0][i]+6656.0/12825*k[2][i]+28561.0/56430*k[3][i]-9.0/50*k[4][i]+2.0/55*k[5][i]);
      }
    } else {
      /* Dormand-Prince 5(4) */
      ok = rhs_at(y, t, k[0]);
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(0.2*k[0][i]); ok=rhs_at(ytmp,t+0.2*hs,k[1]); }
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(3.0/40*k[0][i]+9.0/40*k[1][i]); ok=rhs_at(ytmp,t+0.3*hs,k[2]); }
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(44.0/45*k[0][i]-56.0/15*k[1][i]+32.0/9*k[2][i]); ok=rhs_at(ytmp,t+0.8*hs,k[3]); }
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(19372.0/6561*k[0][i]-25360.0/2187*k[1][i]+64448.0/6561*k[2][i]-212.0/729*k[3][i]); ok=rhs_at(ytmp,t+8.0/9*hs,k[4]); }
      if (ok) { for (size_t i=0;i<dim;++i) ytmp[i]=y[i]+hs*(9017.0/3168*k[0][i]-355.0/33*k[1][i]+46732.0/5247*k[2][i]+49.0/176*k[3][i]-5103.0/18656*k[4][i]); ok=rhs_at(ytmp,t+hs,k[5]); }
      if (ok) for (size_t i=0;i<dim;++i) y5[i]=y[i]+hs*(35.0/384*k[0][i]+500.0/1113*k[2][i]+125.0/192*k[3][i]-2187.0/6784*k[4][i]+11.0/84*k[5][i]);
      if (ok) ok=rhs_at(y5,t+hs,k[6]);
      if (ok) for (size_t i=0;i<dim;++i)
        y4[i]=y[i]+hs*(5179.0/57600*k[0][i]+7571.0/16695*k[2][i]+393.0/640*k[3][i]-92097.0/339200*k[4][i]+187.0/2100*k[5][i]+1.0/40*k[6][i]);
    }
    if (!ok) return false;

    /* error estimate (max norm of y5 - y4, relative) */
    double erot = 0.0;
    for (size_t i = 0; i < dim; ++i) {
      const double sc = 1.0 + std::max(std::fabs(y[i]), std::fabs(y5[i]));
      const double e = std::fabs(y5[i] - y4[i]) / sc;
      if (e > erot) erot = e;
    }
    /* accept if within tol or already at the minimum step */
    if (erot <= tol || h <= hmin * 1.0000001) {
      const std::vector<double> &ynew = dopri ? y5 : y5; /* use higher order */
      for (size_t i = 0; i < dim; ++i) y[i] = ynew[i];
      t += hs;
      remaining -= h;
      app.last_adaptive_dt = h;
      bool finite = true;
      for (double v : y) finite = finite && std::isfinite(v);
      if (!finite) { set_error(err, err_cap, "adaptive step diverged"); return false; }
    }
    /* PI-ish step update: h_new = h * clamp(0.9*(tol/err)^(1/5)) */
    double factor;
    if (erot <= 1e-300) factor = 5.0;
    else factor = 0.9 * std::pow(tol / erot, 0.2);
    factor = std::max(0.2, std::min(5.0, factor));
    h *= factor;
    if (h < hmin) h = hmin;
    if (h > hmax) h = hmax;
  }

  resize_state(*out, dim);
  out->t = t;
  for (size_t i = 0; i < dim; ++i) set_state_at(*out, i, y[i]);
  return true;
}

bool step_ode_state(AppState &app, const State &in, State *out, char *err,
                    size_t err_cap) {
  const size_t dim = app.state_names.size();
  if (app.integrator == Integrator::Euler) {
    if (!eval_rhs(app, in, &app.scratch_k1, err, err_cap)) return false;
    add_scaled_into(in, app.scratch_k1, app.dt, dim, out);
    out->t = in.t + app.dt;
    return true;
  }
  if (app.integrator == Integrator::RK2) {
    if (!eval_rhs(app, in, &app.scratch_k1, err, err_cap)) return false;
    add_scaled_into(in, app.scratch_k1, app.dt * 0.5, dim, &app.scratch_mid);
    if (!eval_rhs(app, app.scratch_mid, &app.scratch_k2, err, err_cap)) return false;
    add_scaled_into(in, app.scratch_k2, app.dt, dim, out);
    out->t = in.t + app.dt;
    return true;
  }
  if (app.integrator == Integrator::Heun) {
    /* explicit trapezoid: predictor Euler, corrector average of slopes */
    if (!eval_rhs(app, in, &app.scratch_k1, err, err_cap)) return false;
    add_scaled_into(in, app.scratch_k1, app.dt, dim, &app.scratch_mid);
    if (!eval_rhs(app, app.scratch_mid, &app.scratch_k2, err, err_cap)) return false;
    resize_state(*out, dim);
    out->t = in.t + app.dt;
    for (size_t i = 0; i < dim; ++i)
      set_state_at(*out, i, state_at(in, i) +
        app.dt * 0.5 * (state_at(app.scratch_k1, i) + state_at(app.scratch_k2, i)));
    return true;
  }
  if (app.integrator == Integrator::RK38) {
    /* classic 3/8 rule (4th order) */
    const double h = app.dt;
    if (!eval_rhs(app, in, &app.scratch_k1, err, err_cap)) return false;
    add_scaled_into(in, app.scratch_k1, h / 3.0, dim, &app.scratch_mid);
    if (!eval_rhs(app, app.scratch_mid, &app.scratch_k2, err, err_cap)) return false;
    /* mid = in + h*(-1/3 k1 + k2) */
    resize_state(app.scratch_mid, dim);
    for (size_t i = 0; i < dim; ++i)
      set_state_at(app.scratch_mid, i, state_at(in, i) + h * (-1.0/3.0 * state_at(app.scratch_k1, i) + state_at(app.scratch_k2, i)));
    if (!eval_rhs(app, app.scratch_mid, &app.scratch_k3, err, err_cap)) return false;
    for (size_t i = 0; i < dim; ++i)
      set_state_at(app.scratch_mid, i, state_at(in, i) + h * (state_at(app.scratch_k1, i) - state_at(app.scratch_k2, i) + state_at(app.scratch_k3, i)));
    if (!eval_rhs(app, app.scratch_mid, &app.scratch_k4, err, err_cap)) return false;
    resize_state(*out, dim);
    out->t = in.t + h;
    for (size_t i = 0; i < dim; ++i)
      set_state_at(*out, i, state_at(in, i) + h / 8.0 *
        (state_at(app.scratch_k1, i) + 3.0 * state_at(app.scratch_k2, i)
         + 3.0 * state_at(app.scratch_k3, i) + state_at(app.scratch_k4, i)));
    return true;
  }
  if (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45) {
    return step_ode_adaptive(app, in, out, err, err_cap);
  }

  /* RK4 — entirely in pre-sized scratch buffers. */
  if (!eval_rhs(app, in, &app.scratch_k1, err, err_cap)) return false;
  add_scaled_into(in, app.scratch_k1, app.dt * 0.5, dim, &app.scratch_mid);
  if (!eval_rhs(app, app.scratch_mid, &app.scratch_k2, err, err_cap)) return false;
  add_scaled_into(in, app.scratch_k2, app.dt * 0.5, dim, &app.scratch_mid);
  if (!eval_rhs(app, app.scratch_mid, &app.scratch_k3, err, err_cap)) return false;
  add_scaled_into(in, app.scratch_k3, app.dt,        dim, &app.scratch_mid);
  if (!eval_rhs(app, app.scratch_mid, &app.scratch_k4, err, err_cap)) return false;

  resize_state(*out, dim);
  out->t = in.t + app.dt;
  for (size_t i = 0; i < dim; ++i) {
    const double next = state_at(in, i) + app.dt *
      (state_at(app.scratch_k1, i)
       + 2.0 * state_at(app.scratch_k2, i)
       + 2.0 * state_at(app.scratch_k3, i)
       + state_at(app.scratch_k4, i)) / 6.0;
    set_state_at(*out, i, next);
  }
  return true;
}

bool step_map_state(AppState &app, const State &in, State *out, char *err,
                    size_t err_cap) {
  const size_t dim = app.state_names.size();
  if (app.next_equation_programs.size() != dim) {
    set_error(err, err_cap, "map equations are not compiled");
    return false;
  }
  resize_state(*out, dim);
  out->t = in.t + 1.0;
  for (size_t i = 0; i < dim; ++i) {
    double value = 0.0;
    if (!eval_program_at(app, app.next_equation_programs[i], in, &value, err, err_cap)) {
      return false;
    }
    set_state_at(*out, i, value);
  }
  return true;
}

bool step_state(AppState &app, const State &in, State *out, char *err,
                size_t err_cap) {
  return app.mode == SystemMode::Map ? step_map_state(app, in, out, err, err_cap)
                                     : step_ode_state(app, in, out, err, err_cap);
}

std::string format_parse_error(const char *label, const char *input,
                               parse_result_t result) {
  char buf[768];
  std::snprintf(buf, sizeof(buf),
                "%s parse error at line %zu, column %zu: %s\n%s\n%*s^",
                label, result.err_span.line, result.err_span.col,
                result.err_msg ? result.err_msg : "unknown parse error", input,
                result.err_span.col > 0 ? static_cast<int>(result.err_span.col - 1) : 0,
                "");
  return std::string(buf);
}

node_t *parse_expression_or_fail(arena_t *arena, const std::string &expr,
                                 const std::string &label, std::string *error) {
  parse_result_t result = parse(expr.c_str(), arena);
  if (!result.ok) {
    if (error != nullptr) {
      *error = format_parse_error(label.c_str(), expr.c_str(), result);
    }
    return nullptr;
  }
  return result.ast;
}

std::string extract_param_expr_and_range(const std::string &rhs_raw,
                                         double *min_out, double *max_out,
                                         bool *has_range) {
  std::string rhs = trim_copy(rhs_raw);
  *has_range = false;
  size_t expr_end = rhs.size();

  const size_t bracket = rhs.find('[');
  if (bracket != std::string::npos) {
    const size_t close = rhs.find(']', bracket + 1);
    if (close != std::string::npos) {
      const std::vector<std::string> nums = split_commas(rhs.substr(bracket + 1, close - bracket - 1));
      if (nums.size() == 2) {
        double lo = 0.0, hi = 0.0;
        if (parse_number_literal(nums[0].c_str(), &lo) && parse_number_literal(nums[1].c_str(), &hi)) {
          *min_out = lo;
          *max_out = hi;
          *has_range = true;
        }
      }
      expr_end = std::min(expr_end, bracket);
    }
  }

  const size_t min_pos = rhs.find(" min=");
  const size_t max_pos = rhs.find(" max=");
  if (min_pos != std::string::npos || max_pos != std::string::npos) {
    expr_end = std::min(expr_end, std::min(min_pos == std::string::npos ? rhs.size() : min_pos,
                                           max_pos == std::string::npos ? rhs.size() : max_pos));
    auto read_kv = [&](const char *key, double *dst) -> bool {
      const std::string needle = std::string(" ") + key + "=";
      size_t pos = rhs.find(needle);
      if (pos == std::string::npos) return false;
      pos += needle.size();
      size_t end = pos;
      while (end < rhs.size() && !std::isspace(static_cast<unsigned char>(rhs[end]))) ++end;
      return parse_number_literal(rhs.substr(pos, end - pos).c_str(), dst);
    };
    double lo = *min_out, hi = *max_out;
    const bool has_min = read_kv("min", &lo);
    const bool has_max = read_kv("max", &hi);
    if (has_min && has_max) {
      *min_out = lo;
      *max_out = hi;
      *has_range = true;
    }
  }

  return trim_copy(rhs.substr(0, expr_end));
}

int effective_dimension(const AppState &app);
int bridge_family_for_system(const AppState &app);

bool compile_system(AppState &app, const char *system_text, std::string *error) {
  /* PHASE0: remember the pre-compile state signature so we can decide,
   * after a successful compile, what analysis to keep vs. invalidate.
   * A located equilibrium is tied to the equations and is cleared on
   * any recompile; user-seeded orbits are kept as long as the state
   * variables (count + names) are unchanged, since their coordinates
   * still mean the same thing. */
  const std::vector<std::string> prev_state_names = app.state_names;
  std::vector<std::string> prev_param_names;
  prev_param_names.reserve(app.params.size());
  for (const auto &pp : app.params) prev_param_names.push_back(pp.name);

  const std::vector<std::string> lines = split_lines(system_text);

  /* PHASE C: IFS systems are parsed separately — they're a list of affine
   * maps run via the chaos game, with no state-vector equations. Detect
   * `mode = ifs` first and handle it here, so we don't run the ODE/map
   * equation machinery (which would reject the lack of x' = ... lines). */
  {
    bool is_ifs = false;
    for (const std::string &raw : lines) {
      const std::string ln = trim_copy(raw);
      if (ln.empty() || ln[0] == '#') continue;
      const size_t eq = ln.find('=');
      if (eq == std::string::npos) continue;
      const std::string lhs = trim_copy(ln.substr(0, eq));
      const std::string rhs = trim_copy(ln.substr(eq + 1));
      if ((lhs == "mode" || lhs == "type") && rhs == "ifs") { is_ifs = true; break; }
    }
    if (is_ifs) {
      /* Parse parameters (so an IFS can be parametrized) and affine maps
       * whose coefficients are EXPRESSIONS (so a coefficient can reference a
       * parameter, e.g. a rotation angle). Coefficients are evaluated to
       * numbers at render time against the current parameter values. */
      if (app.ifs_arena_init) { arena_destroy(&app.ifs_arena); app.ifs_arena_init = false; }
      arena_init(&app.ifs_arena, 256 * 1024);
      app.ifs_arena_init = true;

      std::vector<AppState::Param> ifs_params;
      std::vector<AppState::IFSMapExpr> map_exprs;
      bool ok = true;
      std::string local_err;

      auto parse_field = [&](const std::string &expr, const std::string &label) -> node_t * {
        node_t *n = parse_expression_or_fail(&app.ifs_arena, expr, label, &local_err);
        if (!n) ok = false;
        return n;
      };

      for (size_t i = 0; i < lines.size() && ok; ++i) {
        std::string ln = trim_copy(strip_line_comment(lines[i]));
        if (ln.empty()) continue;
        /* param name = value [min,max] — reuse the same surface syntax */
        if (starts_with(ln, "param ")) {
          const std::string rest = trim_copy(ln.substr(6));
          const size_t eq = rest.find('=');
          if (eq == std::string::npos) { local_err = "line " + std::to_string(i + 1) + ": expected 'param name = value [min,max]'"; ok = false; break; }
          std::string name = trim_copy(rest.substr(0, eq));
          std::string rhs = trim_copy(rest.substr(eq + 1));
          AppState::Param p; p.name = name; p.line = (int)i + 1;
          /* optional [min,max] */
          const size_t lb = rhs.find('[');
          if (lb != std::string::npos) {
            const size_t rb = rhs.find(']', lb);
            if (rb != std::string::npos) {
              std::string range = rhs.substr(lb + 1, rb - lb - 1);
              rhs = trim_copy(rhs.substr(0, lb));
              const size_t comma = range.find(',');
              if (comma != std::string::npos) {
                try { p.min_value = std::stod(trim_copy(range.substr(0, comma)));
                      p.max_value = std::stod(trim_copy(range.substr(comma + 1)));
                      p.has_range = true; } catch (...) {}
              }
            }
          }
          try { p.value = std::stod(trim_copy(rhs)); } catch (...) { local_err = "line " + std::to_string(i + 1) + ": IFS param value must be a number"; ok = false; break; }
          p.default_value = p.value;
          ifs_params.push_back(p);
          continue;
        }
        const size_t eqp = ln.find('=');
        if (eqp == std::string::npos) continue;
        const std::string lhs = trim_copy(ln.substr(0, eqp));
        const std::string rhs = trim_copy(ln.substr(eqp + 1));
        if (lhs == "mode" || lhs == "type") continue;
        if (lhs == "ifs_map" || lhs == "map") {
          /* split into 6 or 7 comma-separated EXPRESSION fields */
          std::vector<std::string> fields; std::string cur;
          for (char ch : rhs + ",") { if (ch == ',') { fields.push_back(trim_copy(cur)); cur.clear(); } else cur += ch; }
          while (!fields.empty() && fields.back().empty()) fields.pop_back();
          if (fields.size() < 6) { local_err = "line " + std::to_string(i + 1) + ": ifs_map needs 6 fields a,b,c,d,e,f (+ optional p)"; ok = false; break; }
          AppState::IFSMapExpr m;
          m.a = parse_field(fields[0], "ifs a"); m.b = parse_field(fields[1], "ifs b");
          m.c = parse_field(fields[2], "ifs c"); m.d = parse_field(fields[3], "ifs d");
          m.e = parse_field(fields[4], "ifs e"); m.f = parse_field(fields[5], "ifs f");
          if (fields.size() >= 7 && !fields[6].empty()) m.p = parse_field(fields[6], "ifs p");
          map_exprs.push_back(m);
        }
      }
      if (!ok) { *error = local_err; arena_destroy(&app.ifs_arena); app.ifs_arena_init = false; return false; }
      if (map_exprs.empty()) { *error = "an IFS needs at least one 'ifs_map = a,b,c,d,e,f' line"; arena_destroy(&app.ifs_arena); app.ifs_arena_init = false; return false; }

      /* commit the IFS system */
      app.mode = SystemMode::IFS;
      app.ifs_map_exprs = std::move(map_exprs);
      app.state_names = {"x", "y"};
      app.params = std::move(ifs_params);
      app.definitions.clear();
      app.observables.clear();
      app.parse_error.clear();
      app.detected_dim = 2;
      sync_param_values(app);             /* make param values available to eval */
      /* A coefficient is directly editable iff it's a compile-time constant
       * — it references no parameter/variable. (A literal like -0.04 parses
       * as SUB(0,0.04), and "1-r" references r, so we can't just check for a
       * bare number node; we recurse looking for any NODE_VAR.) When a
       * coefficient is constant we also fold it to its numeric value so the
       * slider has a clean starting point. System-wide ifs_maps_editable is
       * true only when ALL coefficients are constant (enables add/remove). */
      std::function<bool(node_t *)> is_const_expr = [&](node_t *n) -> bool {
        if (!n) return true;
        switch (n->kind) {
          case NODE_VAR: return false;            /* a parameter/variable ref */
          case NODE_CONST: return true;           /* number or named const */
          case NODE_APP: {
            if (!is_const_expr(n->app.head)) return false;
            for (size_t i = 0; i < n->app.argc; ++i)
              if (!is_const_expr(n->app.args[i])) return false;
            return true;
          }
          default: return false;                  /* binders etc. — treat as non-const */
        }
      };
      app.ifs_coef_literal.clear();
      app.ifs_maps_editable = true;
      for (const auto &me : app.ifs_map_exprs) {
        node_t *src[7] = {me.a,me.b,me.c,me.d,me.e,me.f,me.p};
        std::array<bool,7> lit{};
        for (int k = 0; k < 7; ++k) { lit[(size_t)k] = is_const_expr(src[k]); if (!lit[(size_t)k]) app.ifs_maps_editable = false; }
        app.ifs_coef_literal.push_back(lit);
      }
      app.ifs_dirty = true;
      app.ifs_selected = -1;
      app.ifs_maps.clear();               /* evaluated lazily at render time */
      return true;
    }
  }

  std::vector<std::string> next_state_names = {"x", "y", "z"};
  for (size_t line_no = 0; line_no < lines.size(); ++line_no) {
    const std::string line = trim_copy(strip_line_comment(lines[line_no]));
    if (line.empty()) continue;

    if (starts_with(line, "state ") || starts_with(line, "vars ") || starts_with(line, "variables ")) {
      const size_t space = line.find(' ');
      std::string rest = trim_copy(line.substr(space + 1));
      if (!rest.empty() && rest[0] == '=') rest = trim_copy(rest.substr(1));
      std::string state_error;
      if (!parse_state_name_list(rest, &next_state_names, &state_error)) {
        *error = "line " + std::to_string(line_no + 1) + ": " + state_error;
        return false;
      }
      continue;
    }

    const size_t eq = line.find('=');
    if (eq != std::string::npos) {
      const std::string lhs = trim_copy(line.substr(0, eq));
      if (lhs == "state" || lhs == "vars" || lhs == "variables") {
        std::string state_error;
        if (!parse_state_name_list(trim_copy(line.substr(eq + 1)), &next_state_names, &state_error)) {
          *error = "line " + std::to_string(line_no + 1) + ": " + state_error;
          return false;
        }
      }
    }
  }

  const size_t dim = next_state_names.size();
  arena_t next_arena{};
  arena_init(&next_arena, 1024 * 1024);

  std::vector<AppState::Definition> next_definitions;
  std::vector<AppState::Param> next_params;
  std::vector<AppState::Observable> next_observables;
  std::vector<node_t *> next_equations(dim, nullptr);
  std::vector<node_t *> next_next_equations(dim, nullptr);
  std::vector<bool> seen_derivatives(dim, false);
  std::vector<bool> seen_next(dim, false);
  std::vector<node_t *> initial_bodies(dim, nullptr);
  bool explicit_mode = false;
  SystemMode next_mode = SystemMode::ODE;
  bool next_poincare_enabled = false;
  bool next_view2d_set = false;
  double next_view2d[4] = {0, 0, 0, 0};
  node_t *next_section_body = nullptr;
  node_t *next_section_x_body = nullptr;
  node_t *next_section_y_body = nullptr;
  SectionDirection next_direction = SectionDirection::Positive;
  std::array<node_t *, 3> next_plot3d_bodies = {nullptr, nullptr, nullptr};
  std::array<std::string, 3> next_plot3d_labels = {"x", "y", "z"};

  auto parse_expr = [&](const std::string &expr, const std::string &label) -> node_t * {
    return parse_expression_or_fail(&next_arena, expr, label, error);
  };

  auto default_axis_expr = [&](size_t axis) -> std::string {
    if (axis < next_state_names.size()) return next_state_names[axis];
    return "0";
  };
  for (size_t axis = 0; axis < 3; ++axis) {
    next_plot3d_labels[axis] = default_axis_expr(axis);
    next_plot3d_bodies[axis] = parse_expr(next_plot3d_labels[axis], "default plot3d axis");
    if (next_plot3d_bodies[axis] == nullptr) { arena_destroy(&next_arena); return false; }
  }

  for (size_t line_no = 0; line_no < lines.size(); ++line_no) {
    const std::string line = trim_copy(strip_line_comment(lines[line_no]));
    if (line.empty()) continue;

    if (starts_with(line, "state ") || starts_with(line, "vars ") || starts_with(line, "variables ")) {
      continue;
    }

    if (starts_with(line, "param ")) {
      const std::string rest = trim_copy(line.substr(6));
      const size_t eq = rest.find('=');
      if (eq == std::string::npos) {
        *error = "line " + std::to_string(line_no + 1) + ": expected 'param name = value [min,max]'";
        arena_destroy(&next_arena);
        return false;
      }
      const std::string name = trim_copy(rest.substr(0, eq));
      if (is_reserved_value_name(name)) {
        *error = "line " + std::to_string(line_no + 1) + ": '" + name +
                 "' is a reserved constant (e = 2.71828..., pi = 3.14159..., t = time), "
                 "so it can't be a parameter name. Rename it, e.g. '" + name + name + "'.";
        arena_destroy(&next_arena);
        return false;
      }
      if (!is_identifier(name) || is_state_value_name(next_state_names, name) || is_builtin_function_name(name)) {
        *error = "line " + std::to_string(line_no + 1) + ": invalid parameter name '" + name + "'";
        arena_destroy(&next_arena);
        return false;
      }
      for (const auto &p : next_params) {
        if (p.name == name) {
          *error = "line " + std::to_string(line_no + 1) + ": duplicate parameter '" + name + "'";
          arena_destroy(&next_arena);
          return false;
        }
      }
      double min_value = -10.0;
      double max_value = 10.0;
      bool has_range = false;
      const std::string expr = extract_param_expr_and_range(rest.substr(eq + 1), &min_value, &max_value, &has_range);
      node_t *body = parse_expr(expr, "line " + std::to_string(line_no + 1) + " parameter");
      if (body == nullptr) { arena_destroy(&next_arena); return false; }
      AppState::Param param;
      param.name = name;
      param.default_body = body;
      param.min_value = min_value;
      param.max_value = max_value;
      param.has_range = has_range;
      param.line = static_cast<int>(line_no + 1);
      next_params.push_back(param);
      continue;
    }

    if (starts_with(line, "observe ") || starts_with(line, "observable ")) {
      const bool long_prefix = starts_with(line, "observable ");
      const std::string rest = trim_copy(line.substr(long_prefix ? 11 : 8));
      const size_t eq = rest.find('=');
      if (eq == std::string::npos) {
        *error = "line " + std::to_string(line_no + 1) + ": expected 'observe name = expression'";
        arena_destroy(&next_arena);
        return false;
      }
      const std::string name = trim_copy(rest.substr(0, eq));
      if (!is_identifier(name) || is_state_value_name(next_state_names, name) || is_reserved_value_name(name)) {
        *error = "line " + std::to_string(line_no + 1) + ": invalid observable name '" + name + "'";
        arena_destroy(&next_arena);
        return false;
      }
      const std::string expr = trim_copy(rest.substr(eq + 1));
      node_t *body = parse_expr(expr, "line " + std::to_string(line_no + 1) + " observable");
      if (body == nullptr) { arena_destroy(&next_arena); return false; }
      AppState::Observable obs;
      obs.name = name;
      obs.body = body;
      obs.line = static_cast<int>(line_no + 1);
      next_observables.push_back(obs);
      continue;
    }

    const size_t eq = line.find('=');
    if (eq == std::string::npos) {
      *error = "line " + std::to_string(line_no + 1) + ": expected 'name = expression'";
      arena_destroy(&next_arena);
      return false;
    }

    const std::string lhs = trim_copy(line.substr(0, eq));
    const std::string rhs = trim_copy(line.substr(eq + 1));
    if (rhs.empty()) {
      *error = "line " + std::to_string(line_no + 1) + ": empty RHS";
      arena_destroy(&next_arena);
      return false;
    }

    if (lhs == "state" || lhs == "vars" || lhs == "variables") {
      continue;
    }
    if (lhs == "mode" || lhs == "type") {
      explicit_mode = true;
      if (rhs == "ode" || rhs == "continuous") next_mode = SystemMode::ODE;
      else if (rhs == "map" || rhs == "discrete") next_mode = SystemMode::Map;
      else if (rhs == "ifs") next_mode = SystemMode::IFS;
      else {
        *error = "line " + std::to_string(line_no + 1) + ": mode must be 'ode', 'map', or 'ifs'";
        arena_destroy(&next_arena);
        return false;
      }
      continue;
    }
    if (lhs == "integrator") {
      if (rhs == "euler") app.integrator = Integrator::Euler;
      else if (rhs == "rk2" || rhs == "midpoint") app.integrator = Integrator::RK2;
      else if (rhs == "heun") app.integrator = Integrator::Heun;
      else if (rhs == "rk4") app.integrator = Integrator::RK4;
      else if (rhs == "rk38" || rhs == "rk3/8") app.integrator = Integrator::RK38;
      else if (rhs == "rkf45" || rhs == "fehlberg") app.integrator = Integrator::RKF45;
      else if (rhs == "dopri45" || rhs == "dopri" || rhs == "dormand-prince") app.integrator = Integrator::DOPRI45;
      else {
        *error = "line " + std::to_string(line_no + 1) + ": integrator must be euler, rk2, heun, rk4, rk38, rkf45, or dopri45";
        arena_destroy(&next_arena);
        return false;
      }
      continue;
    }
    if (lhs == "plot3d" || lhs == "render" || lhs == "projection") {
      const auto fields = split_commas(rhs);
      if (fields.size() != 3 || fields[0].empty() || fields[1].empty() || fields[2].empty()) {
        *error = "line " + std::to_string(line_no + 1) + ": plot3d must be three expressions, e.g. plot3d = x, y, z";
        arena_destroy(&next_arena);
        return false;
      }
      for (size_t axis = 0; axis < 3; ++axis) {
        next_plot3d_labels[axis] = fields[axis];
        next_plot3d_bodies[axis] = parse_expr(fields[axis], "line " + std::to_string(line_no + 1) + " plot3d axis");
        if (next_plot3d_bodies[axis] == nullptr) { arena_destroy(&next_arena); return false; }
      }
      continue;
    }
    if (lhs == "view2d" || lhs == "view") {
      /* Optional initial 2D window: view2d = xmin, xmax, ymin, ymax.
       * Gives a preset (or user system) a good fixed starting view
       * instead of auto-growing from the orbit. */
      const auto fields = split_commas(rhs);
      double vals[4];
      bool okv = fields.size() == 4;
      for (size_t k = 0; okv && k < 4; ++k)
        okv = parse_number_literal(fields[k].c_str(), &vals[k]);
      if (!okv) {
        *error = "line " + std::to_string(line_no + 1) +
                 ": view2d must be four numbers: xmin, xmax, ymin, ymax";
        arena_destroy(&next_arena);
        return false;
      }
      next_view2d_set = true;
      next_view2d[0] = vals[0]; next_view2d[1] = vals[1];
      next_view2d[2] = vals[2]; next_view2d[3] = vals[3];
      continue;
    }
    if (lhs == "section" || lhs == "poincare") {
      next_section_body = parse_expr(rhs, "line " + std::to_string(line_no + 1) + " section");
      if (next_section_body == nullptr) { arena_destroy(&next_arena); return false; }
      next_poincare_enabled = true;
      continue;
    }
    if (lhs == "section_direction" || lhs == "direction") {
      if (rhs == "positive" || rhs == "+") next_direction = SectionDirection::Positive;
      else if (rhs == "negative" || rhs == "-") next_direction = SectionDirection::Negative;
      else if (rhs == "any" || rhs == "both") next_direction = SectionDirection::Any;
      else {
        *error = "line " + std::to_string(line_no + 1) + ": direction must be positive, negative, or any";
        arena_destroy(&next_arena);
        return false;
      }
      continue;
    }
    if (lhs == "section_plot" || lhs == "section_projection") {
      const auto fields = split_commas(rhs);
      if (fields.size() != 2 || fields[0].empty() || fields[1].empty()) {
        *error = "line " + std::to_string(line_no + 1) + ": section_plot must be two expressions, e.g. section_plot = x, y";
        arena_destroy(&next_arena);
        return false;
      }
      next_section_x_body = parse_expr(fields[0], "line " + std::to_string(line_no + 1) + " section x");
      next_section_y_body = parse_expr(fields[1], "line " + std::to_string(line_no + 1) + " section y");
      if (next_section_x_body == nullptr || next_section_y_body == nullptr) { arena_destroy(&next_arena); return false; }
      continue;
    }
    if (starts_with(lhs, "initial ") || starts_with(lhs, "start ")) {
      const std::string prefix = starts_with(lhs, "initial ") ? "initial " : "start ";
      const std::string var = trim_copy(lhs.substr(prefix.size()));
      const int idx = state_index_for_name(next_state_names, var);
      if (idx < 0) {
        *error = "line " + std::to_string(line_no + 1) + ": unknown initial-state variable '" + var + "'";
        arena_destroy(&next_arena);
        return false;
      }
      initial_bodies[static_cast<size_t>(idx)] = parse_expr(rhs, "line " + std::to_string(line_no + 1) + " initial value");
      if (initial_bodies[static_cast<size_t>(idx)] == nullptr) { arena_destroy(&next_arena); return false; }
      continue;
    }

    std::string name;
    std::vector<std::string> params;
    std::string lhs_error;
    if (!parse_lhs(lhs, &name, &params, &lhs_error)) {
      *error = "line " + std::to_string(line_no + 1) + ": " + lhs_error;
      arena_destroy(&next_arena);
      return false;
    }

    const int derivative_index = derivative_index_for_name(next_state_names, name);
    const int next_index = next_index_for_name(next_state_names, name);
    if ((derivative_index >= 0 || next_index >= 0) && !params.empty()) {
      *error = "line " + std::to_string(line_no + 1) + ": state equations cannot take parameters";
      arena_destroy(&next_arena);
      return false;
    }
    if (derivative_index < 0 && next_index < 0 && (is_reserved_value_name(name) || is_state_value_name(next_state_names, name))) {
      *error = "line " + std::to_string(line_no + 1) + ": '" + name + "' is a reserved runtime value";
      arena_destroy(&next_arena);
      return false;
    }
    if (derivative_index < 0 && next_index < 0 && is_builtin_function_name(name)) {
      *error = "line " + std::to_string(line_no + 1) + ": '" + name + "' is a builtin function name";
      arena_destroy(&next_arena);
      return false;
    }
    for (const auto &existing : next_definitions) {
      if (existing.name == name && existing.params.size() == params.size()) {
        *error = "line " + std::to_string(line_no + 1) + ": duplicate definition for '" + name + "/" + std::to_string(params.size()) + "'";
        arena_destroy(&next_arena);
        return false;
      }
    }
    for (const auto &p : next_params) {
      if (p.name == name) {
        *error = "line " + std::to_string(line_no + 1) + ": definition conflicts with parameter '" + name + "'";
        arena_destroy(&next_arena);
        return false;
      }
    }

    node_t *body = parse_expr(rhs, "line " + std::to_string(line_no + 1) + " RHS");
    if (body == nullptr) { arena_destroy(&next_arena); return false; }

    AppState::Definition def;
    def.name = name;
    def.params = std::move(params);
    def.body = body;
    def.line = static_cast<int>(line_no + 1);
    def.derivative = derivative_index >= 0;
    def.derivative_index = derivative_index;
    def.next_equation = next_index >= 0;
    def.next_index = next_index;

    if (derivative_index >= 0) {
      seen_derivatives[static_cast<size_t>(derivative_index)] = true;
      next_equations[static_cast<size_t>(derivative_index)] = body;
    }
    if (next_index >= 0) {
      seen_next[static_cast<size_t>(next_index)] = true;
      next_next_equations[static_cast<size_t>(next_index)] = body;
    }
    next_definitions.push_back(std::move(def));
  }

  if (!explicit_mode) {
    const bool has_all_next = std::all_of(seen_next.begin(), seen_next.end(), [](bool v) { return v; });
    const bool has_all_deriv = std::all_of(seen_derivatives.begin(), seen_derivatives.end(), [](bool v) { return v; });
    next_mode = has_all_next && !has_all_deriv ? SystemMode::Map : SystemMode::ODE;
  }

  if (next_mode == SystemMode::ODE) {
    for (size_t i = 0; i < dim; ++i) {
      if (!seen_derivatives[i]) {
        *error = "missing required derivative definition: d" + next_state_names[i] + " = ...";
        arena_destroy(&next_arena);
        return false;
      }
    }
  } else {
    for (size_t i = 0; i < dim; ++i) {
      if (!seen_next[i]) {
        *error = "missing required map definition: " + next_state_names[i] + "_next = ...";
        arena_destroy(&next_arena);
        return false;
      }
    }
  }

  if (next_poincare_enabled && next_section_body == nullptr) {
    next_poincare_enabled = false;
  }
  if (next_poincare_enabled && next_section_x_body == nullptr) {
    next_section_x_body = parse_expr(default_axis_expr(0), "default section x");
    next_section_y_body = parse_expr(default_axis_expr(1), "default section y");
    if (next_section_x_body == nullptr || next_section_y_body == nullptr) { arena_destroy(&next_arena); return false; }
  }

  /* ============================================================
   * IR lowering pass — happens BEFORE the AppState swap so a
   * failure here destroys next_arena cleanly without leaving the
   * live system half-replaced.
   *
   * Two-phase: first build def signatures, then lower every body.
   * This way mutual recursion (def A calls def B which calls A)
   * is detected at run time by the cycle guard, not at lowering.
   * ============================================================ */
  std::vector<std::string> next_param_names;
  next_param_names.reserve(next_params.size());
  for (const auto &p : next_params) next_param_names.push_back(p.name);

  std::vector<dynsys::ir::DefSig> next_def_sigs;
  next_def_sigs.reserve(next_definitions.size());
  for (const auto &d : next_definitions) {
    next_def_sigs.push_back({d.name, d.params.size()});
  }

  auto lower_one = [&](node_t *body,
                       const std::vector<std::string> &locals,
                       dynsys::ir::Program *out,
                       const std::string &what) -> bool {
    out->code.clear();
    out->constants.clear();
    if (body == nullptr) return true;  /* optional/unused — leave empty */
    std::string lerr;
    dynsys::ir::LowerContext ctx{next_state_names, next_param_names,
                                  next_def_sigs, locals};
    if (!dynsys::ir::lower(body, ctx, out, &lerr)) {
      *error = what + ": " + lerr;
      return false;
    }
    return true;
  };

  const std::vector<std::string> no_locals;

  std::vector<dynsys::ir::Program> next_definition_programs(next_definitions.size());
  for (size_t i = 0; i < next_definitions.size(); ++i) {
    const auto &d = next_definitions[i];
    if (!lower_one(d.body, d.params, &next_definition_programs[i],
                   "definition '" + d.name + "'")) {
      arena_destroy(&next_arena); return false;
    }
    next_definition_programs[i].arity = d.params.size();
  }

  std::vector<dynsys::ir::Program> next_equation_programs;
  std::vector<dynsys::ir::Program> next_next_equation_programs;
  if (next_mode == SystemMode::ODE) {
    next_equation_programs.resize(dim);
    for (size_t i = 0; i < dim; ++i) {
      if (!lower_one(next_equations[i], no_locals, &next_equation_programs[i],
                     "equation d" + next_state_names[i])) {
        arena_destroy(&next_arena); return false;
      }
    }
  } else {
    next_next_equation_programs.resize(dim);
    for (size_t i = 0; i < dim; ++i) {
      if (!lower_one(next_next_equations[i], no_locals, &next_next_equation_programs[i],
                     next_state_names[i] + "_next")) {
        arena_destroy(&next_arena); return false;
      }
    }
  }

  std::array<dynsys::ir::Program, 3> next_plot3d_programs;
  for (size_t i = 0; i < 3; ++i) {
    if (!lower_one(next_plot3d_bodies[i], no_locals, &next_plot3d_programs[i],
                   std::string("plot3d axis ") + std::to_string(i))) {
      arena_destroy(&next_arena); return false;
    }
  }

  dynsys::ir::Program next_section_program;
  dynsys::ir::Program next_section_x_program;
  dynsys::ir::Program next_section_y_program;
  if (!lower_one(next_section_body,   no_locals, &next_section_program,   "section") ||
      !lower_one(next_section_x_body, no_locals, &next_section_x_program, "section_plot x") ||
      !lower_one(next_section_y_body, no_locals, &next_section_y_program, "section_plot y")) {
    arena_destroy(&next_arena); return false;
  }

  std::vector<dynsys::ir::Program> next_observable_programs(next_observables.size());
  for (size_t i = 0; i < next_observables.size(); ++i) {
    if (!lower_one(next_observables[i].body, no_locals,
                   &next_observable_programs[i],
                   "observable '" + next_observables[i].name + "'")) {
      arena_destroy(&next_arena); return false;
    }
  }

  State next_start = app.start;
  resize_state(next_start, dim);

  if (app.arena_ready) {
    arena_destroy(&app.system_arena);
  }
  app.system_arena = next_arena;
  app.arena_ready = true;
  app.state_names = std::move(next_state_names);
  app.definitions = std::move(next_definitions);
  app.params = std::move(next_params);
  app.observables = std::move(next_observables);
  app.equations = std::move(next_equations);
  app.next_equations = std::move(next_next_equations);
  app.mode = next_mode;
  /* leaving IFS (or recompiling as ODE/map): drop any IFS maps/expressions
   * so nothing stale lingers (their arena is freed too). */
  app.ifs_map_exprs.clear();
  app.ifs_coef_literal.clear();
  app.ifs_maps.clear();
  if (app.ifs_arena_init) { arena_destroy(&app.ifs_arena); app.ifs_arena_init = false; }
  app.poincare_enabled = next_poincare_enabled;
  app.section_body = next_section_body;
  app.section_x_body = next_section_x_body;
  app.section_y_body = next_section_y_body;
  app.section_direction = next_direction;
  app.plot3d_bodies = next_plot3d_bodies;
  app.plot3d_labels = next_plot3d_labels;
  app.start = next_start;
  resize_state(app.current, dim);

  /* Swap in the lowered IR alongside the AST. */
  app.definition_programs       = std::move(next_definition_programs);
  app.equation_programs         = std::move(next_equation_programs);
  app.next_equation_programs    = std::move(next_next_equation_programs);
  app.plot3d_programs           = std::move(next_plot3d_programs);
  app.section_program           = std::move(next_section_program);
  app.section_x_program         = std::move(next_section_x_program);
  app.section_y_program         = std::move(next_section_y_program);
  app.observable_programs       = std::move(next_observable_programs);

  /* Reusable VM scratch — size it to the def count so cycle/cache
   * arrays are right. Also pre-size the integrator scratch states
   * so the hot loop never reallocates. */
  dynsys::ir::scratch_init(&app.eval_scratch, app.definition_programs.size());
  dynsys::ir::dual_scratch_init(&app.ad_scratch,
                                app.definition_programs.size());
  resize_state(app.scratch_k1,  dim);
  resize_state(app.scratch_k2,  dim);
  resize_state(app.scratch_k3,  dim);
  resize_state(app.scratch_k4,  dim);
  resize_state(app.scratch_mid, dim);

  /* ============================================================
   * IR/AST self-check: for every lowered program, evaluate it via
   * BOTH the new IR and the original eval_expr_at AST walker on a
   * probe state, and require agreement. Fires once per compile so
   * the cost is irrelevant; any divergence aborts before the
   * simulation can pick up a wrong-but-plausible trajectory.
   * ============================================================ */
  {
    /* Reset state for a deterministic probe regardless of how the
     * GUI has been driven previously. */
    State probe = app.start;
    resize_state(probe, dim);
    /* Sync params for both code paths to read the same values. */
    sync_param_values(app);

    auto self_check = [&](const dynsys::ir::Program &prog, node_t *ast,
                          const std::string &what) -> bool {
      if (ast == nullptr || prog.code.empty()) return true;
      double v_ir = 0.0, v_ast = 0.0;
      char e1[256] = {0}, e2[256] = {0};
      bool ok_ir = eval_program_at(app, prog, probe, &v_ir, e1, sizeof e1);
      bool ok_ast = eval_expr_at(app, ast, probe, &v_ast, e2, sizeof e2);
      if (ok_ir != ok_ast) {
        *error = "IR/AST self-check disagreement on " + what + ": ir " +
                 (ok_ir ? "ok" : e1) + " vs ast " + (ok_ast ? "ok" : e2);
        return false;
      }
      if (ok_ir && ok_ast) {
        const double mag = std::fmax(1.0, std::fmax(std::fabs(v_ir), std::fabs(v_ast)));
        if (std::fabs(v_ir - v_ast) > 1e-10 * mag) {
          char buf[256];
          std::snprintf(buf, sizeof buf,
                        "IR=%.17g AST=%.17g (diff=%.3e)",
                        v_ir, v_ast, v_ir - v_ast);
          *error = "IR/AST self-check value mismatch on " + what + ": " + buf;
          return false;
        }
      }
      return true;
    };
    for (size_t i = 0; i < app.equation_programs.size(); ++i) {
      if (!self_check(app.equation_programs[i], app.equations[i],
                      "d" + app.state_names[i])) return false;
    }
    for (size_t i = 0; i < app.next_equation_programs.size(); ++i) {
      if (!self_check(app.next_equation_programs[i], app.next_equations[i],
                      app.state_names[i] + "_next")) return false;
    }
    for (size_t i = 0; i < 3; ++i) {
      if (!self_check(app.plot3d_programs[i], app.plot3d_bodies[i],
                      std::string("plot3d[") + std::to_string(i) + "]"))
        return false;
    }
    if (!self_check(app.section_program,   app.section_body,   "section"))   return false;
    if (!self_check(app.section_x_program, app.section_x_body, "section_x")) return false;
    if (!self_check(app.section_y_program, app.section_y_body, "section_y")) return false;
    for (size_t i = 0; i < app.observable_programs.size(); ++i) {
      if (!self_check(app.observable_programs[i], app.observables[i].body,
                      "observable " + app.observables[i].name))
        return false;
    }
    /* Spot-check user defs at arity 0 (the ones the IR can call as
     * scalar dependencies). N-arity defs are exercised transitively
     * via any caller equation; this catches isolated 0-arity helpers
     * that no equation references but the user typed in. */
    for (size_t i = 0; i < app.definitions.size(); ++i) {
      if (app.definitions[i].params.empty() && i < app.definition_programs.size()) {
        if (!self_check(app.definition_programs[i], app.definitions[i].body,
                        "definition " + app.definitions[i].name))
          return false;
      }
    }
  }

  // Initialize parameter values from their parsed default expressions, allowing earlier params to be referenced.
  char err[256] = {0};
  for (auto &param : app.params) {
    double v = 0.0;
    if (!eval_expr_at(app, param.default_body, app.start, &v, err, sizeof(err))) {
      *error = "could not initialize parameter '" + param.name + "': " + err;
      return false;
    }
    param.value = v;
    param.default_value = v;
    if (!param.has_range) {
      const double span = std::max(1.0, std::fabs(v) * 2.0);
      param.min_value = v - span;
      param.max_value = v + span;
    }
  }

  /* Params are now finalized — publish the flat values array the
   * IR run path reads. The GUI updates this whenever a slider moves
   * (see draw_gui). */
  sync_param_values(app);

  /* PHASE6-UI: detect the effective dimension. Starting from the top
   * state index, a dimension counts as "dummy" (not a real dimension to
   * visualise in 3D) when its dynamics are trivial across several random
   * probes: derivative ~ 0 for an ODE, or next == current for a map.
   * We shrink detected_dim past contiguous dummy top dimensions only, so
   * e.g. a 3-state map with z_next = z detects as 2D, while a genuine
   * 3D system stays 3D. */
  {
    app.detected_dim = static_cast<int>(dim);
    auto trivial_top = [&](size_t k) -> bool {
      const unsigned seeds[4] = {1u, 7u, 13u, 29u};
      for (unsigned s : seeds) {
        State probe = make_state_like(dim, 0.0);
        unsigned r = s;
        for (size_t i = 0; i < dim; ++i) {
          r = r * 1664525u + 1013904223u;
          set_state_at(probe, i, (static_cast<double>(r % 2000u) / 1000.0) - 1.0);
        }
        char e[128] = {0};
        if (app.mode == SystemMode::Map) {
          State nxt{};
          if (!step_map_state(app, probe, &nxt, e, sizeof(e))) return false;
          if (std::fabs(state_at(nxt, k) - state_at(probe, k)) > 1e-9) return false;
        } else {
          State d{};
          if (!eval_rhs(app, probe, &d, e, sizeof(e))) return false;
          if (std::fabs(state_at(d, k)) > 1e-9) return false;
        }
      }
      return true;
    };
    for (size_t k = dim; k-- > 1;) { /* never drop below 1 dim */
      if (trivial_top(k))
        app.detected_dim = static_cast<int>(k);
      else
        break;
    }
  }

  for (size_t i = 0; i < dim; ++i) {
    if (initial_bodies[i] != nullptr) {
      double v = 0.0;
      if (!eval_expr_at(app, initial_bodies[i], app.start, &v, err, sizeof(err))) {
        *error = "could not initialize state variable '" + app.state_names[i] + "': " + err;
        return false;
      }
      set_state_at(app.start, i, v);
    }
  }

  app.parse_error.clear();
  app.runtime_error.clear();
  app.poincare_points.clear();
  app.have_last_section = false;

  /* PHASE6-UI: apply an explicit initial 2D window if the system gave one
   * (view2d = xmin,xmax,ymin,ymax). This makes presets open with a good
   * fixed view instead of auto-growing. */
  if (next_view2d_set) {
    app.phase_auto_bounds = false;
    app.phase_x_min = static_cast<float>(next_view2d[0]);
    app.phase_x_max = static_cast<float>(next_view2d[1]);
    app.phase_y_min = static_cast<float>(next_view2d[2]);
    app.phase_y_max = static_cast<float>(next_view2d[3]);
    app.phase_bounds_valid = false;
    app.home_x_min = app.phase_x_min; app.home_x_max = app.phase_x_max;
    app.home_y_min = app.phase_y_min; app.home_y_max = app.phase_y_max;
    app.home_view_set = true;
  } else {
    /* No explicit view2d: remember a sane default home window centered on the
     * start state so Reset View still has somewhere good to go. */
    const double cx = (!app.start.v.empty() && std::isfinite(app.start.v[0])) ? app.start.v[0] : 0.0;
    const double cy = (app.start.v.size() > 1 && std::isfinite(app.start.v[1])) ? app.start.v[1] : 0.0;
    app.home_x_min = (float)(cx - 4.0); app.home_x_max = (float)(cx + 4.0);
    app.home_y_min = (float)(cy - 4.0); app.home_y_max = (float)(cy + 4.0);
    app.home_view_set = true;
  }

  /* PHASE0: equation-tied analysis is always stale after a recompile —
   * the located equilibrium, its Jacobian/spectrum, and the saddle
   * separatrices were computed for the previous vector field. Clear
   * them. */
  app.fractal_dirty = true; /* PHASE C: recompute the fractal image too */
  app.basin_dirty = true;   /* PHASE B/C: recompute basins for the new system */
  app.scan_dirty = true; app.scan_view_init = false; /* PHASE B: recompute 2-param scan */
  app.fixed_ready = false;
  app.cas_eig_ready = false; /* stale certified spectrum from the previous system */
  app.fixed_jacobian.clear();
  app.fixed_classification.clear();
  app.fixed_eigenvalues.clear();
  app.separatrix_curves.clear();

  /* User-seeded extra orbits keep their meaning as long as the state
   * variables are the same (same count and names). If the state space
   * itself changed, their stored coordinates no longer correspond, so
   * drop them; otherwise keep them on screen. */
  if (app.state_names != prev_state_names) {
    app.phase_trajectories.clear();
  }

  /* PHASE B+: auto-configure the bifurcation sweep for THIS system so the
   * defaults aren't pointing at a previous system's parameter/range (the
   * usual reason "the bifurcation diagram looks wrong"). Only do this when
   * the system actually changed, so we don't stomp a user's tweaks while
   * they iterate on one system. We pick the first parameter and its
   * declared [min,max] range, and the first state variable as observable. */
  if (app.state_names != prev_state_names || [&]{
        std::vector<std::string> cur;
        for (const auto &pp : app.params) cur.push_back(pp.name);
        return cur != prev_param_names;
      }()) {
    if (!app.params.empty()) {
      const auto &p0 = app.params[0];
      std::snprintf(app.bif_param, sizeof(app.bif_param), "%s", p0.name.c_str());
      if (p0.has_range) { app.bif_start = p0.min_value; app.bif_end = p0.max_value; }
      else { app.bif_start = p0.value - 1.0; app.bif_end = p0.value + 1.0; }
      /* sensible cascade defaults for maps vs ODEs */
      if (app.mode == SystemMode::Map) {
        app.bif_discard = 1000; app.bif_keep = 200; app.bif_slices = 800;
      } else {
        app.bif_discard = 4000; app.bif_keep = 2000; app.bif_slices = 400;
      }
    }
    if (!app.state_names.empty())
      std::snprintf(app.bif_observable, sizeof(app.bif_observable), "%s",
                    app.state_names[0].c_str());
    app.bif_view_valid = false;
    app.bif_autorun_done = false; /* let the new system auto-run its sweep on entry */
    app.cont_autorun_done = false; /* and re-trace continuation on entry */
    app.lcc_autorun_done = false; app.lcc_has_data = false; /* re-sweep LC on a new system */
  }
  /* Invalidate stale per-system DATA on EVERY successful (re)compile — not only
   * when the variable/parameter NAMES change. Two different systems can share
   * the same names (e.g. both use x,y and a) but have different equations, in
   * which case the bridge/bifurcation/scan caches would otherwise keep showing
   * the previous system until manually rebuilt. The defaults above (bif range,
   * observable) are name-structural and stay gated; these clears are not. */
  app.bifurcation_points.clear();
  app.bifurcation_lyapunov.clear();
  app.bifurcation_period.clear();
  app.bridge_built = false; app.bridge_cam_init = false;
  app.bridge_point_count = 0;
  /* Tie the projection-solid bridge to the loaded system: pick the family it
   * actually supports (quadratic / cubic / sine), or fall back to the
   * "current system" mode when this system has no valid unified bridge (e.g.
   * Henon, an ODE flow). This is what makes the family match the loaded system
   * instead of always showing Mandelbrot/logistic. */
  {
    const int fam = bridge_family_for_system(app);
    if (fam == 0) app.bridge_family = AppState::BridgeFamily::Quadratic;
    else if (fam == 1) app.bridge_family = AppState::BridgeFamily::Cubic;
    else if (fam == 2) app.bridge_family = AppState::BridgeFamily::Sine;
    if (fam < 0 && app.bridge_mode == AppState::BridgeMode::ProjectionSolid)
      app.bridge_mode = AppState::BridgeMode::CurrentSystem; /* no bridge -> lift loaded system */
    else if (fam >= 0)
      app.bridge_mode = AppState::BridgeMode::ProjectionSolid;
  }
  app.fractal_dirty = true;
  app.basin_dirty = true;
  app.scan_view_init = false;
  app.cas_eig_ready = false;   /* stale certified spectrum from the previous system */
  app.hopf_l1_ready = false;   /* stale Hopf classification */
  app.fold_a_ready = false;    /* stale fold classification */
  app.cont_switched_ready = false; /* stale switched branch */
  app.cas_equi_ready = false;      /* stale exact equilibria */
  app.cyc_ready = false;           /* stale collocation cycle branch */
  app.lpc_ready = false;           /* stale fold-of-cycles curve */
  return true;
}

void copy_system_input(AppState &app, const char *text) {
  std::snprintf(app.system_input, sizeof(app.system_input), "%s", text);
}

struct SystemPreset {
  const char *name;
  const char *category;
  const char *source;
  State start;
  double dt;
  int steps_per_frame;
  Integrator integrator;
  float zoom;
  float center_x;
  float center_y;
  float center_z;
};

const SystemPreset kPresets[] = {
    {"Lorenz", "Continuous 3D", "# Lorenz system\nmode = ode\nintegrator = rk4\nparam sigma = 10 [0,50]\nparam rho = 28 [0,100]\nparam beta = 8 / 3 [0,10]\nobserve r = sqrt(x*x + y*y + z*z)\nsection = z - 27\nsection_direction = positive\nsection_plot = x, y\ndx = sigma * (y - x)\ndy = x * (rho - z) - y\ndz = x * y - beta * z\n", {0.1,0.1,0.1,0.0}, 0.01, 3, Integrator::RK4, 1.0f, 0.0f, 0.0f, 0.0f},
    {"Rossler", "Continuous 3D", "# Rössler system\nmode = ode\nintegrator = rk4\nparam a = 0.2 [0,1]\nparam b = 0.2 [0,1]\nparam c = 5.7 [0,15]\nobserve r = sqrt(x*x + y*y + z*z)\nsection = x\nsection_direction = positive\nsection_plot = y, z\ndx = 0 - (y + z)\ndy = x + a * y\ndz = b + z * (x - c)\n", {0.1,0.1,0.1,0.0}, 0.01, 3, Integrator::RK4, 1.0f,0.0f,0.0f,0.0f},
    {"Thomas", "Continuous 3D", "# Thomas cyclically symmetric attractor\nmode = ode\nintegrator = rk4\nparam b = 0.208186 [0,1]\nwave(u) = sin(u)\nobserve r = sqrt(x*x + y*y + z*z)\nsection = z\nsection_direction = positive\nsection_plot = x, y\ndx = wave(y) - b * x\ndy = wave(z) - b * y\ndz = wave(x) - b * z\n", {0.1,0,0,0}, 0.01, 4, Integrator::RK4, 1.0f,0.0f,0.0f,0.0f},
    {"Langford / Aizawa", "Continuous 3D", "# Langford, also called Aizawa\nmode = ode\nintegrator = rk4\nparam a = 0.95 [0,2]\nparam b = 0.7 [0,2]\nparam c = 0.6 [0,2]\nparam d = 3.5 [0,6]\nparam e = 0.25 [0,1]\nparam f = 0.1 [0,1]\nr2 = pow(x, 2) + pow(y, 2)\ndx = (z - b) * x - d * y\ndy = d * x + (z - b) * y\ndz = c + a * z - pow(z, 3) / 3 - r2 * (1 + e * z) + f * z * pow(x, 3)\n", {0.1,0,0,0}, 0.005, 4, Integrator::RK4, 1.0f,0,0,0},
    {"Dadras", "Continuous 3D", "# Dadras attractor\nmode = ode\nintegrator = rk4\nparam a = 3 [0,8]\nparam b = 2.7 [0,8]\nparam c = 1.7 [0,8]\nparam d = 2 [0,8]\nparam e = 9 [0,15]\ndx = y - a * x + b * y * z\ndy = c * y - x * z + z\ndz = d * x * y - e * z\n", {0.1,0.1,0.1,0}, 0.005, 4, Integrator::RK4,1.0f,0,0,0},
    {"Chen / Chen-Lee", "Continuous 3D", "# Chen / Chen-Lee system\nmode = ode\nintegrator = rk4\nparam alpha = 5 [0,10]\nparam beta = 0 - 10 [-20,0]\nparam delta = 0 - 0.38 [-2,2]\ndx = alpha * x - y * z\ndy = beta * y + x * z\ndz = delta * z + x * y / 3\n", {1,1,1,0}, 0.002, 8, Integrator::RK4,1.0f,0,0,0},
    {"Lorenz83", "Continuous 3D", "# Lorenz83 system\nmode = ode\nintegrator = rk4\nparam a = 0.95 [0,2]\nparam b = 7.91 [0,12]\nparam f = 4.83 [0,8]\nparam g = 4.66 [0,8]\ndx = 0 - a * x - pow(y, 2) - pow(z, 2) + a * f\ndy = 0 - y + x * y - b * x * z + g\ndz = 0 - z + b * x * y + x * z\n", {0.1,0.1,0.1,0}, 0.005, 4, Integrator::RK4,1.0f,0,0,0},
    {"Halvorsen", "Continuous 3D", "# Halvorsen attractor\nmode = ode\nintegrator = rk4\nparam a = 1.89 [0,4]\ndx = 0 - a * x - 4 * y - 4 * z - pow(y, 2)\ndy = 0 - a * y - 4 * z - 4 * x - pow(z, 2)\ndz = 0 - a * z - 4 * x - 4 * y - pow(x, 2)\n", {0.1,0,0,0}, 0.005, 4, Integrator::RK4,1.0f,0,0,0},
    {"Rabinovich-Fabrikant", "Continuous 3D", "# Rabinovich-Fabrikant system\nmode = ode\nintegrator = rk4\nparam alpha = 0.14 [0,1]\nparam gamma = 0.10 [0,1]\ndx = y * (z - 1 + pow(x, 2)) + gamma * x\ndy = x * (3 * z + 1 - pow(x, 2)) + gamma * y\ndz = 0 - 2 * z * (alpha + x * y)\n", {0.1,0.1,0.1,0}, 0.002, 8, Integrator::RK4,1.0f,0,0,0},
    {"Three-Scroll Unified", "Continuous 3D", "# Three-Scroll Unified Chaotic System (param 'ee'; bare 'e' is the constant)\nmode = ode\nintegrator = rk4\nparam a = 32.48 [0,60]\nparam b = 45.84 [0,70]\nparam c = 1.18 [0,5]\nparam d = 0.13 [0,1]\nparam ee = 0.57 [0,2]\nparam f = 14.7 [0,30]\ndx = a * (y - x) + d * x * z\ndy = b * x - x * z + f * y\ndz = c * z + x * y - ee * pow(x, 2)\n", {0.1,0.1,0.1,0}, 0.0005, 16, Integrator::RK4,1.0f,0,0,0},
    {"Sprott", "Continuous 3D", "# Sprott system\nmode = ode\nintegrator = rk4\nparam a = 2.07 [0,4]\nparam b = 1.79 [0,4]\ndx = y + a * x * y + x * z\ndy = 1 - b * pow(x, 2) + y * z\ndz = x - pow(x, 2) - pow(y, 2)\n", {0.1,0.1,0.1,0}, 0.002, 8, Integrator::RK4,1.0f,0,0,0},
    {"Four-Wing", "Continuous 3D", "# Four-Wing attractor\nmode = ode\nintegrator = rk4\nparam a = 0.2 [0,2]\nparam b = 0.01 [0,1]\nparam c = 0 - 0.4 [-2,1]\ndx = a * x + y * z\ndy = b * x + c * y - x * z\ndz = 0 - z - x * y\n", {0.1,0.1,0.1,0}, 0.005, 4, Integrator::RK4,1.0f,0,0,0},
    {"Chua's circuit", "Continuous 3D", R"dyn(# Chua's circuit (double-scroll). f is the piecewise-linear diode,
# written smoothly with abs(): f = m1 x + 0.5 (m0-m1)(|x+1|-|x-1|).
state x, y, z
mode = ode
integrator = rk4
param alpha = 15.6 [0,30]
param beta = 28 [0,60]
param m0 = 0 - 1.143 [-2,0]
param m1 = 0 - 0.714 [-2,0]
dx = alpha * (y - x - (m1 * x + 0.5 * (m0 - m1) * (abs(x + 1) - abs(x - 1))))
dy = x - y + z
dz = 0 - beta * y
)dyn", {0.7,0.0,0.0,0}, 0.01, 4, Integrator::RK4,1.0f,0,0,0},
    {"Sprott B", "Continuous 3D", "# Sprott B attractor (algebraically simple chaos)\nmode = ode\nintegrator = rk4\ndx = y * z\ndy = x - y\ndz = 1 - x * y\n", {0.1,0.1,0.1,0}, 0.01, 4, Integrator::RK4,1.0f,0,0,0},
    {"Nose-Hoover", "Continuous 3D", "# Nose-Hoover oscillator (conservative chaos)\nmode = ode\nintegrator = rk4\nparam a = 1.5 [0,3]\ndx = y\ndy = 0 - x + y * z\ndz = a - y * y\n", {0.1,0.1,0.0,0}, 0.01, 4, Integrator::RK4,1.0f,0,0,0},
    {"Van der Pol", "Phase-plane examples", R"dyn(# Van der Pol oscillator
state x, y
mode = ode
integrator = rk4
param mu = 1.0 [0,5]
plot3d = x, y, 0
view2d = -4, 4, -6, 6
observe energy_like = 0.5 * (x*x + y*y)
initial x = 2
initial y = 0
dx = y
dy = mu * (1 - x*x) * y - x
)dyn", {2.0,0.0,0.0,0}, 0.01, 4, Integrator::RK4,1.0f,0,0,0},
    {"Quasiperiodic torus (2 frequencies)", "Phase-plane examples", R"dyn(# A trajectory on a 2-TORUS built from two independent oscillators at
# frequencies w1 and w2. Each (cos,sin) pair traces a circle; together they
# live on a torus. View this in the 3D SCENE (plot3d embeds it as a donut):
# when w2/w1 is IRRATIONAL the orbit never closes and densely fills the torus
# surface (quasiperiodic). Set w2 to a rational multiple of w1 (try w2 = 2) and
# it closes into a single knotted loop (periodic). The simplest attractor that
# is neither a point nor a plain cycle. State stays bounded (it's two circles).
state u, v
mode = ode
integrator = rk4
param w1 = 1.0 [0,3]
param w2 = 1.6180339887 [0,3]
# u is the angle on the big circle, v the angle on the tube; embed as a torus:
plot3d = (2 + 0.8*cos(v)) * cos(u), (2 + 0.8*cos(v)) * sin(u), 0.8*sin(v)
view2d = -3.2, 3.2, -3.2, 3.2
observe big_angle = u
initial u = 0
initial v = 0
du = w1
dv = w2
)dyn", {0.0,0.0,0.0,0}, 0.01, 4, Integrator::RK4,1.0f,0,0,0},
    {"Lotka-Volterra", "Phase-plane examples", R"dyn(# Lotka-Volterra predator-prey system
state prey, predator
mode = ode
integrator = rk4
view2d = 0, 4, 0, 4
param alpha = 1.1 [0,3]
param beta = 0.4 [0,2]
param delta = 0.1 [0,1]
param gamma = 0.4 [0,2]
plot3d = prey, predator, 0
initial prey = 10
initial predator = 5
dprey = alpha * prey - beta * prey * predator
dpredator = delta * prey * predator - gamma * predator
)dyn", {10.0,5.0,0.0,0}, 0.01, 3, Integrator::RK4,1.0f,0,0,0},
    {"Damped pendulum", "Phase-plane examples", R"dyn(# Autonomous damped pendulum
state theta, omega
view2d = -6.5, 6.5, -5, 5
mode = ode
integrator = rk4
param damping = 0.15 [0,2]
plot3d = theta, omega, 0
initial theta = 2.5
initial omega = 0
dtheta = omega
domega = 0 - sin(theta) - damping * omega
)dyn", {2.5,0.0,0.0,0}, 0.01, 4, Integrator::RK4,1.0f,0,0,0},
    {"Duffing oscillator", "Phase-plane examples", R"dyn(# Unforced Duffing oscillator (double-well)
state x, v
mode = ode
integrator = rk4
param delta = 0.2 [0,2]
param beta = 1.0 [-2,2]
view2d = -2, 2, -2, 2
dx = v
dv = 0 - delta * v + x - beta * x*x*x
)dyn", {1.0,0.0,0.0,0}, 0.01, 4, Integrator::RK4,1.0f,0,0,0},
    {"FitzHugh-Nagumo", "Phase-plane examples", R"dyn(# FitzHugh-Nagumo neuron model
state v, w
mode = ode
integrator = rk4
param I = 0.5 [0,2]
param a = 0.7 [0,1]
param b = 0.8 [0,1]
param tau = 12.5 [1,50]
view2d = -2.5, 2.5, -1, 2
dv = v - v*v*v / 3 - w + I
dw = (v + a - b * w) / tau
)dyn", {-1.0,1.0,0.0,0}, 0.02, 4, Integrator::RK4,1.0f,0,0,0},
    {"Brusselator", "Phase-plane examples", R"dyn(# Brusselator chemical oscillator
state x, y
mode = ode
integrator = rk4
param a = 1.0 [0,5]
param b = 3.0 [0,5]
view2d = 0, 5, 0, 5
dx = a + x*x*y - (b + 1) * x
dy = b * x - x*x*y
)dyn", {1.0,1.0,0.0,0}, 0.01, 4, Integrator::RK4,1.0f,0,0,0},
    {"Selkov glycolysis", "Phase-plane examples", R"dyn(# Sel'kov model of glycolytic oscillations
state x, y
mode = ode
integrator = rk4
param a = 0.08 [0,0.5]
param b = 0.6 [0,2]
view2d = 0, 3, 0, 3
dx = 0 - x + a * y + x*x*y
dy = b - a * y - x*x*y
)dyn", {1.0,1.0,0.0,0}, 0.02, 4, Integrator::RK4,1.0f,0,0,0},
    {"SIR epidemic", "Phase-plane examples", R"dyn(# SIR epidemic model (S-I plane, R = 1 - S - I)
state S, I
mode = ode
integrator = rk4
param beta = 0.6 [0,2]
param gamma = 0.1 [0,1]
view2d = 0, 1, 0, 1
dS = 0 - beta * S * I
dI = beta * S * I - gamma * I
)dyn", {0.99,0.01,0.0,0}, 0.05, 4, Integrator::RK4,1.0f,0,0,0},
    {"Predator-prey (limit cycle)", "Phase-plane examples", R"dyn(# Rosenzweig-MacArthur predator-prey with a stable limit cycle
state x, y
mode = ode
integrator = rk4
param r = 1.0 [0,3]
param K = 10 [1,30]
param a = 1.0 [0,5]
param h = 2.0 [0,5]
param ee = 0.5 [0,2]
param m = 0.3 [0,2]
view2d = 0, 12, 0, 4
dx = r * x * (1 - x / K) - a * x * y / (1 + a * h * x)
dy = ee * a * x * y / (1 + a * h * x) - m * y
)dyn", {5.0,1.0,0.0,0}, 0.02, 4, Integrator::RK4,1.0f,0,0,0},
    {"Simple pendulum (undamped)", "Phase-plane examples", R"dyn(# Undamped pendulum: nested closed orbits + separatrix
state theta, omega
mode = ode
integrator = rk4
view2d = -6.5, 6.5, -3, 3
dtheta = omega
domega = 0 - sin(theta)
)dyn", {0.5,0.0,0.0,0}, 0.02, 4, Integrator::RK4,1.0f,0,0,0},
    {"Henon map", "Discrete maps", "# Hénon map\nmode = map\nparam a = 1.4 [0,2]\nparam b = 0.3 [0,1]\nview2d = -2.0, 2.0, -1.6, 1.6\nx_next = 1 - a * x * x + y\ny_next = b * x\nz_next = z\nobserve radius = sqrt(x*x + y*y)\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Lozi map", "Discrete maps", "# Lozi map (piecewise-linear Henon cousin; chaotic)\nstate x, y\nmode = map\nparam a = 1.7 [0,2]\nparam b = 0.5 [0,1]\nview2d = -2, 2, -1, 1\nx_next = 1 - a * abs(x) + y\ny_next = b * x\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Ikeda map", "Discrete maps", "# Ikeda map (laser cavity; spiral chaotic attractor)\nstate x, y\nmode = map\nparam u = 0.918 [0,1]\nview2d = -0.5, 2, -2.5, 1\nx_next = 1 + u * (x * cos(0.4 - 6 / (1 + x*x + y*y)) - y * sin(0.4 - 6 / (1 + x*x + y*y)))\ny_next = u * (x * sin(0.4 - 6 / (1 + x*x + y*y)) + y * cos(0.4 - 6 / (1 + x*x + y*y)))\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Logistic map", "Discrete maps", "# Logistic map x_{n+1} = r x (1-x)  (genuine 1D map)\n# The classic period-doubling route to chaos.\nstate x\nmode = map\nparam r = 3.9 [0,4]\nview2d = 0, 1, 0, 1\nx_next = r * x * (1 - x)\nobserve x_value = x\n", {0.2,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Logistic: period-3 window (=> chaos)", "Discrete maps", "# Logistic map in the PERIOD-3 WINDOW (r = 3.83).\n# Sharkovskii / Li-Yorke theorem: 'period three implies chaos' -- once a\n# continuous 1-D map has a period-3 orbit, it has periodic orbits of EVERY\n# period and an uncountable chaotic (scrambled) set. Here at r=3.83 the\n# attractor is a clean period-3 cycle (look at the bifurcation diagram: a\n# white band of 3 points inside the chaotic region). Nudge r down to ~3.8284\n# to watch period-3 itself appear by a tangent (saddle-node) bifurcation, and\n# just below it the period-3 cycle period-doubles (3->6->12) into chaos.\nstate x\nmode = map\nparam r = 3.83 [3.7, 3.9]\nview2d = 0, 1, 0, 1\nx_next = r * x * (1 - x)\nobserve x_value = x\n", {0.2,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Tent map", "Discrete maps", "# Tent map: x_next = mu * min(x, 1-x)\n# 1D map; chaotic for mu near 2, x stays in [0,1].\nstate x\nmode = map\nparam mu = 1.9 [0,2]\nview2d = 0, 1, 0, 1\nx_next = mu * min(x, 1 - x)\nobserve x_value = x\n", {0.35,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Sine map", "Discrete maps", "# Sine map x_{n+1} = r sin(pi x)  (1D, period-doubling like logistic)\nstate x\nmode = map\nparam r = 0.9 [0,1]\nview2d = 0, 1, 0, 1\nx_next = r * sin(3.141592653589793 * x)\nobserve x_value = x\n", {0.4,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Gauss / mouse map", "Discrete maps", "# Gauss (mouse) map x_{n+1} = exp(-alpha x^2) + beta  (1D)\nstate x\nmode = map\nparam alpha = 6.2 [1,10]\nparam beta = 0 - 0.5 [-1,1]\nview2d = -1, 1.5, -1, 1.5\nx_next = exp(0 - alpha * x * x) + beta\nobserve x_value = x\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Cubic map", "Discrete maps", "# Cubic map x_{n+1} = a x - x^3  (1D, symmetric period-doubling)\nstate x\nmode = map\nparam a = 2.8 [0,3]\nview2d = -2, 2, -2, 2\nx_next = a * x - x*x*x\nobserve x_value = x\n", {0.2,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Gingerbreadman map", "Discrete maps", "# Gingerbreadman map (piecewise-linear, chaotic)\nstate x, y\nmode = map\nview2d = -8, 8, -8, 8\nx_next = 1 - y + abs(x)\ny_next = x\n", {0.5,3.7,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Standard (Chirikov) map", "Discrete maps", "# Chirikov standard map on the cylinder\nstate theta, p\nmode = map\nparam K = 0.97 [0,5]\nview2d = 0, 6.2832, -3.1416, 3.1416\np_next = p + K * sin(theta)\ntheta_next = theta + p + K * sin(theta)\n", {3.0,0.1,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Tinkerbell map", "Discrete maps", "# Tinkerbell map\nstate x, y\nmode = map\nparam a = 0.9 [0,2]\nparam b = 0 - 0.6013 [-2,2]\nparam c = 2.0 [0,4]\nparam d = 0.5 [0,2]\nview2d = -1.5, 0.5, -1.8, 0.6\nx_next = x*x - y*y + a*x + b*y\ny_next = 2*x*y + c*x + d*y\n", {-0.72,-0.64,0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Newton fractal (z^3 - 1)", "Fractals", "# Newton's method for z^3 = 1, as a 2D map.\n# Three roots -> three basins with fractal boundaries.\n# Switch to the 'basins' view after loading.\nstate x, y\nmode = map\nview2d = -2, 2, -2, 2\n# z^2 = a + i b\nx_next = x - ( ( (x*(x*x-y*y) - y*(2*x*y)) - 1 )*(3*(x*x-y*y)) + ( x*(2*x*y) + y*(x*x-y*y) )*(3*(2*x*y)) ) / ( 9*((x*x-y*y)*(x*x-y*y) + (2*x*y)*(2*x*y)) )\ny_next = y - ( ( x*(2*x*y) + y*(x*x-y*y) )*(3*(x*x-y*y)) - ( (x*(x*x-y*y) - y*(2*x*y)) - 1 )*(3*(2*x*y)) ) / ( 9*((x*x-y*y)*(x*x-y*y) + (2*x*y)*(2*x*y)) )\n", {0.5,0.5,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Complex quadratic (Mandelbrot/Julia)", "Fractals", "# z -> z^2 + c  (x+iy)^2 + (cx+i cy)\n# Fractal view: parameter space = Mandelbrot, state space = Julia.\nstate x, y\nmode = map\nparam cx = 0 [-2,2]\nparam cy = 0 [-2,2]\nview2d = -2, 2, -2, 2\nx_next = x*x - y*y + cx\ny_next = 2*x*y + cy\n", {0.0,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Burning Ship fractal", "Fractals", "# Burning Ship: z -> (|Re z| + i|Im z|)^2 + c.\n# The absolute values break the symmetry and produce the famous 'ship'\n# silhouette with masts/antenna. Open the Fractal view (parameter space).\nstate x, y\nmode = map\nparam cx = 0 [-2.2, 1.2]\nparam cy = 0 [-2, 1]\nview2d = -2.2, 1.2, -2, 1\nx_next = abs(x)*abs(x) - abs(y)*abs(y) + cx\ny_next = 2*abs(x)*abs(y) + cy\n", {0.0,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Multibrot (z^4 + c)", "Fractals", "# Multibrot of degree 4: z -> z^4 + c. Higher powers add lobes to the\n# Mandelbrot set (degree d gives d-1 fold symmetry). Open the Fractal view.\n# z^4 = (x^2-y^2)^2 - (2xy)^2  +  i * 2(x^2-y^2)(2xy)\nstate x, y\nmode = map\nparam cx = 0 [-1.6, 1.6]\nparam cy = 0 [-1.6, 1.6]\nview2d = -1.6, 1.6, -1.6, 1.6\nx_next = (x*x - y*y)*(x*x - y*y) - (2*x*y)*(2*x*y) + cx\ny_next = 2*(x*x - y*y)*(2*x*y) + cy\n", {0.0,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Henon map 2D", "N-dimensional examples", "# Hénon map as a true 2D state\nstate x, y\nmode = map\nview2d = -1.5, 1.5, -0.5, 0.5\nparam a = 1.4 [0,2]\nparam b = 0.3 [0,1]\nplot3d = x, y, 0\nx_next = 1 - a * x * x + y\ny_next = b * x\nobserve radius = sqrt(x*x + y*y)\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"4D coupled oscillator demo", "N-dimensional examples", "# Four-dimensional demo system\nstate x, y, z, w\nmode = ode\nintegrator = rk4\nparam a = 0.1 [0,1]\nplot3d = x, y, w\nobserve r4 = sqrt(x*x + y*y + z*z + w*w)\nsection = w\nsection_direction = positive\nsection_plot = x, z\ninitial x = 0.1\ninitial y = 0\ninitial z = 0\ninitial w = 0.2\ndx = y\ndy = 0 - x - a*y + z\ndz = w\ndw = 0 - z - a*w + x\n", {0.1,0.0,0.0,0}, 0.01, 3, Integrator::RK4,1.0f,0,0,0},
    {"Barnsley fern (IFS)", "Fractals", "# Barnsley fern — an iterated function system.\n# Four affine maps a,b,c,d,e,f with selection probability p.\n# Opens in the IFS view; the chaos game draws the attractor.\nmode = ifs\nifs_map = 0, 0, 0, 0.16, 0, 0, 0.01\nifs_map = 0.85, 0.04, -0.04, 0.85, 0, 1.6, 0.85\nifs_map = 0.2, -0.26, 0.23, 0.22, 0, 1.6, 0.07\nifs_map = -0.15, 0.28, 0.26, 0.24, 0, 0.44, 0.07\n", {0.0,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Sierpinski triangle (IFS)", "Fractals", "# Sierpinski triangle as an IFS (three half-scale maps).\nmode = ifs\nifs_map = 0.5, 0, 0, 0.5, 0, 0, 0.333\nifs_map = 0.5, 0, 0, 0.5, 0.5, 0, 0.333\nifs_map = 0.5, 0, 0, 0.5, 0.25, 0.5, 0.334\n", {0.0,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Spiral IFS (parametrized)", "Fractals", "# A PARAMETRIZED IFS: drag 'theta' (rotation) and 's' (scale) in the\n# Parameters panel and watch the attractor change live.\n# Coefficients are EXPRESSIONS in the parameters.\nmode = ifs\nparam theta = 0.6 [0, 3.14159]\nparam s = 0.62 [0.3, 0.72]\nifs_map = s*cos(theta), 0 - s*sin(theta), s*sin(theta), s*cos(theta), 0, 0, 0.5\nifs_map = s*cos(theta), 0 - s*sin(theta), s*sin(theta), s*cos(theta), 0.5, 0, 0.5\n", {0.0,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
};

constexpr int kPresetCount = static_cast<int>(sizeof(kPresets) / sizeof(kPresets[0]));

void apply_preset(AppState &app, int index) {
  if (index < 0 || index >= kPresetCount) {
    index = 0;
  }
  const SystemPreset &preset = kPresets[index];
  app.selected_preset = index;
  copy_system_input(app, preset.source);
  app.start = preset.start;
  app.dt = preset.dt;
  app.steps_per_frame = preset.steps_per_frame;
  app.integrator = preset.integrator;
  app.zoom = preset.zoom;
  app.center_x = preset.center_x;
  app.center_y = preset.center_y;
  app.center_z = preset.center_z;
}

void set_lorenz(AppState &app) {
  apply_preset(app, 0);
}

/* One-click demo that GUARANTEES the canonical logistic period-doubling
 * tree: find the Logistic preset, load + compile it, configure the sweep
 * over r in [0,4] observing x, run it, and switch to the bifurcation view.
 * This removes every "did I set it up right / am I on the new build" doubt. */
bool compile_system(AppState &app, const char *system_text, std::string *error);
void run_bifurcation(AppState &app);
void reset_simulation(AppState &app);
[[maybe_unused]] void load_logistic_bifurcation_demo(AppState &app) {
  int idx = -1;
  for (int i = 0; i < kPresetCount; ++i)
    if (std::string(kPresets[i].name) == "Logistic map") { idx = i; break; }
  if (idx < 0) return;
  apply_preset(app, idx);
  std::string err;
  if (!compile_system(app, app.system_input, &err)) { app.parse_error = err; return; }
  reset_simulation(app);
  /* configure the sweep explicitly (don't rely on auto-config timing) */
  std::snprintf(app.bif_param, sizeof(app.bif_param), "r");
  std::snprintf(app.bif_observable, sizeof(app.bif_observable), "x");
  app.bif_start = 0.0; app.bif_end = 4.0;
  app.bif_discard = 1000; app.bif_keep = 200; app.bif_slices = 800;
  app.bif_view_valid = false;
  run_bifurcation(app);
  app.active_view = AppState::ActiveView::Bifurcation;
}

void allocate_points(AppState &app) {
  if (app.num_points < 2) app.num_points = 2;
  app.points.assign(static_cast<size_t>(app.num_points), Point{0.0f, 0.0f, 0.0f});
  app.colors.assign(static_cast<size_t>(app.num_points) * 3U, 0.0f);
  app.current_point = 0;
  /* PHASE0: the phase/series views read app.history while the 3D view
   * reads app.points[]. Keeping two independent length caps made the
   * two plots disagree about where the orbit began. Tie the history
   * cap to the 3D buffer so both views show the same span of the same
   * orbit. */
  app.history_limit = app.num_points;
}

void reset_lyapunov(AppState &app) {
  app.lyapunov_ready = false;
  app.lyapunov_shadow = app.start;
  resize_state(app.lyapunov_shadow, app.state_names.size());
  if (!app.lyapunov_shadow.v.empty()) {
    app.lyapunov_shadow.v[0] += app.lyapunov_epsilon;
  }
  app.lyapunov_sum = 0.0;
  app.lyapunov_samples = 0;
  app.lyapunov_estimate = 0.0;
}

void reset_simulation(AppState &app) {
  if (static_cast<int>(app.points.size()) != app.num_points || static_cast<int>(app.colors.size()) != app.num_points * 3) {
    allocate_points(app);
  }
  app.current_point = 0;
  app.current = app.start;
  resize_state(app.current, app.state_names.size());
  app.phase_bounds_valid = false; /* PHASE6: re-seed the eased view */
  for (int i = 0; i < app.num_points; ++i) {
    app.points[static_cast<size_t>(i)] = Point{0.0f, 0.0f, 0.0f};
    app.colors[static_cast<size_t>(i) * 3U + 0U] = 1.0f;
    app.colors[static_cast<size_t>(i) * 3U + 1U] = 0.5f;
    app.colors[static_cast<size_t>(i) * 3U + 2U] = 0.0f;
  }
  app.history.clear();
  app.poincare_points.clear();
  app.have_last_section = false;
  reset_lyapunov(app);
  app.runtime_error.clear();
  /* PHASE0: reset_simulation now resets ONLY the running orbit and its
   * derived buffers. It no longer wipes the analysis layer (extra
   * orbits, separatrices, located equilibria). That layer is cleared
   * explicitly via clear_analysis_objects() at the points where it is
   * genuinely invalidated (see compile_system / equation edits), so a
   * plain reset or a parameter tweak keeps the user's analysis on
   * screen instead of silently forgetting it. */
}

void upload_buffers(AppState &app) {
  glBindVertexArray(app.vao);

  glBindBuffer(GL_ARRAY_BUFFER, app.vbo);
  glBufferData(GL_ARRAY_BUFFER,
               static_cast<GLsizeiptr>(app.points.size() * sizeof(Point)),
               app.points.data(), GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Point), reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, app.cbo);
  glBufferData(GL_ARRAY_BUFFER,
               static_cast<GLsizeiptr>(app.colors.size() * sizeof(float)),
               app.colors.data(), GL_DYNAMIC_DRAW);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(1);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void push_history(AppState &app, const State &s) {
  app.history.push_back(s);
  while (static_cast<int>(app.history.size()) > app.history_limit) {
    app.history.pop_front();
  }
}

uint32_t trajectory_color(size_t index) {
  static const uint32_t colors[] = {
      0xff5fb9ffu, 0xff66df8au, 0xffffcf5au, 0xffff7a7au,
      0xffc783ffu, 0xff5af0dfu, 0xffff9f5au, 0xffc2ff66u,
  };
  return colors[index % (sizeof(colors) / sizeof(colors[0]))];
}

void push_phase_history(AppState &app, AppState::PhaseTrajectory &traj, const State &s) {
  traj.history.push_back(s);
  while (static_cast<int>(traj.history.size()) > app.phase_trace_limit) {
    traj.history.pop_front();
  }
}

void restart_main_from_state(AppState &app, const State &seed) {
  if (static_cast<int>(app.points.size()) != app.num_points || static_cast<int>(app.colors.size()) != app.num_points * 3) {
    allocate_points(app);
  }
  app.current_point = 0;
  app.start = seed;
  resize_state(app.start, app.state_names.size());
  app.current = app.start;
  for (int i = 0; i < app.num_points; ++i) {
    app.points[static_cast<size_t>(i)] = Point{0.0f, 0.0f, 0.0f};
    app.colors[static_cast<size_t>(i) * 3U + 0U] = 1.0f;
    app.colors[static_cast<size_t>(i) * 3U + 1U] = 0.5f;
    app.colors[static_cast<size_t>(i) * 3U + 2U] = 0.0f;
  }
  app.history.clear();
  push_history(app, app.current);
  app.poincare_points.clear();
  app.have_last_section = false;
  reset_lyapunov(app);
  app.runtime_error.clear();
}

/* PHASE6: seed an extra orbit as a COMPLETE curve through the point,
 * integrated both backward and forward in time the way pplane draws an
 * orbit through a click — instead of a forward-only curve that slowly
 * grows one step per frame and looks half-finished. The result is a
 * static trajectory (active=false) so it doesn't keep marching. */
bool step_ode_state(AppState &app, const State &in, State *out, char *err, size_t err_cap);
void run_bifurcation(AppState &app); /* PHASE B+: used by the in-view control strip */
void run_continuation(AppState &app); /* PHASE D: equilibrium branch tracer */
void run_lc_continuation(AppState &app); /* PHASE D step 2: limit-cycle sweep */
void run_homoclinic(AppState &app); /* homoclinic-orbit BVP */
void run_homoclinic_curve(AppState &app); /* 2-param homoclinic locus */
void run_cycle_collocation(AppState &app); /* collocation + arclength + Floquet */

void add_phase_trajectory(AppState &app, const State &seed, const char *label_prefix = "orbit") {
  if (static_cast<int>(app.phase_trajectories.size()) >= app.phase_max_extra_trajectories) {
    app.phase_trajectories.erase(app.phase_trajectories.begin());
  }
  AppState::PhaseTrajectory traj;
  char label[64];
  std::snprintf(label, sizeof(label), "%s %zu", label_prefix, app.phase_trajectories.size() + 1);
  traj.label = label;
  traj.color = trajectory_color(app.phase_trajectories.size());

  State start = seed;
  resize_state(start, app.state_names.size());
  traj.current = start;

  const int half = std::max(50, app.phase_trace_limit / 2);

  /* For ODEs, trace backward then forward by temporarily flipping dt,
   * and assemble one continuous polyline. For maps, only forward makes
   * sense, so we iterate the map forward. */
  std::deque<State> backward, forward;
  char err[256] = {0};

  if (app.mode == SystemMode::ODE) {
    const double dt0 = app.dt;
    /* Stop tracing if the orbit runs far from the seed — this keeps a
     * backward trace that diverges in reverse time from polluting the
     * auto-fit (and from wasting work). Cap is relative to the seed's own
     * magnitude with a generous floor, so it works at any scale without
     * needing the view bounds (which aren't declared yet here). */
    const double sx = state_at(start, (size_t)app.phase_x_index);
    const double sy = state_at(start, (size_t)app.phase_y_index);
    const double seed_mag = std::sqrt(sx * sx + sy * sy);
    const double cap = 1.0e4 + 1.0e3 * seed_mag;
    auto too_far = [&](const State &s) {
      const double x = state_at(s, (size_t)app.phase_x_index) - sx;
      const double y = state_at(s, (size_t)app.phase_y_index) - sy;
      return std::sqrt(x * x + y * y) > cap;
    };
    /* backward */
    app.dt = -std::fabs(dt0);
    State s = start;
    for (int i = 0; i < half; ++i) {
      State nx{};
      if (!step_ode_state(app, s, &nx, err, sizeof(err))) break;
      bool finite = true;
      for (double v : nx.v) finite = finite && std::isfinite(v);
      if (!finite || too_far(nx)) break;
      backward.push_front(nx);
      s = nx;
    }
    /* forward */
    app.dt = std::fabs(dt0);
    s = start;
    for (int i = 0; i < half; ++i) {
      State nx{};
      if (!step_ode_state(app, s, &nx, err, sizeof(err))) break;
      bool finite = true;
      for (double v : nx.v) finite = finite && std::isfinite(v);
      if (!finite || too_far(nx)) break;
      forward.push_back(nx);
      s = nx;
    }
    app.dt = dt0;
  } else {
    State s = start;
    for (int i = 0; i < 2 * half; ++i) {
      State nx{};
      if (!step_state(app, s, &nx, err, sizeof(err))) break;
      forward.push_back(nx);
      s = nx;
    }
  }

  for (const State &s : backward) traj.history.push_back(s);
  traj.history.push_back(start);
  for (const State &s : forward) traj.history.push_back(s);
  /* leave the live cursor at the forward end, but mark inactive so it
   * doesn't keep integrating every frame. */
  traj.current = traj.history.empty() ? start : traj.history.back();
  traj.active = false;
  app.phase_trajectories.push_back(std::move(traj));
}

bool step_phase_trajectory(AppState &app, AppState::PhaseTrajectory &traj) {
  if (!traj.active) return true;
  State next{};
  char err[256] = {0};
  if (!step_state(app, traj.current, &next, err, sizeof(err))) {
    traj.active = false;
    app.runtime_error = err;
    return false;
  }
  traj.current = next;
  push_phase_history(app, traj, traj.current);
  return true;
}

void step_phase_trajectories(AppState &app) {
  for (auto &traj : app.phase_trajectories) {
    step_phase_trajectory(app, traj);
  }
}

[[maybe_unused]] void clear_analysis_objects(AppState &app) {
  app.phase_trajectories.clear();
  app.separatrix_curves.clear();
  app.fixed_ready = false;
  app.fixed_jacobian.clear();
  app.fixed_classification.clear();
  app.fixed_eigenvalues.clear();
}

bool section_crossed(double s0, double s1, SectionDirection dir) {
  if (dir == SectionDirection::Positive) return s0 < 0.0 && s1 >= 0.0;
  if (dir == SectionDirection::Negative) return s0 > 0.0 && s1 <= 0.0;
  return (s0 < 0.0 && s1 >= 0.0) || (s0 > 0.0 && s1 <= 0.0);
}

State lerp_state(const State &a, const State &b, double alpha) {
  const size_t dim = std::max(a.v.size(), b.v.size());
  State out = make_state_like(dim, a.t + (b.t - a.t) * alpha);
  for (size_t i = 0; i < dim; ++i) {
    set_state_at(out, i, state_at(a, i) + (state_at(b, i) - state_at(a, i)) * alpha);
  }
  return out;
}

void maybe_record_poincare(AppState &app, const State &old_s, const State &new_s) {
  if (!app.poincare_enabled || app.section_program.code.empty()) return;
  char err[256] = {0};
  double s0 = 0.0, s1 = 0.0;
  if (!eval_program_at(app, app.section_program, old_s, &s0, err, sizeof(err)) ||
      !eval_program_at(app, app.section_program, new_s, &s1, err, sizeof(err))) {
    app.runtime_error = err;
    return;
  }
  if (!app.have_last_section) {
    app.have_last_section = true;
    app.last_section_state = new_s;
    app.last_section_value = s1;
    return;
  }
  if (section_crossed(s0, s1, app.section_direction)) {
    const double denom = std::fabs(s0) + std::fabs(s1);
    const double alpha = denom > 0.0 ? std::fabs(s0) / denom : 0.5;
    const State cross = lerp_state(old_s, new_s, alpha);
    double px = state_at(cross, 0), py = state_at(cross, 1);
    if (!app.section_x_program.code.empty()) eval_program_at(app, app.section_x_program, cross, &px, err, sizeof(err));
    if (!app.section_y_program.code.empty()) eval_program_at(app, app.section_y_program, cross, &py, err, sizeof(err));
    app.poincare_points.push_back(Point2{px, py});
    if (static_cast<int>(app.poincare_points.size()) > app.poincare_limit) {
      app.poincare_points.erase(app.poincare_points.begin(), app.poincare_points.begin() + (app.poincare_points.size() - app.poincare_limit));
    }
  }
  app.last_section_state = new_s;
  app.last_section_value = s1;
}

void update_lyapunov(AppState &app, const State &old_main, const State &new_main) {
  if (!app.lyapunov_enabled) return;
  char err[256] = {0};
  const size_t dim = app.state_names.size();
  if (!app.lyapunov_ready) {
    app.lyapunov_shadow = old_main;
    resize_state(app.lyapunov_shadow, dim);
    if (!app.lyapunov_shadow.v.empty()) app.lyapunov_shadow.v[0] += app.lyapunov_epsilon;
    app.lyapunov_ready = true;
  }
  State new_shadow{};
  if (!step_state(app, app.lyapunov_shadow, &new_shadow, err, sizeof(err))) {
    return;
  }
  double dist2 = 0.0;
  for (size_t i = 0; i < dim; ++i) {
    const double d = state_at(new_shadow, i) - state_at(new_main, i);
    dist2 += d * d;
  }
  const double dist = std::sqrt(dist2);
  if (dist <= 1e-300 || !std::isfinite(dist)) {
    app.lyapunov_ready = false;
    return;
  }
  app.lyapunov_sum += std::log(dist / app.lyapunov_epsilon);
  app.lyapunov_samples += 1;
  const double scale = app.lyapunov_epsilon / dist;
  app.lyapunov_shadow = make_state_like(dim, new_main.t);
  for (size_t i = 0; i < dim; ++i) {
    const double d = state_at(new_shadow, i) - state_at(new_main, i);
    set_state_at(app.lyapunov_shadow, i, state_at(new_main, i) + d * scale);
  }
  const double elapsed = app.mode == SystemMode::Map ? static_cast<double>(app.lyapunov_samples)
                                                     : std::max(1e-12, new_main.t - app.start.t);
  app.lyapunov_estimate = app.lyapunov_sum / elapsed;
}

bool eval_plot3d(AppState &app, const State &s, Point *out) {
  char err[256] = {0};
  double p[3] = {state_at(s, 0), state_at(s, 1), state_at(s, 2)};
  for (size_t i = 0; i < 3; ++i) {
    if (!app.plot3d_programs[i].code.empty()) {
      if (!eval_program_at(app, app.plot3d_programs[i], s, &p[i], err, sizeof(err))) {
        app.runtime_error = err;
        return false;
      }
    }
  }
  *out = Point{static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2])};
  return true;
}

bool calculate_next_point(AppState &app) {
  State next{};
  char err[256] = {0};
  const State old = app.current;
  if (!step_state(app, old, &next, err, sizeof(err))) {
    app.runtime_error = err;
    app.paused = true;
    return false;
  }

  maybe_record_poincare(app, old, next);
  update_lyapunov(app, old, next);

  app.current = next;
  push_history(app, app.current);

  const auto idx = static_cast<size_t>(app.current_point);
  Point projected{};
  if (!eval_plot3d(app, app.current, &projected)) {
    app.paused = true;
    return false;
  }
  app.points[idx] = projected;
  app.colors[idx * 3U + 0U] = static_cast<float>(app.current_point % 256) / 256.0f;
  app.colors[idx * 3U + 1U] = 0.5f;
  app.colors[idx * 3U + 2U] = static_cast<float>(255 - (app.current_point % 256)) / 256.0f;
  app.current_point = (app.current_point + 1) % app.num_points;
  step_phase_trajectories(app);
  return true;
}

bool key_pressed(AppState &app, GLFWwindow *window, int key, double debounce_time) {
  if (key < 0 || key > GLFW_KEY_LAST) return false;
  const double current_time = glfwGetTime();
  KeyState &state = app.keys[key];
  if (glfwGetKey(window, key) == GLFW_PRESS) {
    if (!state.pressed || (current_time - state.last_pressed_time) > debounce_time) {
      state.pressed = true;
      state.last_pressed_time = current_time;
      return true;
    }
  } else {
    state.pressed = false;
  }
  return false;
}

void process_keyboard(AppState &app, GLFWwindow *window) {
  ImGuiIO &io = ImGui::GetIO();
  if (io.WantCaptureKeyboard) return;
  if (key_pressed(app, window, GLFW_KEY_W, 0.015)) app.angle_x += 5.0f;
  if (key_pressed(app, window, GLFW_KEY_S, 0.015)) app.angle_x -= 5.0f;
  if (key_pressed(app, window, GLFW_KEY_A, 0.015)) app.angle_y -= 5.0f;
  if (key_pressed(app, window, GLFW_KEY_D, 0.015)) app.angle_y += 5.0f;
  /* +/- zoom is handled per-view in render_phase_background /
   * render_scene_background (about cursor or center), so we no longer
   * also handle it here — doing both double-stepped the zoom. */
  if (key_pressed(app, window, GLFW_KEY_SPACE, 0.1)) app.paused = !app.paused;
  if (key_pressed(app, window, GLFW_KEY_X, 0.1)) app.show_axes = !app.show_axes;
  if (key_pressed(app, window, GLFW_KEY_Z, 0.1)) app.show_center_cross = !app.show_center_cross;
  if (key_pressed(app, window, GLFW_KEY_PERIOD, 0.05)) app.steps_per_frame += 1;
  if (key_pressed(app, window, GLFW_KEY_SLASH, 0.05)) app.steps_per_frame *= 2;
  if (key_pressed(app, window, GLFW_KEY_COMMA, 0.05)) app.steps_per_frame = std::max(1, app.steps_per_frame - 1);
  if (key_pressed(app, window, GLFW_KEY_M, 0.05)) app.steps_per_frame = std::max(1, app.steps_per_frame / 2);
  if (key_pressed(app, window, GLFW_KEY_V, 0.05)) app.dt *= 0.5;
  if (key_pressed(app, window, GLFW_KEY_B, 0.05)) app.dt *= 2.0;
  if (key_pressed(app, window, GLFW_KEY_C, 0.015)) reset_simulation(app);
  if (key_pressed(app, window, GLFW_KEY_ESCAPE, 0.01)) glfwSetWindowShouldClose(window, GLFW_TRUE);
}

void check_compile_errors(GLuint shader, const char *type) {
  GLint success = 0;
  GLchar info_log[1024] = {0};
  if (std::strcmp(type, "PROGRAM") == 0) {
    glGetProgramiv(shader, GL_LINK_STATUS, &success);
    if (!success) {
      glGetProgramInfoLog(shader, 1024, nullptr, info_log);
      std::fprintf(stderr, "ERROR::PROGRAM_LINKING_ERROR of type: %s\n%s\n", type, info_log);
    }
  } else {
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
      glGetShaderInfoLog(shader, 1024, nullptr, info_log);
      std::fprintf(stderr, "ERROR::SHADER_COMPILATION_ERROR of type: %s\n%s\n", type, info_log);
    }
  }
}

GLuint create_shader_program(const char *vertex_src, const char *fragment_src) {
  const GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex_shader, 1, &vertex_src, nullptr);
  glCompileShader(vertex_shader);
  check_compile_errors(vertex_shader, "VERTEX");
  const GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment_shader, 1, &fragment_src, nullptr);
  glCompileShader(fragment_shader);
  check_compile_errors(fragment_shader, "FRAGMENT");
  const GLuint shader_program = glCreateProgram();
  glAttachShader(shader_program, vertex_shader);
  glAttachShader(shader_program, fragment_shader);
  glLinkProgram(shader_program);
  check_compile_errors(shader_program, "PROGRAM");
  glDeleteShader(vertex_shader);
  glDeleteShader(fragment_shader);
  return shader_program;
}

void update_projection(AppState &app) {
  const float height = app.window_height > 0 ? static_cast<float>(app.window_height) : 1.0f;
  const float aspect = static_cast<float>(app.window_width) / height;
  if (app.orthographic_3d) {
    const float half_height = 25.0f / std::max(app.zoom, 0.001f);
    glm_ortho(-half_height * aspect, half_height * aspect,
              -half_height, half_height,
              -10000.0f, 10000.0f, app.projection);
  } else {
    glm_perspective(glm_rad(45.0f), aspect, 0.01f, 10000.0f, app.projection);
  }
}

void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
  auto *app = static_cast<AppState *>(glfwGetWindowUserPointer(window));
  if (app == nullptr) return;
  app->window_width = width;
  app->window_height = height;
  glViewport(0, 0, width, height);
  update_projection(*app);
}

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
  (void)mods;
  auto *app = static_cast<AppState *>(glfwGetWindowUserPointer(window));
  if (app == nullptr || ImGui::GetIO().WantCaptureMouse) return;
  if (button == GLFW_MOUSE_BUTTON_LEFT) {
    if (action == GLFW_PRESS) {
      app->mouse_left_pressed = true;
      glfwGetCursorPos(window, &app->last_mouse_x, &app->last_mouse_y);
    } else if (action == GLFW_RELEASE) {
      app->mouse_left_pressed = false;
    }
  }
}

void cursor_position_callback(GLFWwindow *window, double xpos, double ypos) {
  auto *app = static_cast<AppState *>(glfwGetWindowUserPointer(window));
  if (app == nullptr || ImGui::GetIO().WantCaptureMouse) return;
  if (app->mouse_left_pressed) {
    const double dx = xpos - app->last_mouse_x;
    const double dy = ypos - app->last_mouse_y;
    app->angle_y += static_cast<float>(dx) * 0.1f;
    app->angle_x += static_cast<float>(dy) * 0.1f;
    app->last_mouse_x = xpos;
    app->last_mouse_y = ypos;
  }
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
  (void)xoffset;
  auto *app = static_cast<AppState *>(glfwGetWindowUserPointer(window));
  if (app == nullptr || ImGui::GetIO().WantCaptureMouse) return;
  app->zoom *= static_cast<float>(std::pow(1.1, std::max(-1.0, std::min(1.0, yoffset))));
  app->zoom = std::max(0.05f, std::min(app->zoom, 100.0f));
  update_projection(*app);
}

void draw_attractor(AppState &app) {
  glBindVertexArray(app.vao);
  if (app.current_point + 1 < app.num_points) {
    glDrawArrays(GL_LINE_STRIP, app.current_point + 1, app.num_points - app.current_point - 1);
  }
  if (app.current_point > 1) {
    glDrawArrays(GL_LINE_STRIP, 0, app.current_point);
  }
  glBindVertexArray(0);
}

void draw_axes(AppState &app) {
  if (!app.show_axes) return;
  const float axes_vertices[] = {
      -20.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
       20.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        0.0f,-20.0f,0.0f, 0.0f, 1.0f, 0.0f,
        0.0f,20.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f,0.0f,-20.0f, 0.0f, 0.0f, 1.0f,
        0.0f,0.0f,20.0f, 0.0f, 0.0f, 1.0f,
  };
  GLuint axes_vbo = 0, axes_vao = 0;
  glGenVertexArrays(1, &axes_vao);
  glGenBuffers(1, &axes_vbo);
  glBindVertexArray(axes_vao);
  glBindBuffer(GL_ARRAY_BUFFER, axes_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(axes_vertices), axes_vertices, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void *>(3 * sizeof(float)));
  glEnableVertexAttribArray(1);
  glDrawArrays(GL_LINES, 0, 6);
  glDeleteVertexArrays(1, &axes_vao);
  glDeleteBuffers(1, &axes_vbo);
}

/* PHASE6: a system is worth showing in the 3D OpenGL view only when it
 * actually has three independent dimensions to look at. A 2D system
 * (two state variables, projected as x,y,0) would just duplicate the
 * phase plane, so we skip the 3D geometry and show a note instead. */
/* PHASE6-UI: effective dimensionality.
 *   detected_dim is computed at compile time (see compute_detected_dim):
 *   it's the number of state variables whose dynamics are non-trivial.
 *   A map carrying z_next = z, or an ODE with dz = 0, has a dummy 3rd
 *   dimension and detects as 2D. The user can override Auto with an
 *   explicit Force2D / Force3D. */
int effective_dimension(const AppState &app) {
  switch (app.dim_override) {
    case AppState::DimOverride::Force2D: return std::min<int>(2, static_cast<int>(app.state_names.size()));
    case AppState::DimOverride::Force3D: return static_cast<int>(app.state_names.size());
    case AppState::DimOverride::Auto:
    default: return app.detected_dim;
  }
}

bool system_is_3d(const AppState &app) {
  return effective_dimension(app) >= 3 && app.state_names.size() >= 3;
}

/* Which projection-solid bridge family (if any) the CURRENTLY LOADED system
 * supports. The unified bridge is honest only when the loaded map is the real
 * restriction of an iterated complex polynomial (or c*sin z): then one complex
 * iteration produces a single object whose footprint is that family's
 * connectedness set and whose real-axis silhouette is its bifurcation diagram.
 *   - logistic / quadratic / complex-quadratic  -> Quadratic (z^2+c)
 *   - cubic map                                 -> Cubic (z^3+c)
 *   - sine map                                  -> Sine (c*sin z)
 *   - everything else (Henon, Lozi, Ikeda, tent, Gauss, gingerbread, standard,
 *     all ODE flows, phase-plane systems)       -> -1 (NO valid bridge)
 * Detection is by the loaded equations' signatures. */
int bridge_family_for_system(const AppState &app) {
  const char *s = app.system_input;
  const int dim = effective_dimension(app);
  const size_t nstate = app.state_names.size();

  /* The complex-quadratic (Mandelbrot) system specifically: x_next = x^2-y^2+cx,
   * y_next = 2xy + cy and NOTHING else. Other 2-D maps share the quadratic core
   * but add linear state coupling -- Tinkerbell's "+ a*x + b*y", etc. Those are
   * NOT the Mandelbrot system and have no unified bridge, so require the
   * quadratic signatures AND the absence of extra linear state terms. */
  if (app.mode == SystemMode::Map && nstate >= 2 &&
      std::strstr(s, "x*x - y*y") && std::strstr(s, "2*x*y")) {
    /* Linear state coupling like Tinkerbell's "+ a*x + b*y" / "+ c*x + d*y"
     * marks a DIFFERENT 2-D map, not the Mandelbrot system. Look specifically
     * for a parameter letter multiplying x or y (a..d * x|y). The Mandelbrot
     * system only adds a bare parameter (cx, cy), never <param>*x. */
    const bool has_linear_coupling =
        std::strstr(s, "a*x") || std::strstr(s, "b*x") || std::strstr(s, "c*x") || std::strstr(s, "d*x") ||
        std::strstr(s, "a*y") || std::strstr(s, "b*y") || std::strstr(s, "c*y") || std::strstr(s, "d*y") ||
        std::strstr(s, "k*x") || std::strstr(s, "k*y");
    if (!has_linear_coupling)
      return 0; /* Quadratic (complex-quadratic / Mandelbrot) */
  }

  /* 1-D unimodal map families. These MUST be genuinely one-dimensional: a 2-D
   * map that happens to contain sin( or x*x*x (Henon, Lozi, Ikeda, Tinkerbell,
   * the standard map, ...) is NOT one of these and has no unified bridge. */
  if (app.mode == SystemMode::Map && dim == 1 && nstate >= 1) {
    /* logistic and its exact rescalings -> quadratic family */
    if (std::strstr(s, "x * (1 - x)") || std::strstr(s, "x*(1 - x)") ||
        std::strstr(s, "x * (1-x)")  || std::strstr(s, "x*(1-x)"))
      return 0;
    /* cubic map  x_next = a*x - x^3 */
    if (std::strstr(s, "x - x*x*x") || std::strstr(s, "x*x*x"))
      return 1;
    /* sine map  x_next = r*sin(pi*x)  (1-D only) */
    if (std::strstr(s, "sin("))
      return 2;
  }
  return -1; /* no valid unified bridge for this system */
}

/* PHASE-A: which full-area views make sense for the current system.
 *   1D line/cobweb : exactly 1 effective dimension
 *   2D phase       : >= 2 effective dimensions
 *   3D scene       : >= 3 effective dimensions
 *   bifurcation    : any system that has at least one parameter
 * (effective dimension respects the Force2D/Force3D override.) */
/* Human-readable reason a view applies (or what it needs). Paired with
 * view_valid so the UI can explain *why* a view is unavailable instead of
 * silently greying it out — making the consistency model legible. */
const char *view_requirement(AppState::ActiveView v) {
  switch (v) {
    case AppState::ActiveView::Line1D:      return "1-D systems (one state variable)";
    case AppState::ActiveView::Phase2D:     return "2 or more state variables";
    case AppState::ActiveView::Scene3D:     return "3-D systems (three state variables)";
    case AppState::ActiveView::Bifurcation: return "a system with at least one parameter";
    case AppState::ActiveView::Fractal:     return "a map: a parameter (parameter-space) or 2 state vars (state-space Julia)";
    case AppState::ActiveView::Scene3DBridge: return "one 3D object whose two shadows are a fractal set and a bifurcation diagram (z^2+c / z^3+c), or the loaded system lifted into 3D";
    case AppState::ActiveView::Basins:      return "2+ state variables (and ideally several attractors)";
    case AppState::ActiveView::ParamScan2D: return "at least two parameters";
    case AppState::ActiveView::Continuation: return "an ODE with a parameter (equilibrium tracing)";
    case AppState::ActiveView::IFS:         return "always available (a self-contained fractal gallery)";
    case AppState::ActiveView::LimitCycle:  return "an ODE with a parameter (periodic-orbit tracing)";
  }
  return "";
}

bool view_valid(const AppState &app, AppState::ActiveView v) {
  const int d = effective_dimension(app);
  /* IFS systems are a list of affine maps run via the chaos game: they are
   * a planar object, so they display in the 2D phase view (like any other
   * system shows there). No ODE/map-specific views apply (no vector field,
   * parameter, or equilibria). */
  if (app.mode == SystemMode::IFS)
    return v == AppState::ActiveView::Phase2D;
  switch (v) {
    case AppState::ActiveView::Line1D:      return d == 1;
    case AppState::ActiveView::Phase2D:     return d >= 2;
    case AppState::ActiveView::Scene3D:     return system_is_3d(app);
    case AppState::ActiveView::Bifurcation: return !app.params.empty();
    case AppState::ActiveView::Fractal:
      /* parameter-space escape map needs >=1 parameter; state-space (Julia)
       * needs >=2 state variables. A 1-D map (e.g. the logistic map) is valid
       * via parameter space — its escape set IS the classic Mandelbrot-type
       * picture. */
      return app.mode == SystemMode::Map &&
             (!app.params.empty() || app.state_names.size() >= 2);
    case AppState::ActiveView::Scene3DBridge:
      /* Always available: the projection-solid and side-by-side modes are
       * self-contained (they compute z^2+c / z^3+c from scratch and don't need
       * the loaded system), and the current-system mode lifts whatever system
       * is loaded. So the 3D bridge is reachable from anywhere. */
      return true;
    case AppState::ActiveView::Basins:      return app.state_names.size() >= 2;
    case AppState::ActiveView::ParamScan2D: return app.params.size() >= 2;
    case AppState::ActiveView::Continuation: return app.mode == SystemMode::ODE && !app.params.empty();
    case AppState::ActiveView::IFS: return true; /* the standalone gallery is always available */
    case AppState::ActiveView::LimitCycle: return app.mode == SystemMode::ODE && !app.params.empty();
  }
  return false;
}

/* A sensible default view for a freshly compiled system. */
AppState::ActiveView default_view_for(const AppState &app) {
  if (app.mode == SystemMode::IFS) return AppState::ActiveView::Phase2D;
  const int d = effective_dimension(app);
  if (d <= 1) return AppState::ActiveView::Line1D;
  if (d >= 3) return AppState::ActiveView::Scene3D;
  return AppState::ActiveView::Phase2D;
}

/* Snap active_view to something valid if the current choice isn't. */
void ensure_valid_view(AppState &app) {
  if (view_valid(app, app.active_view)) return;
  app.active_view = default_view_for(app);
}

/* PHASE6-UI: (re)create the offscreen color+depth target at a given
 * pixel size. Returns false if FBO setup is incomplete. */
bool ensure_scene_fbo(AppState &app, int w, int h) {
  w = std::max(1, w);
  h = std::max(1, h);
  if (app.scene_fbo != 0 && app.scene_tex_w == w && app.scene_tex_h == h)
    return true;

  if (app.scene_fbo == 0) glGenFramebuffers(1, &app.scene_fbo);
  if (app.scene_tex == 0) glGenTextures(1, &app.scene_tex);
  if (app.scene_depth == 0) glGenRenderbuffers(1, &app.scene_depth);

  glBindTexture(GL_TEXTURE_2D, app.scene_tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE,
               nullptr);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  glBindRenderbuffer(GL_RENDERBUFFER, app.scene_depth);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);

  glBindFramebuffer(GL_FRAMEBUFFER, app.scene_fbo);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                         app.scene_tex, 0);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, app.scene_depth);
  const bool ok =
      glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  if (ok) {
    app.scene_tex_w = w;
    app.scene_tex_h = h;
  }
  return ok;
}

/* Build the model/view/projection for the 3D scene at the given aspect
 * and emit the geometry. Assumes the correct framebuffer + viewport are
 * already bound. */
void render_scene_geometry(AppState &app, float aspect) {
  glEnable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glUseProgram(app.shader_program);

  mat4 projection;
  if (app.orthographic_3d) {
    const float hh = 25.0f / std::max(app.zoom, 0.001f);
    glm_ortho(-hh * aspect, hh * aspect, -hh, hh, -10000.0f, 10000.0f,
              projection);
  } else {
    glm_perspective(glm_rad(45.0f), aspect, 0.01f, 10000.0f, projection);
  }

  mat4 model, view;
  glm_mat4_identity(model);
  glm_mat4_identity(view);
  vec3 translate{}; set_vec3(translate, app.center_x, app.center_y, app.center_z); glm_translate(model, translate);
  const float effective_scale = app.orthographic_3d ? app.scene_scale : app.scene_scale * app.zoom;
  vec3 scale{}; set_vec3(scale, effective_scale, effective_scale, effective_scale); glm_scale(model, scale);
  vec3 axis_x{}; set_vec3(axis_x, 1.0f, 0.0f, 0.0f); glm_rotate(model, glm_rad(app.angle_x), axis_x);
  vec3 axis_y{}; set_vec3(axis_y, 0.0f, 1.0f, 0.0f); glm_rotate(model, glm_rad(app.angle_y), axis_y);
  vec3 view_translation{}; set_vec3(view_translation, 0.0f, 0.0f, -app.camera_distance); glm_translate(view, view_translation);
  glUniformMatrix4fv(glGetUniformLocation(app.shader_program, "model"), 1, GL_FALSE, reinterpret_cast<float *>(model));
  glUniformMatrix4fv(glGetUniformLocation(app.shader_program, "view"), 1, GL_FALSE, reinterpret_cast<float *>(view));
  glUniformMatrix4fv(glGetUniformLocation(app.shader_program, "projection"), 1, GL_FALSE, reinterpret_cast<float *>(projection));
  glBindBuffer(GL_ARRAY_BUFFER, app.vbo);
  glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLsizeiptr>(app.points.size() * sizeof(Point)), app.points.data());
  glBindBuffer(GL_ARRAY_BUFFER, app.cbo);
  glBufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLsizeiptr>(app.colors.size() * sizeof(float)), app.colors.data());
  draw_attractor(app);
  draw_axes(app);
}

/* Render the 3D scene into the offscreen texture at the panel's pixel
 * size. Called from the GUI when the 3D tab is visible. After this the
 * default framebuffer is rebound. */
void render_scene_to_fbo(AppState &app, int w, int h) {
  if (!ensure_scene_fbo(app, w, h)) return;
  glBindFramebuffer(GL_FRAMEBUFFER, app.scene_fbo);
  glViewport(0, 0, w, h);
  const float aspect =
      static_cast<float>(std::max(1, w)) / static_cast<float>(std::max(1, h));
  render_scene_geometry(app, aspect);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

/* PHASE6-UI: render the 3D scene directly to the full window backbuffer
 * (no FBO needed when it's the background), and handle drag-rotate /
 * wheel-zoom when the cursor isn't over a panel. Called from the main
 * loop before the GUI when the active view is the 3D scene. */
void render_scene_background(AppState &app) {
  glViewport(0, 0, app.window_width, app.window_height);
  const float aspect = static_cast<float>(app.window_width) /
                       std::max(1.0f, static_cast<float>(app.window_height));
  render_scene_geometry(app, aspect);

  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse) {
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
      app.angle_y += io.MouseDelta.x * 0.4f;
      app.angle_x += io.MouseDelta.y * 0.4f;
    }
    if (io.MouseWheel != 0.0f) {
      app.zoom *= std::pow(1.15f, (float)clamped_wheel(io.MouseWheel));
      app.zoom = std::max(0.05f, std::min(app.zoom, 100.0f));
      update_projection(app);
    }
    if (ImGui::IsKeyPressed(ImGuiKey_Equal) || ImGui::IsKeyPressed(ImGuiKey_KeypadAdd)) {
      app.zoom = std::min(app.zoom * 1.1f, 100.0f); update_projection(app);
    }
    if (ImGui::IsKeyPressed(ImGuiKey_Minus) || ImGui::IsKeyPressed(ImGuiKey_KeypadSubtract)) {
      app.zoom = std::max(app.zoom / 1.1f, 0.05f); update_projection(app);
    }
  }
}

/* Fallback flat clear (used before the GUI when the 2D phase plane is
 * the active background — the phase plane itself paints an opaque
 * backdrop via the ImGui background draw list). */
void draw_scene(AppState &app) {
  (void)app;
  glDisable(GL_DEPTH_TEST);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

[[maybe_unused]] void draw_center_cross(const AppState &app) {
  if (!app.show_center_cross) return;
  ImDrawList *draw_list = ImGui::GetForegroundDrawList();
  const ImVec2 center(static_cast<float>(app.window_width) * 0.5f, static_cast<float>(app.window_height) * 0.5f);
  const ImU32 color = IM_COL32(180, 180, 180, 180);
  draw_list->AddLine(ImVec2(center.x - 8.0f, center.y), ImVec2(center.x + 8.0f, center.y), color, 1.0f);
  draw_list->AddLine(ImVec2(center.x, center.y - 8.0f), ImVec2(center.x, center.y + 8.0f), color, 1.0f);
}

bool value_by_name(AppState &app, const State &s, const char *name, double *out) {
  char err[128] = {0};
  for (size_t i = 0; i < app.state_names.size(); ++i) {
    if (app.state_names[i] == name) {
      *out = state_at(s, i);
      return true;
    }
  }
  if (std::strcmp(name, "t") == 0) { *out = s.t; return true; }
  for (size_t i = 0; i < app.observables.size(); ++i) {
    if (app.observables[i].name == name) {
      if (i >= app.observable_programs.size()) return false;
      return eval_program_at(app, app.observable_programs[i], s, out, err, sizeof(err));
    }
  }
  if (const auto *param = find_param(app, name)) {
    *out = param->value;
    return true;
  }
  return false;
}


struct PlotBounds {
  double xmin = -1.0;
  double xmax = 1.0;
  double ymin = -1.0;
  double ymax = 1.0;
};

Point2 plot_to_screen(const PlotBounds &b, const ImVec2 &p0, float w, float h, double x, double y) {
  const double dx = std::max(1e-12, b.xmax - b.xmin);
  const double dy = std::max(1e-12, b.ymax - b.ymin);
  const float sx = p0.x + static_cast<float>((x - b.xmin) / dx) * w;
  const float sy = p0.y + h - static_cast<float>((y - b.ymin) / dy) * h;
  return Point2{sx, sy};
}

Point2 screen_to_plot(const PlotBounds &b, const ImVec2 &p0, float w, float h, const ImVec2 &p) {
  const double dx = std::max(1e-12, b.xmax - b.xmin);
  const double dy = std::max(1e-12, b.ymax - b.ymin);
  const double x = b.xmin + (static_cast<double>(p.x - p0.x) / std::max(1.0f, w)) * dx;
  const double y = b.ymin + (1.0 - static_cast<double>(p.y - p0.y) / std::max(1.0f, h)) * dy;
  return Point2{x, y};
}

PlotBounds equalize_bounds_for_canvas(PlotBounds b, float w, float h) {
  const double cx = 0.5 * (b.xmin + b.xmax);
  const double cy = 0.5 * (b.ymin + b.ymax);
  double dx = std::max(1e-12, b.xmax - b.xmin);
  double dy = std::max(1e-12, b.ymax - b.ymin);
  const double canvas_aspect = static_cast<double>(std::max(1.0f, w)) / static_cast<double>(std::max(1.0f, h));
  const double bounds_aspect = dx / dy;
  if (bounds_aspect < canvas_aspect) {
    dx = dy * canvas_aspect;
  } else {
    dy = dx / canvas_aspect;
  }
  b.xmin = cx - 0.5 * dx;
  b.xmax = cx + 0.5 * dx;
  b.ymin = cy - 0.5 * dy;
  b.ymax = cy + 0.5 * dy;
  return b;
}

/* PHASE6-UI: guarantee finite, correctly-ordered bounds with a sane
 * minimum span. This is the single choke point that prevents a bad
 * zoom/pan (or a stray NaN) from collapsing the view to nothing or
 * poisoning every subsequent frame with NaN coordinates. */
PlotBounds sanitize_bounds(PlotBounds b) {
  auto fix = [](double &lo, double &hi) {
    if (!std::isfinite(lo) || !std::isfinite(hi)) { lo = -10.0; hi = 10.0; }
    if (hi < lo) std::swap(lo, hi);
    double span = hi - lo;
    if (!(span > 1e-9)) { /* collapsed or tiny */
      const double c = std::isfinite(0.5 * (lo + hi)) ? 0.5 * (lo + hi) : 0.0;
      lo = c - 0.5; hi = c + 0.5;
    }
    /* clamp absurd magnitudes that would also break rendering */
    const double LIM = 1e12;
    if (lo < -LIM) lo = -LIM;
    if (hi > LIM) hi = LIM;
  };
  fix(b.xmin, b.xmax);
  fix(b.ymin, b.ymax);
  return b;
}

void set_manual_phase_bounds(AppState &app, const PlotBounds &b_in) {
  const PlotBounds b = sanitize_bounds(b_in);
  app.phase_auto_bounds = false;
  app.phase_x_min = static_cast<float>(b.xmin);
  app.phase_x_max = static_cast<float>(b.xmax);
  app.phase_y_min = static_cast<float>(b.ymin);
  app.phase_y_max = static_cast<float>(b.ymax);
}

void pan_phase_bounds(AppState &app, PlotBounds b, double dx, double dy) {
  b.xmin += dx; b.xmax += dx;
  b.ymin += dy; b.ymax += dy;
  set_manual_phase_bounds(app, b);
}

void zoom_phase_bounds(AppState &app, PlotBounds b, const Point2 &anchor_in, double factor) {
  /* Guard against a non-finite anchor (which could come from a prior bad
   * frame) -- fall back to the view center so zoom can always recover. */
  Point2 anchor = anchor_in;
  if (!std::isfinite(anchor.x) || !std::isfinite(anchor.y)) {
    anchor.x = 0.5 * (b.xmin + b.xmax);
    anchor.y = 0.5 * (b.ymin + b.ymax);
  }
  factor = std::max(0.02, std::min(factor, 50.0));
  const double new_w = (b.xmax - b.xmin) * factor;
  const double new_h = (b.ymax - b.ymin) * factor;
  /* keep anchor's fractional position within the view fixed */
  const double fx = (b.xmax - b.xmin) > 1e-300
                        ? (anchor.x - b.xmin) / (b.xmax - b.xmin)
                        : 0.5;
  const double fy = (b.ymax - b.ymin) > 1e-300
                        ? (anchor.y - b.ymin) / (b.ymax - b.ymin)
                        : 0.5;
  b.xmin = anchor.x - fx * new_w;
  b.xmax = b.xmin + new_w;
  b.ymin = anchor.y - fy * new_h;
  b.ymax = b.ymin + new_h;
  set_manual_phase_bounds(app, b);
}

/* Reset the zoom/pan of whatever view is active back to a sane framing. The
 * 2D plot views (phase, basins) return to the preset's "home" window; the
 * fractal returns to the full Mandelbrot/Julia view; the 3D views reset the
 * camera. This is the reliable escape hatch from any runaway zoom/pan. */
void reset_current_view(AppState &app) {
  switch (app.active_view) {
    case AppState::ActiveView::Phase2D:
    case AppState::ActiveView::Basins:
      if (app.home_view_set) {
        PlotBounds h{app.home_x_min, app.home_x_max, app.home_y_min, app.home_y_max};
        set_manual_phase_bounds(app, h);
      } else {
        app.phase_auto_bounds = true; /* fall back to auto-fit */
      }
      app.phase_bounds_valid = false;
      app.basin_dirty = true;
      break;
    case AppState::ActiveView::Fractal:
      if (app.params.empty()) {
        /* state-space fractal (e.g. Newton): the plane is initial (x,y); use a
         * symmetric window around the origin rather than the Mandelbrot box. */
        app.fractal_xmin = -2.0; app.fractal_xmax = 2.0;
        app.fractal_ymin = -2.0; app.fractal_ymax = 2.0;
      } else {
        app.fractal_xmin = -2.5; app.fractal_xmax = 1.0;
        app.fractal_ymin = -1.5; app.fractal_ymax = 1.5;
      }
      app.fractal_dirty = true;
      break;
    case AppState::ActiveView::Scene3D:
    case AppState::ActiveView::Scene3DBridge:
      app.angle_x = 22.0f; app.angle_y = -37.0f; app.zoom = 1.0f;
      update_projection(app);
      break;
    case AppState::ActiveView::ParamScan2D:
      app.scan_view_init = false; app.scan_dirty = true;
      break;
    default:
      /* line/bifurcation/continuation/limit-cycle views refit themselves */
      app.phase_auto_bounds = true; app.phase_bounds_valid = false;
      app.bif_view_valid = false; app.cont_view_valid = false;
      break;
  }
}

uint32_t fade_color(uint32_t color, int alpha) {
  return (color & 0x00ffffffu) | (static_cast<uint32_t>(std::max(0, std::min(alpha, 255))) << 24);
}

void draw_plot_grid(ImDrawList *draw, const ImVec2 &p0, float w, float h, const PlotBounds &b) {
  const ImVec2 p1(p0.x + w, p0.y + h);
  draw->AddRectFilled(p0, p1, IM_COL32(20, 20, 24, 235));
  draw->AddRect(p0, p1, IM_COL32(120, 120, 120, 200));
  for (int i = 1; i < 4; ++i) {
    const float x = p0.x + w * static_cast<float>(i) / 4.0f;
    const float y = p0.y + h * static_cast<float>(i) / 4.0f;
    draw->AddLine(ImVec2(x, p0.y), ImVec2(x, p1.y), IM_COL32(70, 70, 80, 120), 1.0f);
    draw->AddLine(ImVec2(p0.x, y), ImVec2(p1.x, y), IM_COL32(70, 70, 80, 120), 1.0f);
    char xb[32], yb[32];
    std::snprintf(xb, sizeof(xb), "%.3g", b.xmin + (b.xmax - b.xmin) * static_cast<double>(i) / 4.0);
    std::snprintf(yb, sizeof(yb), "%.3g", b.ymax - (b.ymax - b.ymin) * static_cast<double>(i) / 4.0);
    draw->AddText(ImVec2(x + 3.0f, p1.y - 18.0f), IM_COL32(150, 150, 160, 170), xb);
    draw->AddText(ImVec2(p0.x + 4.0f, y + 2.0f), IM_COL32(150, 150, 160, 170), yb);
  }
  if (b.xmin < 0.0 && b.xmax > 0.0) {
    const Point2 a = plot_to_screen(b, p0, w, h, 0.0, b.ymin);
    const Point2 c = plot_to_screen(b, p0, w, h, 0.0, b.ymax);
    draw->AddLine(ImVec2(static_cast<float>(a.x), static_cast<float>(a.y)), ImVec2(static_cast<float>(c.x), static_cast<float>(c.y)), IM_COL32(130, 130, 140, 160), 1.0f);
  }
  if (b.ymin < 0.0 && b.ymax > 0.0) {
    const Point2 a = plot_to_screen(b, p0, w, h, b.xmin, 0.0);
    const Point2 c = plot_to_screen(b, p0, w, h, b.xmax, 0.0);
    draw->AddLine(ImVec2(static_cast<float>(a.x), static_cast<float>(a.y)), ImVec2(static_cast<float>(c.x), static_cast<float>(c.y)), IM_COL32(130, 130, 140, 160), 1.0f);
  }
}

bool eval_phase_vector(AppState &app, const State &s, size_t ix, size_t iy, double *vx, double *vy, char *err, size_t err_cap) {
  if (app.mode == SystemMode::Map) {
    State next{};
    if (!step_state(app, s, &next, err, err_cap)) return false;
    *vx = state_at(next, ix) - state_at(s, ix);
    *vy = state_at(next, iy) - state_at(s, iy);
    return true;
  }
  State rhs{};
  if (!eval_rhs(app, s, &rhs, err, err_cap)) return false;
  *vx = state_at(rhs, ix);
  *vy = state_at(rhs, iy);
  return true;
}

bool eval_phase_component(AppState &app, const State &s, size_t component, double *out, char *err, size_t err_cap) {
  if (app.mode == SystemMode::Map) {
    State next{};
    if (!step_state(app, s, &next, err, err_cap)) return false;
    *out = state_at(next, component) - state_at(s, component);
    return true;
  }
  State rhs{};
  if (!eval_rhs(app, s, &rhs, err, err_cap)) return false;
  *out = state_at(rhs, component);
  return true;
}

PlotBounds current_phase_bounds(AppState &app, size_t ix, size_t iy) {
  PlotBounds b{};
  if (!app.phase_auto_bounds) {
    b.xmin = app.phase_x_min;
    b.xmax = app.phase_x_max;
    b.ymin = app.phase_y_min;
    b.ymax = app.phase_y_max;
  } else {
    /* PHASE B+: robust auto-fit. Collect the data coordinates and frame
     * to a percentile range rather than the raw min/max, so a few extreme
     * points — e.g. from a backward ODE trace that diverges in reverse
     * time — don't blow the view up. This fixes the "first trajectory
     * appears then zooms way out" problem. We still include all points,
     * but the initial window is set from the bulk of the data. */
    std::vector<double> xs, ys;
    xs.reserve(4096); ys.reserve(4096);
    auto collect = [&](const State &s) {
      const double x = state_at(s, ix);
      const double y = state_at(s, iy);
      if (std::isfinite(x) && std::isfinite(y)) { xs.push_back(x); ys.push_back(y); }
    };
    for (const State &s : app.history) collect(s);
    collect(app.current);
    for (const auto &traj : app.phase_trajectories) {
      if (!traj.visible) continue;
      for (const State &s : traj.history) collect(s);
    }
    for (const auto &curve : app.separatrix_curves) {
      if (!curve.visible) continue;
      for (const State &s : curve.history) collect(s);
    }
    if (app.fixed_ready) collect(app.fixed_point);

    if (xs.empty()) {
      const double x = state_at(app.current, ix);
      const double y = state_at(app.current, iy);
      b.xmin = x - 10.0; b.xmax = x + 10.0;
      b.ymin = y - 10.0; b.ymax = y + 10.0;
    } else {
      auto pct = [](std::vector<double> v, double lo, double hi, double *out_lo, double *out_hi) {
        std::sort(v.begin(), v.end());
        const size_t n = v.size();
        auto at = [&](double q) {
          double idx = q * (n - 1);
          size_t i = (size_t)idx;
          double f = idx - i;
          if (i + 1 < n) return v[i] * (1 - f) + v[i + 1] * f;
          return v[i];
        };
        *out_lo = at(lo); *out_hi = at(hi);
      };
      /* if the data is small, just use full extent; for larger sets clip
       * to 1st..99th percentile to drop divergent outliers */
      if (xs.size() < 64) {
        b.xmin = *std::min_element(xs.begin(), xs.end());
        b.xmax = *std::max_element(xs.begin(), xs.end());
        b.ymin = *std::min_element(ys.begin(), ys.end());
        b.ymax = *std::max_element(ys.begin(), ys.end());
      } else {
        pct(xs, 0.01, 0.99, &b.xmin, &b.xmax);
        pct(ys, 0.01, 0.99, &b.ymin, &b.ymax);
      }
    }
  }
  if (b.xmax <= b.xmin) b.xmax = b.xmin + 1.0;
  if (b.ymax <= b.ymin) b.ymax = b.ymin + 1.0;
  const double padx = 0.06 * (b.xmax - b.xmin);
  const double pady = 0.06 * (b.ymax - b.ymin);
  b.xmin -= padx; b.xmax += padx; b.ymin -= pady; b.ymax += pady;

  /* PHASE6: stable auto-fit. Rather than snapping to the exact data
   * extent every frame (which made the view twitch as the orbit moved)
   * or easing (which lagged), we hold a fixed view and only grow it
   * when the data actually leaves it, with a small margin so we don't
   * grow on every step. The view never shrinks on its own — use
   * double-click to refit. This matches how pplane shows a steady
   * window you can watch the orbit fill. */
  if (app.phase_auto_bounds) {
    if (!app.phase_bounds_valid) {
      app.phase_view_xmin = b.xmin; app.phase_view_xmax = b.xmax;
      app.phase_view_ymin = b.ymin; app.phase_view_ymax = b.ymax;
      app.phase_bounds_valid = true;
    } else {
      /* Grow to include new data, with a 10% margin added when growing
       * so we make room in chunks instead of nudging every frame. */
      auto grow_lo = [](double cur, double want) {
        if (want < cur) return want - 0.10 * std::fabs(cur - want) - 1e-9;
        return cur;
      };
      auto grow_hi = [](double cur, double want) {
        if (want > cur) return want + 0.10 * std::fabs(want - cur) + 1e-9;
        return cur;
      };
      app.phase_view_xmin = grow_lo(app.phase_view_xmin, b.xmin);
      app.phase_view_xmax = grow_hi(app.phase_view_xmax, b.xmax);
      app.phase_view_ymin = grow_lo(app.phase_view_ymin, b.ymin);
      app.phase_view_ymax = grow_hi(app.phase_view_ymax, b.ymax);
    }
    b.xmin = app.phase_view_xmin; b.xmax = app.phase_view_xmax;
    b.ymin = app.phase_view_ymin; b.ymax = app.phase_view_ymax;
  }
  return sanitize_bounds(b);
}

void draw_arrow(ImDrawList *draw, ImVec2 a, ImVec2 b, ImU32 color) {
  draw->AddLine(a, b, color, 1.0f);
  const float dx = b.x - a.x;
  const float dy = b.y - a.y;
  const float len = std::sqrt(dx * dx + dy * dy);
  if (len < 1e-3f) return;
  const float ux = dx / len;
  const float uy = dy / len;
  const float size = 4.0f;
  ImVec2 p1(b.x - ux * size - uy * size * 0.55f, b.y - uy * size + ux * size * 0.55f);
  ImVec2 p2(b.x - ux * size + uy * size * 0.55f, b.y - uy * size - ux * size * 0.55f);
  draw->AddTriangleFilled(b, p1, p2, color);
}

/* PHASE6: pplane-style nullclines via proper marching squares.
 *
 * We sample the chosen vector-field component on a regular grid of
 * (grid+1) x (grid+1) NODES once (so each node is evaluated a single
 * time and adjacent cells share node values — no seams, ~4x fewer
 * evals than the old per-cell sampling). For each cell we then emit
 * the marching-squares contour of f == 0 using the full 16-case
 * table, including both diagonals for the ambiguous saddle cases, and
 * draw it as a smooth thick polyline so it reads like pplane's
 * nullclines rather than a scatter of disconnected ticks. */
void draw_nullcline_for_component(AppState &app, ImDrawList *draw,
                                  const PlotBounds &bounds, const ImVec2 &p0,
                                  float w, float h, size_t ix, size_t iy,
                                  size_t component, ImU32 color) {
  const int grid = std::max(8, std::min(app.phase_grid, 120));
  const int nx = grid + 1;
  const int ny = grid + 1;
  const double dx = (bounds.xmax - bounds.xmin) / static_cast<double>(grid);
  const double dy = (bounds.ymax - bounds.ymin) / static_cast<double>(grid);
  char err[256] = {0};

  /* Sample nodes once. */
  std::vector<double> f(static_cast<size_t>(nx) * ny, 0.0);
  std::vector<uint8_t> ok(static_cast<size_t>(nx) * ny, 0);
  State st = app.current;
  resize_state(st, app.state_names.size());
  for (int j = 0; j < ny; ++j) {
    const double y = bounds.ymin + j * dy;
    for (int i = 0; i < nx; ++i) {
      const double x = bounds.xmin + i * dx;
      set_state_at(st, ix, x);
      set_state_at(st, iy, y);
      double val = 0.0;
      const bool good =
          eval_phase_component(app, st, component, &val, err, sizeof(err)) &&
          std::isfinite(val);
      const size_t idx = static_cast<size_t>(j) * nx + i;
      f[idx] = val;
      ok[idx] = good ? 1u : 0u;
    }
  }

  auto node = [&](int i, int j) -> double { return f[static_cast<size_t>(j) * nx + i]; };
  auto good = [&](int i, int j) -> bool { return ok[static_cast<size_t>(j) * nx + i] != 0; };
  auto screen = [&](double x, double y) {
    const Point2 q = plot_to_screen(bounds, p0, w, h, x, y);
    return ImVec2(static_cast<float>(q.x), static_cast<float>(q.y));
  };
  /* Zero-crossing point along an edge between values fa, fb. */
  auto interp = [](double a, double fa, double b, double fb) -> double {
    const double denom = fb - fa;
    if (std::fabs(denom) < 1e-300) return 0.5 * (a + b);
    return a + (b - a) * (-fa / denom);
  };

  const float thickness = 2.2f;
  std::vector<ImVec2> seg_mid; /* midpoints for direction arrows */
  for (int j = 0; j < grid; ++j) {
    const double y0 = bounds.ymin + j * dy;
    const double y1 = y0 + dy;
    for (int i = 0; i < grid; ++i) {
      if (!good(i, j) || !good(i + 1, j) || !good(i + 1, j + 1) ||
          !good(i, j + 1))
        continue;
      const double x0 = bounds.xmin + i * dx;
      const double x1 = x0 + dx;
      const double f00 = node(i, j);       /* bottom-left  */
      const double f10 = node(i + 1, j);   /* bottom-right */
      const double f11 = node(i + 1, j + 1); /* top-right  */
      const double f01 = node(i, j + 1);   /* top-left     */

      /* Marching-squares case index (corner sign bits). */
      int c = 0;
      if (f00 > 0.0) c |= 1;
      if (f10 > 0.0) c |= 2;
      if (f11 > 0.0) c |= 4;
      if (f01 > 0.0) c |= 8;
      if (c == 0 || c == 15) continue; /* no crossing */

      /* Edge crossing points (only computed when needed):
       *   bottom: (x0..x1, y0)   right: (x1, y0..y1)
       *   top:    (x0..x1, y1)   left:  (x0, y0..y1) */
      const ImVec2 eb = screen(interp(x0, f00, x1, f10), y0);
      const ImVec2 er = screen(x1, interp(y0, f10, y1, f11));
      const ImVec2 et = screen(interp(x0, f01, x1, f11), y1);
      const ImVec2 el = screen(x0, interp(y0, f00, y1, f01));

      auto seg = [&](const ImVec2 &a, const ImVec2 &b) {
        draw->AddLine(a, b, color, thickness);
        seg_mid.push_back(ImVec2(0.5f * (a.x + b.x), 0.5f * (a.y + b.y)));
      };

      switch (c) {
        case 1: case 14: seg(el, eb); break;
        case 2: case 13: seg(eb, er); break;
        case 3: case 12: seg(el, er); break;
        case 4: case 11: seg(er, et); break;
        case 6: case 9:  seg(eb, et); break;
        case 7: case 8:  seg(el, et); break;
        /* Ambiguous saddle cells: connect both diagonals. */
        case 5:  seg(el, eb); seg(er, et); break;
        case 10: seg(eb, er); seg(el, et); break;
        default: break;
      }
    }
  }

  /* PHASE6-UI: pplane-style direction arrows ALONG the nullcline. On a
   * nullcline one component of the field is zero, so the flow is parallel
   * to one axis there; pplane marks this with little arrows. We sample
   * the full field at a thinned set of segment midpoints and draw a short
   * arrow in the flow direction, in the nullcline's color. */
  if (app.phase_show_nullcline_arrows && !seg_mid.empty()) {
    const double bx = bounds.xmax - bounds.xmin;
    const double by = bounds.ymax - bounds.ymin;
    const float arrow_px = 12.0f;
    const size_t stride = std::max<size_t>(1, seg_mid.size() / 24);
    for (size_t k = 0; k < seg_mid.size(); k += stride) {
      /* convert midpoint back to data coords */
      const double xd = bounds.xmin + (seg_mid[k].x - p0.x) / std::max(1.0f, w) * bx;
      const double yd = bounds.ymax - (seg_mid[k].y - p0.y) / std::max(1.0f, h) * by;
      State s2 = app.current;
      resize_state(s2, app.state_names.size());
      set_state_at(s2, ix, xd);
      set_state_at(s2, iy, yd);
      double vx = 0.0, vy = 0.0;
      char e2[128] = {0};
      if (!eval_phase_vector(app, s2, ix, iy, &vx, &vy, e2, sizeof(e2))) continue;
      const double mag = std::sqrt(vx * vx + vy * vy);
      if (mag < 1e-12 || !std::isfinite(mag)) continue;
      /* screen-space direction (y flips); normalize to a fixed length */
      const double sxv = vx / mag, syv = -vy / mag;
      const ImVec2 tip(seg_mid[k].x + static_cast<float>(sxv) * arrow_px,
                       seg_mid[k].y + static_cast<float>(syv) * arrow_px);
      draw_arrow(draw, seg_mid[k], tip, color);
    }
  }
}

void draw_trajectory_curve(ImDrawList *draw, const PlotBounds &bounds, const ImVec2 &p0,
                           float w, float h, const std::deque<State> &history,
                           size_t ix, size_t iy, uint32_t color, float thickness) {
  if (history.size() < 2) return;
  bool have = false;
  ImVec2 last{};
  for (const State &st : history) {
    const double x = state_at(st, ix);
    const double y = state_at(st, iy);
    if (!std::isfinite(x) || !std::isfinite(y)) { have = false; continue; }
    const Point2 p = plot_to_screen(bounds, p0, w, h, x, y);
    const ImVec2 now(static_cast<float>(p.x), static_cast<float>(p.y));
    if (have) draw->AddLine(last, now, static_cast<ImU32>(color), thickness);
    last = now;
    have = true;
  }
}

State phase_seed_state(const AppState &app, size_t ix, size_t iy, double x, double y) {
  State seed = app.current;
  resize_state(seed, app.state_names.size());
  set_state_at(seed, ix, x);
  set_state_at(seed, iy, y);
  seed.t = app.current.t;
  return seed;
}

void seed_phase_grid(AppState &app, const PlotBounds &bounds, size_t ix, size_t iy, int nx, int ny) {
  nx = std::max(2, std::min(nx, 12));
  ny = std::max(2, std::min(ny, 12));
  for (int gy = 0; gy < ny; ++gy) {
    for (int gx = 0; gx < nx; ++gx) {
      const double x = bounds.xmin + (gx + 0.5) * (bounds.xmax - bounds.xmin) / nx;
      const double y = bounds.ymin + (gy + 0.5) * (bounds.ymax - bounds.ymin) / ny;
      add_phase_trajectory(app, phase_seed_state(app, ix, iy, x, y), "grid");
    }
  }
}

void seed_phase_circle(AppState &app, const PlotBounds &bounds, size_t ix, size_t iy, int count) {
  count = std::max(4, std::min(count, 64));
  const double cx = 0.5 * (bounds.xmin + bounds.xmax);
  const double cy = 0.5 * (bounds.ymin + bounds.ymax);
  const double r = 0.25 * std::min(bounds.xmax - bounds.xmin, bounds.ymax - bounds.ymin);
  for (int i = 0; i < count; ++i) {
    const double a = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(count);
    add_phase_trajectory(app, phase_seed_state(app, ix, iy, cx + r * std::cos(a), cy + r * std::sin(a)), "circle");
  }
}

Point2 eigenvector_2x2(double a, double b, double c, double d, double lambda);

void draw_fixed_point_marker(AppState &app, ImDrawList *draw, const PlotBounds &bounds, const ImVec2 &p0,
                             float w, float h, size_t ix, size_t iy) {
  if (!app.fixed_ready || !app.phase_show_fixed_points) return;
  const Point2 p = plot_to_screen(bounds, p0, w, h, state_at(app.fixed_point, ix), state_at(app.fixed_point, iy));
  const ImVec2 q(static_cast<float>(p.x), static_cast<float>(p.y));
  ImU32 color = IM_COL32(255, 255, 255, 255);
  if (app.fixed_classification.find("saddle") != std::string::npos) color = IM_COL32(255, 95, 95, 255);
  else if (app.fixed_classification.find("stable") != std::string::npos) color = IM_COL32(120, 220, 135, 255);
  else if (app.fixed_classification.find("unstable") != std::string::npos) color = IM_COL32(255, 190, 80, 255);
  else if (app.fixed_classification.find("center") != std::string::npos) color = IM_COL32(120, 210, 255, 255);
  draw->AddCircle(q, 7.0f, color, 24, 2.0f);
  draw->AddLine(ImVec2(q.x - 5, q.y), ImVec2(q.x + 5, q.y), color, 1.3f);
  draw->AddLine(ImVec2(q.x, q.y - 5), ImVec2(q.x, q.y + 5), color, 1.3f);

  /* PHASE6: draw the linearized eigendirections at a 2D fixed point.
   * Each real eigenvalue's eigenvector becomes a short ray, colored by
   * stability (green = stable, orange = unstable). This turns the
   * classification into something you can see in the plane, the way a
   * hand-drawn phase portrait shows the slow/fast/saddle directions. */
  if (ix < 2 && iy < 2 && app.fixed_jacobian.size() == 4 &&
      app.fixed_eigenvalues.size() == 2) {
    const double a = app.fixed_jacobian[0], b = app.fixed_jacobian[1];
    const double c = app.fixed_jacobian[2], d = app.fixed_jacobian[3];
    const double ray_px = 46.0; /* on-screen ray half-length */
    for (int k = 0; k < 2; ++k) {
      const auto &ev = app.fixed_eigenvalues[static_cast<size_t>(k)];
      if (std::fabs(ev.second) > 1e-9) continue; /* complex: no real dir */
      const Point2 v = eigenvector_2x2(a, b, c, d, ev.first);
      /* map the data-space direction to screen direction (y flips) */
      const double sx = v.x / std::max(1e-12, bounds.xmax - bounds.xmin) * w;
      const double sy = -v.y / std::max(1e-12, bounds.ymax - bounds.ymin) * h;
      const double len = std::sqrt(sx * sx + sy * sy);
      if (len < 1e-9) continue;
      const double ux = sx / len, uy = sy / len;
      const ImU32 dir_col = ev.first < 0.0 ? IM_COL32(120, 220, 135, 230)
                                           : IM_COL32(255, 190, 80, 230);
      const ImVec2 a1(q.x + static_cast<float>(ux * ray_px),
                      q.y + static_cast<float>(uy * ray_px));
      const ImVec2 a2(q.x - static_cast<float>(ux * ray_px),
                      q.y - static_cast<float>(uy * ray_px));
      draw->AddLine(q, a1, dir_col, 2.0f);
      draw->AddLine(q, a2, dir_col, 2.0f);
      /* small arrowheads indicating flow direction (in for stable) */
      const double sgn = ev.first < 0.0 ? -1.0 : 1.0;
      const ImVec2 tip(q.x + static_cast<float>(sgn * ux * ray_px),
                       q.y + static_cast<float>(sgn * uy * ray_px));
      const double ah = 6.0;
      const ImVec2 b1(tip.x - static_cast<float>(sgn * ux * ah + uy * ah * 0.6),
                      tip.y - static_cast<float>(sgn * uy * ah - ux * ah * 0.6));
      const ImVec2 b2(tip.x - static_cast<float>(sgn * ux * ah - uy * ah * 0.6),
                      tip.y - static_cast<float>(sgn * uy * ah + ux * ah * 0.6));
      draw->AddLine(tip, b1, dir_col, 2.0f);
      draw->AddLine(tip, b2, dir_col, 2.0f);
    }
  }

  if (!app.fixed_classification.empty()) {
    draw->AddText(ImVec2(q.x + 9, q.y - 10), color, app.fixed_classification.c_str());
  }
}

/* PHASE6-UI: build a planar field for the current axis pair (other state
 * variables held at the app's current values) and rescan all fixed
 * points when the view, axis pair, or parameters have changed. Throttled
 * by signature comparison so it isn't recomputed every frame. */
void maybe_rescan_fixed_points(AppState &app, const PlotBounds &b, size_t ix,
                               size_t iy) {
  if (!app.phase_auto_fixed_points || app.mode != SystemMode::ODE ||
      app.state_names.size() < 2) {
    app.phase_fixed_points.clear();
    return;
  }
  double psig = 0.0;
  for (double v : app.param_values) psig += v;

  const bool same =
      app.phase_fp_scan_px == static_cast<int>(ix) &&
      app.phase_fp_scan_py == static_cast<int>(iy) &&
      std::fabs(app.phase_fp_scan_xmin - b.xmin) < 1e-9 * (1 + std::fabs(b.xmin)) &&
      std::fabs(app.phase_fp_scan_xmax - b.xmax) < 1e-9 * (1 + std::fabs(b.xmax)) &&
      std::fabs(app.phase_fp_scan_ymin - b.ymin) < 1e-9 * (1 + std::fabs(b.ymin)) &&
      std::fabs(app.phase_fp_scan_ymax - b.ymax) < 1e-9 * (1 + std::fabs(b.ymax)) &&
      std::fabs(app.phase_fp_scan_param_sig - psig) < 1e-12 * (1 + std::fabs(psig));
  if (same) return;

  dynsys::analysis::PlanarField field;
  field.eval = [&app, ix, iy](double x, double y, double *u, double *v) -> bool {
    State s = app.current;
    resize_state(s, app.state_names.size());
    set_state_at(s, ix, x);
    set_state_at(s, iy, y);
    double vx = 0.0, vy = 0.0;
    char err[128] = {0};
    if (!eval_phase_vector(app, s, ix, iy, &vx, &vy, err, sizeof(err)))
      return false;
    *u = vx;
    *v = vy;
    return true;
  };
  app.phase_fixed_points = dynsys::analysis::scan_fixed_points_2d(
      field, b.xmin, b.xmax, b.ymin, b.ymax, 13);
  app.phase_fp_scan_px = static_cast<int>(ix);
  app.phase_fp_scan_py = static_cast<int>(iy);
  app.phase_fp_scan_xmin = b.xmin; app.phase_fp_scan_xmax = b.xmax;
  app.phase_fp_scan_ymin = b.ymin; app.phase_fp_scan_ymax = b.ymax;
  app.phase_fp_scan_param_sig = psig;
}

/* Draw all auto-scanned fixed points with stability color, label, and
 * (for real eigenvalues) the stable/unstable eigendirection manifolds. */
void draw_auto_fixed_points(AppState &app, ImDrawList *draw,
                            const PlotBounds &bounds, const ImVec2 &p0, float w,
                            float h) {
  for (const auto &fp : app.phase_fixed_points) {
    const Point2 ps = plot_to_screen(bounds, p0, w, h, fp.x, fp.y);
    const ImVec2 q(static_cast<float>(ps.x), static_cast<float>(ps.y));
    ImU32 color = IM_COL32(230, 230, 235, 255);
    if (fp.is_saddle) color = IM_COL32(255, 95, 95, 255);
    else if (fp.label.find("stable") != std::string::npos &&
             fp.label.find("unstable") == std::string::npos)
      color = IM_COL32(120, 220, 135, 255);
    else if (fp.label.find("unstable") != std::string::npos)
      color = IM_COL32(255, 190, 80, 255);
    else if (fp.label.find("center") != std::string::npos)
      color = IM_COL32(120, 210, 255, 255);

    /* eigendirection manifolds */
    if (app.phase_show_manifolds && fp.directions.size() == 2 &&
        fp.eigenvalues.size() == 2) {
      const double ray = 60.0;
      for (size_t k = 0; k < 2; ++k) {
        const double lam = fp.eigenvalues[k].real();
        const auto &d = fp.directions[k];
        const double sx = d.first / std::max(1e-12, bounds.xmax - bounds.xmin) * w;
        const double sy = -d.second / std::max(1e-12, bounds.ymax - bounds.ymin) * h;
        const double n = std::hypot(sx, sy);
        if (n < 1e-9) continue;
        const double ux = sx / n, uy = sy / n;
        const ImU32 mc = lam < 0.0 ? IM_COL32(120, 220, 135, 200)
                                   : IM_COL32(255, 150, 80, 200);
        draw->AddLine(ImVec2(q.x + static_cast<float>(ux * ray), q.y + static_cast<float>(uy * ray)),
                      ImVec2(q.x - static_cast<float>(ux * ray), q.y - static_cast<float>(uy * ray)),
                      mc, 1.6f);
      }
    }

    draw->AddCircleFilled(q, 5.0f, color);
    draw->AddCircle(q, 5.0f, IM_COL32(20, 20, 24, 255), 16, 1.5f);
    draw->AddText(ImVec2(q.x + 8, q.y - 6), color, fp.label.c_str());
  }
}

[[maybe_unused]] void draw_phase_plane(AppState &app, ImVec2 size) {
  if (app.state_names.empty()) {
    ImGui::TextDisabled("No state variables declared.");
    return;
  }
  app.phase_x_index = std::max(0, std::min(app.phase_x_index, static_cast<int>(app.state_names.size() - 1)));
  app.phase_y_index = std::max(0, std::min(app.phase_y_index, static_cast<int>(app.state_names.size() - 1)));
  if (app.state_names.size() > 1 && app.phase_y_index == app.phase_x_index) {
    app.phase_y_index = (app.phase_x_index + 1) % static_cast<int>(app.state_names.size());
  }

  if (ImGui::BeginCombo("x axis", app.state_names[static_cast<size_t>(app.phase_x_index)].c_str())) {
    for (int i = 0; i < static_cast<int>(app.state_names.size()); ++i) {
      if (ImGui::Selectable(app.state_names[static_cast<size_t>(i)].c_str(), i == app.phase_x_index)) app.phase_x_index = i;
    }
    ImGui::EndCombo();
  }
  ImGui::SameLine();
  if (ImGui::BeginCombo("y axis", app.state_names[static_cast<size_t>(app.phase_y_index)].c_str())) {
    for (int i = 0; i < static_cast<int>(app.state_names.size()); ++i) {
      if (ImGui::Selectable(app.state_names[static_cast<size_t>(i)].c_str(), i == app.phase_y_index)) app.phase_y_index = i;
    }
    ImGui::EndCombo();
  }
  ImGui::Checkbox("trajectory", &app.phase_show_trajectory); ImGui::SameLine();
  ImGui::Checkbox("vector field", &app.phase_show_vector_field); ImGui::SameLine();
  ImGui::Checkbox("nullclines", &app.phase_show_nullclines);
  ImGui::Checkbox("auto fixed points", &app.phase_auto_fixed_points); ImGui::SameLine();
  ImGui::Checkbox("manifolds", &app.phase_show_manifolds); ImGui::SameLine();
  ImGui::Checkbox("separatrices", &app.phase_show_separatrices);
  ImGui::Checkbox("equal aspect", &app.phase_equal_aspect);
  ImGui::Checkbox("auto bounds", &app.phase_auto_bounds); ImGui::SameLine();
  ImGui::Checkbox("normalized vectors", &app.phase_normalize_vectors); ImGui::SameLine();
  ImGui::SliderInt("grid", &app.phase_grid, 8, 40);
  if (app.phase_auto_fixed_points && app.mode == SystemMode::ODE) {
    ImGui::TextDisabled("equilibria in view: %zu (auto-found & classified)",
                        app.phase_fixed_points.size());
  }
  if (!app.phase_auto_bounds) {
    ImGui::InputFloat("x min", &app.phase_x_min, 0.1f, 1.0f, "%.4g"); ImGui::SameLine();
    ImGui::InputFloat("x max", &app.phase_x_max, 0.1f, 1.0f, "%.4g");
    ImGui::InputFloat("y min", &app.phase_y_min, 0.1f, 1.0f, "%.4g"); ImGui::SameLine();
    ImGui::InputFloat("y max", &app.phase_y_max, 0.1f, 1.0f, "%.4g");
  }

  /* CRASH FIX (same class as the basins one): clamp the plane indices to the
   * current state dimension before using them to index state_names / states.
   * A stale phase_y_index left over from a higher-dimensional system would
   * otherwise read out of bounds after switching to a smaller system. */
  if (!app.state_names.empty()) {
    app.phase_x_index = std::max(0, std::min(app.phase_x_index, (int)app.state_names.size() - 1));
    app.phase_y_index = std::max(0, std::min(app.phase_y_index, (int)app.state_names.size() - 1));
  }
  const size_t ix = static_cast<size_t>(app.phase_x_index);
  const size_t iy = static_cast<size_t>(app.phase_y_index);
  const float w = size.x <= 0.0f ? ImGui::GetContentRegionAvail().x : size.x;
  const float h = size.y;
  PlotBounds bounds = current_phase_bounds(app, ix, iy);
  if (app.phase_equal_aspect) bounds = equalize_bounds_for_canvas(bounds, w, h);

  /* PHASE6-UI: auto-find every equilibrium in the current view. */
  maybe_rescan_fixed_points(app, bounds, ix, iy);

  ImGui::TextDisabled("L-click: restart orbit | shift+L: add orbit | R-drag: pan | wheel or +/-: zoom | dbl-click: auto-fit");
  ImGui::InvisibleButton("phase_canvas", ImVec2(w, h), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
  const bool hovered = ImGui::IsItemHovered();
  const bool active = ImGui::IsItemActive();
  ImDrawList *draw = ImGui::GetWindowDrawList();
  const ImVec2 p0 = ImGui::GetItemRectMin();
  const ImVec2 p1 = ImGui::GetItemRectMax();
  const ImVec2 mouse = ImGui::GetIO().MousePos;
  const Point2 mouse_data = screen_to_plot(bounds, p0, w, h, mouse);

  if (hovered) {
    ImGuiIO &io = ImGui::GetIO();
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
      app.phase_auto_bounds = true;
      app.phase_bounds_valid = false; /* snap to data on explicit auto-fit */
    } else if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      State seed = phase_seed_state(app, ix, iy, mouse_data.x, mouse_data.y);
      if (io.KeyShift) {
        add_phase_trajectory(app, seed);
      } else {
        restart_main_from_state(app, seed);
      }
    }
    if (io.MouseWheel != 0.0f) {
      /* wheel up -> zoom in (shrink span). 0.8 per notch is responsive
       * without being twitchy. */
      const double factor = std::pow(0.8, clamped_wheel(static_cast<double>(io.MouseWheel)));
      zoom_phase_bounds(app, bounds, mouse_data, factor);
    }
    /* PHASE6-UI: keyboard zoom for the phase plane, about the view
     * center. +/= zoom in, -/_ zoom out. (The arrow/WASD camera keys
     * only affect the 3D view; the 2D plane has its own bounds.) */
    const Point2 view_center{0.5 * (bounds.xmin + bounds.xmax),
                             0.5 * (bounds.ymin + bounds.ymax)};
    if (ImGui::IsKeyPressed(ImGuiKey_Equal) || ImGui::IsKeyPressed(ImGuiKey_KeypadAdd))
      zoom_phase_bounds(app, bounds, view_center, 0.83);
    if (ImGui::IsKeyPressed(ImGuiKey_Minus) || ImGui::IsKeyPressed(ImGuiKey_KeypadSubtract))
      zoom_phase_bounds(app, bounds, view_center, 1.20);
  }
  if (active && ImGui::IsMouseDragging(ImGuiMouseButton_Right)) {
    ImGuiIO &io = ImGui::GetIO();
    const double dx = -static_cast<double>(io.MouseDelta.x) / std::max(1.0f, w) * (bounds.xmax - bounds.xmin);
    const double dy = static_cast<double>(io.MouseDelta.y) / std::max(1.0f, h) * (bounds.ymax - bounds.ymin);
    pan_phase_bounds(app, bounds, dx, dy);
  }

  draw_plot_grid(draw, p0, w, h, bounds);
  draw->PushClipRect(p0, p1, true);

  if (app.phase_show_vector_field) {
    const int grid = std::max(4, std::min(app.phase_grid, 40));
    char err[256] = {0};
    const double bx = bounds.xmax - bounds.xmin;
    const double by = bounds.ymax - bounds.ymin;
    for (int gy = 0; gy < grid; ++gy) {
      for (int gx = 0; gx < grid; ++gx) {
        const double x = bounds.xmin + (gx + 0.5) * bx / grid;
        const double y = bounds.ymin + (gy + 0.5) * by / grid;
        State st = app.current;
        resize_state(st, app.state_names.size());
        set_state_at(st, ix, x);
        set_state_at(st, iy, y);
        double vx = 0.0, vy = 0.0;
        if (!eval_phase_vector(app, st, ix, iy, &vx, &vy, err, sizeof(err))) continue;
        const double len = std::sqrt(vx * vx + vy * vy);
        if (len <= 1e-14 || !std::isfinite(len)) continue;
        const double max_len_data = 0.32 * std::min(bx, by) / grid;
        if (app.phase_normalize_vectors) {
          vx = vx / len * max_len_data;
          vy = vy / len * max_len_data;
        } else {
          const double scale = max_len_data / (1.0 + len);
          vx *= scale;
          vy *= scale;
        }
        const Point2 a = plot_to_screen(bounds, p0, w, h, x - 0.5 * vx, y - 0.5 * vy);
        const Point2 b = plot_to_screen(bounds, p0, w, h, x + 0.5 * vx, y + 0.5 * vy);
        draw_arrow(draw,
                   ImVec2(static_cast<float>(a.x), static_cast<float>(a.y)),
                   ImVec2(static_cast<float>(b.x), static_cast<float>(b.y)),
                   IM_COL32(95, 130, 165, 140));
      }
    }
  }

  if (app.phase_show_nullclines && app.mode == SystemMode::ODE) {
    /* pplane convention: one nullcline per state variable, where that
     * variable's rate is zero. Solid, saturated colors; a legend that
     * names the actual variables. */
    const ImU32 col_x = IM_COL32(255, 110, 80, 235);
    const ImU32 col_y = IM_COL32(90, 210, 130, 235);
    draw_nullcline_for_component(app, draw, bounds, p0, w, h, ix, iy, ix, col_x);
    if (iy != ix)
      draw_nullcline_for_component(app, draw, bounds, p0, w, h, ix, iy, iy, col_y);
    char legend[160];
    std::snprintf(legend, sizeof(legend),
                  "nullclines —  %s' = 0 (red)   %s' = 0 (green)",
                  app.state_names[ix].c_str(), app.state_names[iy].c_str());
    draw->AddText(ImVec2(p0.x + 8, p0.y + 8), IM_COL32(235, 235, 240, 230),
                  legend);
  }

  if (app.phase_show_trajectory) {
    draw_trajectory_curve(draw, bounds, p0, w, h, app.history, ix, iy, 0xff5fcdffu, 1.6f);
    for (const auto &traj : app.phase_trajectories) {
      if (traj.visible) draw_trajectory_curve(draw, bounds, p0, w, h, traj.history, ix, iy, fade_color(traj.color, 230), 1.2f);
    }
  }
  if (app.phase_show_separatrices) {
    for (const auto &curve : app.separatrix_curves) {
      if (curve.visible) draw_trajectory_curve(draw, bounds, p0, w, h, curve.history, ix, iy, fade_color(curve.color, 245), 1.8f);
    }
  }

  /* PHASE6-UI: all auto-scanned equilibria with manifolds (pplane-style).
   * Falls back to the single Newton-found marker if auto-scan is off. */
  if (app.phase_auto_fixed_points)
    draw_auto_fixed_points(app, draw, bounds, p0, w, h);
  else
    draw_fixed_point_marker(app, draw, bounds, p0, w, h, ix, iy);
  {
    const double cx = state_at(app.current, ix), cy = state_at(app.current, iy);
    if (std::isfinite(cx) && std::isfinite(cy)) {
      const Point2 cp = plot_to_screen(bounds, p0, w, h, cx, cy);
      if (std::isfinite(cp.x) && std::isfinite(cp.y))
        draw->AddCircleFilled(ImVec2(static_cast<float>(cp.x), static_cast<float>(cp.y)), 4.0f, IM_COL32(255, 255, 255, 255));
    }
  }

  if (hovered) {
    char coord[160];
    std::snprintf(coord, sizeof(coord), "%s=%.6g, %s=%.6g", app.state_names[ix].c_str(), mouse_data.x, app.state_names[iy].c_str(), mouse_data.y);
    draw->AddText(ImVec2(mouse.x + 12, mouse.y + 12), IM_COL32(230, 230, 235, 240), coord);
  }
  draw->PopClipRect();

  char label[256];
  std::snprintf(label, sizeof(label), "%s=[%.5g, %.5g], %s=[%.5g, %.5g] | extra orbits: %zu | separatrix curves: %zu",
                app.state_names[ix].c_str(), bounds.xmin, bounds.xmax,
                app.state_names[iy].c_str(), bounds.ymin, bounds.ymax,
                app.phase_trajectories.size(), app.separatrix_curves.size());
  draw->AddText(ImVec2(p0.x + 8, p0.y + h - 22), IM_COL32(200, 200, 210, 220), label);

  if (ImGui::Button("Add current as extra orbit")) add_phase_trajectory(app, app.current);
  ImGui::SameLine();
  if (ImGui::Button("Seed 5x5 grid")) seed_phase_grid(app, bounds, ix, iy, 5, 5);
  ImGui::SameLine();
  if (ImGui::Button("Seed circle")) seed_phase_circle(app, bounds, ix, iy, 16);
  ImGui::SameLine();
  if (ImGui::Button("Clear extra orbits")) app.phase_trajectories.clear();

  if (!app.phase_trajectories.empty() && ImGui::TreeNode("extra orbits")) {
    for (size_t i = 0; i < app.phase_trajectories.size(); ++i) {
      auto &traj = app.phase_trajectories[i];
      ImGui::PushID(static_cast<int>(i));
      ImGui::Checkbox("visible", &traj.visible); ImGui::SameLine();
      ImGui::Checkbox("active", &traj.active); ImGui::SameLine();
      ImGui::Text("%s (%zu points)", traj.label.c_str(), traj.history.size());
      ImGui::PopID();
    }
    ImGui::TreePop();
  }
}

/* ============================================================
 * PHASE6-UI: fullscreen background phase plane.
 *
 * Renders the 2D phase portrait across the whole window using the
 * background draw list (behind all panels), and handles interaction by
 * reading the mouse directly — but only when no ImGui panel wants it
 * (io.WantCaptureMouse). This is why dragging/zooming the plot no longer
 * moves any window: the plot isn't inside a window at all.
 *
 * Standard pplane-style mouse scheme:
 *   left-click           : add an orbit through the clicked point
 *   left-drag            : pan the view
 *   wheel / + / -        : zoom (about cursor for wheel, center for keys)
 *   double-click         : auto-fit to data
 * Orbits accumulate; they are only removed by an explicit Clear.
 * ============================================================ */

/* ============================================================
 * PHASE-A: 1D view, for systems with a single effective dimension
 * (tent map, logistic-as-1D, 1D ODEs).
 *
 * For MAPS: a cobweb (Verhulst) diagram — plots the identity line y=x,
 * the map curve y=f(x), and the staircase x0 -> f(x0) -> f(f(x0)) ...
 * that shows convergence to a fixed point, a period-n cycle, or chaos.
 *
 * For 1D ODEs: a flow-on-a-line — plots y = f(x) (the velocity), the
 * zero line, and marks each equilibrium (f=0) as stable (filled) or
 * unstable (hollow) with little flow arrows on the x-axis.
 *
 * The mouse sets the seed: click an x to start the cobweb / drop the
 * state there. Wheel / +/- zoom the x-range; double-click auto-fits.
 * ============================================================ */
void render_line1d_background(AppState &app) {
  if (app.state_names.empty()) return;
  const size_t ix = 0; /* the single state variable */

  const ImVec2 p0(0.0f, 0.0f);
  const float w = static_cast<float>(app.window_width);
  const float h = static_cast<float>(app.window_height);
  const ImVec2 p1(w, h);

  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  ImGuiIO &io = ImGui::GetIO();
  const bool over_plot = !io.WantCaptureMouse;

  /* x-range: reuse phase x bounds machinery on the single axis; map the
   * y-axis to the same range for the cobweb so the y=x line is diagonal. */
  PlotBounds b = current_phase_bounds(app, ix, ix);
  /* current_phase_bounds with ix==iy collapses y to x's range, which is
   * exactly what we want for a cobweb (square domain). */
  b = sanitize_bounds(b);

  const Point2 mouse_data = screen_to_plot(b, p0, w, h, io.MousePos);

  /* interaction */
  if (over_plot) {
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
      app.phase_auto_bounds = true; app.phase_bounds_valid = false;
    } else if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      /* set the live state's x to the clicked value and restart */
      State seed = app.start;
      resize_state(seed, app.state_names.size());
      set_state_at(seed, ix, mouse_data.x);
      restart_main_from_state(app, seed);
    }
    if (io.MouseWheel != 0.0f)
      zoom_phase_bounds(app, b, mouse_data, std::pow(0.8, clamped_wheel((double)io.MouseWheel)));
    const Point2 c{0.5 * (b.xmin + b.xmax), 0.5 * (b.ymin + b.ymax)};
    if (ImGui::IsKeyPressed(ImGuiKey_Equal) || ImGui::IsKeyPressed(ImGuiKey_KeypadAdd))
      zoom_phase_bounds(app, b, c, 0.83);
    if (ImGui::IsKeyPressed(ImGuiKey_Minus) || ImGui::IsKeyPressed(ImGuiKey_KeypadSubtract))
      zoom_phase_bounds(app, b, c, 1.20);
  }

  /* backdrop + axes */
  draw->AddRectFilled(p0, p1, IM_COL32(16, 17, 22, 255));
  draw_plot_grid(draw, p0, w, h, b);
  draw->PushClipRect(p0, p1, true);

  auto S = [&](double x, double y) {
    const Point2 q = plot_to_screen(b, p0, w, h, x, y);
    return ImVec2((float)q.x, (float)q.y);
  };
  char err[256] = {0};

  if (app.mode == SystemMode::Map) {
    /* identity line y = x */
    draw->AddLine(S(b.xmin, b.xmin), S(b.xmax, b.xmax), IM_COL32(120, 120, 130, 200), 1.4f);
    /* map curve y = f(x) (here f is the next-value of the single state) */
    const int N = std::max(64, std::min((int)w, 1200));
    ImVec2 prev; bool have_prev = false;
    for (int i = 0; i <= N; ++i) {
      const double x = b.xmin + (b.xmax - b.xmin) * (double)i / N;
      State s = app.current; resize_state(s, app.state_names.size());
      set_state_at(s, ix, x);
      State nx{};
      if (!step_state(app, s, &nx, err, sizeof(err))) { have_prev = false; continue; }
      const double y = state_at(nx, ix);
      if (!std::isfinite(y)) { have_prev = false; continue; }
      const ImVec2 cur = S(x, y);
      if (have_prev) draw->AddLine(prev, cur, IM_COL32(255, 170, 70, 235), 2.0f);
      prev = cur; have_prev = true;
    }
    /* cobweb staircase from the current/seed x */
    double x = state_at(app.current, ix);
    if (std::isfinite(x)) {
      const int steps = 200;
      for (int i = 0; i < steps; ++i) {
        State s = app.current; resize_state(s, app.state_names.size());
        set_state_at(s, ix, x);
        State nx{};
        if (!step_state(app, s, &nx, err, sizeof(err))) break;
        const double y = state_at(nx, ix);
        if (!std::isfinite(y)) break;
        /* vertical to the curve, then horizontal to the diagonal */
        const ImU32 col = IM_COL32(120, 210, 255, 150);
        draw->AddLine(S(x, x), S(x, y), col, 1.2f);
        draw->AddLine(S(x, y), S(y, y), col, 1.2f);
        x = y;
      }
    }
    draw->AddText(ImVec2(p0.x + 12, p0.y + 12), IM_COL32(235, 235, 240, 230),
                  "cobweb: orange y=f(x), grey y=x, blue staircase. Click to set start x.");
  } else {
    /* 1D ODE: flow on a line. Plot velocity y = f(x). */
    const double y0 = 0.0;
    draw->AddLine(S(b.xmin, y0), S(b.xmax, y0), IM_COL32(120, 120, 130, 200), 1.4f);
    const int N = std::max(64, std::min((int)w, 1200));
    ImVec2 prev; bool have_prev = false;
    double prevf = 0.0; bool prevok = false;
    for (int i = 0; i <= N; ++i) {
      const double x = b.xmin + (b.xmax - b.xmin) * (double)i / N;
      double f = 0.0;
      State s = app.current; resize_state(s, app.state_names.size());
      set_state_at(s, ix, x);
      const bool ok = eval_phase_component(app, s, ix, &f, err, sizeof(err)) && std::isfinite(f);
      if (ok) {
        const ImVec2 cur = S(x, f);
        if (have_prev) draw->AddLine(prev, cur, IM_COL32(255, 170, 70, 235), 2.0f);
        prev = cur; have_prev = true;
        /* equilibrium where f changes sign */
        if (prevok && ((prevf <= 0 && f >= 0) || (prevf >= 0 && f <= 0))) {
          const double xr = x; /* approximate root location */
          const bool stable = (prevf > 0 && f < 0); /* f decreasing through 0 */
          const ImVec2 q = S(xr, 0.0);
          if (stable) draw->AddCircleFilled(q, 5.0f, IM_COL32(120, 220, 135, 255));
          else { draw->AddCircle(q, 5.0f, IM_COL32(255, 150, 80, 255), 16, 2.0f); }
        }
      } else have_prev = false;
      prevf = f; prevok = ok;
    }
    /* flow arrows along the x-axis: sign of f pushes left/right */
    const int A = 24;
    for (int i = 0; i < A; ++i) {
      const double x = b.xmin + (b.xmax - b.xmin) * (i + 0.5) / A;
      double f = 0.0;
      State s = app.current; resize_state(s, app.state_names.size());
      set_state_at(s, ix, x);
      if (!eval_phase_component(app, s, ix, &f, err, sizeof(err)) || !std::isfinite(f) || f == 0) continue;
      const ImVec2 base = S(x, 0.0);
      const float dir = f > 0 ? 10.0f : -10.0f;
      draw_arrow(draw, base, ImVec2(base.x + dir, base.y), IM_COL32(150, 180, 210, 200));
    }
    /* live state marker */
    const double xs = state_at(app.current, ix);
    if (std::isfinite(xs)) draw->AddCircleFilled(S(xs, 0.0), 4.0f, IM_COL32(255, 255, 255, 255));
    draw->AddText(ImVec2(p0.x + 12, p0.y + 12), IM_COL32(235, 235, 240, 230),
                  "flow on a line: orange y=f(x); filled=stable, hollow=unstable eq. Click to place state.");
  }

  if (over_plot) {
    char coord[96];
    std::snprintf(coord, sizeof(coord), "%s=%.5g", app.state_names[ix].c_str(), mouse_data.x);
    draw->AddText(ImVec2(io.MousePos.x + 14, io.MousePos.y + 14), IM_COL32(230, 230, 235, 240), coord);
  }
  draw->PopClipRect();
}

void render_phase_background(AppState &app) {
  if (app.state_names.size() < 2) return;
  app.phase_x_index = std::max(0, std::min(app.phase_x_index, static_cast<int>(app.state_names.size() - 1)));
  app.phase_y_index = std::max(0, std::min(app.phase_y_index, static_cast<int>(app.state_names.size() - 1)));
  if (app.phase_y_index == app.phase_x_index)
    app.phase_y_index = (app.phase_x_index + 1) % static_cast<int>(app.state_names.size());
  const size_t ix = static_cast<size_t>(app.phase_x_index);
  const size_t iy = static_cast<size_t>(app.phase_y_index);

  const ImVec2 p0(0.0f, 0.0f);
  const float w = static_cast<float>(app.window_width);
  const float h = static_cast<float>(app.window_height);
  const ImVec2 p1(w, h);

  PlotBounds bounds = current_phase_bounds(app, ix, iy);
  if (app.phase_equal_aspect) bounds = equalize_bounds_for_canvas(bounds, w, h);
  maybe_rescan_fixed_points(app, bounds, ix, iy);

  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  ImGuiIO &io = ImGui::GetIO();
  const ImVec2 mouse = io.MousePos;
  const Point2 mouse_data = screen_to_plot(bounds, p0, w, h, mouse);
  const bool over_plot = !io.WantCaptureMouse; /* cursor not over a panel */

  /* --- interaction (only when the cursor is over the plot, not UI) --- */
  if (over_plot) {
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
      app.phase_auto_bounds = true;
      app.phase_bounds_valid = false;
    } else if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !io.KeyShift) {
      /* standard pplane: click drops an orbit through the point.
       * (Shift+click is reserved for the live "main" orbit restart.) */
      add_phase_trajectory(app, phase_seed_state(app, ix, iy, mouse_data.x, mouse_data.y));
    } else if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && io.KeyShift) {
      restart_main_from_state(app, phase_seed_state(app, ix, iy, mouse_data.x, mouse_data.y));
    }
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) &&
        !ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      const double dx = -static_cast<double>(io.MouseDelta.x) / std::max(1.0f, w) * (bounds.xmax - bounds.xmin);
      const double dy = static_cast<double>(io.MouseDelta.y) / std::max(1.0f, h) * (bounds.ymax - bounds.ymin);
      pan_phase_bounds(app, bounds, dx, dy);
    }
    if (io.MouseWheel != 0.0f)
      zoom_phase_bounds(app, bounds, mouse_data, std::pow(0.8, clamped_wheel(static_cast<double>(io.MouseWheel))));
    const Point2 center{0.5 * (bounds.xmin + bounds.xmax), 0.5 * (bounds.ymin + bounds.ymax)};
    if (ImGui::IsKeyPressed(ImGuiKey_Equal) || ImGui::IsKeyPressed(ImGuiKey_KeypadAdd))
      zoom_phase_bounds(app, bounds, center, 0.83);
    if (ImGui::IsKeyPressed(ImGuiKey_Minus) || ImGui::IsKeyPressed(ImGuiKey_KeypadSubtract))
      zoom_phase_bounds(app, bounds, center, 1.20);
  }

  /* --- drawing --- */
  draw->AddRectFilled(p0, p1, IM_COL32(16, 17, 22, 255)); /* opaque backdrop */
  draw_plot_grid(draw, p0, w, h, bounds);
  draw->PushClipRect(p0, p1, true);

  if (app.phase_show_vector_field) {
    const int grid = std::max(6, std::min(app.phase_grid, 48));
    char err[256] = {0};
    const double bx = bounds.xmax - bounds.xmin;
    const double by = bounds.ymax - bounds.ymin;
    /* first pass: sample the field and record the max speed, so the
     * speed->color mapping is scaled to THIS view instead of a hard-coded
     * constant (the old mapping collapsed the speed term and tinted nearly
     * everything the same). */
    const double max_len_data = 0.34 * std::min(bx, by) / grid;
    std::vector<double> fvx((size_t)grid * grid, 0.0), fvy((size_t)grid * grid, 0.0);
    std::vector<unsigned char> fok((size_t)grid * grid, 0);
    double speed_max = 1e-12;
    for (int gy = 0; gy < grid; ++gy) {
      for (int gx = 0; gx < grid; ++gx) {
        const double x = bounds.xmin + (gx + 0.5) * bx / grid;
        const double y = bounds.ymin + (gy + 0.5) * by / grid;
        State st = app.current; resize_state(st, app.state_names.size());
        set_state_at(st, ix, x); set_state_at(st, iy, y);
        double vx = 0.0, vy = 0.0;
        if (!eval_phase_vector(app, st, ix, iy, &vx, &vy, err, sizeof(err))) continue;
        const double len = std::sqrt(vx * vx + vy * vy);
        if (len <= 1e-14 || !std::isfinite(len)) continue;
        const size_t k = (size_t)gy * grid + gx;
        fvx[k] = vx; fvy[k] = vy; fok[k] = 1;
        if (len > speed_max) speed_max = len;
      }
    }
    const double inv_speed = 1.0 / speed_max;
    for (int gy = 0; gy < grid; ++gy) {
      for (int gx = 0; gx < grid; ++gx) {
        const size_t k = (size_t)gy * grid + gx;
        if (!fok[k]) continue;
        const double x = bounds.xmin + (gx + 0.5) * bx / grid;
        const double y = bounds.ymin + (gy + 0.5) * by / grid;
        const double vx = fvx[k], vy = fvy[k];
        const double len = std::sqrt(vx * vx + vy * vy);
        double nvx = vx, nvy = vy;
        if (app.phase_normalize_vectors) { nvx = vx / len * max_len_data; nvy = vy / len * max_len_data; }
        else { const double s = max_len_data / (1.0 + len); nvx *= s; nvy *= s; }
        /* speed-colored: slow = blue, fast = warm; sqrt for visual balance */
        const double tcol = std::min(1.0, std::sqrt(len * inv_speed));
        const int rr = static_cast<int>(60 + 180 * tcol);
        const int gg = static_cast<int>(140 - 30 * tcol);
        const int bb = static_cast<int>(210 - 140 * tcol);
        const Point2 a = plot_to_screen(bounds, p0, w, h, x - 0.5 * nvx, y - 0.5 * nvy);
        const Point2 b = plot_to_screen(bounds, p0, w, h, x + 0.5 * nvx, y + 0.5 * nvy);
        draw_arrow(draw, ImVec2((float)a.x, (float)a.y), ImVec2((float)b.x, (float)b.y),
                   IM_COL32(rr, gg, bb, 160));
      }
    }
  }

  if (app.phase_show_nullclines && app.mode == SystemMode::ODE) {
    draw_nullcline_for_component(app, draw, bounds, p0, w, h, ix, iy, ix, IM_COL32(255, 110, 80, 235));
    if (iy != ix)
      draw_nullcline_for_component(app, draw, bounds, p0, w, h, ix, iy, iy, IM_COL32(90, 210, 130, 235));
    char legend[160];
    std::snprintf(legend, sizeof(legend), "nullclines:  %s' = 0 (red)   %s' = 0 (green)",
                  app.state_names[ix].c_str(), app.state_names[iy].c_str());
    draw->AddText(ImVec2(p0.x + 12, p0.y + 40), IM_COL32(235, 235, 240, 230), legend);
  }
  if (app.mode == SystemMode::Map && app.phase_show_nullclines) {
    /* Maps have no vector field, so nullclines / equilibria in the
     * continuous sense don't apply — make that explicit instead of
     * silently drawing nothing. */
    draw->AddText(ImVec2(p0.x + 12, p0.y + 40), IM_COL32(220, 200, 140, 230),
                  "nullclines/fixed-point fields apply to ODEs; this is a map "
                  "(showing iterated points)");
  }

  if (app.phase_show_trajectory) {
    draw_trajectory_curve(draw, bounds, p0, w, h, app.history, ix, iy, 0xff5fcdffu, 1.8f);
    for (const auto &traj : app.phase_trajectories)
      if (traj.visible) draw_trajectory_curve(draw, bounds, p0, w, h, traj.history, ix, iy, fade_color(traj.color, 235), 1.4f);
  }
  if (app.phase_show_separatrices)
    for (const auto &curve : app.separatrix_curves)
      if (curve.visible) draw_trajectory_curve(draw, bounds, p0, w, h, curve.history, ix, iy, fade_color(curve.color, 245), 1.8f);

  if (app.phase_auto_fixed_points) draw_auto_fixed_points(app, draw, bounds, p0, w, h);
  else draw_fixed_point_marker(app, draw, bounds, p0, w, h, ix, iy);

  const double cx = state_at(app.current, ix), cy = state_at(app.current, iy);
  if (std::isfinite(cx) && std::isfinite(cy)) {
    const Point2 cp = plot_to_screen(bounds, p0, w, h, cx, cy);
    if (std::isfinite(cp.x) && std::isfinite(cp.y))
      draw->AddCircleFilled(ImVec2((float)cp.x, (float)cp.y), 4.0f, IM_COL32(255, 255, 255, 255));
  }

  if (over_plot) {
    char coord[160];
    std::snprintf(coord, sizeof(coord), "%s=%.5g, %s=%.5g", app.state_names[ix].c_str(), mouse_data.x, app.state_names[iy].c_str(), mouse_data.y);
    draw->AddText(ImVec2(mouse.x + 14, mouse.y + 14), IM_COL32(230, 230, 235, 240), coord);
  }
  draw->PopClipRect();
}

/* PHASE-A: bifurcation diagram as a full-window background view (the same
 * data the side-panel diagram shows, but filling the plot area). */
/* ============================================================
 * PHASE C: escape-time fractal view — the dynamical-systems <-> fractal
 * bridge.
 *
 * Iterates the system's 2D map over a region of the plane. In
 * ParameterSpace mode each pixel is a value of two chosen parameters
 * (c on the plane) with the initial state fixed → a Mandelbrot-type set
 * (the locus where the orbit stays bounded). In StateSpace mode each
 * pixel is an initial condition with the parameters fixed → a Julia-type
 * set. For the complex-quadratic preset z->z^2+c these are exactly the
 * Mandelbrot and Julia sets; for any other 2D map they are the natural
 * generalizations, computed by the same engine.
 *
 * Coloring: escape time (how many iterations until |state| exceeds the
 * escape radius), with optional smooth (continuous) coloring. Points
 * that never escape (in the set) are drawn dark.
 * ============================================================ */
ImU32 fractal_palette(double t) {
  /* t in [0,1): a smooth cyclic palette (à la Ultra Fractal). */
  const double tau = 6.2831853;
  const int r = (int)(std::round(127.5 * (1.0 + std::sin(tau * (t + 0.00)))));
  const int g = (int)(std::round(127.5 * (1.0 + std::sin(tau * (t + 0.33)))));
  const int b = (int)(std::round(127.5 * (1.0 + std::sin(tau * (t + 0.67)))));
  return IM_COL32(r, g, b, 255);
}

/* Buddhabrot: instead of colouring each c by its escape time, we accumulate the
 * TRAJECTORIES of escaping points into a density histogram. Points that escape
 * leave a ghostly trace of where their orbit wandered; summed over many samples
 * this produces the famous Buddha-like figure. This is meaningful for the
 * complex-quadratic map z->z^2+c (state-space orbit in the (Re,Im) plane); for
 * other systems we fall back to the normal escape-time fractal.
 *
 * We render progressively: each call adds another batch of random samples to a
 * persistent accumulation buffer (app.buddha_accum) and re-maps it to colour,
 * so the image converges and sharpens the longer you watch. */
void compute_buddhabrot(AppState &app, int W, int H, std::vector<uint32_t> &out, bool reset) {
  if (W < 2) W = 2;
  if (H < 2) H = 2;
  const size_t npix = (size_t)W * H;
  out.assign(npix, 0xff000000u);
  const size_t n = app.state_names.size();
  if (n < 2 || app.mode != SystemMode::Map) return;

  if (reset || app.buddha_accum.size() != npix ||
      app.buddha_w != W || app.buddha_h != H) {
    app.buddha_accum.assign(npix, 0u);
    app.buddha_samples = 0;
    app.buddha_w = W; app.buddha_h = H;
    app.buddha_max = 1u;
  }

  const double x0 = app.fractal_xmin, x1 = app.fractal_xmax;
  const double y0 = app.fractal_ymin, y1 = app.fractal_ymax;
  const double R2 = app.fractal_escape_r * app.fractal_escape_r;
  const int maxit = std::max(20, app.fractal_max_iter);

  State cur = make_state_like(n, 0.0), nx = make_state_like(n, 0.0);
  std::vector<double> saved = app.param_values;
  char err[128] = {0};

  /* a deterministic-ish RNG so the image is reproducible per session */
  uint32_t &rng = app.buddha_rng;
  auto urand = [&]() { rng = rng * 1664525u + 1013904223u; return (double)(rng >> 8) * (1.0 / 16777216.0); };

  /* trajectory buffer to record visited (re,im) before knowing if it escapes */
  static std::vector<std::pair<float,float>> traj;
  traj.clear();

  /* sample budget per call: enough to make progress, capped to stay responsive */
  const int batch = 20000;
  /* sample c over a slightly padded box so orbits that wander in from outside
   * the view are represented too. */
  const double sx0 = -2.2, sx1 = 0.8, sy0 = -1.4, sy1 = 1.4;
  for (int sIdx = 0; sIdx < batch; ++sIdx) {
    const double cre = sx0 + (sx1 - sx0) * urand();
    const double cim = sy0 + (sy1 - sy0) * urand();
    /* iterate z=0 under z^2+c, recording the orbit; if it escapes, splat it */
    set_state_at(cur, 0, 0.0); set_state_at(cur, 1, 0.0);
    for (size_t i = 2; i < n; ++i) set_state_at(cur, i, 0.0);
    /* the complex-quadratic preset reads c from the FIRST two params (cx,cy) */
    if (app.param_values.size() >= 1) app.param_values[0] = cre;
    if (app.param_values.size() >= 2) app.param_values[1] = cim;
    traj.clear();
    bool escaped = false;
    for (int it = 0; it < maxit; ++it) {
      if (!step_map_state(app, cur, &nx, err, sizeof(err))) { escaped = false; break; }
      const double xx = state_at(nx, 0), yy = state_at(nx, 1);
      traj.emplace_back((float)xx, (float)yy);
      if (!std::isfinite(xx) || xx * xx + yy * yy > R2) { escaped = true; break; }
      std::swap(cur.v, nx.v); cur.t = nx.t;
    }
    if (escaped) {
      for (const auto &pt : traj) {
        if (pt.first < x0 || pt.first > x1 || pt.second < y0 || pt.second > y1) continue;
        const int px = (int)((pt.first - x0) / (x1 - x0) * (W - 1) + 0.5);
        const int py = (int)((pt.second - y0) / (y1 - y0) * (H - 1) + 0.5);
        if (px < 0 || px >= W || py < 0 || py >= H) continue;
        uint32_t &cell = app.buddha_accum[(size_t)py * W + px];
        cell++;
        if (cell > app.buddha_max) app.buddha_max = cell;
      }
    }
  }
  app.buddha_samples += batch;
  app.param_values = saved; sync_param_values(app);

  /* map accumulated density -> colour. Use a robust high-percentile as the
   * white point (a few hot pixels shouldn't wash everything out) and a steep
   * curve so the structure glows on a DARK background rather than a bright haze
   * swamping it. */
  uint32_t hi = 1;
  {
    /* find ~99.5th percentile of nonzero counts as the white point */
    std::vector<uint32_t> nz;
    nz.reserve(npix / 4);
    for (size_t i = 0; i < npix; ++i) if (app.buddha_accum[i]) nz.push_back(app.buddha_accum[i]);
    if (!nz.empty()) {
      size_t k = (size_t)(nz.size() * 0.995);
      if (k >= nz.size()) k = nz.size() - 1;
      std::nth_element(nz.begin(), nz.begin() + k, nz.end());
      hi = std::max<uint32_t>(1, nz[k]);
    }
  }
  const double inv_hi = 1.0 / (double)hi;
  for (size_t i = 0; i < npix; ++i) {
    const uint32_t d = app.buddha_accum[i];
    double t = (double)d * inv_hi;       /* 0..1 (clamped) against the white point */
    if (t > 1.0) t = 1.0;
    t = std::pow(t, 0.45);               /* lift faint filaments */
    const uint8_t v = (uint8_t)(t * 255.0);
    const uint8_t r = v;                 /* warm parchment */
    const uint8_t g = (uint8_t)(v * 0.82);
    const uint8_t b = (uint8_t)(v * 0.55);
    out[i] = 0xff000000u | r | ((uint32_t)g << 8) | ((uint32_t)b << 16);
  }
}

/* Thread-safe stepper for parallel grid sweeps (basins). Owns its OWN eval
 * scratch so multiple instances can run concurrently without touching the
 * shared app.eval_scratch / app.scratch_k* buffers. Evaluates the IR programs
 * directly (the AST fallback path is not thread-safe, so callers must check
 * app.use_ast_fallback is false before using this). Replicates exactly the two
 * integrators the basin sweep uses: a single map iteration, or one RK4 step. */
struct ThreadStepper {
  const AppState *app;
  dynsys::ir::Scratch scratch;
  std::vector<double> params;
  size_t dim = 0;
  bool is_map = false;
  double dt = 0.01;

  void init(const AppState &a) {
    app = &a;
    dim = a.state_names.size();
    is_map = (a.mode == SystemMode::Map);
    dt = a.dt;
    params = a.param_values;
    dynsys::ir::scratch_init(&scratch, a.definition_programs.size());
  }
  bool eval_prog(const dynsys::ir::Program &prog, const double *state, double t, double *out) {
    dynsys::ir::scratch_reset_eval(&scratch);
    dynsys::ir::RunContext rc;
    rc.state = state; rc.n_state = dim; rc.t = t;
    rc.params = params.data(); rc.n_params = params.size();
    rc.defs = app->definition_programs.data(); rc.n_defs = app->definition_programs.size();
    char e[8]; return dynsys::ir::run(prog, rc, scratch, out, e, sizeof(e));
  }
  /* override one parameter in this thread's PRIVATE snapshot (param-space
   * fractal: each thread sweeps its own rows with its own param values). */
  void set_param(int idx, double v) { if (idx >= 0 && (size_t)idx < params.size()) params[(size_t)idx] = v; }
  /* one map iteration: x_next[i] = prog_i(x) */
  bool map_step(const double *x, double *xn) {
    if (app->next_equation_programs.size() != dim) return false;
    for (size_t i = 0; i < dim; ++i)
      if (!eval_prog(app->next_equation_programs[i], x, 0.0, &xn[i])) return false;
    return true;
  }
  /* RHS f(x) into k */
  bool rhs(const double *x, double *k) {
    if (app->equation_programs.size() != dim) return false;
    for (size_t i = 0; i < dim; ++i)
      if (!eval_prog(app->equation_programs[i], x, 0.0, &k[i])) return false;
    return true;
  }
  /* one RK4 step of size dt */
  bool rk4_step(const double *x, double *xn) {
    std::vector<double> k1(dim), k2(dim), k3(dim), k4(dim), tmp(dim);
    if (!rhs(x, k1.data())) return false;
    for (size_t i = 0; i < dim; ++i) tmp[i] = x[i] + 0.5 * dt * k1[i];
    if (!rhs(tmp.data(), k2.data())) return false;
    for (size_t i = 0; i < dim; ++i) tmp[i] = x[i] + 0.5 * dt * k2[i];
    if (!rhs(tmp.data(), k3.data())) return false;
    for (size_t i = 0; i < dim; ++i) tmp[i] = x[i] + dt * k3[i];
    if (!rhs(tmp.data(), k4.data())) return false;
    for (size_t i = 0; i < dim; ++i)
      xn[i] = x[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    return true;
  }
};

void compute_fractal_image(AppState &app, int W, int H, std::vector<uint32_t> &out, int step = 1) {
  if (W < 2) W = 2;
  if (H < 2) H = 2;
  if (step < 1) step = 1;
  out.assign((size_t)W * H, 0xff000000u);
  const size_t n = app.state_names.size();
  if (n < 2) return; /* needs a 2D plane */
  if (app.mode != SystemMode::Map) return; /* escape-time is for maps */

  /* Honor the performance governor: when throttling, cap the per-pixel
   * iteration budget (down to ~25% at full throttle, floor 10) so a deep zoom
   * on a heavy fractal can't lock up the frame. The progressive refine still
   * sharpens later once load drops and the throttle eases. */
  int maxit = std::max(10, app.fractal_max_iter);
  if (app.perf_throttle > 0.02) {
    maxit = (int)(maxit * (1.0 - 0.75 * app.perf_throttle));
    if (maxit < 10) maxit = 10;
  }
  const double R2 = app.fractal_escape_r * app.fractal_escape_r;
  const double x0 = app.fractal_xmin, x1 = app.fractal_xmax;
  const double y0 = app.fractal_ymin, y1 = app.fractal_ymax;

  bool param_mode = (app.fractal_mode == AppState::FractalMode::ParameterSpace);
  /* Parameter-space mode needs at least one parameter to vary. A map with no
   * params (e.g. the Newton fractal) must use state-space: the plane is the
   * initial (x,y). Otherwise we'd vary nothing and get a meaningless image. */
  if (param_mode && app.params.empty()) param_mode = false;
  const int cxi = std::max(0, std::min(app.fractal_param_cx_index, (int)app.params.size() - 1));
  const int cyi = std::max(0, std::min(app.fractal_param_cy_index, (int)app.params.size() - 1));

  /* Save params we will temporarily override (param mode, serial path only). */
  std::vector<double> saved = app.param_values;

  /* Per-pixel escape/convergence work for one row py, using a stepper `st` that
   * owns private eval scratch + a private parameter snapshot (so multiple
   * threads can run disjoint rows safely). Mirrors the serial logic exactly. */
  auto do_row = [&](ThreadStepper &st, int py) {
    const double b_im = y0 + (y1 - y0) * (double)py / (H - 1);
    std::vector<double> sx(n), cur(n), nx(n);
    for (int px = 0; px < W; px += step) {
      const double a_re = x0 + (x1 - x0) * (double)px / (W - 1);
      for (size_t i = 0; i < n; ++i) sx[i] = state_at(app.start, i);
      if (param_mode) {
        st.set_param(cxi, a_re);
        st.set_param(cyi, b_im);
      } else {
        sx[0] = a_re; if (n > 1) sx[1] = b_im;
      }
      int it = 0; double r2 = 0.0;
      for (size_t i = 0; i < n; ++i) cur[i] = sx[i];
      bool escaped = false, converged = false; double conv_x = 0.0, conv_y = 0.0;
      const bool allow_converge = app.params.empty();
      for (; it < maxit; ++it) {
        if (!st.map_step(cur.data(), nx.data())) { it = maxit; break; }
        const double xx = nx[0], yy = (n > 1) ? nx[1] : 0.0;
        r2 = xx * xx + yy * yy;
        if (!std::isfinite(r2) || r2 > R2) { ++it; escaped = true; break; }
        if (allow_converge) {
          const double dx = xx - cur[0], dy = (n > 1) ? yy - cur[1] : 0.0;
          if (dx * dx + dy * dy < 1e-20) { converged = true; conv_x = xx; conv_y = yy; ++it; break; }
        }
        std::swap(cur, nx);
      }
      uint32_t color = 0xff101014u;
      if (escaped && it < maxit && r2 > R2) {
        double mu = it;
        if (app.fractal_smooth && r2 > 1.0) mu = it + 1.0 - std::log(std::log(std::sqrt(r2))) / std::log(2.0);
        const double t = std::fmod(mu * 0.04, 1.0);
        color = fractal_palette(t < 0 ? t + 1.0 : t);
      } else if (converged) {
        double ang = std::atan2(conv_y, conv_x) / (2.0 * 3.14159265358979) + 0.5;
        const double shade = 1.0 - 0.6 * std::min(1.0, (double)it / 40.0);
        uint32_t base = fractal_palette(std::fmod(ang, 1.0));
        const uint8_t r = (uint8_t)(((base >> 0) & 0xff) * shade);
        const uint8_t g = (uint8_t)(((base >> 8) & 0xff) * shade);
        const uint8_t b = (uint8_t)(((base >> 16) & 0xff) * shade);
        color = 0xff000000u | r | ((uint32_t)g << 8) | ((uint32_t)b << 16);
      }
      for (int by = py; by < py + step && by < H; ++by)
        for (int bx = px; bx < px + step && bx < W; ++bx)
          out[(size_t)by * W + bx] = color;
    }
  };

  /* Use the PARALLEL path when the thread-safe IR evaluator is active (each
   * thread gets a private ThreadStepper); else the serial fallback. Identical
   * pixels either way. */
  std::vector<int> rows;
  for (int py = 0; py < H; py += step) rows.push_back(py);
  const bool can_parallel = !app.use_ast_fallback &&
                            std::thread::hardware_concurrency() > 1 && (int)rows.size() >= 8;
  if (can_parallel) {
    unsigned hw = std::thread::hardware_concurrency();
    unsigned nth = std::min<unsigned>(hw, (unsigned)rows.size());
    std::atomic<size_t> next{0};
    auto worker = [&]() {
      ThreadStepper st; st.init(app);
      for (;;) {
        size_t k = next.fetch_add(1, std::memory_order_relaxed);
        if (k >= rows.size()) break;
        do_row(st, rows[k]);
      }
    };
    std::vector<std::thread> pool; pool.reserve(nth);
    for (unsigned t = 0; t < nth; ++t) pool.emplace_back(worker);
    for (auto &th : pool) th.join();
  } else {
    ThreadStepper st; st.init(app);
    for (int py : rows) do_row(st, py);
  }
  app.param_values = saved; /* restore (untouched on the parallel path) */
  sync_param_values(app);
}

void render_fractal_background(AppState &app) {
  const float w = (float)app.window_width, h = (float)app.window_height;
  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  draw->AddRectFilled(ImVec2(0, 0), ImVec2(w, h), IM_COL32(10, 10, 14, 255));

  if (app.state_names.size() < 2 || app.mode != SystemMode::Map) {
    draw->AddText(ImVec2(20, 48), IM_COL32(220, 220, 230, 230),
                  "Fractal view needs a 2D MAP. Load the 'Complex quadratic' preset");
    draw->AddText(ImVec2(20, 66), IM_COL32(220, 220, 230, 230),
                  "(Mandelbrot/Julia), or any 2D map in map mode.");
    return;
  }

  /* init view once from the configured window */
  if (!app.fractal_view_init) {
    app.fractal_view_init = true;
    if (app.params.empty()) {
      /* No-parameter maps (Newton) are state-space fractals; default to a
       * symmetric window around the origin so the basins are framed. */
      app.fractal_mode = AppState::FractalMode::StateSpace;
      app.fractal_xmin = -2.0; app.fractal_xmax = 2.0;
      app.fractal_ymin = -2.0; app.fractal_ymax = 2.0;
    } else if (app.home_view_set) {
      /* Use the preset's view2d window so fractals like Burning Ship / Multibrot
       * open framed on their actual structure rather than the Mandelbrot box. */
      app.fractal_xmin = app.home_x_min; app.fractal_xmax = app.home_x_max;
      app.fractal_ymin = app.home_y_min; app.fractal_ymax = app.home_y_max;
    }
    app.fractal_dirty = true;
  }

  /* interaction: pan/zoom the complex window; mark dirty to recompute */
  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse) {
    auto invx = [&](float pxf){ return app.fractal_xmin + (double)pxf / w * (app.fractal_xmax - app.fractal_xmin); };
    auto invy = [&](float pyf){ return app.fractal_ymin + (double)pyf / h * (app.fractal_ymax - app.fractal_ymin); };
    if (io.MouseWheel != 0.0f) {
      const double mx = invx(io.MousePos.x), my = invy(io.MousePos.y);
      const double f = std::pow(0.8, clamped_wheel((double)io.MouseWheel));
      app.fractal_xmin = mx + (app.fractal_xmin - mx) * f;
      app.fractal_xmax = mx + (app.fractal_xmax - mx) * f;
      app.fractal_ymin = my + (app.fractal_ymin - my) * f;
      app.fractal_ymax = my + (app.fractal_ymax - my) * f;
      app.fractal_dirty = true;
    }
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) && !ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      const double ddx = (double)io.MouseDelta.x / w * (app.fractal_xmax - app.fractal_xmin);
      const double ddy = (double)io.MouseDelta.y / h * (app.fractal_ymax - app.fractal_ymin);
      app.fractal_xmin -= ddx; app.fractal_xmax -= ddx;
      app.fractal_ymin -= ddy; app.fractal_ymax -= ddy;
      app.fractal_dirty = true;
    }
    /* In Julia mode, right-drag picks the constant c from the plane and
     * recomputes — a lovely way to explore the Mandelbrot/Julia link. */
  }

  /* (re)compute the image at a capped resolution. Two refinements for a
   * smooth experience:
   *  - Debounce: a view change waits a few still frames before recompute,
   *    so dragging/zooming stays responsive.
   *  - Progressive: after the wait, render a coarse pass (step 8) for
   *    instant feedback, then refine 8->4->2->1 on successive frames. This
   *    removes the long stall on the first paint / after a big zoom. */
  const int CW = std::min(app.window_width, 900);
  const int CH = std::min(app.window_height, 700);
  const bool force = (app.fractal_tex == 0 || app.fractal_tex_w != CW || app.fractal_tex_h != CH);

  if (app.fractal_mode == AppState::FractalMode::Buddhabrot) {
    /* Buddhabrot accumulates over frames: reset when the view changed, then add
     * a batch of samples each frame so the figure sharpens as you watch. */
    const bool reset = force || app.fractal_dirty;
    app.fractal_dirty = false;
    std::vector<uint32_t> img;
    compute_buddhabrot(app, CW, CH, img, reset);
    if (app.fractal_tex == 0) glGenTextures(1, &app.fractal_tex);
    glBindTexture(GL_TEXTURE_2D, app.fractal_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, CW, CH, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    app.fractal_tex_w = CW; app.fractal_tex_h = CH;
    if (app.fractal_tex != 0)
      draw->AddImage((ImTextureID)(uintptr_t)app.fractal_tex, ImVec2(0, 0), ImVec2(w, h));
    char bhud[256];
    std::snprintf(bhud, sizeof(bhud),
                  "Buddhabrot (trajectory density of escaping orbits)  |  %ld samples  |  re [%.3g, %.3g]  im [%.3g, %.3g]",
                  app.buddha_samples, app.fractal_xmin, app.fractal_xmax, app.fractal_ymin, app.fractal_ymax);
    draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), bhud);
    draw->AddText(ImVec2(14, app.window_toolbar_h + 26.0f), IM_COL32(170, 170, 180, 220),
                  "accumulates as you watch (longer = sharper).  drag: pan   wheel: zoom");
    return;
  }

  if (app.fractal_dirty) { app.fractal_settle = 5; app.fractal_dirty = false; }
  bool start_progress = force;
  if (app.fractal_settle > 0) { if (--app.fractal_settle == 0) start_progress = true; }
  if (start_progress) app.fractal_prog_level = 8; /* begin coarse */

  if (app.fractal_prog_level > 0) {
    const int step = app.fractal_prog_level;
    std::vector<uint32_t> img;
    compute_fractal_image(app, CW, CH, img, step);
    if (app.fractal_tex == 0) glGenTextures(1, &app.fractal_tex);
    glBindTexture(GL_TEXTURE_2D, app.fractal_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, CW, CH, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    app.fractal_tex_w = CW; app.fractal_tex_h = CH;
    /* refine on the next frame: 8 -> 4 -> 2 -> 1 -> done */
    app.fractal_prog_level = (step > 1) ? step / 2 : 0;
  }

  if (app.fractal_tex != 0)
    draw->AddImage((ImTextureID)(uintptr_t)app.fractal_tex, ImVec2(0, 0), ImVec2(w, h));

  const char *mode = app.fractal_mode == AppState::FractalMode::ParameterSpace
                         ? "parameter space (Mandelbrot-type: orbit of the start point)"
                         : "state space (Julia-type: sweeping the initial condition)";
  char hud[256];
  std::snprintf(hud, sizeof(hud), "Fractal — %s   |  re [%.4g, %.4g]  im [%.4g, %.4g]  iters %d",
                mode, app.fractal_xmin, app.fractal_xmax, app.fractal_ymin, app.fractal_ymax, app.fractal_max_iter);
  draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), hud);
  if (!io.WantCaptureMouse)
    draw->AddText(ImVec2(14, h - 24), IM_COL32(150, 150, 160, 200),
                  "drag: pan   wheel: zoom   (set mode/params/iterations in the Setup tab)");
}

/* ============================================================
 * PHASE C: IFS / chaos game.
 * A built-in gallery of iterated function systems; the chaos game renders
 * each attractor into a density texture. The fractal dimension of the
 * result is measured with box_counting_dimension (cross-checking it).
 * ============================================================ */
struct IFSPreset {
  const char *name;
  std::vector<dynsys::analysis::AffineMap> maps;
};

const std::vector<IFSPreset> &ifs_gallery() {
  static const std::vector<IFSPreset> g = {
    { "Barnsley fern", {
        {0.0,  0.0,  0.0,  0.16, 0.0, 0.0,  0.01},
        {0.85, 0.04, -0.04,0.85, 0.0, 1.6,  0.85},
        {0.2, -0.26, 0.23, 0.22, 0.0, 1.6,  0.07},
        {-0.15,0.28, 0.26, 0.24, 0.0, 0.44, 0.07} } },
    { "Sierpinski triangle", {
        {0.5,0,0,0.5, 0.0,  0.0, 1.0/3},
        {0.5,0,0,0.5, 0.5,  0.0, 1.0/3},
        {0.5,0,0,0.5, 0.25, 0.5, 1.0/3} } },
    { "Dragon curve", {
        {0.824074, 0.281482, -0.212346, 0.864198, -1.882290, -0.110607, 0.787473},
        {0.088272, 0.520988, -0.463889, -0.377778, 0.785360, 8.095795, 0.212527} } },
    { "Fractal tree", {
        {0.0,  0.0,  0.0,  0.5,  0.0, 0.0,  0.05},
        {0.42,-0.42, 0.42, 0.42, 0.0, 0.2,  0.40},
        {0.42, 0.42,-0.42, 0.42, 0.0, 0.2,  0.40},
        {0.1,  0.0,  0.0,  0.1,  0.0, 0.2,  0.15} } },
    { "Maple leaf", {
        {0.14, 0.01, 0.0,  0.51, -0.08, -1.31, 0.10},
        {0.43, 0.52, -0.45,0.5,  1.49,  -0.75, 0.35},
        {0.45,-0.49, 0.47, 0.47, -1.62, -0.74, 0.35},
        {0.49, 0.0,  0.0,  0.51, 0.02,  1.62,  0.20} } },
  };
  return g;
}

/* Evaluate the IFS coefficient expressions against the current parameter
 * values into numeric affine maps. On the initial build (ifs_maps empty)
 * every coefficient is evaluated. On a refresh (e.g. a parameter changed),
 * pass preserve_literals=true to update ONLY the parameter-driven
 * coefficients and keep any literal coefficients the user edited via their
 * sliders. Returns false on an eval error. */
bool evaluate_ifs_maps(AppState &app, bool preserve_literals = false) {
  if (app.ifs_map_exprs.empty()) { app.ifs_maps.clear(); return false; }
  State dummy = app.start; resize_state(dummy, app.state_names.size());
  char err[128] = {0};
  auto ev = [&](node_t *n, double fallback) -> double {
    if (!n) return fallback;
    double v = fallback;
    if (!eval_expr_at(app, n, dummy, &v, err, sizeof(err))) return fallback;
    return std::isfinite(v) ? v : fallback;
  };
  const bool have_old = preserve_literals && app.ifs_maps.size() == app.ifs_map_exprs.size()
                        && app.ifs_coef_literal.size() == app.ifs_map_exprs.size();
  std::vector<dynsys::analysis::AffineMap> out;
  for (size_t i = 0; i < app.ifs_map_exprs.size(); ++i) {
    const auto &me = app.ifs_map_exprs[i];
    dynsys::analysis::AffineMap m;
    double *dst[7] = {&m.a,&m.b,&m.c,&m.d,&m.e,&m.f,&m.p};
    node_t *src[7] = {me.a,me.b,me.c,me.d,me.e,me.f,me.p};
    double oldv[7];
    if (have_old) {
      const auto &o = app.ifs_maps[i];
      oldv[0]=o.a; oldv[1]=o.b; oldv[2]=o.c; oldv[3]=o.d; oldv[4]=o.e; oldv[5]=o.f; oldv[6]=o.p;
    }
    for (int k = 0; k < 7; ++k) {
      const bool lit = (i < app.ifs_coef_literal.size()) ? app.ifs_coef_literal[i][(size_t)k] : true;
      if (have_old && lit) *dst[k] = oldv[k];  /* keep the user's edited literal */
      else *dst[k] = ev(src[k], 0);            /* (re)evaluate param-driven coef */
    }
    out.push_back(m);
  }
  app.ifs_maps = std::move(out);
  return true;
}

void render_ifs_background(AppState &app) {
  const float w = (float)app.window_width, h = (float)app.window_height;
  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  draw->AddRectFilled(ImVec2(0, 0), ImVec2(w, h), IM_COL32(10, 11, 15, 255));

  const auto &gallery = ifs_gallery();
  /* For a loaded IFS model, (re)evaluate the coefficient expressions against
   * the current parameter values whenever the maps are dirty — this is what
   * lets an IFS respond to its parameters. */
  const bool is_ifs_model = (app.mode == SystemMode::IFS && !app.ifs_map_exprs.empty());
  /* Initial build evaluates everything. On a dirty refresh, preserve any
   * literal coefficients the user edited via sliders and only re-evaluate
   * the parameter-driven ones. (A fully-literal IFS never needs refresh.) */
  if (is_ifs_model) {
    if (app.ifs_maps.empty()) evaluate_ifs_maps(app, false);
    else if (app.ifs_dirty && !app.ifs_maps_editable) evaluate_ifs_maps(app, true);
  }

  const bool use_model = (app.mode == SystemMode::IFS && !app.ifs_maps.empty());
  std::vector<dynsys::analysis::AffineMap> maps_to_use;
  std::string ifs_name;
  if (use_model) {
    maps_to_use = app.ifs_maps;
    ifs_name = app.params.empty() ? "loaded IFS model"
                                  : "loaded IFS model (parametrized)";
  } else {
    if (app.ifs_selected < 0) app.ifs_selected = 0;
    if (app.ifs_selected >= (int)gallery.size()) app.ifs_selected = (int)gallery.size() - 1;
    maps_to_use = gallery[(size_t)app.ifs_selected].maps;
    ifs_name = gallery[(size_t)app.ifs_selected].name;
  }

  const int CW = std::max(64, (int)w);
  const int CH = std::max(64, (int)h);

  if (app.ifs_dirty || app.ifs_tex == 0 || app.ifs_tex_w != CW || app.ifs_tex_h != CH) {
    /* run the chaos game */
    dynsys::analysis::IFSResult R =
        dynsys::analysis::chaos_game(maps_to_use, app.ifs_iterations, 12345u);
    if (R.ok && !R.xs.empty()) {
      /* measure the dimension (cross-check) */
      std::vector<double> xs(R.xs.begin(), R.xs.end()), ys(R.ys.begin(), R.ys.end());
      auto D = dynsys::analysis::box_counting_dimension(xs, ys, 12);
      app.ifs_box_dim = D.dimension; app.ifs_box_dim_ready = D.ok;

      /* fit the attractor into the window with a margin, preserving aspect */
      const double spanx = std::max(1e-9, R.xmax - R.xmin);
      const double spany = std::max(1e-9, R.ymax - R.ymin);
      const double margin = 0.06;
      const double sx = (CW * (1 - 2 * margin)) / spanx;
      const double sy = (CH * (1 - 2 * margin)) / spany;
      const double s = std::min(sx, sy);
      const double cx = 0.5 * (R.xmin + R.xmax), cy = 0.5 * (R.ymin + R.ymax);

      /* density accumulation */
      std::vector<float> dens((size_t)CW * CH, 0.0f);
      float dmax = 0.0f;
      for (size_t i = 0; i < R.xs.size(); ++i) {
        const int px = (int)(CW * 0.5 + (R.xs[i] - cx) * s);
        const int py = (int)(CH * 0.5 - (R.ys[i] - cy) * s); /* flip y */
        if (px < 0 || px >= CW || py < 0 || py >= CH) continue;
        float &c = dens[(size_t)py * CW + px];
        c += 1.0f; if (c > dmax) dmax = c;
      }
      /* colorize: green-ish for botanical sets, gamma-lifted density */
      std::vector<uint32_t> img((size_t)CW * CH, 0xff0a0b0fu);
      const float inv = dmax > 0 ? 1.0f / std::log(1.0f + dmax) : 0.0f;
      for (size_t i = 0; i < img.size(); ++i) {
        const float c = dens[i];
        if (c <= 0) continue;
        float t = std::log(1.0f + c) * inv;
        t = std::pow(t, 0.6f);
        const int r = (int)(40 + 120 * t);
        const int g = (int)(120 + 135 * t);
        const int b = (int)(70 + 80 * t);
        img[i] = 0xff000000u | ((uint32_t)b << 16) | ((uint32_t)g << 8) | (uint32_t)r;
      }
      if (app.ifs_tex == 0) glGenTextures(1, &app.ifs_tex);
      glBindTexture(GL_TEXTURE_2D, app.ifs_tex);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, CW, CH, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      app.ifs_tex_w = CW; app.ifs_tex_h = CH;
    }
    app.ifs_dirty = false;
  }

  if (app.ifs_tex != 0)
    draw->AddImage((ImTextureID)(uintptr_t)app.ifs_tex, ImVec2(0, 0), ImVec2(w, h));

  char hud[256];
  if (app.ifs_box_dim_ready)
    std::snprintf(hud, sizeof(hud), "IFS chaos game — %s   |  %ld points  |  measured box-dimension D = %.3f",
                  ifs_name.c_str(), app.ifs_iterations, app.ifs_box_dim);
  else
    std::snprintf(hud, sizeof(hud), "IFS chaos game — %s   |  %ld points", ifs_name.c_str(), app.ifs_iterations);
  draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), hud);
  draw->AddText(ImVec2(14, h - 24), IM_COL32(150, 150, 160, 200),
                use_model ? "IFS attractor (mode = ifs) — rendered via the chaos game"
                          : "IFS attractor — load a 'mode = ifs' model (e.g. the Barnsley fern preset)");
}

/* ============================================================
 * PHASE B/C: basins of attraction.
 * Integrate from a grid of initial conditions over the current 2D view;
 * color each pixel by which attractor its trajectory settles onto
 * (attractors discovered automatically by clustering endpoints). The
 * basin boundaries are the interesting object — often fractal (Newton's
 * method, forced oscillators). Works for both ODEs (one integrator step
 * via step_state) and maps (one iteration).
 * ============================================================ */
ImU32 basin_color(int label, int nlabels, float speed, bool shade) {
  if (label == -1) return IM_COL32(8, 8, 12, 255);     /* diverged to infinity (near-black) */
  if (label == -2) {
    /* Bounded but did not settle to a fixed point — i.e. a chaotic (or
     * long-period) attractor, as in the Henon map. Colour it a clear teal so
     * the bounded basin stands out sharply against the near-black escape set,
     * shaded by how quickly it was confirmed bounded. Previously this was a
     * flat dark grey almost indistinguishable from "escaped", which made such
     * basins look uniformly blank. */
    float br = shade ? (0.45f + 0.55f * std::max(0.0f, std::min(1.0f, speed))) : 1.0f;
    return IM_COL32((int)(40 * br), (int)(150 * br), (int)(165 * br), 255);
  }
  /* distinct hues around the wheel */
  const double h = nlabels > 0 ? (double)label / std::max(1, nlabels) : 0.0;
  const double H = h * 6.0;
  const int i = (int)H; const double f = H - i;
  double r, g, b;
  const double q = 1 - f, t = f;
  switch (i % 6) {
    case 0: r=1; g=t; b=0; break;
    case 1: r=q; g=1; b=0; break;
    case 2: r=0; g=1; b=t; break;
    case 3: r=0; g=q; b=1; break;
    case 4: r=t; g=0; b=1; break;
    default: r=1; g=0; b=q; break;
  }
  float br = shade ? (0.35f + 0.65f * std::max(0.0f, std::min(1.0f, speed))) : 1.0f;
  return IM_COL32((int)(r * 255 * br), (int)(g * 255 * br), (int)(b * 255 * br), 255);
}


void compute_basin_image(AppState &app, int W, int H, std::vector<uint32_t> &out, int step = 1) {
  if (W < 2) W = 2;
  if (H < 2) H = 2;
  if (step < 1) step = 1;
  out.assign((size_t)W * H, 0xff000000u);
  /* Progressive refinement: when step>1 we classify a coarse (W/step x H/step)
   * grid and block-upscale it into the full image, so the first frames are
   * cheap and the picture sharpens 8 -> 4 -> 2 -> 1 over subsequent frames
   * instead of one multi-hundred-million-step stall. The compute grid is the
   * reduced size; cw/ch below are that grid's dimensions. */
  const int cw = std::max(2, W / step);
  const int ch = std::max(2, H / step);
  const size_t n = app.state_names.size();
  if (n < 2) return;

  app.phase_x_index = std::max(0, std::min(app.phase_x_index, (int)n - 1));
  app.phase_y_index = std::max(0, std::min(app.phase_y_index, (int)n - 1));
  const size_t ix = (size_t)app.phase_x_index, iy = (size_t)app.phase_y_index;
  PlotBounds b = sanitize_bounds(current_phase_bounds(app, ix, iy));

  /* advance(x,y): set the two plane coords on a copy of the start state,
   * step once, read them back. Other state dims are held at their start
   * values (the basin is a 2D slice of the full state space).
   *
   * IMPORTANT: force a cheap fixed-step method during the grid sweep. The
   * adaptive integrators (RKF45/DOPRI45) can take many internal substeps
   * per call, which — multiplied by thousands of steps over hundreds of
   * thousands of cells — would hang the UI. RK4 with the current dt is
   * accurate enough for basin classification. Saved/restored around the
   * sweep. */
  const Integrator saved_integrator = app.integrator;
  if (app.mode == SystemMode::ODE &&
      (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45))
    app.integrator = Integrator::RK4;

  char err[128] = {0};
  bool eval_err = false;
  auto advance = [&](double x, double y, double *nx, double *ny) -> bool {
    State s = app.start;
    resize_state(s, n);
    set_state_at(s, ix, x);
    set_state_at(s, iy, y);
    State out_s{};
    if (!step_state(app, s, &out_s, err, sizeof(err))) { eval_err = true; return false; }
    *nx = state_at(out_s, ix);
    *ny = state_at(out_s, iy);
    return true;
  };

  dynsys::analysis::BasinOptions opt;
  opt.xmin = b.xmin; opt.xmax = b.xmax; opt.ymin = b.ymin; opt.ymax = b.ymax;
  opt.width = cw; opt.height = ch;            /* compute on the reduced grid */
  /* Per-cell iteration budget. Maps settle or escape FAST (a Henon point is
   * clearly bounded-or-divergent within a couple hundred iterations), so the
   * 1500-step default — fine for slowly-settling ODE flows — would do ~half a
   * billion map evaluations on the full grid and freeze the view. Cap maps
   * much lower. */
  opt.max_steps = (app.mode == SystemMode::Map)
                      ? std::min(std::max(50, app.basin_steps), 250)
                      : std::max(50, app.basin_steps);
  opt.cluster_tol = app.basin_cluster_tol;
  opt.settle_tol = (app.mode == SystemMode::Map) ? 1e-6 : 1e-5;

  /* Use the PARALLEL basin solver when the (thread-safe) IR eval path is in
   * use: each worker thread gets its own ThreadStepper (private eval scratch),
   * giving a big speedup on multi-core machines with identical results. Fall
   * back to the serial path when the AST-fallback evaluator is active (it is
   * not thread-safe). The grid coordinates baked into the stepper's advance
   * match compute_basins' own (x0,y0) mapping exactly. */
  dynsys::analysis::BasinResult R;
  const bool can_parallel = !app.use_ast_fallback &&
                            std::thread::hardware_concurrency() > 1 && ch >= 8;
  if (can_parallel) {
    const bool is_map = (app.mode == SystemMode::Map);
    auto make_advance = [&app, n, ix, iy, is_map](int /*tid*/) {
      auto stepper = std::make_shared<ThreadStepper>();
      stepper->init(app);
      std::vector<double> s(n), sn(n);
      for (size_t i = 0; i < n; ++i) s[i] = state_at(app.start, i);
      return [stepper, n, ix, iy, is_map, s, sn](double x, double y, double *nx, double *ny) mutable -> bool {
        for (size_t i = 0; i < n; ++i) s[i] = state_at(stepper->app->start, i);
        s[ix] = x; s[iy] = y;
        bool ok = is_map ? stepper->map_step(s.data(), sn.data())
                         : stepper->rk4_step(s.data(), sn.data());
        if (!ok) return false;
        *nx = sn[ix]; *ny = sn[iy];
        return true;
      };
    };
    R = dynsys::analysis::compute_basins_mt(make_advance, opt);
  } else {
    R = dynsys::analysis::compute_basins(advance, opt);
  }
  app.integrator = saved_integrator; /* restore the user's choice */
  app.basin_attractor_count = (int)R.attractors.size();
  app.basin_n_converged = R.n_converged;
  app.basin_n_diverged = R.n_diverged;
  app.basin_n_nonconvergent = R.n_nonconvergent;
  if (!R.ok) return;
  const int nl = (int)R.attractors.size();
  /* block-upscale the coarse cw x ch classification into the full W x H image
   * (nearest-neighbor). When step==1, cw==W and ch==H so this is 1:1. */
  for (int j = 0; j < H; ++j) {
    int cj = (j * ch) / H; if (cj >= ch) cj = ch - 1;
    for (int i = 0; i < W; ++i) {
      int ci = (i * cw) / W; if (ci >= cw) ci = cw - 1;
      const size_t cidx = (size_t)cj * cw + ci;
      out[(size_t)j * W + i] =
          basin_color(R.cell_attractor[cidx], nl, R.cell_speed[cidx], app.basin_shade_speed);
    }
  }
}

void render_basin_background(AppState &app) {
  const float w = (float)app.window_width, h = (float)app.window_height;
  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  draw->AddRectFilled(ImVec2(0, 0), ImVec2(w, h), IM_COL32(10, 10, 14, 255));
  if (app.state_names.size() < 2) {
    draw->AddText(ImVec2(20, 48), IM_COL32(220, 220, 230, 230),
                  "Basins need a system with >= 2 state variables.");
    return;
  }

  /* pan/zoom reuse the phase bounds so basins and the phase portrait share
   * a coordinate frame; any view change marks the image dirty.
   * CRASH FIX: clamp the plane indices to the CURRENT state dimension before
   * using them. When a system with fewer variables is loaded (e.g. switching
   * from a 3-var to a 2-var system) the stale phase_y_index could index past
   * state_names, reading out of bounds in the HUD below. Clamp (and persist)
   * here exactly as compute_basin_image does. */
  const size_t nstate = app.state_names.size();
  app.phase_x_index = std::max(0, std::min(app.phase_x_index, (int)nstate - 1));
  app.phase_y_index = std::max(0, std::min(app.phase_y_index, (int)nstate - 1));
  const size_t ix = (size_t)app.phase_x_index, iy = (size_t)app.phase_y_index;
  PlotBounds b = sanitize_bounds(current_phase_bounds(app, ix, iy));
  /* Safety clamp for basins: if we're auto-fitting (e.g. after a stray
   * double-click) the fit can be dragged enormously wide by diverging seed
   * orbits in a map basin. Cap the window to a few times the home window so the
   * basin can never "zoom way out" to a useless scale. */
  if (app.phase_auto_bounds && app.home_view_set) {
    const double homw = (double)app.home_x_max - app.home_x_min;
    const double homh = (double)app.home_y_max - app.home_y_min;
    const double cx = 0.5 * (b.xmin + b.xmax), cy = 0.5 * (b.ymin + b.ymax);
    double bw = b.xmax - b.xmin, bh = b.ymax - b.ymin;
    const double maxw = homw * 6.0, maxh = homh * 6.0;
    bool clamped = false;
    if (bw > maxw) { b.xmin = cx - maxw * 0.5; b.xmax = cx + maxw * 0.5; clamped = true; }
    if (bh > maxh) { b.ymin = cy - maxh * 0.5; b.ymax = cy + maxh * 0.5; clamped = true; }
    if (clamped) { set_manual_phase_bounds(app, b); b = sanitize_bounds(b); }
  }
  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse) {
    const ImVec2 p0(0, 0);
    const Point2 m = screen_to_plot(b, p0, w, h, io.MousePos);
    if (io.MouseWheel != 0.0f) { zoom_phase_bounds(app, b, m, std::pow(0.8, clamped_wheel((double)io.MouseWheel))); app.basin_dirty = true; }
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) && !ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      /* clamp the per-frame delta: a spurious large jump (e.g. focus regained,
       * or the very first drag frame reporting accumulated motion) must not be
       * able to fling the view across the plane. Cap at a quarter window. */
      float mdx = io.MouseDelta.x, mdy = io.MouseDelta.y;
      mdx = std::max(-w * 0.25f, std::min(mdx, w * 0.25f));
      mdy = std::max(-h * 0.25f, std::min(mdy, h * 0.25f));
      const double ddx = (double)mdx / w * (b.xmax - b.xmin);
      const double ddy = (double)mdy / h * (b.ymax - b.ymin);
      pan_phase_bounds(app, b, -ddx, ddy);
      app.basin_dirty = true;
    }
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
      /* Reset the basin view to the preset's home window (the framing it loaded
       * with), NOT an auto-fit. Auto-fitting a basin is dangerous: a map basin
       * contains diverging seed orbits, so fitting to the data would zoom the
       * view enormously far out -- the "double-click / stray click zooms me way
       * out" bug. The home window is always sane. */
      if (app.home_view_set) {
        PlotBounds hb{app.home_x_min, app.home_x_max, app.home_y_min, app.home_y_max};
        set_manual_phase_bounds(app, hb);
      } else {
        app.phase_auto_bounds = true; /* only if we have no home to go to */
      }
      app.basin_dirty = true;
    }
  }

  /* ODE basins integrate a trajectory PER CELL, far costlier than a map's
   * single step, so cap their grid lower; maps can afford the full grid. */
  const bool basin_is_ode = (app.mode == SystemMode::ODE);
  const int CW = std::min(app.window_width, basin_is_ode ? 380 : 480);
  const int CH = std::min(app.window_height, basin_is_ode ? 300 : 384);
  const bool bforce = (app.basin_tex == 0 || app.basin_tex_w != CW || app.basin_tex_h != CH);
  if (app.basin_dirty) { app.basin_settle = 6; app.basin_dirty = false; }
  bool start_progress = bforce;
  if (app.basin_settle > 0) { if (--app.basin_settle == 0) start_progress = true; }
  /* begin a coarse->fine pass: 16 -> 8 -> 4 -> 2 -> 1 (the coarse first frame
   * appears almost instantly). Under the performance governor's throttle, stop
   * refining early (a coarser final level) so a heavy basin can't lock the UI. */
  if (start_progress) app.basin_prog_level = 16;
  if (app.basin_prog_level > 0) {
    int step = app.basin_prog_level;
    std::vector<uint32_t> img;
    compute_basin_image(app, CW, CH, img, step);
    if (app.basin_tex == 0) glGenTextures(1, &app.basin_tex);
    glBindTexture(GL_TEXTURE_2D, app.basin_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, CW, CH, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    app.basin_tex_w = CW; app.basin_tex_h = CH;
    /* next finer level; if heavily throttled, stop at a coarser floor */
    int floor_level = 1;
    if (app.perf_throttle > 0.6) floor_level = 4;
    else if (app.perf_throttle > 0.3) floor_level = 2;
    app.basin_prog_level = (step > floor_level) ? step / 2 : 0;
  }
  if (app.basin_tex != 0)
    draw->AddImage((ImTextureID)(uintptr_t)app.basin_tex, ImVec2(0, 0), ImVec2(w, h));

  char hud[360];
  std::snprintf(hud, sizeof(hud),
                "Basins of attraction — %d basin(s)  |  %s vs %s  |  converged %ld, escaped %ld, non-convergent %ld%s",
                app.basin_attractor_count, app.state_names[ix].c_str(), app.state_names[iy].c_str(),
                app.basin_n_converged, app.basin_n_diverged, app.basin_n_nonconvergent,
                app.basin_prog_level > 1 ? "   [refining…]" : "");
  draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), hud);

  /* If almost nothing converged, the system likely has a single chaotic
   * attractor (e.g. Lorenz) — basins aren't meaningful there. Say so and
   * point to a multistable example, instead of showing a confusing wash. */
  const long total_cells = app.basin_n_converged + app.basin_n_diverged + app.basin_n_nonconvergent;
  if (total_cells > 0 && app.basin_n_nonconvergent > 0.6 * total_cells) {
    draw->AddText(ImVec2(14, 34), IM_COL32(255, 210, 120, 240),
                  "Most orbits don't settle to a point/cycle — this system has a chaotic or single global attractor.");
    draw->AddText(ImVec2(14, 52), IM_COL32(255, 210, 120, 240),
                  "Basins are meaningful for MULTISTABLE systems. Try the Duffing oscillator or the Newton fractal preset.");
  }
  if (!io.WantCaptureMouse)
    draw->AddText(ImVec2(14, h - 24), IM_COL32(150, 150, 160, 200),
                  "grey = didn't settle · black = escaped · colors = distinct attractors    drag: pan  wheel: zoom  double-click: reset");
}

/* ============================================================
 * PHASE B: 2-parameter scan ("shrimp" map).
 * Sweep two parameters over a grid; at each (p1,p2) estimate the largest
 * Lyapunov exponent (renormalized shadow orbit) and color by it. Periodic
 * windows (lambda < 0) form the characteristic shrimp shapes embedded in
 * the chaotic sea (lambda > 0). Works for ODEs (RK4 forced, like basins)
 * and maps.
 * ============================================================ */
ImU32 scan_color(double lyap, double lo, double hi) {
  if (!std::isfinite(lyap)) return IM_COL32(0, 0, 0, 255);
  if (lyap < 0) {
    /* periodic: cool, darker for more strongly stable */
    double t = (lo < 0) ? (lyap / lo) : 0.0; /* 0..1, 1 = most negative */
    t = std::max(0.0, std::min(1.0, t));
    int v = (int)(40 + 150 * (1.0 - t));
    return IM_COL32(v / 3, v / 2, v, 255); /* blue-ish, dark when stable */
  }
  /* chaotic: warm, brighter for larger exponent */
  double t = (hi > 0) ? (lyap / hi) : 0.0;
  t = std::max(0.0, std::min(1.0, t));
  int r = (int)(160 + 95 * t), g = (int)(90 + 120 * t), b = (int)(30 * (1.0 - t));
  return IM_COL32(r, g, b, 255);
}

void compute_scan_image(AppState &app, int W, int H, std::vector<uint32_t> &out, int step = 1) {
  if (W < 2) W = 2;
  if (H < 2) H = 2;
  if (step < 1) step = 1;
  out.assign((size_t)W * H, 0xff000000u);
  const size_t dim = app.state_names.size();
  if (dim == 0 || app.params.size() < 2) return;

  const int pxi = std::max(0, std::min(app.scan_px_index, (int)app.params.size() - 1));
  const int pyi = std::max(0, std::min(app.scan_py_index, (int)app.params.size() - 1));
  const double leps = app.lyapunov_epsilon > 0 ? app.lyapunov_epsilon : 1e-8;
  const int transient = std::max(50, app.scan_transient);
  const int iters = std::max(100, app.scan_iterations);

  /* force fixed-step for ODEs (adaptive would hang per cell) */
  const Integrator saved_integrator = app.integrator;
  if (app.mode == SystemMode::ODE &&
      (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45))
    app.integrator = Integrator::RK4;
  const std::vector<double> saved_params = app.param_values;

  std::vector<float> ly((size_t)W * H, std::numeric_limits<float>::quiet_NaN());
  double lo = 1e300, hi = -1e300;
  char err[128] = {0};

  for (int j = 0; j < H; j += step) {
    const double py = app.scan_ymin + (app.scan_ymax - app.scan_ymin) * (double)j / (H - 1);
    for (int i = 0; i < W; i += step) {
      const double px = app.scan_xmin + (app.scan_xmax - app.scan_xmin) * (double)i / (W - 1);
      if ((size_t)pxi < app.param_values.size()) app.param_values[(size_t)pxi] = px;
      if ((size_t)pyi < app.param_values.size()) app.param_values[(size_t)pyi] = py;

      State s = app.start; resize_state(s, dim);
      bool bad = false;
      for (int k = 0; k < transient; ++k) {
        State nx{};
        if (!step_state(app, s, &nx, err, sizeof(err))) { bad = true; break; }
        s = nx;
      }
      if (bad) { ly[(size_t)j * W + i] = std::numeric_limits<float>::quiet_NaN(); continue; }
      State shadow = s; resize_state(shadow, dim);
      shadow.v[0] += leps;
      double ly_sum = 0.0; long ly_n = 0;
      for (int k = 0; k < iters; ++k) {
        State n1{}, n2{};
        if (!step_state(app, s, &n1, err, sizeof(err))) { bad = true; break; }
        if (!step_state(app, shadow, &n2, err, sizeof(err))) { bad = true; break; }
        double d2 = 0.0;
        for (size_t q = 0; q < dim; ++q) { const double d = state_at(n2, q) - state_at(n1, q); d2 += d * d; }
        const double dist = std::sqrt(d2);
        if (dist > 1e-300 && std::isfinite(dist)) {
          ly_sum += std::log(dist / leps); ly_n++;
          const double sc = leps / dist;
          shadow = make_state_like(dim, n1.t);
          for (size_t q = 0; q < dim; ++q) {
            const double d = state_at(n2, q) - state_at(n1, q);
            set_state_at(shadow, q, state_at(n1, q) + d * sc);
          }
        }
        s = n1;
      }
      double val = std::numeric_limits<float>::quiet_NaN();
      if (!bad && ly_n > 0) {
        const double dt_factor = (app.mode == SystemMode::Map) ? 1.0 : app.dt;
        val = ly_sum / (ly_n * (dt_factor > 0 ? dt_factor : 1.0));
        if (std::isfinite(val)) { lo = std::min(lo, val); hi = std::max(hi, val); }
      }
      /* fill the whole step x step block with this sample (progressive) */
      for (int jj = j; jj < std::min(j + step, H); ++jj)
        for (int ii = i; ii < std::min(i + step, W); ++ii)
          ly[(size_t)jj * W + ii] = (float)val;
    }
  }
  app.param_values = saved_params; sync_param_values(app);
  app.integrator = saved_integrator;
  if (lo > hi) { lo = -1; hi = 1; }
  app.scan_lyap_min = lo; app.scan_lyap_max = hi;
  for (int idx = 0; idx < W * H; ++idx)
    out[(size_t)idx] = scan_color(ly[(size_t)idx], lo, hi);
}

void render_scan_background(AppState &app) {
  const float w = (float)app.window_width, h = (float)app.window_height;
  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  draw->AddRectFilled(ImVec2(0, 0), ImVec2(w, h), IM_COL32(10, 10, 14, 255));
  if (app.params.size() < 2) {
    draw->AddText(ImVec2(20, 48), IM_COL32(220, 220, 230, 230),
                  "2-parameter scan needs a system with >= 2 parameters.");
    draw->AddText(ImVec2(20, 66), IM_COL32(180, 180, 190, 220),
                  "Try the Henon or Tinkerbell map, or add a second parameter.");
    return;
  }

  if (!app.scan_view_init) {
    /* default the ranges to the two params' declared [min,max] */
    const int pxi = std::max(0, std::min(app.scan_px_index, (int)app.params.size() - 1));
    const int pyi = std::max(0, std::min(app.scan_py_index, (int)app.params.size() - 1));
    if (app.params[pxi].has_range) { app.scan_xmin = app.params[pxi].min_value; app.scan_xmax = app.params[pxi].max_value; }
    if (app.params[pyi].has_range) { app.scan_ymin = app.params[pyi].min_value; app.scan_ymax = app.params[pyi].max_value; }
    app.scan_view_init = true; app.scan_dirty = true;
  }

  /* pan/zoom the parameter window */
  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse) {
    auto invx = [&](float pxf){ return app.scan_xmin + (double)pxf / w * (app.scan_xmax - app.scan_xmin); };
    auto invy = [&](float pyf){ return app.scan_ymin + (double)(h - pyf) / h * (app.scan_ymax - app.scan_ymin); };
    if (io.MouseWheel != 0.0f) {
      const double mx = invx(io.MousePos.x), my = invy(io.MousePos.y);
      const double f = std::pow(0.8, clamped_wheel((double)io.MouseWheel));
      app.scan_xmin = mx + (app.scan_xmin - mx) * f; app.scan_xmax = mx + (app.scan_xmax - mx) * f;
      app.scan_ymin = my + (app.scan_ymin - my) * f; app.scan_ymax = my + (app.scan_ymax - my) * f;
      app.scan_dirty = true;
    }
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) && !ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      const double ddx = (double)io.MouseDelta.x / w * (app.scan_xmax - app.scan_xmin);
      const double ddy = (double)io.MouseDelta.y / h * (app.scan_ymax - app.scan_ymin);
      app.scan_xmin -= ddx; app.scan_xmax -= ddx;
      app.scan_ymin += ddy; app.scan_ymax += ddy;
      app.scan_dirty = true;
    }
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) { app.scan_view_init = false; app.scan_dirty = true; }
  }

  /* scans are expensive (each cell integrates a trajectory), so render
   * PROGRESSIVELY: a coarse pass fills the screen immediately, then it refines
   * over the next frames. This is what keeps ODE scans from freezing the UI. */
  const int CW = std::min(app.window_width, 400);
  const int CH = std::min(app.window_height, 320);
  const bool sforce = (app.scan_tex == 0 || app.scan_tex_w != CW || app.scan_tex_h != CH);
  if (app.scan_dirty) { app.scan_settle = 6; app.scan_dirty = false; }
  bool sstart = sforce;
  if (app.scan_settle > 0) { if (--app.scan_settle == 0) sstart = true; }
  if (sstart) app.scan_prog_level = 16;          /* (re)start the coarse->fine ladder */
  if (app.scan_prog_level > 0) {
    const int step = app.scan_prog_level;
    std::vector<uint32_t> img;
    compute_scan_image(app, CW, CH, img, step);
    if (app.scan_tex == 0) glGenTextures(1, &app.scan_tex);
    glBindTexture(GL_TEXTURE_2D, app.scan_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, CW, CH, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    app.scan_tex_w = CW; app.scan_tex_h = CH;
    /* next, finer level: 16 -> 8 -> 4 -> 2 -> 1 -> 0 (done) */
    app.scan_prog_level = (step > 1) ? step / 2 : 0;
  }
  if (app.scan_tex != 0)
    draw->AddImage((ImTextureID)(uintptr_t)app.scan_tex, ImVec2(0, 0), ImVec2(w, h));

  const int pxi = std::max(0, std::min(app.scan_px_index, (int)app.params.size() - 1));
  const int pyi = std::max(0, std::min(app.scan_py_index, (int)app.params.size() - 1));
  char hud[320];
  std::snprintf(hud, sizeof(hud),
                "2-parameter scan (largest Lyapunov)  |  x=%s [%.3g,%.3g]  y=%s [%.3g,%.3g]  |  lambda in [%.3f, %.3f]",
                app.params[pxi].name.c_str(), app.scan_xmin, app.scan_xmax,
                app.params[pyi].name.c_str(), app.scan_ymin, app.scan_ymax,
                app.scan_lyap_min, app.scan_lyap_max);
  draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), hud);
  draw->AddText(ImVec2(14, app.window_toolbar_h + 26.0f), IM_COL32(170, 170, 180, 220),
                "blue/dark = periodic (shrimps);  warm/bright = chaotic.  drag: pan  wheel: zoom");
  /* For the complex-quadratic system the two parameters ARE Re(c), Im(c), so
   * the Lyapunov scan of the c-plane genuinely reproduces the Mandelbrot set
   * (bounded orbits = negative exponent). Say so, since it understandably
   * looks just like the fractal view. */
  if (app.state_names.size() >= 2 &&
      (std::strstr(app.system_input, "x*x - y*y") != nullptr ||
       std::strstr(app.system_input, "z^2 + c") != nullptr)) {
    draw->AddText(ImVec2(14, app.window_toolbar_h + 44.0f), IM_COL32(150, 175, 150, 220),
                  "(this is the c-plane: the Lyapunov scan here IS the Mandelbrot set -- bounded = periodic = dark)");
  }
}

/* ============================================================
 * PHASE C+: the "3D bridge" — Mandelbrot set + bifurcation diagram in one
 * 3D scene.
 *
 * The deep fact being visualized: the logistic map x -> r x(1-x) and the
 * complex map z -> z^2 + c are conjugate, so the logistic period-doubling
 * cascade corresponds exactly to the real-axis ("needle") of the
 * Mandelbrot set. We draw the Mandelbrot set lying flat in the world x-z
 * plane (re -> x, im -> z) at height y=0, and the logistic bifurcation
 * diagram rising in world-y above the real axis (im=0 line). Where a
 * Mandelbrot bulb sits on the real axis, the bifurcation diagram shows
 * the matching period — the same cascade seen two ways.
 *
 * Geometry is built once into a dedicated point-cloud VBO and rendered
 * with the existing camera (drag-rotate / wheel-zoom).
 *
 * Coordinate mapping (world units, roughly in [-25,25] like the attractor
 * scenes): the Mandelbrot c-plane re in [-2.2, 0.7] -> x, im in
 * [-1.2,1.2] -> z; logistic r in [2.6, 4.0] mapped onto the same real-axis
 * span as the Mandelbrot needle so the cascades align; attractor x in
 * [0,1] -> world-y * bridge_height.
 * ============================================================ */
namespace bridge_detail {
/* Mandelbrot escape time from z0 = 0; returns iters (maxit if bounded). */
inline int mandel_escape(double cre, double cim, int maxit, double R2) {
  double x = 0, y = 0;
  int it = 0;
  for (; it < maxit; ++it) {
    const double nx = x * x - y * y + cre;
    const double ny = 2 * x * y + cim;
    x = nx; y = ny;
    if (x * x + y * y > R2) { ++it; break; }
  }
  return it;
}

/* For a point c INSIDE the Mandelbrot set (bounded orbit), detect the period
 * of the attracting cycle: settle the orbit, then find the smallest p such
 * that z_{n+p} ~ z_n. Returns 0 if it escapes, or the detected period (capped
 * at maxp; 0/uncertain -> returns maxp+1 meaning "high/chaotic"). This is what
 * makes the cascade<->bulb correspondence visible: period-2 window under the
 * period-2 bulb, period-3 under the period-3 bulb, and so on. */
inline int mandel_period(double cre, double cim, int settle, int maxp, double R2) {
  double x = 0, y = 0;
  for (int k = 0; k < settle; ++k) {
    const double nx = x * x - y * y + cre, ny = 2 * x * y + cim;
    x = nx; y = ny;
    if (x * x + y * y > R2) return 0; /* escaped */
  }
  const double rx = x, ry = y;            /* reference point on the cycle */
  for (int p = 1; p <= maxp; ++p) {
    const double nx = x * x - y * y + cre, ny = 2 * x * y + cim;
    x = nx; y = ny;
    if (x * x + y * y > R2) return 0;
    const double dx = x - rx, dy = y - ry;
    if (dx * dx + dy * dy < 1e-12) return p;
  }
  return maxp + 1; /* long/!detected period within budget */
}

/* a distinct hue per period, so bulbs of different period read at a glance.
 * Writes r,g,b in [0,1]. period 1 (main cardioid) = warm; 2,3,4,... cycle hues. */
inline void period_color(int period, float *r, float *g, float *b) {
  if (period <= 0) { *r = 0.05f; *g = 0.06f; *b = 0.12f; return; } /* exterior */
  static const float pal[8][3] = {
    {0.95f, 0.75f, 0.25f}, /* 1  amber  */
    {0.30f, 0.70f, 1.00f}, /* 2  blue   */
    {0.45f, 0.90f, 0.45f}, /* 3  green  */
    {1.00f, 0.45f, 0.55f}, /* 4  rose   */
    {0.75f, 0.55f, 1.00f}, /* 5  violet */
    {1.00f, 0.65f, 0.30f}, /* 6  orange */
    {0.40f, 0.85f, 0.85f}, /* 7  teal   */
    {0.85f, 0.85f, 0.40f}, /* 8+ olive  */
  };
  const int idx = (period >= 1 && period <= 8) ? period - 1 : 7;
  *r = pal[idx][0]; *g = pal[idx][1]; *b = pal[idx][2];
}
}  // namespace bridge_detail

void build_bridge_geometry(AppState &app) {
  std::vector<Point> pts;
  std::vector<float> cols;
  pts.reserve(200000);
  cols.reserve(600000);

  /* ---- CURRENT-SYSTEM mode: lift THIS system's bifurcation diagram into 3D.
   * x = primary parameter (bif_param), y = attractor value (height), z = a
   * second parameter (bridge_param2) swept over a few stacked slices. This
   * generalizes the bridge to ANY 1-D observable of ANY system, beyond the
   * classical logistic/Mandelbrot correspondence. The sweep is self-contained
   * (its own step loop on a state copy) so it doesn't disturb the 2D
   * bifurcation view's data or the integrator setting. ---- */
  if (app.bridge_mode == AppState::BridgeMode::CurrentSystem) {
    AppState::Param *p1 = find_param(app, app.bif_param);
    AppState::Param *p2 = find_param(app, app.bridge_param2);
    const size_t dim = app.state_names.size();
    if (p1 == nullptr || dim == 0) { app.bridge_point_count = 0; app.bridge_built = true; return; }

    /* observable: a state index if bif_observable names one, else state 0. */
    size_t obs_idx = 0;
    for (size_t i = 0; i < app.state_names.size(); ++i)
      if (app.state_names[i] == app.bif_observable) { obs_idx = i; break; }

    /* world extents (shared by the floor and the lifted diagram so they line
     * up along the parameter axis — this is the whole point of the bridge). */
    const double WX = 44.0, WZ = 38.0;
    auto wx = [&](double u) { return (float)(u * WX - WX * 0.5); };       /* u in [0,1] */
    auto wz = [&](double v) { return (float)(v * WZ - WZ * 0.5); };       /* v in [0,1] */

    /* ---- the FRACTAL FLOOR (z=0 plane) -----------------------------------
     * The bridge's real idea: show the object that GENERATES the cascade with
     * the cascade rising out of it, sharing the parameter axis. For a MAP that
     * has an escape-to-infinity criterion, the natural floor is the
     * parameter-space escape set: x = the bifurcation parameter (same axis as
     * the diagram above), and the perpendicular axis sweeps the system's first
     * state variable's initial value. A point is "in the set" when the orbit
     * stays bounded — exactly the Mandelbrot construction, but for THIS map.
     * (For ODEs there's no escape fractal; the floor is skipped and only the
     * lifted diagram is drawn.) The floor is laid in the plane y=0, with the
     * bifurcation parameter mapped to the SAME wx() the diagram uses. */
    if (app.bridge_show_mandelbrot && app.mode == SystemMode::Map) {
      const int RES = std::max(80, app.bridge_mandel_res);
      const int fmaxit = 120;
      const double FR2 = 100.0;            /* escape radius^2 for the floor */
      const double p1lo = app.bif_start, p1hi = app.bif_end;
      /* perpendicular axis: first state initial value over a sensible window */
      const double s0 = state_at(app.start, 0);
      const double ic_lo = s0 - 1.0, ic_hi = s0 + 1.0;
      const double pin = (p1 ? p1->value : 0.0);
      State cur = app.start, nxt = app.start;
      resize_state(cur, dim); resize_state(nxt, dim);
      char ferr[128] = {0};
      for (int jz = 0; jz < RES; ++jz) {
        const double v = (double)jz / (RES - 1);
        const double ic = ic_lo + (ic_hi - ic_lo) * v;
        for (int ix2 = 0; ix2 < RES; ++ix2) {
          const double u = (double)ix2 / (RES - 1);
          if (p1) p1->value = p1lo + (p1hi - p1lo) * u;
          sync_param_values(app);
          for (size_t i = 0; i < dim; ++i) set_state_at(cur, i, state_at(app.start, i));
          set_state_at(cur, 0, ic);
          bool escaped = false; int it = 0;
          for (; it < fmaxit; ++it) {
            if (!step_map_state(app, cur, &nxt, ferr, sizeof(ferr))) { escaped = true; break; }
            const double xx = state_at(nxt, 0), yy = (dim > 1) ? state_at(nxt, 1) : 0.0;
            const double r2 = xx * xx + yy * yy;
            if (!std::isfinite(r2) || r2 > FR2) { escaped = true; break; }
            for (size_t i = 0; i < dim; ++i) set_state_at(cur, i, state_at(nxt, i));
          }
          if (!escaped) { /* bounded: part of the floor "set" */
            pts.push_back(Point{wx(u), 0.0f, wz(v)});
            cols.push_back(0.10f); cols.push_back(0.18f); cols.push_back(0.42f);
          } else if (it > 3 && (it % 2 == 0)) { /* sparse exterior bands for context */
            const double t = std::fmod(it * 0.06, 1.0);
            pts.push_back(Point{wx(u), 0.0f, wz(v)});
            cols.push_back(0.15f + 0.25f * (float)t);
            cols.push_back(0.10f + 0.20f * (float)t);
            cols.push_back(0.25f + 0.15f * (float)t);
          }
        }
      }
      if (p1) p1->value = pin;
      sync_param_values(app);
    }

    /* force fixed-step during the sweep (restore after) */
    const Integrator saved_integrator = app.integrator;
    if (app.mode == SystemMode::ODE &&
        (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45))
      app.integrator = Integrator::RK4;
    const double p1_old = p1->value;
    const double p2_old = p2 ? p2->value : 0.0;

    const int xslices = std::max(60, app.bridge_bif_slices);
    const int zslices = p2 ? std::max(1, app.bridge_p2_slices) : 1;
    const int discard = (app.mode == SystemMode::Map) ? 600 : 1500;
    const int keep = std::max(20, app.bridge_bif_keep);

    /* estimate the observable's vertical range for a sane height scale */
    double obs_lo = 1e300, obs_hi = -1e300;
    char err[256] = {0};
    /* points pushed so far are the FLOOR (y=0); the diagram points come next.
     * Only the diagram points get height-rescaled/recolored below. */
    const size_t floor_count = pts.size();

    for (int zi = 0; zi < zslices; ++zi) {
      /* When there's no 2nd parameter we draw ONE diagram; place it at the
       * floor's center line (v=0.5) so it visibly rises out of the fractal
       * floor, mirroring how the logistic cascade sits on the Mandelbrot
       * real axis. With a 2nd parameter we stack diagrams across z. */
      const double v = (zslices > 1) ? (double)zi / (zslices - 1) : 0.5;
      if (p2) { p2->value = app.bridge_p2_min + v * (app.bridge_p2_max - app.bridge_p2_min); }
      for (int xi = 0; xi < xslices; ++xi) {
        const double u = (double)xi / (xslices - 1);
        p1->value = app.bif_start + u * (app.bif_end - app.bif_start);
        sync_param_values(app); /* IR reads param_values, not ->value */
        State s = app.start; resize_state(s, dim);
        bool bad = false;
        for (int k = 0; k < discard; ++k) {
          State o{};
          if (!step_state(app, s, &o, err, sizeof(err))) { bad = true; break; }
          s = o;
        }
        if (bad) continue;
        /* Sampling differs by system kind, exactly like the 2-D bifurcation
         * view: for a MAP every iterate is a diagram point; for an ODE we take
         * a Poincare-style section — the LOCAL MAXIMA of the observable — via
         * a 3-sample window, otherwise an ODE just smears into a band. */
        if (app.mode == SystemMode::Map) {
          for (int k = 0; k < keep; ++k) {
            State o{};
            if (!step_state(app, s, &o, err, sizeof(err))) break;
            s = o;
            const double val = state_at(s, obs_idx);
            if (!std::isfinite(val)) break;
            if (val < obs_lo) obs_lo = val;
            if (val > obs_hi) obs_hi = val;
            pts.push_back(Point{wx(u), (float)val, wz(v)});
            cols.push_back(0); cols.push_back(0); cols.push_back(0); /* recolored below */
          }
        } else {
          /* ODE: prefer local maxima of the observable (a Poincare-style
           * section), but COLLECT enough trajectory to find them — short
           * windows on slow oscillators yield no maxima and a black slice.
           * If a slice still finds too few maxima, fall back to raw subsampled
           * values so the slice is never empty. */
          const int odekeep = std::max(keep, 1200);
          std::vector<double> raw; raw.reserve((size_t)odekeep);
          double w0 = 0, w1 = 0, w2 = 0; int filled = 0; int maxima = 0;
          for (int k = 0; k < odekeep; ++k) {
            State o{};
            if (!step_state(app, s, &o, err, sizeof(err))) break;
            s = o;
            const double val = state_at(s, obs_idx);
            if (!std::isfinite(val)) break;
            raw.push_back(val);
            w0 = w1; w1 = w2; w2 = val; if (filled < 3) ++filled;
            if (filled == 3 && w1 > w0 && w1 >= w2) { /* local maximum */
              if (w1 < obs_lo) obs_lo = w1;
              if (w1 > obs_hi) obs_hi = w1;
              pts.push_back(Point{wx(u), (float)w1, wz(v)});
              cols.push_back(0); cols.push_back(0); cols.push_back(0); /* recolored below */
              ++maxima;
            }
          }
          /* fallback: an equilibrium or slow drift produced (almost) no maxima
           * — subsample the raw trajectory so the slice still shows something. */
          if (maxima < 3 && raw.size() >= 8) {
            const int want = std::min((int)raw.size(), std::max(20, keep));
            const int stride = std::max(1, (int)raw.size() / want);
            for (size_t k = 0; k < raw.size(); k += (size_t)stride) {
              const double val = raw[k];
              if (val < obs_lo) obs_lo = val;
              if (val > obs_hi) obs_hi = val;
              pts.push_back(Point{wx(u), (float)val, wz(v)});
              cols.push_back(0); cols.push_back(0); cols.push_back(0);
            }
          }
        }
      }
    }
    p1->value = p1_old; if (p2) p2->value = p2_old; sync_param_values(app);
    app.integrator = saved_integrator;

    /* rescale heights into a pleasing band and recolor by normalized height —
     * for the DIAGRAM points only (indices >= floor_count); the floor keeps
     * its y=0 and its blue/exterior coloring. */
    const double span = (obs_hi > obs_lo) ? (obs_hi - obs_lo) : 1.0;
    for (size_t i = floor_count; i < pts.size(); ++i) {
      const double tnorm = (pts[i].y - obs_lo) / span; /* 0..1 */
      pts[i].y = (float)(tnorm * app.bridge_height);
      const float t = (float)tnorm;
      cols[i * 3 + 0] = 0.2f + 0.8f * t;
      cols[i * 3 + 1] = 0.7f + 0.3f * t;
      cols[i * 3 + 2] = 0.9f - 0.6f * t;
    }

    app.bridge_point_count = (int)pts.size();
    if (pts.empty()) { app.bridge_built = true; return; }
    if (app.bridge_vao == 0) glGenVertexArrays(1, &app.bridge_vao);
    if (app.bridge_vbo == 0) glGenBuffers(1, &app.bridge_vbo);
    if (app.bridge_cbo == 0) glGenBuffers(1, &app.bridge_cbo);
    glBindVertexArray(app.bridge_vao);
    glBindBuffer(GL_ARRAY_BUFFER, app.bridge_vbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(pts.size() * sizeof(Point)), pts.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Point), (void *)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, app.bridge_cbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(cols.size() * sizeof(float)), cols.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    app.bridge_built = true;
    return;
  }

  /* =====================================================================
   * PROJECTION-SOLID bridge.
   *
   * A SINGLE 3D point cloud whose two perpendicular 2D shadows are exactly
   * the Mandelbrot set and the logistic bifurcation diagram.
   *
   *   X = c          the real/parameter axis, shared by both images
   *                  (logistic r maps in by the conjugacy c = r/2 - r^2/4)
   *   Z = im         imaginary part  -> the X-Z shadow is the Mandelbrot set
   *   Y = x          attractor value -> the X-Y shadow is the bifurcation diagram
   *
   * For each slice at fixed c the cross-section is the Cartesian product
   *     (Mandelbrot imaginary membership in Z)  x  (logistic attractor fiber in Y).
   * Projecting that product onto X-Z collapses Y and recovers the set; onto
   * X-Y collapses Z and recovers the cascade. So the one solid casts both
   * pictures, and you rotate it freely to see how they belong to each other.
   * ===================================================================== */
  if (app.bridge_mode == AppState::BridgeMode::ProjectionSolid) {
    /* The unified construction (after thejoggeli/mandelbrot-logistic):
     * iterate z -> z^2 + c over the WHOLE complex plane, X = Re(c), Z = Im(c).
     * For each c IN the Mandelbrot set, take the orbit's periodic tail and plot
     * every value it visits as a point at height  Y = -Re(z) * scale.
     * This is ONE object from ONE computation:
     *   - its shadow on X-Z (look down Y) is the Mandelbrot set (only in-set
     *     points exist), and
     *   - its shadow on X-Y (look down Z) is the bifurcation diagram: along the
     *     real axis z^2+c is conjugate to the logistic map, so the orbit tail
     *     traces the period-doubling cascade, and the off-axis orbits fill the
     *     body that connects the set to the cascade.
     * The Mandelbrot set and the logistic bifurcation are literally two
     * projections of the same 3-D structure. */
    const bool cubic = (app.bridge_family == AppState::BridgeFamily::Cubic);
    const bool sine  = (app.bridge_family == AppState::BridgeFamily::Sine);
    /* Box framing per family: quadratic set sits around c_re in [-2,0.25];
     * cubic is disc-shaped about the origin (|c|<~1.5); the sine family
     * c*sin(z) has its interesting parameter range along positive real c with
     * a tall imaginary extent. */
    const double pos_x = sine ? 1.9 : (cubic ? 0.0 : -1.0); /* center of the c-plane box */
    const double pos_z = 0.0;
    const double box_x = sine ? 3.6 : (cubic ? 3.2 : 4.0);
    const double box_z = sine ? 5.0 : (cubic ? 3.2 : 2.0);
    const double scale_y = sine ? 0.5 : (2.0 / 3.0);
    /* world scale: map the c-box to ~50 world units (matches camera) */
    const double WSCALE = sine ? 11.0 : (cubic ? 16.0 : 13.0);
    auto W = [&](double v) { return (float)(v * WSCALE); };

    const int density = std::max(60, app.bridge_mandel_res / 4); /* voxels per unit */
    int nx = (int)(density * box_x); if (nx % 2 == 0) ++nx;
    int nz = (int)(density * box_z); if (nz % 2 == 0) ++nz;
    const double step_x = box_x / nx, step_z = box_z / nz;
    const double off_x = pos_x - step_x * (nx - 1) * 0.5;
    const double off_z = pos_z - step_z * (nz - 1) * 0.5;

    const int max_iter = 4096;
    /* The sine family is chaotic over much of its range and its bounded orbits
     * run long, which would emit millions of points. Cap its tail length hard
     * and subsample below so the cloud stays comparable to the polynomial
     * families (~a few hundred k points). */
    const int max_tail = sine ? std::min(96, std::max(48, app.bridge_bif_keep))
                              : std::min(max_iter, std::max(64, app.bridge_bif_keep * 3));
    const int emit_stride = sine ? 3 : 1;   /* keep every Nth voxel-orbit point */
    /* Boundedness test. Quadratic/cubic: |z|>R diverges (escape radius). Sine:
     * c*sin(z) sends orbits to i*infinity when they escape, so the right test
     * is on |Im z| growing without bound, not |z|. */
    const double esc2 = cubic ? 9.0 : 4.0;
    const double sine_im_bail = 30.0;
    std::vector<double> zre((size_t)max_iter);

    for (int vx = 0; vx < nx; ++vx) {
      for (int vz = 0; vz < nz; ++vz) {
        const double c_re = off_x + vx * step_x;
        const double c_im = off_z + vz * step_z;
        /* sine orbits start from z=0.5 (a generic point); polynomial families
         * start from the critical point z=0. */
        double z_re = sine ? 0.5 : 0.0, z_im = 0.0;
        bool in_set = true;
        int iters = 0;
        for (int i = 0; i < max_iter; ++i) {
          const double a = z_re, b = z_im;
          if (sine) {
            /* c * sin(z), z=a+bi: sin(a+bi) = sin a cosh b + i cos a sinh b,
             * then multiply by c = c_re + i c_im. */
            const double sr = std::sin(a) * std::cosh(b);
            const double si = std::cos(a) * std::sinh(b);
            z_re = c_re * sr - c_im * si;
            z_im = c_re * si + c_im * sr;
            if (!std::isfinite(z_re) || std::fabs(z_im) > sine_im_bail) { in_set = false; break; }
          } else if (cubic) {
            /* z^3 + c, with z = a+bi:  z^3 = (a^3 - 3ab^2) + (3a^2 b - b^3) i */
            z_re = a * a * a - 3.0 * a * b * b + c_re;
            z_im = 3.0 * a * a * b - b * b * b + c_im;
            if (z_re * z_re + z_im * z_im > esc2) { in_set = false; break; }
          } else {
            z_re = a * a - b * b + c_re;
            z_im = 2 * a * b + c_im;
            if (z_re * z_re + z_im * z_im > esc2) { in_set = false; break; }
          }
          zre[(size_t)i] = z_re;
          iters = i + 1;
        }
        if (!in_set || iters < 2) continue;

        /* find the periodic tail: walk back from the last value to the most
         * recent exact repeat (the cycle), like the reference. Bounded by
         * max_tail so dense periodic regions stay affordable. */
        const int end = iters;
        const int lo = std::max(0, end - max_tail);
        int start = lo;
        const double last = zre[(size_t)end - 1];
        for (int j = end - 2; j >= lo; --j) {
          if (zre[(size_t)j] == last) { start = j + 1; break; }
        }
        int tail = end - start;
        if (tail < 1) { start = lo; tail = end - start; }

        /* colour: by cycle length (period), brighter for short periods; the
         * real axis (the logistic cascade) is emphasised slightly. */
        float pr, pg, pb;
        if (app.bridge_color_by_period) {
          int per = tail; if (per < 1) per = 1; if (per > 64) per = 0;
          bridge_detail::period_color(per, &pr, &pg, &pb);
        } else {
          const float g = std::min(1.0f, (float)tail / (float)max_tail);
          pr = 0.5f + g; pg = g; pb = g;
        }
        const float fx = W(c_re), fz = W(c_im);
        /* dim each point so additive blending accumulates to the right
         * brightness where many orbit points overlap (rather than clipping). */
        const float dim = 0.32f;
        for (int j = start; j < end; j += emit_stride) {
          const float fy = (float)(-zre[(size_t)j] * scale_y) * (float)WSCALE;
          pts.push_back(Point{fx, fy, fz});
          cols.push_back(pr * dim); cols.push_back(pg * dim); cols.push_back(pb * dim);
        }
      }
    }

    app.bridge_point_count = (int)pts.size();
    if (pts.empty()) { app.bridge_built = true; return; }
    if (app.bridge_vao == 0) glGenVertexArrays(1, &app.bridge_vao);
    if (app.bridge_vbo == 0) glGenBuffers(1, &app.bridge_vbo);
    if (app.bridge_cbo == 0) glGenBuffers(1, &app.bridge_cbo);
    glBindVertexArray(app.bridge_vao);
    glBindBuffer(GL_ARRAY_BUFFER, app.bridge_vbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(pts.size() * sizeof(Point)), pts.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Point), (void *)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, app.bridge_cbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(cols.size() * sizeof(float)), cols.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    app.bridge_built = true;
    return;
  }

  /* world mapping helpers */
  const double re_lo = -2.2, re_hi = 0.7;
  const double im_lo = -1.2, im_hi = 1.2;
  const double WX = 44.0; /* world width across re range */
  const double WZ = 38.0; /* world depth across im range */
  auto wx = [&](double re) { return (float)((re - re_lo) / (re_hi - re_lo) * WX - WX * 0.5); };
  auto wz = [&](double im) { return (float)((im - im_lo) / (im_hi - im_lo) * WZ - WZ * 0.5); };

  /* ---- Mandelbrot set in the y=0 plane ---- */
  if (app.bridge_show_mandelbrot) {
    const int RES = std::max(80, app.bridge_mandel_res);
    const int maxit = 200;
    const double R2 = 16.0;
    for (int j = 0; j < RES; ++j) {
      const double im = im_lo + (im_hi - im_lo) * j / (RES - 1);
      for (int i = 0; i < RES; ++i) {
        const double re = re_lo + (re_hi - re_lo) * i / (RES - 1);
        const int it = bridge_detail::mandel_escape(re, im, maxit, R2);
        const bool inside = (it >= maxit);
        if (inside) {
          /* INTERIOR: color by the period of the attracting cycle, so the
           * period-n bulbs are visually distinct and line up with the
           * period-n windows of the cascade rising above the real axis. */
          float r, g, b;
          if (app.bridge_color_by_period) {
            const int per = bridge_detail::mandel_period(re, im, 250, 64, R2);
            bridge_detail::period_color(per, &r, &g, &b);
            /* dim slightly so the rising cascade stays the bright foreground */
            r *= 0.55f; g *= 0.55f; b *= 0.6f;
          } else {
            r = 0.10f; g = 0.18f; b = 0.42f;
          }
          pts.push_back(Point{wx(re), 0.0f, wz(im)});
          cols.push_back(r); cols.push_back(g); cols.push_back(b);
        } else if (it > 4 && (it % 2 == 0)) {
          /* sparse colored exterior bands for context, near the set only */
          const double t = std::fmod(it * 0.06, 1.0);
          pts.push_back(Point{wx(re), 0.0f, wz(im)});
          cols.push_back(0.15f + 0.25f * (float)t);
          cols.push_back(0.10f + 0.20f * (float)t);
          cols.push_back(0.25f + 0.15f * (float)t);
        }
      }
    }
  }

  /* ---- logistic bifurcation diagram rising in +y above the real axis ---- *
   * Map logistic r in [1,4] onto the Mandelbrot real-axis span via the
   * conjugacy c = r/2 - r^2/4: r=1->c=0.25 (the cusp of the cardioid) down to
   * r=4->c=-2 (the antenna tip). Spanning the FULL real-axis range makes the
   * cascade's footprint coincide with the Mandelbrot set's real-axis slice, so
   * the two diagrams share the same c-axis edge and read as one structure. The
   * diagram lives in the vertical plane z=0 (the real axis), the set in the
   * horizontal floor — rotate to see each face and their shared spine. */
  if (app.bridge_show_bifurcation) {
    const int slices = std::max(100, app.bridge_bif_slices);
    const int discard = 800;
    const int keep = std::max(20, app.bridge_bif_keep);
    auto c_of_r = [](double r) { return r / 2.0 - r * r / 4.0; };
    /* detect the logistic period at parameter r (post-transient): smallest p
     * with x_{n+p} ~ x_n. Lets us color the cascade to MATCH the bulb colors
     * and draw period-keyed droplines, making the correspondence explicit. */
    auto logistic_period = [](double r, int settle, int maxp) -> int {
      double x = 0.5;
      for (int k = 0; k < settle; ++k) x = r * x * (1.0 - x);
      const double rx = x;
      for (int p = 1; p <= maxp; ++p) {
        x = r * x * (1.0 - x);
        if (std::fabs(x - rx) < 1e-9) return p;
      }
      return maxp + 1;
    };
    for (int s = 0; s < slices; ++s) {
      const double r = 1.0 + (4.0 - 1.0) * s / (slices - 1);
      double x = 0.5;
      for (int k = 0; k < discard; ++k) x = r * x * (1.0 - x);
      const double cre = c_of_r(r);
      const int per = app.bridge_color_by_period ? logistic_period(r, 2000, 64) : 0;
      float pr = 0.5f, pg = 0.5f, pb = 0.5f;
      if (app.bridge_color_by_period) bridge_detail::period_color(per > 64 ? 0 : per, &pr, &pg, &pb);
      for (int k = 0; k < keep; ++k) {
        x = r * x * (1.0 - x);
        if (!std::isfinite(x)) break;
        const float wy = (float)x * app.bridge_height;
        pts.push_back(Point{wx(cre), wy, wz(0.0)});
        if (app.bridge_color_by_period && per <= 8) {
          /* match the bulb hue, brightened toward the top of the cascade */
          const float t = (float)x;
          cols.push_back(std::min(1.0f, pr + 0.25f * t));
          cols.push_back(std::min(1.0f, pg + 0.15f * t));
          cols.push_back(std::min(1.0f, pb));
        } else {
          const float t = (float)x;
          cols.push_back(0.2f + 0.8f * t);
          cols.push_back(0.7f + 0.3f * t);
          cols.push_back(0.9f - 0.6f * t);
        }
      }
      /* periodic-window droplines: in a low-period window, drop a faint line
       * from the cascade down to the real-axis floor at the same c, visually
       * tying the window to the bulb it sits above. */
      if (app.bridge_color_by_period && per >= 1 && per <= 8 && (s % 3 == 0)) {
        const int LN = 10;
        for (int t = 0; t <= LN; ++t) {
          const float wy = (float)((double)t / LN) * app.bridge_height * 0.5f;
          pts.push_back(Point{wx(cre), wy, wz(0.0)});
          cols.push_back(pr * 0.5f); cols.push_back(pg * 0.5f); cols.push_back(pb * 0.5f);
        }
      }
    }

    /* The SHARED SPINE: the real axis c in [-2, 0.25] at y=0, z=0 is the common
     * edge where the Mandelbrot floor and the vertical bifurcation diagram meet.
     * Drawing it as a bright line makes the join explicit — as you orbit, this
     * is the hinge between the set (seen from the side) and the cascade (seen
     * from the front). */
    {
      const int SN = std::max(200, slices);
      for (int s = 0; s <= SN; ++s) {
        const double c = -2.0 + (0.25 - (-2.0)) * (double)s / SN;
        pts.push_back(Point{wx(c), 0.0f, wz(0.0)});
        cols.push_back(1.0f); cols.push_back(0.92f); cols.push_back(0.55f);
      }
    }
  }

  app.bridge_point_count = (int)pts.size();
  if (pts.empty()) { app.bridge_built = true; return; }

  if (app.bridge_vao == 0) glGenVertexArrays(1, &app.bridge_vao);
  if (app.bridge_vbo == 0) glGenBuffers(1, &app.bridge_vbo);
  if (app.bridge_cbo == 0) glGenBuffers(1, &app.bridge_cbo);
  glBindVertexArray(app.bridge_vao);
  glBindBuffer(GL_ARRAY_BUFFER, app.bridge_vbo);
  glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(pts.size() * sizeof(Point)), pts.data(), GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Point), (void *)0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, app.bridge_cbo);
  glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(cols.size() * sizeof(float)), cols.data(), GL_STATIC_DRAW);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(1);
  glBindVertexArray(0);
  app.bridge_built = true;
}

void render_bridge_scene(AppState &app) {
  glViewport(0, 0, app.window_width, app.window_height);
  const float aspect = (float)app.window_width / std::max(1.0f, (float)app.window_height);

  if (!app.bridge_built) build_bridge_geometry(app);

  /* Frame the geometry on first entry (or when it was just rebuilt for a new
   * mode). The point clouds span roughly +/-23 in X, 0..26 in Y, +/-17 in Z,
   * far larger than the default scene_scale expects, so without this the whole
   * structure renders as a dot in the corner. A 3/4 view shows both shadow
   * faces of the projection solid at once. */
  if (!app.bridge_cam_init) {
    app.bridge_cam_init = true;
    const bool proj = (app.bridge_mode == AppState::BridgeMode::ProjectionSolid);
    const bool classic = (app.bridge_mode == AppState::BridgeMode::Classic);
    const bool cubic = (app.bridge_family == AppState::BridgeFamily::Cubic);
    const bool sine  = (app.bridge_family == AppState::BridgeFamily::Sine);
    /* Projection-solid uses an orthographic camera so that looking down an axis
     * gives a TRUE 2D projection (the shadows ARE the two images). The other
     * modes use the normal perspective camera. */
    app.orthographic_3d = proj;
    app.zoom = 1.15f;
    app.camera_distance = 60.0f;
    app.center_z = 0.0f;
    app.angle_x = 26.0f;
    app.angle_y = -48.0f;
    if (proj) {
      app.scene_scale = sine ? 1.05f : (cubic ? 0.55f : 0.62f);
      app.center_x = sine ? -19.0f : (cubic ? 0.0f : 13.0f);
      app.center_y = sine ? 1.0f : -1.0f;
    } else if (classic) {
      /* side-by-side: Mandelbrot floor + vertical cascade, ~±22 world units. */
      app.scene_scale = 0.42f;
      app.center_x = 0.0f;
      app.center_y = -4.0f;
      app.angle_x = 32.0f;
      app.angle_y = -46.0f;
    } else {
      /* current-system: bifurcation lifted into 3D, normalized geometry. */
      app.scene_scale = 0.42f;
      app.center_x = 0.0f;
      app.center_y = 0.0f;
      app.angle_x = 24.0f;
      app.angle_y = -40.0f;
    }
    if (const char *ax = std::getenv("DYNSYS_ANGLE_X")) app.angle_x = (float)std::atof(ax);
    if (const char *ay = std::getenv("DYNSYS_ANGLE_Y")) app.angle_y = (float)std::atof(ay);
    if (const char *sc = std::getenv("DYNSYS_SCALE"))   app.scene_scale = (float)std::atof(sc);
    if (const char *cy = std::getenv("DYNSYS_CENTER_Y")) app.center_y = (float)std::atof(cy);
    if (const char *cx = std::getenv("DYNSYS_CENTER_X")) app.center_x = (float)std::atof(cx);
    update_projection(app);
  }

  const bool proj_solid = (app.bridge_mode == AppState::BridgeMode::ProjectionSolid);
  if (proj_solid) {
    /* Translucent additive look (as in the reference): overlapping orbit
     * points accumulate so the dense periodic blobs become see-through and the
     * real-axis cascade shows through the body. Depth test off so nothing is
     * hidden; additive blend brightens where many orbits coincide. */
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  } else {
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
  }
  glClearColor(0.03f, 0.03f, 0.05f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glUseProgram(app.shader_program);

  mat4 projection;
  if (app.orthographic_3d) {
    const float hh = 30.0f / std::max(app.zoom, 0.001f);
    glm_ortho(-hh * aspect, hh * aspect, -hh, hh, -10000.0f, 10000.0f, projection);
  } else {
    glm_perspective(glm_rad(45.0f), aspect, 0.01f, 10000.0f, projection);
  }
  mat4 model, view;
  glm_mat4_identity(model);
  glm_mat4_identity(view);
  vec3 tr{}; set_vec3(tr, app.center_x, app.center_y, app.center_z); glm_translate(model, tr);
  const float es = app.orthographic_3d ? app.scene_scale : app.scene_scale * app.zoom;
  vec3 sc{}; set_vec3(sc, es, es, es); glm_scale(model, sc);
  vec3 ax{}; set_vec3(ax, 1, 0, 0); glm_rotate(model, glm_rad(app.angle_x), ax);
  vec3 ay{}; set_vec3(ay, 0, 1, 0); glm_rotate(model, glm_rad(app.angle_y), ay);
  vec3 vt{}; set_vec3(vt, 0, 0, -app.camera_distance); glm_translate(view, vt);
  glUniformMatrix4fv(glGetUniformLocation(app.shader_program, "model"), 1, GL_FALSE, (float *)model);
  glUniformMatrix4fv(glGetUniformLocation(app.shader_program, "view"), 1, GL_FALSE, (float *)view);
  glUniformMatrix4fv(glGetUniformLocation(app.shader_program, "projection"), 1, GL_FALSE, (float *)projection);

  if (app.bridge_point_count > 0 && app.bridge_vao != 0) {
    glBindVertexArray(app.bridge_vao);
    glDrawArrays(GL_POINTS, 0, app.bridge_point_count);
    glBindVertexArray(0);
  }
  /* restore default GL state so the additive/no-depth bridge settings don't
   * leak into other views drawn later in the frame */
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
  /* camera interaction (same scheme as the attractor scene) */
  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse) {
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
      app.angle_y += io.MouseDelta.x * 0.4f;
      app.angle_x += io.MouseDelta.y * 0.4f;
    }
    /* PAN: right- or middle-drag shifts the view center so the structure moves
     * around the screen. center_x/y are in WORLD units (applied after the model
     * scale), so to track the cursor 1:1 we convert pixels->world using the
     * actual on-screen extent: orthographic spans +-hh*aspect horizontally over
     * the window width; perspective approximates the same at the structure's
     * depth. Previously this used 0.09/scene_scale, which was ~30x too small in
     * ortho -- the structure barely budged, so panning looked broken. */
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Right) ||
        ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
      const float winw = std::max(1.0f, (float)app.window_width);
      const float winh = std::max(1.0f, (float)app.window_height);
      float world_per_px_x, world_per_px_y;
      if (app.orthographic_3d) {
        const float hh = 30.0f / std::max(app.zoom, 0.001f);
        world_per_px_x = (2.0f * hh * aspect) / winw;
        world_per_px_y = (2.0f * hh) / winh;
      } else {
        /* perspective: world extent at the camera distance for a 45deg fovy */
        const float hh = std::tan(glm_rad(22.5f)) * app.camera_distance;
        world_per_px_x = (2.0f * hh * aspect) / winw;
        world_per_px_y = (2.0f * hh) / winh;
      }
      app.center_x += io.MouseDelta.x * world_per_px_x;
      app.center_y -= io.MouseDelta.y * world_per_px_y;
    }
    if (io.MouseWheel != 0.0f) {
      app.zoom *= std::pow(1.15f, (float)clamped_wheel(io.MouseWheel));
      app.zoom = std::max(0.05f, std::min(app.zoom, 100.0f));
    }
  }

  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  if (app.bridge_mode == AppState::BridgeMode::CurrentSystem) {
    char hud[420];
    const char *obs = (app.bif_observable[0] != '\0') ? app.bif_observable : "x";
    if (app.bridge_param2[0] != '\0')
      std::snprintf(hud, sizeof(hud),
                    "Bifurcation diagram in 3D — x: %s in [%.3g, %.3g],  height: %s,  z: %s in [%.3g, %.3g]  |  %d pts",
                    app.bif_param, app.bif_start, app.bif_end, obs,
                    app.bridge_param2, app.bridge_p2_min, app.bridge_p2_max, app.bridge_point_count);
    else
      std::snprintf(hud, sizeof(hud),
                    "Bifurcation diagram in 3D — x: %s in [%.3g, %.3g],  height: %s  (set a 2nd param for a stacked surface)  |  %d pts",
                    app.bif_param, app.bif_start, app.bif_end, obs, app.bridge_point_count);
    draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), hud);
    draw->AddText(ImVec2(14, app.window_toolbar_h + 26.0f), IM_COL32(170, 170, 180, 220),
                  "this system has NO unified fractal<->bifurcation bridge (only logistic/cubic/sine/Mandelbrot do), "
                  "so this lifts its own bifurcation diagram into 3D instead.  drag: rotate   wheel: zoom");
    if (app.bridge_point_count == 0)
      draw->AddText(ImVec2(14, 52), IM_COL32(255, 210, 120, 240),
                    "no points — check the x param name/range and observable in the Setup tab, then Rebuild bridge.");
  } else if (app.bridge_mode == AppState::BridgeMode::ProjectionSolid) {
    const bool cubic = (app.bridge_family == AppState::BridgeFamily::Cubic);
    const bool sine  = (app.bridge_family == AppState::BridgeFamily::Sine);
    char hud[384];
    std::snprintf(hud, sizeof(hud),
                  "3D bridge — ONE structure (iterate %s): footprint = %s, real-axis silhouette = %s  |  %d pts",
                  sine ? "c*sin(z)" : (cubic ? "z^3+c" : "z^2+c"),
                  sine ? "sine connectedness set" : (cubic ? "cubic Mandelbrot set" : "Mandelbrot set"),
                  sine ? "sine-map bifurcation" : (cubic ? "cubic-map bifurcation" : "logistic bifurcation"),
                  app.bridge_point_count);
    draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), hud);
    draw->AddText(ImVec2(14, app.window_toolbar_h + 26.0f), IM_COL32(170, 170, 180, 220),
                  "look down the Y axis -> the set; down the Z axis -> the bifurcation diagram. self-contained (independent of the loaded system).  left-drag: rotate   right-drag: move   wheel: zoom");
  } else {
    draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235),
                  "3D bridge (side-by-side) — the Mandelbrot set is the floor, the logistic bifurcation is the vertical face, joined along the real axis (gold spine)");
    draw->AddText(ImVec2(14, app.window_toolbar_h + 26.0f), IM_COL32(170, 170, 180, 220),
                  "ROTATE to read each face: from the side the set, from the front the cascade — they share the c-axis because z^2+c and the logistic map are conjugate.  drag: rotate   wheel: zoom");
  }
}

/* PHASE B+: high-quality bifurcation view.
 * Instead of plotting each point as a faint dot (which looks sparse and
 * loses density information), we accumulate hits into a per-pixel-column
 * density buffer and color by log-density, so dense bands glow and the
 * fine structure of the period-doubling cascade is visible. The view
 * supports pan (left-drag) and zoom (wheel / +/-) over both axes, and
 * draws tick labels. The Lyapunov curve is overlaid with its own scale. */
void render_bifurcation_background(AppState &app) {
  const ImVec2 p0(0.0f, 0.0f);
  const float w = static_cast<float>(app.window_width);
  const float h = static_cast<float>(app.window_height);
  const ImVec2 p1(w, h);
  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  draw->AddRectFilled(p0, p1, IM_COL32(12, 13, 17, 255));

  /* PHASE B+: a small always-visible control strip inside the view, so the
   * bifurcation diagram is one click to (re)compute without hunting in the
   * side panel. This is also what makes "enter the view -> see something"
   * work even before the user finds the Analysis tab. */
  const auto &pts = app.bifurcation_points;
  if (pts.empty()) {
    /* Auto-run once on entry if we have a parameter to sweep — so the view
     * shows the diagram immediately instead of a blank panel the user has
     * to hunt to populate. (Only auto-runs when truly empty.) */
    if (!app.params.empty() && !app.bif_autorun_done) {
      app.bif_autorun_done = true;
      run_bifurcation(app);
    }
    if (app.bifurcation_points.empty()) {
      draw->AddText(ImVec2(p0.x + 18, p0.y + 44), IM_COL32(225, 225, 235, 235),
                    "Bifurcation view. Pick a parameter & range in the top toolbar, then 'Run bifurcation'.");
      draw->AddText(ImVec2(p0.x + 18, p0.y + 64), IM_COL32(190, 190, 200, 220),
                    "For the classic period-doubling tree, load the 'Logistic map' preset, then Run.");
      return;
    }
  }

  const float pad_l = 64.0f, pad_b = 40.0f, pad_t = 36.0f, pad_r = 18.0f;
  const float plot_w = std::max(16.0f, w - pad_l - pad_r);
  const float plot_h = std::max(16.0f, h - pad_b - pad_t);

  /* data extents (full sweep) */
  double dx_min = std::min(app.bif_start, app.bif_end);
  double dx_max = std::max(app.bif_start, app.bif_end);
  if (dx_max <= dx_min) dx_max = dx_min + 1.0;
  double dy_min = DBL_MAX, dy_max = -DBL_MAX;
  for (const auto &p : pts) if (std::isfinite(p.y)) { dy_min = std::min(dy_min, p.y); dy_max = std::max(dy_max, p.y); }
  if (!(dy_max > dy_min)) { dy_min -= 1.0; dy_max += 1.0; }

  /* initialize the interactive view window to the data extents once */
  if (!app.bif_view_valid) {
    const double mx = 0.02 * (dx_max - dx_min), my = 0.04 * (dy_max - dy_min);
    app.bif_view_xmin = dx_min - mx; app.bif_view_xmax = dx_max + mx;
    app.bif_view_ymin = dy_min - my; app.bif_view_ymax = dy_max + my;
    app.bif_view_valid = true;
  }
  double vx0 = app.bif_view_xmin, vx1 = app.bif_view_xmax;
  double vy0 = app.bif_view_ymin, vy1 = app.bif_view_ymax;
  if (vx1 <= vx0) vx1 = vx0 + 1e-9;
  if (vy1 <= vy0) vy1 = vy0 + 1e-9;

  auto sx = [&](double x){ return pad_l + (float)((x - vx0) / (vx1 - vx0)) * plot_w; };
  auto sy = [&](double y){ return pad_t + plot_h - (float)((y - vy0) / (vy1 - vy0)) * plot_h; };
  auto inv_x = [&](float px){ return vx0 + (double)(px - pad_l) / plot_w * (vx1 - vx0); };
  auto inv_y = [&](float py){ return vy0 + (double)(pad_t + plot_h - py) / plot_h * (vy1 - vy0); };

  /* ---- interaction (when not over a panel) ---- */
  ImGuiIO &io = ImGui::GetIO();
  const bool over = !io.WantCaptureMouse;
  if (over) {
    const double mx = inv_x(io.MousePos.x), my = inv_y(io.MousePos.y);
    if (io.MouseWheel != 0.0f) {
      const double f = std::pow(0.85, (double)io.MouseWheel);
      vx0 = mx + (vx0 - mx) * f; vx1 = mx + (vx1 - mx) * f;
      vy0 = my + (vy0 - my) * f; vy1 = my + (vy1 - my) * f;
    }
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) && !ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      const double ddx = (double)io.MouseDelta.x / plot_w * (vx1 - vx0);
      const double ddy = -(double)io.MouseDelta.y / plot_h * (vy1 - vy0);
      vx0 -= ddx; vx1 -= ddx; vy0 += ddy; vy1 += ddy;
    }
    if (ImGui::IsKeyPressed(ImGuiKey_Equal) || ImGui::IsKeyPressed(ImGuiKey_KeypadAdd)) {
      const double cx=0.5*(vx0+vx1), cy=0.5*(vy0+vy1), f=0.83;
      vx0=cx+(vx0-cx)*f; vx1=cx+(vx1-cx)*f; vy0=cy+(vy0-cy)*f; vy1=cy+(vy1-cy)*f;
    }
    if (ImGui::IsKeyPressed(ImGuiKey_Minus) || ImGui::IsKeyPressed(ImGuiKey_KeypadSubtract)) {
      const double cx=0.5*(vx0+vx1), cy=0.5*(vy0+vy1), f=1.20;
      vx0=cx+(vx0-cx)*f; vx1=cx+(vx1-cx)*f; vy0=cy+(vy0-cy)*f; vy1=cy+(vy1-cy)*f;
    }
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) app.bif_view_valid = false;
    app.bif_view_xmin=vx0; app.bif_view_xmax=vx1; app.bif_view_ymin=vy0; app.bif_view_ymax=vy1;
  }

  /* ---- density accumulation into an integer grid ---- */
  const int GW = std::max(1, (int)plot_w);
  const int GH = std::max(1, (int)plot_h);
  static std::vector<float> dens; /* reused buffer */
  dens.assign((size_t)GW * GH, 0.0f);
  float dmax = 0.0f;
  for (const auto &p : pts) {
    if (!std::isfinite(p.x) || !std::isfinite(p.y)) continue;
    if (p.x < vx0 || p.x > vx1 || p.y < vy0 || p.y > vy1) continue;
    const int gx = (int)((p.x - vx0) / (vx1 - vx0) * (GW - 1) + 0.5);
    const int gy = (int)((p.y - vy0) / (vy1 - vy0) * (GH - 1) + 0.5);
    if (gx < 0 || gx >= GW || gy < 0 || gy >= GH) continue;
    float &c = dens[(size_t)gy * GW + gx];
    c += 1.0f;
    if (c > dmax) dmax = c;
  }

  /* draw the density as short vertical runs per (gx,gy) cell, colored by
   * log-density. A gamma lift makes the thin low-density branches (the
   * single-period curve, the early forks) clearly visible instead of
   * near-black — previously they could read as faint scattered dots. */
  if (dmax > 0.0f) {
    const float inv_log = 1.0f / std::log(1.0f + dmax);
    for (int gy = 0; gy < GH; ++gy) {
      for (int gx = 0; gx < GW; ++gx) {
        const float c = dens[(size_t)gy * GW + gx];
        if (c <= 0.0f) continue;
        float t = std::log(1.0f + c) * inv_log; /* 0..1 */
        t = std::pow(t, 0.45f);                  /* gamma lift for faint branches */
        /* ramp with a bright floor: even a single hit is clearly visible
         * (light blue), dense bands go cyan->white. */
        int r, g, b;
        if (t < 0.5f) { const float u = t / 0.5f; r = (int)(70 + 30*u); g = (int)(150 + 90*u); b = (int)(210 + 45*u); }
        else { const float u = (t - 0.5f) / 0.5f; r = (int)(100 + 155*u); g = (int)(240 + 15*u); b = (int)(255); }
        const float px = pad_l + gx;
        const float py = pad_t + gy;
        draw->AddRectFilled(ImVec2(px, py), ImVec2(px + 1.6f, py + 1.6f),
                            IM_COL32(r, g, b, 255));
      }
    }
  }

  /* ---- axes + tick labels ---- */
  const ImU32 axcol = IM_COL32(140, 140, 150, 200);
  draw->AddRect(ImVec2(pad_l, pad_t), ImVec2(pad_l + plot_w, pad_t + plot_h), axcol);
  auto nice_step = [](double range, int target){
    double raw = range / std::max(1, target);
    double mag = std::pow(10.0, std::floor(std::log10(raw)));
    double n = raw / mag;
    double s = n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10;
    return s * mag;
  };
  const double xs = nice_step(vx1 - vx0, 8);
  for (double xv = std::ceil(vx0 / xs) * xs; xv <= vx1; xv += xs) {
    const float X = sx(xv);
    draw->AddLine(ImVec2(X, pad_t + plot_h), ImVec2(X, pad_t + plot_h + 4), axcol);
    char b[32]; std::snprintf(b, sizeof(b), "%.4g", xv);
    draw->AddText(ImVec2(X - 12, pad_t + plot_h + 6), IM_COL32(190,190,200,220), b);
  }
  const double ys = nice_step(vy1 - vy0, 6);
  for (double yv = std::ceil(vy0 / ys) * ys; yv <= vy1; yv += ys) {
    const float Y = sy(yv);
    draw->AddLine(ImVec2(pad_l - 4, Y), ImVec2(pad_l, Y), axcol);
    char b[32]; std::snprintf(b, sizeof(b), "%.4g", yv);
    draw->AddText(ImVec2(6, Y - 7), IM_COL32(190,190,200,220), b);
  }

  /* ---- Lyapunov overlay (own scale, zero line) ---- */
  if (app.bif_compute_lyapunov && !app.bifurcation_lyapunov.empty()) {
    double lo=DBL_MAX, hi=-DBL_MAX;
    for (const auto &p : app.bifurcation_lyapunov) if (std::isfinite(p.y)) { lo=std::min(lo,p.y); hi=std::max(hi,p.y); }
    if (std::isfinite(lo) && hi>lo) {
      hi=std::max(hi,0.05); lo=std::min(lo,-0.05);
      auto lyy=[&](double y){ return pad_t + plot_h - (float)((y-lo)/(hi-lo))*plot_h; };
      const float zl = lyy(0.0);
      draw->AddLine(ImVec2(pad_l, zl), ImVec2(pad_l + plot_w, zl), IM_COL32(160,160,160,110), 1.0f);
      bool have=false; ImVec2 prev;
      for (const auto &pt : app.bifurcation_lyapunov) {
        if (!std::isfinite(pt.y)) { have=false; continue; }
        if (pt.x < vx0 || pt.x > vx1) { have=false; continue; }
        const ImVec2 cur(sx(pt.x), lyy(pt.y));
        const ImU32 col = pt.y > 0 ? IM_COL32(255,110,90,235) : IM_COL32(120,230,140,235);
        if (have) draw->AddLine(prev, cur, col, 1.6f);
        prev = cur; have = true;
      }
    }
  }

  /* ---- period strip (maps): a thin band at the top encoding the detected
   * period at each parameter value. Makes the period-doubling cascade and
   * the periodic windows (e.g. the period-3 window) legible at a glance. */
  if (app.mode == SystemMode::Map && app.bif_show_period && !app.bifurcation_period.empty()) {
    const float band_y0 = pad_t + 2.0f, band_h = 7.0f;
    auto period_color = [](int per) -> ImU32 {
      if (per <= 0) return IM_COL32(220, 90, 70, 255);         /* chaotic: red */
      switch (per) {
        case 1: return IM_COL32(60, 120, 220, 255);            /* blue */
        case 2: return IM_COL32(60, 200, 200, 255);            /* cyan */
        case 3: return IM_COL32(255, 215, 80, 255);            /* gold (period-3!) */
        case 4: return IM_COL32(120, 220, 120, 255);           /* green */
        case 8: return IM_COL32(200, 160, 240, 255);           /* violet */
        default: return IM_COL32(180, 180, 190, 255);          /* other periodic: grey */
      }
    };
    for (const auto &pp : app.bifurcation_period) {
      if (pp.x < vx0 || pp.x > vx1) continue;
      const float X = sx(pp.x);
      draw->AddRectFilled(ImVec2(X, band_y0), ImVec2(X + 2.0f, band_y0 + band_h), period_color((int)pp.y));
    }
    draw->AddText(ImVec2(pad_l + plot_w - 220, band_y0 + band_h + 1),
                  IM_COL32(170,170,180,210), "period band: blue=1 cyan=2 gold=3 green=4 red=chaos");
  }

  /* ---- labels ---- */
  char title[360];
  std::snprintf(title, sizeof(title), "Bifurcation diagram   %s (vertical)  vs  %s (horizontal)   |  %zu points",
                app.bif_observable, app.bif_param, pts.size());
  draw->AddText(ImVec2(pad_l, 10), IM_COL32(235, 235, 240, 235), title);
  if (over)
    draw->AddText(ImVec2(pad_l, pad_t + plot_h + 22),
                  IM_COL32(150,150,160,200),
                  "drag: pan   wheel/+/-: zoom   double-click: reset   (blue->white = density;  red Lyapunov>0 = chaos)");
}

/* ============================================================
 * PHASE D: equilibrium continuation diagram (MatCont-style).
 * x-axis = continuation parameter, y-axis = a chosen state component.
 * Stable branch solid & bright, unstable branch dashed & dim. Fold (o)
 * and Hopf (<>) points marked and labeled.
 * ============================================================ */
void render_continuation_background(AppState &app) {
  const ImVec2 p0(0.0f, 0.0f);
  const float w = (float)app.window_width, h = (float)app.window_height;
  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  draw->AddRectFilled(p0, ImVec2(w, h), IM_COL32(12, 13, 17, 255));

  if (app.mode != SystemMode::ODE) {
    draw->AddText(ImVec2(20, 48), IM_COL32(225, 225, 235, 235),
                  "Continuation is for ODE systems (equilibria vs. a parameter).");
    return;
  }
  if (app.params.empty()) {
    draw->AddText(ImVec2(20, 48), IM_COL32(225, 225, 235, 235),
                  "This system has no parameter to continue.");
    return;
  }
  if (!app.cont_has_branch) {
    /* auto-run once on entry so the diagram appears without hunting */
    if (!app.params.empty() && !app.cont_autorun_done) {
      app.cont_autorun_done = true;
      run_continuation(app);
    }
    if (!app.cont_has_branch) {
      draw->AddText(ImVec2(20, 48), IM_COL32(225, 225, 235, 235),
                    "Equilibrium continuation. Pick a parameter & range in the top toolbar, then 'Continue'.");
      draw->AddText(ImVec2(20, 68), IM_COL32(185, 185, 195, 220),
                    "It traces the equilibrium both ways from the current state and marks fold/Hopf points.");
      if (!app.cont_message.empty())
        draw->AddText(ImVec2(20, 92), IM_COL32(255, 200, 120, 230), app.cont_message.c_str());
      return;
    }
  }

  const float pad_l = 64.0f, pad_b = 40.0f, pad_t = 36.0f, pad_r = 18.0f;
  const float plot_w = std::max(16.0f, w - pad_l - pad_r);
  const float plot_h = std::max(16.0f, h - pad_t - pad_b);

  /* fit view to the branch on (re)compute */
  if (!app.cont_view_valid) {
    double xlo = 1e300, xhi = -1e300, ylo = 1e300, yhi = -1e300;
    for (size_t i = 0; i < app.cont_pp.size(); ++i) {
      xlo = std::min(xlo, app.cont_pp[i]); xhi = std::max(xhi, app.cont_pp[i]);
      ylo = std::min(ylo, app.cont_yy[i]); yhi = std::max(yhi, app.cont_yy[i]);
    }
    if (xlo > xhi) { xlo = -1; xhi = 1; }
    if (ylo > yhi) { ylo = -1; yhi = 1; }
    const double mx = 0.06 * (xhi - xlo) + 1e-9, my = 0.08 * (yhi - ylo) + 1e-9;
    app.cont_view_xmin = xlo - mx; app.cont_view_xmax = xhi + mx;
    app.cont_view_ymin = ylo - my; app.cont_view_ymax = yhi + my;
    app.cont_view_valid = true;
  }
  double vx0 = app.cont_view_xmin, vx1 = app.cont_view_xmax;
  double vy0 = app.cont_view_ymin, vy1 = app.cont_view_ymax;
  if (vx1 <= vx0) vx1 = vx0 + 1e-9;
  if (vy1 <= vy0) vy1 = vy0 + 1e-9;

  auto sx = [&](double x){ return pad_l + (float)((x - vx0) / (vx1 - vx0)) * plot_w; };
  auto sy = [&](double y){ return pad_t + plot_h - (float)((y - vy0) / (vy1 - vy0)) * plot_h; };
  auto inv_x = [&](float px){ return vx0 + (double)(px - pad_l) / plot_w * (vx1 - vx0); };
  auto inv_y = [&](float py){ return vy0 + (double)(pad_t + plot_h - py) / plot_h * (vy1 - vy0); };

  ImGuiIO &io = ImGui::GetIO();
  const bool over = !io.WantCaptureMouse;
  if (over) {
    const double mx = inv_x(io.MousePos.x), my = inv_y(io.MousePos.y);
    if (io.MouseWheel != 0.0f) {
      const double f = std::pow(0.85, (double)io.MouseWheel);
      vx0 = mx + (vx0 - mx) * f; vx1 = mx + (vx1 - mx) * f;
      vy0 = my + (vy0 - my) * f; vy1 = my + (vy1 - my) * f;
    }
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) && !ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      const double ddx = (double)io.MouseDelta.x / plot_w * (vx1 - vx0);
      const double ddy = -(double)io.MouseDelta.y / plot_h * (vy1 - vy0);
      vx0 -= ddx; vx1 -= ddx; vy0 += ddy; vy1 += ddy;
    }
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) app.cont_view_valid = false;
    app.cont_view_xmin=vx0; app.cont_view_xmax=vx1; app.cont_view_ymin=vy0; app.cont_view_ymax=vy1;
  }

  /* axes box */
  const ImU32 axcol = IM_COL32(140, 140, 150, 200);
  draw->AddRect(ImVec2(pad_l, pad_t), ImVec2(pad_l + plot_w, pad_t + plot_h), axcol);
  auto nice_step = [](double range, int target){
    double raw = range / std::max(1, target);
    double mag = std::pow(10.0, std::floor(std::log10(raw)));
    double nn = raw / mag; double st = nn < 1.5 ? 1 : nn < 3 ? 2 : nn < 7 ? 5 : 10;
    return st * mag;
  };
  const double xs = nice_step(vx1 - vx0, 8);
  for (double xv = std::ceil(vx0 / xs) * xs; xv <= vx1; xv += xs) {
    const float X = sx(xv);
    draw->AddLine(ImVec2(X, pad_t + plot_h), ImVec2(X, pad_t + plot_h + 4), axcol);
    char bb[32]; std::snprintf(bb, sizeof(bb), "%.4g", xv);
    draw->AddText(ImVec2(X - 12, pad_t + plot_h + 6), IM_COL32(190,190,200,220), bb);
  }
  const double ys = nice_step(vy1 - vy0, 6);
  for (double yv = std::ceil(vy0 / ys) * ys; yv <= vy1; yv += ys) {
    const float Y = sy(yv);
    draw->AddLine(ImVec2(pad_l - 4, Y), ImVec2(pad_l, Y), axcol);
    char bb[32]; std::snprintf(bb, sizeof(bb), "%.4g", yv);
    draw->AddText(ImVec2(8, Y - 7), IM_COL32(190,190,200,220), bb);
  }

  /* draw the branch as connected segments. A segment is drawn only when
   * consecutive points are not separated by a polyline break; stable
   * segments are solid bright green, unstable are dim red and dashed. */
  auto is_break = [&](size_t i){
    for (size_t b : app.cont_break) if (b == i) return true;
    return false;
  };
  const ImU32 col_stable = IM_COL32(90, 230, 120, 255);
  const ImU32 col_unstable = IM_COL32(225, 110, 90, 200);
  for (size_t i = 1; i < app.cont_pp.size(); ++i) {
    if (is_break(i)) continue;
    const ImVec2 a(sx(app.cont_pp[i-1]), sy(app.cont_yy[i-1]));
    const ImVec2 b(sx(app.cont_pp[i]),   sy(app.cont_yy[i]));
    const bool stable = app.cont_stable[i] && app.cont_stable[i-1];
    if (stable) {
      draw->AddLine(a, b, col_stable, 2.0f);
    } else {
      /* dashed: draw short pieces */
      const float dx = b.x - a.x, dy = b.y - a.y;
      const float len = std::sqrt(dx*dx + dy*dy);
      const int seg = std::max(1, (int)(len / 6.0f));
      for (int k = 0; k < seg; k += 2) {
        const float t0 = (float)k / seg, t1 = std::min(1.0f, (float)(k+1) / seg);
        draw->AddLine(ImVec2(a.x + dx*t0, a.y + dy*t0), ImVec2(a.x + dx*t1, a.y + dy*t1), col_unstable, 1.6f);
      }
    }
  }

  /* special-point markers: fold (circle), Hopf (diamond), and codim-2 / branch
   * points (squares, distinctly colored and labelled). */
  for (size_t i = 0; i < app.cont_special.size(); ++i) {
    const int sp = app.cont_special[i];
    if (sp == 0) continue;
    const ImVec2 c(sx(app.cont_pp[i]), sy(app.cont_yy[i]));
    const float r = 6.0f;
    if (sp == 1) {
      draw->AddCircle(c, r, IM_COL32(255, 220, 90, 255), 16, 2.0f);
      draw->AddText(ImVec2(c.x + 8, c.y - 6), IM_COL32(255, 220, 90, 255), "Fold (LP)");
    } else if (sp == 2) {
      draw->AddQuad(ImVec2(c.x, c.y-r), ImVec2(c.x+r, c.y), ImVec2(c.x, c.y+r), ImVec2(c.x-r, c.y),
                    IM_COL32(120, 200, 255, 255), 2.0f);
      draw->AddText(ImVec2(c.x + 8, c.y - 6), IM_COL32(120, 200, 255, 255), "Hopf (H)");
    } else {
      /* codim-2 / branch point: filled square + label */
      ImU32 col; const char *lab;
      switch (sp) {
        case 3: col = IM_COL32(255, 120, 120, 255); lab = "BT (Bogdanov-Takens)"; break;
        case 4: col = IM_COL32(255, 160, 80, 255);  lab = "CP (cusp)";            break;
        case 5: col = IM_COL32(180, 140, 255, 255); lab = "GH (gen. Hopf)";       break;
        case 6: col = IM_COL32(120, 255, 160, 255); lab = "BP (branch point)";    break;
        default: col = IM_COL32(230,230,230,255);   lab = "special";              break;
      }
      draw->AddRectFilled(ImVec2(c.x-r, c.y-r), ImVec2(c.x+r, c.y+r), col);
      draw->AddRect(ImVec2(c.x-r, c.y-r), ImVec2(c.x+r, c.y+r), IM_COL32(20,20,25,255), 0, 0, 1.5f);
      draw->AddText(ImVec2(c.x + 8, c.y - 6), col, lab);
    }
  }

  /* axis labels + legend */
  const int yi = std::max(0, std::min(app.cont_y_index, (int)app.state_names.size() - 1));
  char xlabel[96]; std::snprintf(xlabel, sizeof(xlabel), "parameter %s", app.cont_param);
  draw->AddText(ImVec2(pad_l + plot_w * 0.5f - 40, pad_t + plot_h + 22), IM_COL32(210,210,220,230), xlabel);
  draw->AddText(ImVec2(12, pad_t - 22), IM_COL32(210,210,220,230), app.state_names[yi].c_str());

  char hud[320];
  std::snprintf(hud, sizeof(hud),
                "Equilibrium continuation — %zu points  |  %d fold, %d Hopf  |  solid green = stable, dashed red = unstable",
                app.cont_pp.size(), app.cont_n_fold, app.cont_n_hopf);
  draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), hud);
}

/* ============================================================
 * PHASE D step 2: limit-cycle continuation diagram.
 * x-axis = parameter; left axis (green) = oscillation amplitude;
 * right axis (orange) = period. Amplitude rising from zero marks a Hopf
 * bifurcation; gaps mean "no cycle here" (settled to a fixed point).
 * ============================================================ */
void render_lc_continuation_background(AppState &app) {
  const float w = (float)app.window_width, h = (float)app.window_height;
  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  draw->AddRectFilled(ImVec2(0, 0), ImVec2(w, h), IM_COL32(12, 13, 17, 255));

  if (app.mode != SystemMode::ODE) {
    draw->AddText(ImVec2(20, 48), IM_COL32(225, 225, 235, 235),
                  "Limit-cycle continuation is for ODE systems.");
    return;
  }
  if (app.params.empty()) {
    draw->AddText(ImVec2(20, 48), IM_COL32(225, 225, 235, 235),
                  "This system has no parameter to sweep.");
    return;
  }
  if (!app.lcc_has_data) {
    if (!app.lcc_autorun_done) {
      app.lcc_autorun_done = true;
      if (app.lcc_param[0] == '\0' && !app.params.empty())
        std::snprintf(app.lcc_param, sizeof(app.lcc_param), "%s", app.params[0].name.c_str());
      run_lc_continuation(app);
      /* Also run the rigorous collocation + pseudo-arclength continuation in the
       * SAME parameter, so the main view shows the branch that follows folds of
       * cycles with Floquet stability -- not just the simulate-and-measure sweep
       * (which only finds stable cycles). */
      std::snprintf(app.cont_param, sizeof(app.cont_param), "%s", app.lcc_param);
      run_cycle_collocation(app);
    }
    if (!app.lcc_has_data) {
      draw->AddText(ImVec2(20, 48), IM_COL32(225, 225, 235, 235),
                    "Limit-cycle continuation. Pick a parameter & range in the toolbar, then 'Sweep'.");
      draw->AddText(ImVec2(20, 68), IM_COL32(185, 185, 195, 220),
                    "Integrates at each parameter, measures the cycle's period & amplitude.");
      if (!app.lcc_msg.empty() && app.lcc_msg != "ok")
        draw->AddText(ImVec2(20, 92), IM_COL32(255, 200, 120, 230), app.lcc_msg.c_str());
      return;
    }
  }

  const float pad_l = 64.0f, pad_r = 64.0f, pad_b = 42.0f, pad_t = 38.0f;
  const float plot_w = std::max(16.0f, w - pad_l - pad_r);
  const float plot_h = std::max(16.0f, h - pad_t - pad_b);

  double pmin = app.lcc_pp.front(), pmax = app.lcc_pp.back();
  if (pmax <= pmin) pmax = pmin + 1e-9;
  double amp_hi = app.lcc_amp_max * 1.08 + 1e-9;
  double per_hi = app.lcc_period_max * 1.08 + 1e-9;
  /* widen the ranges so the overlaid collocation branch fits on the same axes */
  const bool have_coll = app.cyc_ready && app.cyc_branch.samples.size() > 1;
  if (have_coll) {
    for (const auto &s : app.cyc_branch.samples) {
      pmin = std::min(pmin, s.p); pmax = std::max(pmax, s.p);
      amp_hi = std::max(amp_hi, s.amplitude * 1.08);
      per_hi = std::max(per_hi, s.period * 1.08);
    }
    if (pmax <= pmin) pmax = pmin + 1e-9;
  }

  auto sx = [&](double p){ return pad_l + (float)((p - pmin) / (pmax - pmin)) * plot_w; };
  auto sy_amp = [&](double a){ return pad_t + plot_h - (float)((a / amp_hi)) * plot_h; };
  auto sy_per = [&](double t){ return pad_t + plot_h - (float)((t / per_hi)) * plot_h; };

  /* axes box */
  const ImU32 axcol = IM_COL32(140, 140, 150, 200);
  draw->AddRect(ImVec2(pad_l, pad_t), ImVec2(pad_l + plot_w, pad_t + plot_h), axcol);

  auto nice_step = [](double range, int target){
    double raw = range / std::max(1, target);
    double mag = std::pow(10.0, std::floor(std::log10(raw)));
    double nn = raw / mag; double st = nn < 1.5 ? 1 : nn < 3 ? 2 : nn < 7 ? 5 : 10;
    return st * mag;
  };
  /* x ticks (parameter) */
  const double xs = nice_step(pmax - pmin, 8);
  for (double xv = std::ceil(pmin / xs) * xs; xv <= pmax; xv += xs) {
    const float X = sx(xv);
    draw->AddLine(ImVec2(X, pad_t + plot_h), ImVec2(X, pad_t + plot_h + 4), axcol);
    char b[32]; std::snprintf(b, sizeof(b), "%.4g", xv);
    draw->AddText(ImVec2(X - 12, pad_t + plot_h + 6), IM_COL32(190,190,200,220), b);
  }
  /* left y ticks (amplitude, green) */
  const double as = nice_step(amp_hi, 5);
  for (double av = 0; av <= amp_hi; av += as) {
    const float Y = sy_amp(av);
    draw->AddLine(ImVec2(pad_l - 4, Y), ImVec2(pad_l, Y), IM_COL32(90,210,120,200));
    char b[32]; std::snprintf(b, sizeof(b), "%.3g", av);
    draw->AddText(ImVec2(6, Y - 7), IM_COL32(120,220,150,220), b);
  }
  /* right y ticks (period, orange) */
  const double ps = nice_step(per_hi, 5);
  for (double pv = 0; pv <= per_hi; pv += ps) {
    const float Y = sy_per(pv);
    draw->AddLine(ImVec2(pad_l + plot_w, Y), ImVec2(pad_l + plot_w + 4, Y), IM_COL32(240,170,80,200));
    char b[32]; std::snprintf(b, sizeof(b), "%.3g", pv);
    draw->AddText(ImVec2(pad_l + plot_w + 8, Y - 7), IM_COL32(245,180,90,220), b);
  }

  /* amplitude curve (green): connect consecutive points that both have a
   * cycle; leave gaps where there's no oscillation. Mark cycle points. */
  const ImU32 col_amp = IM_COL32(90, 230, 120, 255);
  for (size_t i = 1; i < app.lcc_pp.size(); ++i) {
    const double a0 = app.lcc_amp[i-1], a1 = app.lcc_amp[i];
    const bool c0 = std::isfinite(app.lcc_period[i-1]) && a0 > 0;
    const bool c1 = std::isfinite(app.lcc_period[i]) && a1 > 0;
    if (c0 && c1)
      draw->AddLine(ImVec2(sx(app.lcc_pp[i-1]), sy_amp(a0)), ImVec2(sx(app.lcc_pp[i]), sy_amp(a1)), col_amp, 2.0f);
  }
  /* period curve (orange) */
  const ImU32 col_per = IM_COL32(240, 170, 80, 255);
  for (size_t i = 1; i < app.lcc_pp.size(); ++i) {
    const double t0 = app.lcc_period[i-1], t1 = app.lcc_period[i];
    if (std::isfinite(t0) && std::isfinite(t1))
      draw->AddLine(ImVec2(sx(app.lcc_pp[i-1]), sy_per(t0)), ImVec2(sx(app.lcc_pp[i]), sy_per(t1)), col_per, 1.8f);
  }
  /* dots on cycle points (amplitude) */
  int n_cycle = 0;
  for (size_t i = 0; i < app.lcc_pp.size(); ++i) {
    if (std::isfinite(app.lcc_period[i]) && app.lcc_amp[i] > 0) {
      draw->AddCircleFilled(ImVec2(sx(app.lcc_pp[i]), sy_amp(app.lcc_amp[i])), 2.2f, col_amp);
      ++n_cycle;
    }
  }

  /* OVERLAY: the rigorous collocation + pseudo-arclength branch. Drawn as a
   * brighter amplitude line coloured by Floquet stability (green stable / red
   * unstable -- the sweep can't show unstable cycles at all), with the period
   * as a brighter orange, plus LPC / PD / NS markers. This is the branch that
   * follows folds of cycles. */
  if (have_coll) {
    const auto &S = app.cyc_branch.samples;
    for (size_t i = 1; i < S.size(); ++i) {
      draw->AddLine(ImVec2(sx(S[i-1].p), sy_per(S[i-1].period)),
                    ImVec2(sx(S[i].p),   sy_per(S[i].period)),
                    IM_COL32(255, 200, 110, 200), 1.4f);
      const ImU32 col = S[i].stable ? IM_COL32(150, 255, 165, 255) : IM_COL32(255, 120, 95, 255);
      draw->AddLine(ImVec2(sx(S[i-1].p), sy_amp(S[i-1].amplitude)),
                    ImVec2(sx(S[i].p),   sy_amp(S[i].amplitude)), col, 2.6f);
    }
    for (const auto &s : S) {
      const ImVec2 c(sx(s.p), sy_amp(s.amplitude));
      if (s.is_fold) { draw->AddCircleFilled(c, 5.0f, IM_COL32(255,210,80,255)); draw->AddText(ImVec2(c.x+6,c.y-6), IM_COL32(255,210,80,255), "LPC"); }
      if (s.is_pd)   { draw->AddCircleFilled(c, 4.5f, IM_COL32(255,120,200,255)); draw->AddText(ImVec2(c.x+6,c.y+4), IM_COL32(255,120,200,255), "PD"); }
      if (s.is_ns)   { draw->AddCircleFilled(c, 4.5f, IM_COL32(150,200,255,255)); draw->AddText(ImVec2(c.x+6,c.y+4), IM_COL32(150,200,255,255), "NS"); }
    }
  }

  /* labels */
  char xlabel[96]; std::snprintf(xlabel, sizeof(xlabel), "parameter %s", app.lcc_param);
  draw->AddText(ImVec2(pad_l + plot_w * 0.5f - 40, pad_t + plot_h + 22), IM_COL32(210,210,220,230), xlabel);
  draw->AddText(ImVec2(8, pad_t - 22), IM_COL32(120,220,150,235), "amplitude");
  draw->AddText(ImVec2(pad_l + plot_w - 28, pad_t - 22), IM_COL32(245,180,90,235), "period");

  char hud[320];
  if (have_coll)
    std::snprintf(hud, sizeof(hud),
                  "Limit-cycle continuation — sweep: %zu slices (%d cyclic, faint dots) + collocation branch (bold: green=stable, red=unstable; LPC/PD/NS marked)%s",
                  app.lcc_pp.size(), n_cycle, app.cyc_branch.turned ? "  [turned a fold]" : "");
  else
    std::snprintf(hud, sizeof(hud),
                  "Limit-cycle continuation — %zu slices, %d with a cycle  |  green = amplitude (left), orange = period (right)",
                  app.lcc_pp.size(), n_cycle);
  draw->AddText(ImVec2(14, app.window_toolbar_h + 8.0f), IM_COL32(235, 235, 240, 235), hud);
}

void draw_series_plot(AppState &app, const char *title, const char *value_name, ImVec2 size) {
  ImGui::TextUnformatted(title);
  ImDrawList *draw = ImGui::GetWindowDrawList();
  ImVec2 p0 = ImGui::GetCursorScreenPos();
  ImVec2 p1(p0.x + size.x, p0.y + size.y);
  draw->AddRectFilled(p0, p1, IM_COL32(20, 20, 24, 220));
  draw->AddRect(p0, p1, IM_COL32(120, 120, 120, 200));
  if (app.history.size() > 1) {
    double mn = DBL_MAX, mx = -DBL_MAX;
    std::vector<double> vals;
    vals.reserve(app.history.size());
    for (const State &s : app.history) {
      double v = 0.0;
      if (value_by_name(app, s, value_name, &v) && std::isfinite(v)) {
        vals.push_back(v);
        mn = std::min(mn, v);
        mx = std::max(mx, v);
      } else {
        vals.push_back(0.0);
      }
    }
    if (mx <= mn) { mx = mn + 1.0; }
    for (size_t i = 1; i < vals.size(); ++i) {
      const float x0 = p0.x + static_cast<float>(i - 1) / static_cast<float>(vals.size() - 1) * size.x;
      const float x1 = p0.x + static_cast<float>(i) / static_cast<float>(vals.size() - 1) * size.x;
      const float y0 = p1.y - static_cast<float>((vals[i - 1] - mn) / (mx - mn)) * size.y;
      const float y1 = p1.y - static_cast<float>((vals[i] - mn) / (mx - mn)) * size.y;
      draw->AddLine(ImVec2(x0, y0), ImVec2(x1, y1), IM_COL32(120, 210, 255, 255), 1.0f);
    }
    char label[128];
    std::snprintf(label, sizeof(label), "%.4g .. %.4g", mn, mx);
    draw->AddText(ImVec2(p0.x + 4, p0.y + 4), IM_COL32(180, 180, 180, 220), label);
  }
  ImGui::Dummy(size);
}

void draw_scatter_plot(const char *title, const std::vector<Point2> &points, ImVec2 size) {
  ImGui::Text("%s (%zu points)", title, points.size());
  ImDrawList *draw = ImGui::GetWindowDrawList();
  ImVec2 p0 = ImGui::GetCursorScreenPos();
  ImVec2 p1(p0.x + size.x, p0.y + size.y);
  draw->AddRectFilled(p0, p1, IM_COL32(20, 20, 24, 220));
  draw->AddRect(p0, p1, IM_COL32(120, 120, 120, 200));
  if (!points.empty()) {
    double min_x = DBL_MAX, max_x = -DBL_MAX, min_y = DBL_MAX, max_y = -DBL_MAX;
    for (const auto &p : points) {
      if (!std::isfinite(p.x) || !std::isfinite(p.y)) continue;
      min_x = std::min(min_x, p.x); max_x = std::max(max_x, p.x);
      min_y = std::min(min_y, p.y); max_y = std::max(max_y, p.y);
    }
    if (max_x <= min_x) max_x = min_x + 1.0;
    if (max_y <= min_y) max_y = min_y + 1.0;
    for (const auto &p : points) {
      const float sx = p0.x + static_cast<float>((p.x - min_x) / (max_x - min_x)) * size.x;
      const float sy = p1.y - static_cast<float>((p.y - min_y) / (max_y - min_y)) * size.y;
      draw->AddCircleFilled(ImVec2(sx, sy), 1.2f, IM_COL32(255, 190, 80, 230));
    }
  }
  ImGui::Dummy(size);
}

/* PHASE7: render the bifurcation (orbit) diagram, optionally with the
 * Lyapunov-exponent curve overlaid on the same parameter axis. */
[[maybe_unused]] void draw_bifurcation_diagram(AppState &app, ImVec2 size) {
  ImGui::Text("bifurcation: %s vs %s  (%zu points)", app.bif_observable,
              app.bif_param, app.bifurcation_points.size());
  ImDrawList *draw = ImGui::GetWindowDrawList();
  ImVec2 p0 = ImGui::GetCursorScreenPos();
  ImVec2 p1(p0.x + size.x, p0.y + size.y);
  draw->AddRectFilled(p0, p1, IM_COL32(16, 17, 22, 255));
  draw->AddRect(p0, p1, IM_COL32(120, 120, 120, 200));

  const auto &pts = app.bifurcation_points;
  if (pts.empty()) {
    draw->AddText(ImVec2(p0.x + 10, p0.y + 10), IM_COL32(200, 200, 210, 220),
                  "press 'Run bifurcation' to compute");
    ImGui::Dummy(size);
    return;
  }

  double px_min = std::min(app.bif_start, app.bif_end);
  double px_max = std::max(app.bif_start, app.bif_end);
  if (px_max <= px_min) px_max = px_min + 1.0;
  double py_min = DBL_MAX, py_max = -DBL_MAX;
  for (const auto &p : pts) {
    if (!std::isfinite(p.y)) continue;
    py_min = std::min(py_min, p.y);
    py_max = std::max(py_max, p.y);
  }
  if (!(py_max > py_min)) { py_min -= 1.0; py_max += 1.0; }

  auto sx = [&](double x) {
    return p0.x + static_cast<float>((x - px_min) / (px_max - px_min)) * size.x;
  };
  auto sy = [&](double y) {
    return p1.y - static_cast<float>((y - py_min) / (py_max - py_min)) * size.y;
  };

  for (const auto &p : pts) {
    if (!std::isfinite(p.x) || !std::isfinite(p.y)) continue;
    draw->AddRectFilled(ImVec2(sx(p.x), sy(p.y)),
                        ImVec2(sx(p.x) + 1.0f, sy(p.y) + 1.0f),
                        IM_COL32(120, 200, 255, 90));
  }

  if (app.bif_compute_lyapunov && !app.bifurcation_lyapunov.empty()) {
    double ly_min = DBL_MAX, ly_max = -DBL_MAX;
    for (const auto &p : app.bifurcation_lyapunov) {
      if (!std::isfinite(p.y)) continue;
      ly_min = std::min(ly_min, p.y);
      ly_max = std::max(ly_max, p.y);
    }
    if (std::isfinite(ly_min) && ly_max > ly_min) {
      ly_max = std::max(ly_max, 0.05);
      ly_min = std::min(ly_min, -0.05);
      auto lyy = [&](double y) {
        return p1.y - static_cast<float>((y - ly_min) / (ly_max - ly_min)) * size.y;
      };
      const float zline = lyy(0.0);
      draw->AddLine(ImVec2(p0.x, zline), ImVec2(p1.x, zline),
                    IM_COL32(150, 150, 150, 120), 1.0f);
      for (size_t i = 1; i < app.bifurcation_lyapunov.size(); ++i) {
        const auto &a = app.bifurcation_lyapunov[i - 1];
        const auto &b = app.bifurcation_lyapunov[i];
        if (!std::isfinite(a.y) || !std::isfinite(b.y)) continue;
        const ImU32 col = b.y > 0 ? IM_COL32(255, 110, 90, 230)
                                  : IM_COL32(120, 230, 140, 230);
        draw->AddLine(ImVec2(sx(a.x), lyy(a.y)), ImVec2(sx(b.x), lyy(b.y)), col, 1.6f);
      }
      draw->AddText(ImVec2(p0.x + 8, p0.y + 8), IM_COL32(235, 235, 240, 230),
                    "red Lyapunov>0 (chaos), green<0 (stable); blue: attractor");
    }
  }
  ImGui::Dummy(size);
}

void export_trajectory_csv(AppState &app) {
  FILE *f = std::fopen(app.export_path, "w");
  if (!f) {
    app.analysis_message = std::string("could not open export path: ") + app.export_path;
    return;
  }
  std::fprintf(f, "t");
  for (const std::string &name : app.state_names) {
    std::fprintf(f, ",%s", name.c_str());
  }
  for (const auto &obs : app.observables) {
    std::fprintf(f, ",%s", obs.name.c_str());
  }
  std::fprintf(f, "\n");
  char err[256] = {0};
  for (const State &s : app.history) {
    std::fprintf(f, "%.17g", s.t);
    for (size_t i = 0; i < app.state_names.size(); ++i) {
      std::fprintf(f, ",%.17g", state_at(s, i));
    }
    for (const auto &obs : app.observables) {
      double value = 0.0;
      if (!eval_expr_at(app, obs.body, s, &value, err, sizeof(err))) value = NAN;
      std::fprintf(f, ",%.17g", value);
    }
    std::fprintf(f, "\n");
  }
  std::fclose(f);
  app.analysis_message = std::string("exported trajectory CSV: ") + app.export_path;
}

void export_poincare_csv(AppState &app) {
  FILE *f = std::fopen(app.export_path, "w");
  if (!f) {
    app.analysis_message = std::string("could not open export path: ") + app.export_path;
    return;
  }
  std::fprintf(f, "u,v\n");
  for (const auto &p : app.poincare_points) {
    std::fprintf(f, "%.17g,%.17g\n", p.x, p.y);
  }
  std::fclose(f);
  app.analysis_message = std::string("exported Poincare CSV: ") + app.export_path;
}

bool solve_linear_system(std::vector<double> A, std::vector<double> b, std::vector<double> *x) {
  const size_t n = b.size();
  if (A.size() != n * n) return false;
  std::vector<double> m(n * (n + 1), 0.0);
  for (size_t r = 0; r < n; ++r) {
    for (size_t c = 0; c < n; ++c) m[r * (n + 1) + c] = A[r * n + c];
    m[r * (n + 1) + n] = b[r];
  }
  for (size_t col = 0; col < n; ++col) {
    size_t pivot = col;
    for (size_t r = col + 1; r < n; ++r) {
      if (std::fabs(m[r * (n + 1) + col]) > std::fabs(m[pivot * (n + 1) + col])) pivot = r;
    }
    if (std::fabs(m[pivot * (n + 1) + col]) < 1e-12) return false;
    if (pivot != col) {
      for (size_t c = col; c <= n; ++c) std::swap(m[pivot * (n + 1) + c], m[col * (n + 1) + c]);
    }
    const double div = m[col * (n + 1) + col];
    for (size_t c = col; c <= n; ++c) m[col * (n + 1) + c] /= div;
    for (size_t r = 0; r < n; ++r) if (r != col) {
      const double factor = m[r * (n + 1) + col];
      for (size_t c = col; c <= n; ++c) m[r * (n + 1) + c] -= factor * m[col * (n + 1) + c];
    }
  }
  x->assign(n, 0.0);
  for (size_t r = 0; r < n; ++r) (*x)[r] = m[r * (n + 1) + n];
  return true;
}

bool compute_jacobian(AppState &app, const State &s, std::vector<double> *J, char *err, size_t err_cap) {
  const size_t n = app.state_names.size();
  const double eps = 1e-5;
  J->assign(n * n, 0.0);
  for (size_t col = 0; col < n; ++col) {
    State a = s, b = s;
    set_state_at(a, col, state_at(a, col) + eps);
    set_state_at(b, col, state_at(b, col) - eps);
    State fa{}, fb{};
    if (!eval_rhs(app, a, &fa, err, err_cap) || !eval_rhs(app, b, &fb, err, err_cap)) return false;
    for (size_t row = 0; row < n; ++row) {
      (*J)[row * n + col] = (state_at(fa, row) - state_at(fb, row)) / (2 * eps);
    }
  }
  return true;
}

std::string classify_fixed_jacobian_2d(AppState &app) {
  app.fixed_eigenvalues.clear();
  const size_t n = app.state_names.size();
  if (n != 2 || app.fixed_jacobian.size() != 4) {
    return "classification available for 2D ODEs";
  }
  const double a = app.fixed_jacobian[0];
  const double b = app.fixed_jacobian[1];
  const double c = app.fixed_jacobian[2];
  const double d = app.fixed_jacobian[3];
  const double tr = a + d;
  const double det = a * d - b * c;
  const double disc = tr * tr - 4.0 * det;
  constexpr double eps = 1e-9;
  if (disc >= 0.0) {
    const double root = std::sqrt(std::max(0.0, disc));
    app.fixed_eigenvalues.push_back({0.5 * (tr + root), 0.0});
    app.fixed_eigenvalues.push_back({0.5 * (tr - root), 0.0});
  } else {
    app.fixed_eigenvalues.push_back({0.5 * tr, 0.5 * std::sqrt(-disc)});
    app.fixed_eigenvalues.push_back({0.5 * tr, -0.5 * std::sqrt(-disc)});
  }
  if (det < -eps) return "saddle";
  if (std::fabs(det) <= eps) return "degenerate";
  if (disc > eps) return tr < -eps ? "stable node" : (tr > eps ? "unstable node" : "degenerate node");
  if (disc < -eps) {
    if (tr < -eps) return "stable spiral";
    if (tr > eps) return "unstable spiral";
    return "center";
  }
  if (tr < -eps) return "stable improper node";
  if (tr > eps) return "unstable improper node";
  return "degenerate";
}

Point2 eigenvector_2x2(double a, double b, double c, double d, double lambda) {
  double vx = b;
  double vy = lambda - a;
  if (std::fabs(vx) + std::fabs(vy) < 1e-12) {
    vx = lambda - d;
    vy = c;
  }
  if (std::fabs(vx) + std::fabs(vy) < 1e-12) {
    vx = 1.0;
    vy = 0.0;
  }
  const double len = std::sqrt(vx * vx + vy * vy);
  return Point2{vx / len, vy / len};
}

void append_separatrix_curve(AppState &app, const State &seed, double signed_dt,
                             uint32_t color, const char *label) {
  AppState::PhaseTrajectory curve;
  curve.label = label;
  curve.current = seed;
  curve.color = color;
  curve.active = false;
  curve.visible = true;
  push_phase_history(app, curve, curve.current);

  const double old_dt = app.dt;
  app.dt = signed_dt;
  char err[256] = {0};
  for (int i = 0; i < app.separatrix_steps; ++i) {
    State next{};
    if (!step_ode_state(app, curve.current, &next, err, sizeof(err))) {
      app.runtime_error = err;
      break;
    }
    curve.current = next;
    push_phase_history(app, curve, curve.current);
    bool finite = true;
    for (double v : curve.current.v) finite = finite && std::isfinite(v);
    if (!finite) break;
  }
  app.dt = old_dt;
  app.separatrix_curves.push_back(std::move(curve));
}

void trace_separatrices(AppState &app) {
  if (app.mode != SystemMode::ODE || app.state_names.size() != 2) {
    app.analysis_message = "separatrix tracing is currently implemented for 2D ODEs only";
    return;
  }
  if (!app.fixed_ready || app.fixed_jacobian.size() != 4) {
    app.analysis_message = "find a 2D fixed point before tracing separatrices";
    return;
  }
  if (app.fixed_classification.find("saddle") == std::string::npos) {
    app.analysis_message = "separatrix tracing needs a saddle fixed point";
    return;
  }
  const double a = app.fixed_jacobian[0];
  const double b = app.fixed_jacobian[1];
  const double c = app.fixed_jacobian[2];
  const double d = app.fixed_jacobian[3];
  const double tr = a + d;
  const double det = a * d - b * c;
  const double disc = tr * tr - 4.0 * det;
  if (disc <= 0.0) {
    app.analysis_message = "saddle separatrix tracing needs real eigenvalues";
    return;
  }
  const double root = std::sqrt(disc);
  const double l1 = 0.5 * (tr + root);
  const double l2 = 0.5 * (tr - root);
  const double lambdas[2] = {l1, l2};
  app.separatrix_curves.clear();
  const double base_dt = std::max(1e-6, std::fabs(app.dt));
  for (int k = 0; k < 2; ++k) {
    const double lambda = lambdas[k];
    const Point2 v = eigenvector_2x2(a, b, c, d, lambda);
    const bool unstable = lambda > 0.0;
    const double signed_dt = unstable ? base_dt : -base_dt;
    for (int sign = -1; sign <= 1; sign += 2) {
      State seed = app.fixed_point;
      resize_state(seed, 2);
      set_state_at(seed, 0, state_at(seed, 0) + sign * app.separatrix_epsilon * v.x);
      set_state_at(seed, 1, state_at(seed, 1) + sign * app.separatrix_epsilon * v.y);
      append_separatrix_curve(app, seed, signed_dt,
                              unstable ? 0xffff7070u : 0xff70dfffu,
                              unstable ? "unstable manifold" : "stable manifold");
    }
  }
  app.analysis_message = "traced local saddle separatrices";
}

void find_fixed_point(AppState &app) {
  if (app.mode != SystemMode::ODE) {
    app.analysis_message = "fixed-point Newton solver currently applies to ODE mode only";
    return;
  }
  const size_t n = app.state_names.size();
  State s = app.current;
  resize_state(s, n);
  char err[256] = {0};
  double residual = 0.0;
  app.fixed_jacobian.clear();
  for (int iter = 0; iter < 40; ++iter) {
    State f{};
    if (!eval_rhs(app, s, &f, err, sizeof(err))) { app.analysis_message = err; return; }
    double residual2 = 0.0;
    for (size_t i = 0; i < n; ++i) residual2 += state_at(f, i) * state_at(f, i);
    residual = std::sqrt(residual2);
    /* Keep the latest Jacobian — also serves as the final one if we
     * exit via convergence or step-size below. */
    if (!compute_jacobian(app, s, &app.fixed_jacobian, err, sizeof(err))) {
      app.analysis_message = err; return;
    }
    if (residual < 1e-9) break;
    std::vector<double> rhs(n, 0.0);
    for (size_t i = 0; i < n; ++i) rhs[i] = -state_at(f, i);
    std::vector<double> delta;
    if (!solve_linear_system(app.fixed_jacobian, rhs, &delta)) {
      app.analysis_message = "Newton solve failed: singular Jacobian"; return;
    }
    double delta2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
      set_state_at(s, i, state_at(s, i) + delta[i]);
      delta2 += delta[i] * delta[i];
    }
    if (std::sqrt(delta2) < 1e-10) {
      /* Re-evaluate Jacobian at the converged point so it
       * matches the reported `fixed_point` exactly. */
      if (!compute_jacobian(app, s, &app.fixed_jacobian, err, sizeof(err))) {
        app.analysis_message = err; return;
      }
      break;
    }
  }
  app.fixed_point = s;
  app.fixed_residual = residual;
  app.fixed_ready = true;
  app.separatrix_curves.clear();
  app.fixed_classification = classify_fixed_jacobian_2d(app);

  /* PHASE1: general N-D classification via the from-scratch eigensolver.
   * Works in any dimension and supersedes the 2D-only trace/det path
   * for the spectrum and the saddle/Hopf/fold flags. */
  app.fixed_general_ready = false;
  app.cas_eig_ready = false; /* invalidate any cached CAS spectrum */
  app.hopf_l1_ready = false; /* invalidate cached Hopf classification */
  if (!app.fixed_jacobian.empty()) {
    app.fixed_general =
        dynsys::analysis::classify_equilibrium(app.fixed_jacobian, n);
    app.fixed_general_ready = true;
  }

  if (app.fixed_general_ready) {
    app.analysis_message =
        std::string("fixed point search completed: ") + app.fixed_general.label;
  } else if (app.state_names.size() == 2) {
    app.analysis_message = std::string("fixed point search completed: ") + app.fixed_classification;
  } else {
    app.analysis_message = "fixed point search completed";
  }
}

/* ============================================================
 * PHASE1/2/3: general analysis bridge.
 *
 * build_model wraps the compiled system as a dynsys::analysis::Model.
 * The vector field temporarily installs the continuation parameter
 * value, evaluates the RHS through the existing IR path, and restores
 * the parameter. The Jacobian is supplied by the analysis core's
 * finite-difference path for now; the forward-mode-AD-over-IR
 * implementation (Phase 3) will set model.jacobian_x without changing
 * any caller. n_state is fixed; the parameter is z[n].
 * ============================================================ */
dynsys::analysis::Model build_model(AppState &app, AppState::Param *param) {
  dynsys::analysis::Model model;
  model.n = app.state_names.size();

  /* Capture by value of the pointer; the closures run synchronously
   * during continuation, while `app` and `param` stay alive. */
  model.vector_field = [&app, param](const double *x, double p, double *f_out,
                                     std::string *err) -> bool {
    const size_t n = app.state_names.size();
    const double saved = param->value;
    param->value = p;
    sync_param_values(app);
    State s = make_state_like(n, app.current.t);
    for (size_t i = 0; i < n; ++i) set_state_at(s, i, x[i]);
    State deriv{};
    char buf[256] = {0};
    const bool ok = eval_rhs(app, s, &deriv, buf, sizeof(buf));
    param->value = saved;
    sync_param_values(app);
    if (!ok) {
      if (err) *err = buf;
      return false;
    }
    for (size_t i = 0; i < n; ++i) f_out[i] = state_at(deriv, i);
    return true;
  };
  /* PHASE3: exact Jacobian and df/dp via forward-mode AD over the IR.
   * For column j we evaluate every equation program once with the dual
   * seed on state variable j; the derivative component is J[row][j].
   * df/dp seeds the active continuation parameter instead. This makes
   * the continuation engine's linear algebra exact rather than
   * finite-difference-noisy, with no change to any caller. */
  size_t param_index = 0;
  for (size_t i = 0; i < app.params.size(); ++i)
    if (&app.params[i] == param) param_index = i;

  model.jacobian_x = [&app, param, param_index](
                         const double *x, double p, double *jac_out,
                         std::string *err) -> bool {
    const size_t n = app.state_names.size();
    if (app.equation_programs.size() != n) {
      if (err) *err = "equations are not compiled";
      return false;
    }
    const double saved = param->value;
    param->value = p;
    sync_param_values(app);

    State s = make_state_like(n, app.current.t);
    for (size_t i = 0; i < n; ++i) set_state_at(s, i, x[i]);

    dynsys::ir::RunContext rc;
    rc.state = s.v.data();
    rc.n_state = s.v.size();
    rc.t = s.t;
    rc.params = app.param_values.data();
    rc.n_params = app.param_values.size();
    rc.defs = app.definition_programs.data();
    rc.n_defs = app.definition_programs.size();

    char buf[256] = {0};
    bool ok = true;
    for (size_t col = 0; col < n && ok; ++col) {
      dynsys::ir::DualSeed seed{dynsys::ir::DualSeed::Kind::State, col};
      for (size_t row = 0; row < n; ++row) {
        double v = 0.0, d = 0.0;
        if (!dynsys::ir::run_dual(app.equation_programs[row], rc, seed,
                                  app.ad_scratch, &v, &d, buf, sizeof(buf))) {
          ok = false;
          break;
        }
        jac_out[row * n + col] = d;
      }
    }
    param->value = saved;
    sync_param_values(app);
    (void)param_index;
    if (!ok && err) *err = buf;
    return ok;
  };

  model.dfdp = [&app, param, param_index](const double *x, double p,
                                          double *dfdp_out,
                                          std::string *err) -> bool {
    const size_t n = app.state_names.size();
    if (app.equation_programs.size() != n) {
      if (err) *err = "equations are not compiled";
      return false;
    }
    const double saved = param->value;
    param->value = p;
    sync_param_values(app);

    State s = make_state_like(n, app.current.t);
    for (size_t i = 0; i < n; ++i) set_state_at(s, i, x[i]);

    dynsys::ir::RunContext rc;
    rc.state = s.v.data();
    rc.n_state = s.v.size();
    rc.t = s.t;
    rc.params = app.param_values.data();
    rc.n_params = app.param_values.size();
    rc.defs = app.definition_programs.data();
    rc.n_defs = app.definition_programs.size();

    char buf[256] = {0};
    bool ok = true;
    dynsys::ir::DualSeed seed{dynsys::ir::DualSeed::Kind::Param, param_index};
    for (size_t row = 0; row < n; ++row) {
      double v = 0.0, d = 0.0;
      if (!dynsys::ir::run_dual(app.equation_programs[row], rc, seed,
                                app.ad_scratch, &v, &d, buf, sizeof(buf))) {
        ok = false;
        break;
      }
      dfdp_out[row] = d;
    }
    param->value = saved;
    sync_param_values(app);
    if (!ok && err) *err = buf;
    return ok;
  };

  return model;
}

/* Two-parameter model for fold/Hopf-curve continuation: the vector field as a
 * function of (x, p, q) where p and q are two chosen parameters. Finite-diff
 * Jacobians inside two_param_curve do the rest, so we only need the field. */
dynsys::analysis::Model2 build_model2(AppState &app, AppState::Param *pp, AppState::Param *qp) {
  dynsys::analysis::Model2 model;
  model.n = app.state_names.size();
  model.vector_field = [&app, pp, qp](const double *x, double p, double q,
                                      double *f_out, std::string *err) -> bool {
    const size_t n = app.state_names.size();
    const double sp = pp->value, sq = qp->value;
    pp->value = p; qp->value = q;
    sync_param_values(app);
    State s = make_state_like(n, app.current.t);
    for (size_t i = 0; i < n; ++i) set_state_at(s, i, x[i]);
    State deriv{};
    char buf[256] = {0};
    const bool ok = eval_rhs(app, s, &deriv, buf, sizeof(buf));
    pp->value = sp; qp->value = sq;
    sync_param_values(app);
    if (!ok) { if (err) *err = buf; return false; }
    for (size_t i = 0; i < n; ++i) f_out[i] = state_at(deriv, i);
    return true;
  };
  return model;
}

/* PHASE E: trace a fold or Hopf CURVE in a two-parameter plane. The first
 * parameter is the continuation parameter (cont_param), the second is
 * twopar_p2; the start point is the current located fixed point at the current
 * (p,q). Reuses build_model2 + analysis::two_param_curve. */
void run_two_param_curve(AppState &app, dynsys::analysis::TwoParamKind kind) {
  app.twopar_ready = false;
  if (app.mode != SystemMode::ODE) { app.twopar_msg = "two-parameter curves are for ODE systems"; return; }
  AppState::Param *pp = find_param(app, app.cont_param);
  AppState::Param *qp = find_param(app, app.twopar_p2);
  if (pp == nullptr || qp == nullptr) { app.twopar_msg = "pick two distinct parameters (cont param + 2nd param)"; return; }
  if (pp == qp) { app.twopar_msg = "the two parameters must be different"; return; }
  const size_t n = app.state_names.size();
  if (n == 0) { app.twopar_msg = "no system"; return; }

  std::vector<double> x0(n);
  for (size_t i = 0; i < n; ++i) x0[i] = app.fixed_ready ? state_at(app.fixed_point, i) : state_at(app.current, i);

  dynsys::analysis::Model2 model = build_model2(app, pp, qp);
  dynsys::analysis::TwoParamSettings s;
  if (pp->has_range) { s.p_min = pp->min_value; s.p_max = pp->max_value; }
  else { s.p_min = pp->value - 5.0; s.p_max = pp->value + 5.0; }
  if (qp->has_range) { s.q_min = qp->min_value; s.q_max = qp->max_value; }
  else { s.q_min = qp->value - 5.0; s.q_max = qp->value + 5.0; }
  s.h0 = 0.05; s.max_points = 800;

  app.twopar_curve = dynsys::analysis::two_param_curve(model, kind, x0, pp->value, qp->value, s);
  app.twopar_ready = app.twopar_curve.ok;
  sync_param_values(app);
  app.twopar_msg = app.twopar_curve.ok
      ? (std::string(kind == dynsys::analysis::TwoParamKind::Fold ? "fold" : "Hopf") +
         " curve: " + app.twopar_curve.message)
      : std::string("two-parameter curve: ") + app.twopar_curve.message;
  app.analysis_message = app.twopar_msg;
}

/* Homoclinic-orbit search. Pick a SADDLE (from the auto-scanned planar fixed
 * points, else correct one from the current state), nudge a point off it along
 * the unstable eigendirection, integrate forward to trace one excursion that
 * (for a near-homoclinic system) loops back near the saddle, and hand that as
 * the seed to analysis::solve_homoclinic (truncated BVP + projection BCs). */
void run_homoclinic(AppState &app) {
  app.homoclinic_ready = false;
  if (app.mode != SystemMode::ODE) { app.homoclinic_msg = "homoclinic orbits are for ODE systems"; return; }
  const size_t n = app.state_names.size();
  if (n < 2) { app.homoclinic_msg = "need a 2+ dimensional ODE"; return; }
  AppState::Param *param = app.params.empty() ? nullptr : find_param(app, app.cont_param);
  if (param == nullptr && !app.params.empty()) param = &app.params[0];
  const double pval = param ? param->value : 0.0;

  dynsys::analysis::Model model = build_model(app, param);

  /* 1. choose a saddle. Prefer a planar saddle from the auto-scan; otherwise
   * correct an equilibrium from the current state and require it be a saddle. */
  std::vector<double> x0(n, 0.0);
  for (size_t i = 0; i < n; ++i) x0[i] = state_at(app.current, i);
  bool have_saddle = false;
  if (n == 2 && !app.phase_fixed_points.empty()) {
    for (const auto &fp : app.phase_fixed_points)
      if (fp.is_saddle) { x0 = {fp.x, fp.y}; have_saddle = true; break; }
  }
  /* Newton-correct the equilibrium from x0 using the model's vector field +
   * finite-difference Jacobian (a few iterations; local to keep this driver
   * self-contained). */
  {
    std::vector<double> xc = x0, f(n), Jl, delta;
    std::string er; bool ok = false;
    for (int it = 0; it < 80; ++it) {
      if (!model.vector_field(xc.data(), pval, f.data(), &er)) break;
      double fn = 0; for (double v : f) fn += v*v; fn = std::sqrt(fn);
      if (fn < 1e-11) { ok = true; break; }
      if (!dynsys::analysis::finite_diff_jacobian(model, xc.data(), pval, &Jl, &er, 1e-7)) break;
      std::vector<double> rhs(n); for (size_t i = 0; i < n; ++i) rhs[i] = -f[i];
      if (!dynsys::analysis::solve_linear(Jl, rhs, &delta)) break;
      double dn = 0; for (double v : delta) dn += v*v; dn = std::sqrt(dn);
      double damp = dn > 1.0 ? 1.0/dn : 1.0;
      for (size_t i = 0; i < n; ++i) xc[i] += damp*delta[i];
    }
    if (ok) { x0 = xc; have_saddle = true; }
  }
  if (!have_saddle) { app.homoclinic_msg = "could not locate a saddle equilibrium to connect"; return; }

  /* Jacobian + unstable eigendirection for the nudge. */
  std::vector<double> J; std::string e;
  if (!dynsys::analysis::finite_diff_jacobian(model, x0.data(), pval, &J, &e, 1e-7)) {
    app.homoclinic_msg = "Jacobian failed at the saddle"; return; }
  std::vector<dynsys::analysis::Complex> ev;
  dynsys::analysis::eigenvalues(J, n, &ev);
  int n_pos = 0; for (auto &z : ev) if (z.real() > 1e-7) ++n_pos;
  if (n_pos == 0 || n_pos == (int)n) { app.homoclinic_msg = "equilibrium is not a saddle (no homoclinic possible)"; return; }

  /* 2. seed by integrating the unstable manifold (robust helper: tries both
   * orientations, keeps the excursion that returns nearest the saddle). */
  const double dt = (app.dt > 0 ? app.dt : 0.01);
  std::vector<std::vector<double>> seed;
  if (!dynsys::analysis::seed_homoclinic_by_integration(model, x0, pval, dt, 2000.0*dt, &seed, &e)) {
    app.homoclinic_msg = std::string("homoclinic seeding failed: ") + e; return;
  }
  if (seed.size() < 20) { app.homoclinic_msg = "could not trace an excursion off the saddle"; return; }

  /* 3. solve the truncated BVP. The truncation half-time should be a MODERATE
   * multiple of the saddle's slowest timescale (1/min|Re lambda|) -- NOT the
   * raw simulation length, which can be huge when the saddle quantity is
   * nonzero and the seed orbit spirals slowly. Too-large T wrecks the uniform
   * mesh resolution. */
  double lam_min = 1e9;
  for (auto &z : ev) { double a = std::fabs(z.real()); if (a > 1e-7) lam_min = std::min(lam_min, a); }
  if (!std::isfinite(lam_min) || lam_min <= 0) lam_min = 1.0;
  double Tsane = 8.0 / lam_min;                     /* ~8 e-foldings to the saddle */
  double Tseed = 0.5 * dt * (double)seed.size();
  dynsys::analysis::HomoclinicSettings hs;
  hs.mesh = std::min(300, std::max(80, (int)seed.size()/8));
  hs.T = std::min(Tseed, std::max(8.0, Tsane));     /* moderate, capped         */
  hs.free_T = false; hs.newton_iters = 80;
  app.homoclinic = dynsys::analysis::solve_homoclinic(model, x0, pval, seed, hs);
  app.homoclinic_ready = app.homoclinic.ok;
  app.homoclinic_msg = "homoclinic: " + app.homoclinic.message;
  app.analysis_message = app.homoclinic_msg;
}

/* Two-parameter continuation of the homoclinic locus. Requires a homoclinic to
 * have been found first (run_homoclinic): its converged orbit seeds the curve.
 * Continues in the secondary parameter twopar_p2, finding the primary parameter
 * cont_param at which the connection closes at each step. */
void run_homoclinic_curve(AppState &app) {
  app.homoclinic_curve_ready = false;
  if (!app.homoclinic_ready || app.homoclinic.orbit.size() < 4) {
    app.homoclinic_msg = "find a homoclinic orbit first (it seeds the curve)"; return;
  }
  AppState::Param *pp = app.params.empty() ? nullptr : find_param(app, app.cont_param);
  AppState::Param *qp = find_param(app, app.twopar_p2);
  if (pp == nullptr) pp = app.params.empty() ? nullptr : &app.params[0];
  if (pp == nullptr || qp == nullptr || pp == qp) {
    app.homoclinic_msg = "set two distinct parameters (cont_param + twopar_p2) for the curve"; return;
  }
  dynsys::analysis::Model2 m2 = build_model2(app, pp, qp);
  dynsys::analysis::HomoclinicContSettings cs;
  cs.dq = std::max(1e-3, (qp->has_range ? (qp->max_value-qp->min_value)/60.0 : 0.05));
  cs.max_steps = 40; cs.both_directions = true;
  if (qp->has_range) { cs.q_min = qp->min_value; cs.q_max = qp->max_value; }
  cs.store_orbits = false;
  cs.bvp.mesh = (int)app.homoclinic.orbit.size() - 1;
  cs.bvp.T = app.homoclinic.T; cs.bvp.free_T = false; cs.bvp.newton_iters = 60;
  app.homoclinic_curve = dynsys::analysis::continue_homoclinic(
      m2, app.homoclinic.saddle, pp->value, qp->value, app.homoclinic.orbit, cs);
  app.homoclinic_curve_ready = app.homoclinic_curve.ok;
  app.homoclinic_msg = "homoclinic curve: " + app.homoclinic_curve.message;
  app.analysis_message = app.homoclinic_msg;
}

/* PHASE E: periodic-orbit continuation by COLLOCATION. Simulate at the current
 * continuation parameter to land on a cycle, extract one loop as the initial
 * guess (mesh points + period), then hand it to analysis::continue_limit_cycle,
 * which solves the periodic BVP exactly (cycle + period) and follows it — so it
 * tracks unstable cycles and folds of cycles, unlike the lcc_* sweep. */
void run_cycle_collocation(AppState &app) {
  app.cyc_ready = false;
  if (app.mode != SystemMode::ODE) { app.cyc_msg = "periodic-orbit continuation is for ODE systems"; return; }
  AppState::Param *param = app.params.empty() ? nullptr : find_param(app, app.cont_param);
  if (param == nullptr && !app.params.empty()) param = &app.params[0];
  if (param == nullptr) { app.cyc_msg = "declare a parameter to continue in"; return; }
  const size_t n = app.state_names.size();
  if (n < 2) { app.cyc_msg = "need at least 2 state variables for a periodic orbit"; return; }

  const Integrator saved = app.integrator;
  if (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45)
    app.integrator = Integrator::RK4;
  const double dt = app.dt > 0 ? app.dt : 0.01;

  /* settle onto the attractor, then record until the first state returns near
   * its starting value with matching direction (one loop). */
  char err[256] = {0};
  State s = app.start; resize_state(s, n);
  bool ok = true;
  for (int j = 0; j < 40000; ++j) { State nx{}; if (!step_ode_state(app, s, &nx, err, sizeof(err))) { ok = false; break; } s = nx; }
  if (!ok) { app.integrator = saved; app.cyc_msg = "integration failed while settling onto a cycle"; return; }

  /* record a long sample, capture the section crossing of x0 through its mean */
  std::vector<std::vector<double>> traj;
  std::vector<double> times;
  const int rec = 60000;
  traj.reserve(rec);
  double mean0 = 0;
  { State t = s; for (int j = 0; j < rec; ++j) { State nx{}; if (!step_ode_state(app, t, &nx, err, sizeof(err))) { ok=false; break; }
      std::vector<double> row(n); for (size_t i=0;i<n;i++) row[i]=state_at(nx,i); traj.push_back(row); times.push_back((j+1)*dt); mean0+=row[0]; t=nx; } }
  if (!ok || traj.size() < 100) { app.integrator = saved; app.cyc_msg = "could not record a cycle"; return; }
  mean0 /= (double)traj.size();

  /* find two successive upward crossings of x0 == mean0 -> one period */
  int c1 = -1, c2 = -1;
  for (size_t j = 1; j < traj.size(); ++j) {
    if (traj[j-1][0] < mean0 && traj[j][0] >= mean0) {
      if (c1 < 0) c1 = (int)j; else { c2 = (int)j; break; }
    }
  }
  if (c1 < 0 || c2 < 0 || c2 <= c1) { app.integrator = saved; app.cyc_msg = "no clean periodic crossing found (is there a stable cycle here?)"; return; }
  const double period_guess = (c2 - c1) * dt;

  /* resample the loop [c1,c2) into `mesh` evenly spaced points */
  const int mesh = 60;
  std::vector<std::vector<double>> guess(mesh, std::vector<double>(n, 0.0));
  for (int k = 0; k < mesh; ++k) {
    const double frac = (double)k / mesh;
    const int idx = c1 + (int)std::floor(frac * (c2 - c1));
    guess[(size_t)k] = traj[(size_t)std::min(idx, (int)traj.size()-1)];
  }
  app.integrator = saved;

  dynsys::analysis::Model model = build_model(app, param);
  dynsys::analysis::CycleSettings cs;
  if (param->has_range) { cs.p_min = param->min_value; cs.p_max = param->max_value; }
  else { cs.p_min = param->value - 3.0; cs.p_max = param->value + 3.0; }
  cs.mesh = mesh; cs.dp = std::max(1e-3, (cs.p_max - cs.p_min) / 120.0); cs.max_steps = 400;
  /* pseudo-arclength so the branch follows folds of cycles (LPC), with Floquet
   * multipliers for stability / period-doubling / torus detection. The
   * arclength step is sized to the parameter window so we cover it in a sane
   * number of steps but can still turn around folds. */
  cs.arclength = true;
  cs.compute_floquet = true;
  cs.ds = std::max(1e-3, (cs.p_max - cs.p_min) / 80.0);

  app.cyc_branch = dynsys::analysis::continue_limit_cycle(model, guess, period_guess, param->value, cs);
  app.cyc_ready = app.cyc_branch.ok;
  sync_param_values(app);
  app.cyc_msg = app.cyc_branch.ok
      ? "collocation: " + app.cyc_branch.message
      : std::string("collocation: ") + app.cyc_branch.message;
  app.analysis_message = app.cyc_msg;
}

/* PHASE E (nicety): trace a FOLD-OF-CYCLES (LPC) curve in the
 * (cont_param, twopar_p2) plane. Seeds a cycle by simulation at the current
 * (p,q), then calls analysis::lpc_curve which, for each q, continues the cycle
 * in p and locates the saddle-node of cycles. */
void run_lpc_curve(AppState &app) {
  app.lpc_ready = false;
  if (app.mode != SystemMode::ODE) { app.lpc_msg = "fold-of-cycles curves are for ODE systems"; return; }
  AppState::Param *pp = find_param(app, app.cont_param);
  AppState::Param *qp = find_param(app, app.twopar_p2);
  if (pp == nullptr || qp == nullptr || pp == qp) { app.lpc_msg = "pick two distinct parameters (cont param + 2nd param)"; return; }
  const size_t n = app.state_names.size();
  if (n < 2) { app.lpc_msg = "need at least 2 state variables for a periodic orbit"; return; }

  /* seed a cycle at the current parameters by simulation (same approach as the
   * collocation runner) */
  const Integrator saved = app.integrator;
  if (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45) app.integrator = Integrator::RK4;
  const double dt = app.dt > 0 ? app.dt : 0.01;
  char err[256] = {0};
  State s = app.start; resize_state(s, n); bool ok = true;
  for (int j = 0; j < 40000; ++j) { State nx{}; if (!step_ode_state(app, s, &nx, err, sizeof(err))) { ok=false; break; } s = nx; }
  std::vector<std::vector<double>> traj; double mean0 = 0;
  if (ok) { State t = s; for (int j = 0; j < 60000; ++j) { State nx{}; if (!step_ode_state(app, t, &nx, err, sizeof(err))) { ok=false; break; }
      std::vector<double> row(n); for (size_t i=0;i<n;i++) row[i]=state_at(nx,i); traj.push_back(row); mean0+=row[0]; t=nx; } }
  if (!ok || traj.size() < 100) { app.integrator = saved; app.lpc_msg = "could not seed a cycle at the current parameters"; return; }
  mean0 /= (double)traj.size();
  int c1=-1,c2=-1;
  for (size_t j=1;j<traj.size();++j) if (traj[j-1][0]<mean0 && traj[j][0]>=mean0) { if(c1<0)c1=(int)j; else {c2=(int)j; break;} }
  if (c1<0||c2<0||c2<=c1) { app.integrator = saved; app.lpc_msg = "no clean cycle to seed from at the current parameters"; return; }
  const double period_guess = (c2-c1)*dt;
  const int mesh = 50;
  std::vector<std::vector<double>> guess(mesh, std::vector<double>(n,0.0));
  for (int k=0;k<mesh;k++){ const int idx=c1+(int)std::floor((double)k/mesh*(c2-c1)); guess[(size_t)k]=traj[(size_t)std::min(idx,(int)traj.size()-1)]; }
  app.integrator = saved;

  dynsys::analysis::Model2 model = build_model2(app, pp, qp);
  dynsys::analysis::TwoParamSettings s2;
  if (pp->has_range) { s2.p_min=pp->min_value; s2.p_max=pp->max_value; } else { s2.p_min=pp->value-3.0; s2.p_max=pp->value+3.0; }
  if (qp->has_range) { s2.q_min=qp->min_value; s2.q_max=qp->max_value; } else { s2.q_min=qp->value-3.0; s2.q_max=qp->value+3.0; }
  s2.max_points = 40;
  dynsys::analysis::CycleSettings cs; cs.mesh=mesh; cs.dp=std::max(1e-3,(s2.p_max-s2.p_min)/120.0); cs.max_steps=200;

  app.lpc_curve_data = dynsys::analysis::lpc_curve(model, guess, period_guess, pp->value, qp->value, s2, cs);
  app.lpc_ready = app.lpc_curve_data.ok;
  sync_param_values(app);
  app.lpc_msg = app.lpc_curve_data.ok ? ("fold-of-cycles: " + app.lpc_curve_data.message)
                                      : std::string("fold-of-cycles: ") + app.lpc_curve_data.message;
  app.analysis_message = app.lpc_msg;
}

/* PHASE E: exact equilibria via the CAS. Extract each state equation's RHS
 * text from the system source, substitute current parameter values as exact
 * rationals, and ask Sangaku to solve f(x)=0 exactly (solve-poly in 1-D,
 * Groebner in N-D). Polynomial fields only; anything transcendental is
 * reported as such. Results land in app.cas_equi_* for the UI. */
void run_exact_equilibria(AppState &app) {
  app.cas_equi_ready = false;
  app.cas_equi_lines.clear();
  if (!dynsys::cas::is_available()) {
    app.cas_equi_msg = "CAS unavailable (set LIZARD and SANGAKU_ROOT, or use the Nix shell)";
    return;
  }
  if (app.mode == SystemMode::IFS) { app.cas_equi_msg = "exact equilibria apply to ODE/map systems"; return; }
  const size_t n = app.state_names.size();
  if (n == 0) { app.cas_equi_msg = "no system"; return; }

  /* pull RHS text per state variable from the source */
  std::vector<std::string> rhs(n);
  std::vector<bool> got(n, false);
  {
    std::string src(app.system_input);
    std::istringstream is(src);
    std::string line;
    while (std::getline(is, line)) {
      /* strip comments (after '#') and skip blanks */
      auto hash = line.find('#'); if (hash != std::string::npos) line = line.substr(0, hash);
      auto eq = line.find('=');
      if (eq == std::string::npos) continue;
      std::string lhs = line.substr(0, eq), r = line.substr(eq + 1);
      std::string name, perr; std::vector<std::string> ps;
      if (!parse_lhs(lhs, &name, &ps, &perr)) continue;
      int di = derivative_index_for_name(app.state_names, name);
      int ni = next_index_for_name(app.state_names, name);
      int idx = di >= 0 ? di : ni;
      if (idx < 0 || idx >= (int)n) continue;
      /* for a MAP, the equilibrium condition is x_next = x, i.e. RHS - x = 0 */
      if (ni >= 0) r = "(" + r + ") - (" + app.state_names[(size_t)idx] + ")";
      rhs[(size_t)idx] = r; got[(size_t)idx] = true;
    }
  }
  for (size_t i = 0; i < n; ++i)
    if (!got[i]) { app.cas_equi_msg = "could not find the equation for " + app.state_names[i]; return; }

  /* current parameter values as exact rationals */
  std::vector<std::pair<std::string, dynsys::cas::Rational>> subst;
  for (const auto &p : app.params) {
    bool ex = false;
    dynsys::cas::Rational rr = dynsys::cas::rationalize(p.value, &ex, 1e-9);
    subst.emplace_back(p.name, rr);
  }

  dynsys::cas::EquilibriaReport rep =
      dynsys::cas::solve_equilibria(rhs, app.state_names, subst);
  app.cas_equi_ready = rep.ok;
  app.cas_equi_lines = rep.solutions;
  app.cas_equi_is_poly = rep.polynomial;
  app.cas_equi_msg = rep.ok ? std::string("exact (") + std::to_string(rep.solutions.size()) + " result line(s))"
                            : rep.message;
}

void run_equilibrium_continuation(AppState &app) {
  if (app.mode != SystemMode::ODE) {
    app.analysis_message = "continuation currently applies to ODE mode only";
    return;
  }
  if (!app.fixed_ready) {
    app.analysis_message = "find a fixed point before starting continuation";
    return;
  }
  AppState::Param *param = find_param(app, app.cont_param);
  if (param == nullptr) {
    app.analysis_message =
        std::string("unknown continuation parameter: ") + app.cont_param;
    return;
  }

  dynsys::analysis::Model model = build_model(app, param);
  std::vector<double> x0(model.n, 0.0);
  for (size_t i = 0; i < model.n; ++i) x0[i] = state_at(app.fixed_point, i);

  dynsys::analysis::ContinuationSettings settings;
  settings.h0 = app.cont_h0;
  settings.max_points = app.cont_max_points;
  settings.p_min = app.cont_p_min;
  settings.p_max = app.cont_p_max;
  settings.direction = app.cont_direction;
  settings.detect_fold = app.cont_detect_fold;
  settings.detect_hopf = app.cont_detect_hopf;

  const double p0 = param->value;
  app.cont_branch =
      dynsys::analysis::continue_equilibrium(model, x0, p0, settings);
  app.cont_ready = app.cont_branch.ok;
  app.cont_switched_ready = false; /* the previous switched branch is stale */

  /* Restore the live parameter — the model closures already restore it
   * per call, but be explicit. */
  sync_param_values(app);

  size_t n_fold = 0, n_hopf = 0, n_bt = 0, n_cusp = 0, n_gh = 0, n_bp = 0;
  for (size_t i : app.cont_branch.special_indices) {
    using K = dynsys::analysis::SpecialPointKind;
    switch (app.cont_branch.points[i].special) {
      case K::Fold: ++n_fold; break;
      case K::Hopf: ++n_hopf; break;
      case K::BogdanovTakens: ++n_bt; break;
      case K::Cusp: ++n_cusp; break;
      case K::GeneralizedHopf: ++n_gh; break;
      case K::BranchPoint: ++n_bp; break;
      default: break;
    }
  }
  char msg[384];
  std::snprintf(msg, sizeof(msg),
                "continuation: %zu points, %zu fold(s), %zu Hopf; codim-2: %zu BT, %zu cusp, %zu gen-Hopf; %zu branch point(s) — %s",
                app.cont_branch.points.size(), n_fold, n_hopf, n_bt, n_cusp, n_gh, n_bp,
                app.cont_branch.message.c_str());
  app.analysis_message = msg;
}

/* PHASE D: trace an equilibrium branch vs. a parameter, both directions
 * from the current state, and flatten it for drawing. Marks stability and
 * fold/Hopf points. Reuses build_model (AD Jacobian) and the tested
 * continue_equilibrium engine. */
void run_continuation(AppState &app) {
  app.cont_pp.clear(); app.cont_yy.clear(); app.cont_stable.clear();
  app.cont_special.clear(); app.cont_break.clear();
  app.cont_has_branch = false; app.cont_n_fold = 0; app.cont_n_hopf = 0;

  if (app.mode != SystemMode::ODE) { app.cont_message = "continuation is for ODE systems"; return; }
  const size_t n = app.state_names.size();
  if (n == 0) { app.cont_message = "no system"; return; }
  if (app.params.empty()) { app.cont_message = "system has no parameter to continue"; return; }

  AppState::Param *param = find_param(app, app.cont_param);
  if (!param) { param = &app.params[0]; std::snprintf(app.cont_param, sizeof(app.cont_param), "%s", param->name.c_str()); }

  const int yi = std::max(0, std::min(app.cont_y_index, (int)n - 1));

  dynsys::analysis::Model model = build_model(app, param);

  /* seed from the current integrator state (engine corrects it to an
   * equilibrium first). */
  std::vector<double> x0(n);
  for (size_t i = 0; i < n; ++i) x0[i] = state_at(app.current, i);
  const double p0 = param->value;

  dynsys::analysis::ContinuationSettings s;
  s.p_min = std::min(app.cont_p_min, app.cont_p_max);
  s.p_max = std::max(app.cont_p_min, app.cont_p_max);
  s.max_points = std::max(50, app.cont_max_points);
  s.detect_fold = true; s.detect_hopf = true;

  auto flatten = [&](const dynsys::analysis::Branch &b, bool reversed) {
    if (b.points.empty()) return;
    /* mark a polyline break at the join between the two directions */
    if (!app.cont_pp.empty()) app.cont_break.push_back(app.cont_pp.size());
    const size_t cnt = b.points.size();
    for (size_t k = 0; k < cnt; ++k) {
      const size_t idx = reversed ? (cnt - 1 - k) : k;
      const auto &pt = b.points[idx];
      app.cont_pp.push_back(pt.p);
      app.cont_yy.push_back(idx < b.points.size() && yi < (int)pt.x.size() ? pt.x[(size_t)yi] : 0.0);
      app.cont_stable.push_back(pt.stable ? 1 : 0);
      int sp = 0;
      using K = dynsys::analysis::SpecialPointKind;
      switch (pt.special) {
        case K::Fold:            sp = 1; ++app.cont_n_fold; break;
        case K::Hopf:            sp = 2; ++app.cont_n_hopf; break;
        case K::BogdanovTakens:  sp = 3; break;
        case K::Cusp:            sp = 4; break;
        case K::GeneralizedHopf: sp = 5; break;
        case K::BranchPoint:     sp = 6; break;
        default: sp = 0; break;
      }
      app.cont_special.push_back(sp);
    }
  };

  /* direction -1 first (reversed so points run min->max), then +1 */
  s.direction = -1;
  dynsys::analysis::Branch bneg = dynsys::analysis::continue_equilibrium(model, x0, p0, s);
  s.direction = +1;
  dynsys::analysis::Branch bpos = dynsys::analysis::continue_equilibrium(model, x0, p0, s);

  flatten(bneg, true);
  flatten(bpos, false);

  /* fold/hopf counts double-count the shared seed region across the two
   * sweeps; halve to a sensible display count (min 0). */
  app.cont_has_branch = !app.cont_pp.empty();
  app.cont_view_valid = false;
  if (!app.cont_has_branch)
    app.cont_message = bneg.message.empty() ? bpos.message : bneg.message;
  else
    app.cont_message = "branch traced";
}

void run_bifurcation(AppState &app) {
  /* Robustness: if the configured parameter/observable don't match the
   * current system (e.g. switched systems via the preset dropdown and the
   * old names linger), snap them to sane defaults so the sweep is always
   * meaningful instead of silently producing a flat/degenerate diagram. */
  if (!app.params.empty() && find_param(app, app.bif_param) == nullptr) {
    std::snprintf(app.bif_param, sizeof(app.bif_param), "%s", app.params[0].name.c_str());
    if (app.params[0].has_range) { app.bif_start = app.params[0].min_value; app.bif_end = app.params[0].max_value; }
  }
  {
    /* observable must be a known state name (or keep a user expression that
     * parses); default to the first state if it's empty or not a state. */
    bool obs_is_state = false;
    for (const auto &nm : app.state_names) if (nm == app.bif_observable) { obs_is_state = true; break; }
    if ((app.bif_observable[0] == '\0' || !obs_is_state) && !app.state_names.empty()) {
      /* only override if it doesn't parse as an expression either */
      arena_t a{}; arena_init(&a, 4096); std::string e;
      node_t *probe = parse_expression_or_fail(&a, app.bif_observable, "obs", &e);
      if (probe == nullptr) std::snprintf(app.bif_observable, sizeof(app.bif_observable), "%s", app.state_names[0].c_str());
      arena_destroy(&a);
    }
  }

  AppState::Param *param = find_param(app, app.bif_param);
  if (param == nullptr) {
    app.analysis_message = std::string("unknown parameter for bifurcation: ") + app.bif_param;
    return;
  }
  arena_t tmp_arena{};
  arena_init(&tmp_arena, 65536);
  std::string err_msg;
  node_t *obs = parse_expression_or_fail(&tmp_arena, app.bif_observable, "bifurcation observable", &err_msg);
  if (obs == nullptr) { app.analysis_message = err_msg; arena_destroy(&tmp_arena); return; }
  app.bifurcation_points.clear();
  app.bifurcation_lyapunov.clear();
  app.bifurcation_period.clear();
  const double old_value = param->value;
  /* Force a cheap fixed-step method for the sweep (adaptive integrators
   * would take many substeps per step over hundreds of slices and hang). */
  const Integrator saved_integrator = app.integrator;
  if (app.mode == SystemMode::ODE &&
      (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45))
    app.integrator = Integrator::RK4;
  char err[256] = {0};
  /* Ensure enough samples to actually fill the diagram. For MAPS the
   * orbit-diagram needs many iterates per slice and many slices, or it
   * renders as a sparse single curve that breaks into scattered dots in
   * the chaotic region (the "horizontal dots" symptom). Enforce sensible
   * minimums here so the picture is right regardless of how the sweep was
   * configured (preset dropdown, leftover ODE values, etc.). Larger
   * user-set values are kept. */
  int eff_slices = std::max(2, app.bif_slices);
  int eff_keep = app.bif_keep;
  int eff_discard = app.bif_discard;
  if (app.mode == SystemMode::Map) {
    if (eff_slices < 600) eff_slices = 600;
    if (eff_keep < 200) eff_keep = 200;
    if (eff_discard < 500) eff_discard = 500;
  } else {
    if (eff_keep < 1000) eff_keep = 1000;     /* ODEs: need many steps to catch maxima */
    if (eff_discard < 2000) eff_discard = 2000;
    if (eff_slices < 300) eff_slices = 300;
  }
  const int slices = eff_slices;
  const size_t dim = app.state_names.size();
  const double leps = app.lyapunov_epsilon > 0 ? app.lyapunov_epsilon : 1e-8;
  for (int i = 0; i < slices; ++i) {
    const double u = static_cast<double>(i) / static_cast<double>(slices - 1);
    param->value = app.bif_start + u * (app.bif_end - app.bif_start);
    /* CRITICAL: the IR evaluator reads app.param_values, not param->value.
     * Without this sync the swept parameter never reaches the stepping
     * function, so every slice produces the SAME orbit -> the diagram comes
     * out as horizontal bands instead of the bifurcation tree. */
    sync_param_values(app);
    State s = app.start;
    for (int j = 0; j < eff_discard; ++j) {
      State next{};
      if (!step_state(app, s, &next, err, sizeof(err))) { app.analysis_message = err; param->value = old_value; app.integrator = saved_integrator; arena_destroy(&tmp_arena); return; }
      s = next;
    }
    /* Lyapunov accumulation via a renormalized shadow orbit, run over the
     * same kept iterations that build the orbit diagram. */
    double ly_sum = 0.0;
    long ly_n = 0;
    /* ODE local-maxima detection window (for the bifurcation section) */
    double obs_w0 = 0, obs_w1 = 0, obs_w2 = 0;
    int obs_have = 0;
    std::vector<double> kept_vals; /* for map period detection */
    kept_vals.reserve(eff_keep);
    State shadow = s;
    if (app.bif_compute_lyapunov && dim > 0) {
      resize_state(shadow, dim);
      shadow.v[0] += leps;
    }
    for (int j = 0; j < eff_keep; ++j) {
      State next{};
      if (!step_state(app, s, &next, err, sizeof(err))) { app.analysis_message = err; param->value = old_value; app.integrator = saved_integrator; arena_destroy(&tmp_arena); return; }
      if (app.bif_compute_lyapunov && dim > 0) {
        State nshadow{};
        if (step_state(app, shadow, &nshadow, err, sizeof(err))) {
          double d2 = 0.0;
          for (size_t k = 0; k < dim; ++k) {
            const double d = state_at(nshadow, k) - state_at(next, k);
            d2 += d * d;
          }
          const double dist = std::sqrt(d2);
          if (dist > 1e-300 && std::isfinite(dist)) {
            ly_sum += std::log(dist / leps);
            ly_n += 1;
            /* renormalize shadow back to leps along the separation dir */
            const double scale = leps / dist;
            shadow = make_state_like(dim, next.t);
            for (size_t k = 0; k < dim; ++k) {
              const double d = state_at(nshadow, k) - state_at(next, k);
              set_state_at(shadow, k, state_at(next, k) + d * scale);
            }
          }
        }
      }
      s = next;
      double val = 0.0;
      if (!eval_expr_at(app, obs, s, &val, err, sizeof(err))) continue;
      if (app.mode == SystemMode::Map) {
        /* maps: every iterate is an attractor sample (the classic tree) */
        app.bifurcation_points.push_back(Point2{param->value, val});
        kept_vals.push_back(val);
      } else {
        /* ODEs: raw time-steps trace the whole orbit and smear into a band
         * (or a single horizontal line at a fixed point). The standard
         * bifurcation diagram for a flow records a Poincare-style section
         * — here, the LOCAL MAXIMA of the observable. We keep a 3-sample
         * window and emit the middle when it's a peak. */
        obs_w0 = obs_w1; obs_w1 = obs_w2; obs_w2 = val;
        obs_have = std::min(obs_have + 1, 3);
        if (obs_have == 3 && obs_w1 > obs_w0 && obs_w1 >= obs_w2)
          app.bifurcation_points.push_back(Point2{param->value, obs_w1});
      }
    }
    if (app.bif_compute_lyapunov && ly_n > 0)
      app.bifurcation_lyapunov.push_back(Point2{param->value, ly_sum / static_cast<double>(ly_n)});

    /* period detection (maps): the period is the number of distinct values
     * in the settled orbit, capped; 0 means chaotic/high-period. */
    if (app.mode == SystemMode::Map && !kept_vals.empty()) {
      const double tol = 1e-4;
      int period = 0;
      const int maxP = 16;
      const int m = (int)kept_vals.size();
      /* compare the last value against earlier ones to find the cycle length */
      for (int pgap = 1; pgap <= maxP && pgap < m; ++pgap) {
        bool match = true;
        for (int q = 0; q < pgap && (m - 1 - q - pgap) >= 0; ++q) {
          if (std::fabs(kept_vals[m - 1 - q] - kept_vals[m - 1 - q - pgap]) > tol) { match = false; break; }
        }
        if (match) { period = pgap; break; }
      }
      app.bifurcation_period.push_back(Point2{param->value, (double)period});
    }
  }
  param->value = old_value;
  app.integrator = saved_integrator; /* restore the user's solver choice */
  arena_destroy(&tmp_arena);
  app.bif_view_valid = false; /* refit the diagram view to the new sweep */
  app.analysis_message = "bifurcation scan completed";
}

/* PHASE D step 2 (foundation): measure the limit-cycle period & amplitude
 * of the current trajectory. Pulls the y-axis state component from the
 * (settled) history and runs the validated detector. */
void run_limit_cycle(AppState &app) {
  app.lc_ready = false;
  if (app.mode != SystemMode::ODE) { app.lc_msg = "limit cycles are an ODE feature"; return; }
  const size_t n = app.state_names.size();
  if (n < 1) { app.lc_msg = "no system"; return; }
  const size_t yi = (size_t)std::max(0, std::min(app.phase_y_index, (int)n - 1));

  /* gather the observable signal; use the second half (more settled) */
  std::vector<double> sig;
  sig.reserve(app.history.size());
  for (const State &st : app.history) {
    const double v = state_at(st, yi);
    if (std::isfinite(v)) sig.push_back(v);
  }
  if (sig.size() < 200) { app.lc_msg = "let the trajectory run longer, then retry"; return; }
  /* drop the first third as transient */
  const size_t drop = sig.size() / 3;
  std::vector<double> tail(sig.begin() + (long)drop, sig.end());

  const double dt = app.dt > 0 ? app.dt : 0.01;
  dynsys::analysis::LimitCycleResult R = dynsys::analysis::limit_cycle_period_amplitude(tail, dt);
  if (!R.ok) { app.lc_msg = R.message; return; }
  app.lc_ready = true;
  app.lc_period = R.period;
  app.lc_amplitude = R.amplitude;
  app.lc_msg = "ok";
}

/* PHASE D step 2: sweep a parameter and measure the limit cycle's period &
 * amplitude at each value, producing the continuation curves. Mirrors the
 * bifurcation sweep's safety (force fixed-step RK4, save/restore integrator
 * and parameter, sync param_values each slice). */
void run_lc_continuation(AppState &app) {
  app.lcc_pp.clear(); app.lcc_period.clear(); app.lcc_amp.clear();
  app.lcc_has_data = false; app.lcc_amp_max = 1.0; app.lcc_period_max = 1.0;

  if (app.mode != SystemMode::ODE) { app.lcc_msg = "limit-cycle continuation is for ODE systems"; return; }
  const size_t n = app.state_names.size();
  if (n < 1) { app.lcc_msg = "no system"; return; }
  if (app.params.empty()) { app.lcc_msg = "system has no parameter to sweep"; return; }

  AppState::Param *param = find_param(app, app.lcc_param);
  if (!param) { param = &app.params[0]; std::snprintf(app.lcc_param, sizeof(app.lcc_param), "%s", param->name.c_str()); }

  const size_t yi = (size_t)std::max(0, std::min(app.phase_y_index, (int)n - 1));
  const double pmin = std::min(app.lcc_p_min, app.lcc_p_max);
  const double pmax = std::max(app.lcc_p_min, app.lcc_p_max);
  const int slices = std::max(8, app.lcc_slices);

  /* force fixed-step RK4 for the sweep; restore after */
  const Integrator saved_integrator = app.integrator;
  if (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45)
    app.integrator = Integrator::RK4;
  const double old_value = param->value;
  const double dt = app.dt > 0 ? app.dt : 0.01;

  char err[128] = {0};
  const int settle_steps = 30000;   /* integrate to settle onto the attractor */
  const int sample_steps = 40000;   /* then sample the orbit */

  double amax = 0.0, tmax = 0.0;
  for (int i = 0; i < slices; ++i) {
    const double u = (double)i / (double)(slices - 1);
    param->value = pmin + u * (pmax - pmin);
    sync_param_values(app); /* CRITICAL: evaluator reads param_values */

    State s = app.start; resize_state(s, n);
    bool ok = true;
    for (int j = 0; j < settle_steps; ++j) {
      State nx{};
      if (!step_ode_state(app, s, &nx, err, sizeof(err))) { ok = false; break; }
      s = nx;
    }
    double period = std::nan(""), amp = 0.0;
    if (ok) {
      std::vector<double> sig; sig.reserve(sample_steps);
      for (int j = 0; j < sample_steps; ++j) {
        State nx{};
        if (!step_ode_state(app, s, &nx, err, sizeof(err))) { ok = false; break; }
        s = nx;
        const double v = state_at(s, yi);
        if (std::isfinite(v)) sig.push_back(v);
      }
      if (ok && sig.size() > 100) {
        dynsys::analysis::LimitCycleResult R = dynsys::analysis::limit_cycle_period_amplitude(sig, dt);
        if (R.ok) { period = R.period; amp = R.amplitude; }
        else { amp = 0.0; period = std::nan(""); } /* settled to a fixed point */
      }
    }
    app.lcc_pp.push_back(param->value);
    app.lcc_period.push_back(period);
    app.lcc_amp.push_back(amp);
    if (std::isfinite(amp) && amp > amax) amax = amp;
    if (std::isfinite(period) && period > tmax) tmax = period;
  }

  param->value = old_value;
  sync_param_values(app);
  app.integrator = saved_integrator;

  app.lcc_amp_max = amax > 0 ? amax : 1.0;
  app.lcc_period_max = tmax > 0 ? tmax : 1.0;
  app.lcc_has_data = !app.lcc_pp.empty();
  app.lcc_msg = app.lcc_has_data ? "ok" : "no data";
}

/* PHASE B/C: estimate the box-counting fractal dimension of the on-screen
 * set in the current 2D plane. For a MAP we iterate from the start state
 * (after a transient) to sample the attractor densely; for an ODE we use
 * the live trajectory history (already the settled attractor). */
void run_box_dimension(AppState &app) {
  app.boxdim_ready = false;
  const size_t n = app.state_names.size();
  if (n < 1) { app.boxdim_msg = "no system"; return; }
  const size_t ix = (size_t)std::max(0, std::min(app.phase_x_index, (int)n - 1));
  const size_t iy = (size_t)std::max(0, std::min(app.phase_y_index, (int)n - 1));

  std::vector<double> xs, ys;
  char err[128] = {0};

  if (app.mode == SystemMode::Map) {
    /* iterate the map to sample the attractor */
    State s = app.start; resize_state(s, n);
    const int transient = 1000;
    const int samples = 200000;
    for (int i = 0; i < transient; ++i) {
      State nx{}; if (!step_map_state(app, s, &nx, err, sizeof(err))) { app.boxdim_msg = err; return; }
      s = nx;
    }
    xs.reserve(samples); ys.reserve(samples);
    for (int i = 0; i < samples; ++i) {
      State nx{}; if (!step_map_state(app, s, &nx, err, sizeof(err))) break;
      s = nx;
      const double x = state_at(s, ix), y = state_at(s, iy);
      if (std::isfinite(x) && std::isfinite(y)) { xs.push_back(x); ys.push_back(y); }
    }
  } else {
    /* ODE: use the trajectory history (the visible orbit/attractor) */
    for (const State &st : app.history) {
      const double x = state_at(st, ix), y = state_at(st, iy);
      if (std::isfinite(x) && std::isfinite(y)) { xs.push_back(x); ys.push_back(y); }
    }
    if (xs.size() < 500) {
      app.boxdim_msg = "need a longer trajectory — let the ODE run, then retry";
      return;
    }
  }

  dynsys::analysis::BoxCountResult R = dynsys::analysis::box_counting_dimension(xs, ys, 12);
  app.boxdim_n_points = (long)xs.size();
  if (!R.ok) { app.boxdim_msg = R.message; return; }
  app.boxdim_ready = true;
  app.boxdim_value = R.dimension;
  app.boxdim_r2 = R.r_squared;
  app.boxdim_msg = "ok";
}

/* PHASE B: compute the full Lyapunov spectrum + Kaplan-Yorke dimension
 * for the current system at its current parameters, starting from the
 * current state. Builds a Model whose vector_field is the ODE rhs (for
 * ODEs) or the map's next-state (for maps), with the Jacobian supplied by
 * forward-mode AD over the appropriate program set. */
void run_lyapunov_spectrum(AppState &app) {
  const size_t n = app.state_names.size();
  if (n == 0) { app.lyap_spectrum_msg = "no state variables"; app.lyap_spectrum_ready = false; return; }
  const bool is_map = (app.mode == SystemMode::Map);
  const auto &progs = is_map ? app.next_equation_programs : app.equation_programs;
  if (progs.size() != n) { app.lyap_spectrum_msg = "equations not compiled"; app.lyap_spectrum_ready = false; return; }

  dynsys::analysis::Model model;
  model.n = n;
  model.vector_field = [&app, n, is_map](const double *x, double, double *f_out, std::string *err) -> bool {
    State s = make_state_like(n, app.current.t);
    for (size_t i = 0; i < n; ++i) set_state_at(s, i, x[i]);
    char buf[256] = {0};
    if (is_map) {
      State nx{};
      if (!step_state(app, s, &nx, buf, sizeof(buf))) { if (err) *err = buf; return false; }
      for (size_t i = 0; i < n; ++i) f_out[i] = state_at(nx, i);
    } else {
      State d{};
      if (!eval_rhs(app, s, &d, buf, sizeof(buf))) { if (err) *err = buf; return false; }
      for (size_t i = 0; i < n; ++i) f_out[i] = state_at(d, i);
    }
    return true;
  };
  model.jacobian_x = [&app, n, &progs](const double *x, double, double *jac_out, std::string *err) -> bool {
    State s = make_state_like(n, app.current.t);
    for (size_t i = 0; i < n; ++i) set_state_at(s, i, x[i]);
    dynsys::ir::RunContext rc;
    rc.state = s.v.data(); rc.n_state = s.v.size(); rc.t = s.t;
    rc.params = app.param_values.data(); rc.n_params = app.param_values.size();
    rc.defs = app.definition_programs.data(); rc.n_defs = app.definition_programs.size();
    char buf[256] = {0};
    for (size_t col = 0; col < n; ++col) {
      dynsys::ir::DualSeed seed{dynsys::ir::DualSeed::Kind::State, col};
      for (size_t row = 0; row < n; ++row) {
        double v = 0.0, d = 0.0;
        if (!dynsys::ir::run_dual(progs[row], rc, seed, app.ad_scratch, &v, &d, buf, sizeof(buf))) {
          if (err) *err = buf;
          return false;
        }
        jac_out[row * n + col] = d;
      }
    }
    return true;
  };

  std::vector<double> x0(n);
  for (size_t i = 0; i < n; ++i) x0[i] = state_at(app.current, i);

  dynsys::analysis::LyapunovOptions opt;
  opt.is_map = is_map;
  opt.dt = app.dt > 0 ? app.dt : 0.01;
  opt.transient = std::max(0, app.lyap_spec_transient);
  opt.steps = std::max(100, app.lyap_spec_steps);
  opt.reorth_every = 1;

  dynsys::analysis::LyapunovResult r =
      dynsys::analysis::lyapunov_spectrum(model, x0, 0.0, opt);
  app.lyap_spectrum = r.exponents;
  app.lyap_kaplan_yorke = r.kaplan_yorke;
  app.lyap_spectrum_sum = r.sum;
  app.lyap_spectrum_ready = r.ok;
  app.lyap_spectrum_msg = r.ok ? "spectrum computed" : r.message;
}

/* PHASE6-UI: the 3D view as an ImGui panel. Renders the scene to the
 * offscreen texture at the panel's pixel size, shows it as an image,
 * and handles drag-to-rotate and wheel-to-zoom while the image is
 * hovered — so the 3D controls live with the 3D view instead of being
 * global window-background interactions. */
[[maybe_unused]] void draw_scene_panel(AppState &app) {
  ImGui::Checkbox("axes", &app.show_axes);
  ImGui::SameLine();
  if (ImGui::Checkbox("orthographic", &app.orthographic_3d)) update_projection(app);
  ImGui::SameLine();
  if (ImGui::SmallButton("reset view")) {
    app.angle_x = 0.0f; app.angle_y = 0.0f; app.zoom = 1.0f;
    app.scene_scale = 0.05f; app.camera_distance = 50.0f;
    app.center_x = app.center_y = app.center_z = 0.0f;
    update_projection(app);
  }

  ImVec2 avail = ImGui::GetContentRegionAvail();
  const int pw = std::max(64, static_cast<int>(avail.x));
  const int ph = std::max(64, static_cast<int>(avail.y));
  render_scene_to_fbo(app, pw, ph);

  const ImVec2 img_pos = ImGui::GetCursorScreenPos();
  /* Flip V so the GL texture (origin bottom-left) shows upright. */
  /* ImTextureID is an integer in this ImGui build and a pointer in
   * others. A C-style cast through uintptr_t compiles for both (it
   * selects the valid integer- or pointer-conversion); static_cast and
   * reinterpret_cast each only work for one of the two configurations.
   * This is the idiom ImGui's own examples use for GL texture names. */
  ImGui::Image((ImTextureID)(uintptr_t)app.scene_tex,
               ImVec2(static_cast<float>(pw), static_cast<float>(ph)),
               ImVec2(0, 1), ImVec2(1, 0));

  const bool hovered = ImGui::IsItemHovered();
  app.scene_hovered = hovered;
  ImGuiIO &io = ImGui::GetIO();
  if (hovered) {
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
      app.angle_y += io.MouseDelta.x * 0.4f;
      app.angle_x += io.MouseDelta.y * 0.4f;
    }
    if (io.MouseWheel != 0.0f) {
      app.zoom *= std::pow(1.15f, (float)clamped_wheel(io.MouseWheel));
      app.zoom = std::max(0.05f, std::min(app.zoom, 100.0f));
      update_projection(app);
    }
  }
  (void)img_pos;
  ImGui::TextDisabled(
      "drag: rotate   wheel: zoom   (this is the real 3D view)");
}

/* PHASE C+: save the current framebuffer (the full-window plot) to a PNG.
 * Reads the GL backbuffer, flips it (GL origin is bottom-left), and writes
 * a timestamped file next to the binary. Returns the path or "". */
std::string capture_screenshot_png(AppState &app) {
  const int w = app.window_width, h = app.window_height;
  if (w <= 0 || h <= 0) return "";
  std::vector<unsigned char> buf((size_t)w * h * 4);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  /* read the buffer we just rendered into (called post-RenderDrawData,
   * pre-SwapBuffers, so the back buffer holds the complete frame) */
  glReadBuffer(0x0405 /*GL_BACK*/);
  glReadPixels(0, 0, w, h, GL_RGBA, /*GL_UNSIGNED_BYTE*/ 0x1401, buf.data());

  /* When capturing a clean (chrome-suppressed) frame, save the FULL window so
   * nothing is hidden and the image is full-size. Otherwise fall back to
   * cropping the toolbar/panel out of a normal frame. */
  int cx0 = 0, cy0_top = 0, cx1 = w, cy1_top = h;
  if (!app.capturing_clean) {
    cy0_top = (int)std::lround(app.window_toolbar_h);
    const bool panel_shown = app.show_side_panel && app.panel_x1 > app.panel_x0 + 1.0f;
    if (panel_shown && app.panel_x0 <= 2.0f) cx0 = std::max(cx0, (int)std::lround(app.panel_x1));
  }
  cx0 = std::max(0, std::min(cx0, w - 1));
  cx1 = std::max(cx0 + 1, std::min(cx1, w));
  cy0_top = std::max(0, std::min(cy0_top, h - 1));
  cy1_top = std::max(cy0_top + 1, std::min(cy1_top, h));
  const int cw = cx1 - cx0, ch = cy1_top - cy0_top;

  std::vector<unsigned char> rgb((size_t)cw * ch * 3);
  for (int y = 0; y < ch; ++y) {
    /* destination row y corresponds to window row (cy0_top + y) from the top;
     * in the bottom-left framebuffer that is row (h - 1 - (cy0_top + y)). */
    const int srcrow = h - 1 - (cy0_top + y);
    const unsigned char *src = &buf[((size_t)srcrow * w + cx0) * 4];
    unsigned char *dst = &rgb[(size_t)y * cw * 3];
    for (int x = 0; x < cw; ++x) {
      dst[x * 3 + 0] = src[x * 4 + 0];
      dst[x * 3 + 1] = src[x * 4 + 1];
      dst[x * 3 + 2] = src[x * 4 + 2];
    }
  }
  const int outw = cw, outh = ch;
  char name[128];
  std::time_t t = std::time(nullptr);
  std::tm *lt = std::localtime(&t);
  int yr = lt ? lt->tm_year + 1900 : 0;
  if (yr < 2000 || yr > 3000) {
    /* the system clock is implausible; fall back to a counter so we still
     * produce a uniquely-named file rather than a misleading date */
    static int seq = 0; ++seq;
    std::snprintf(name, sizeof(name), "dynsys_capture_%03d.png", seq);
  } else {
    std::snprintf(name, sizeof(name), "dynsys_%04d%02d%02d_%02d%02d%02d.png",
                  yr, lt->tm_mon + 1, lt->tm_mday, lt->tm_hour, lt->tm_min, lt->tm_sec);
  }
  if (!png_write_rgb(name, rgb.data(), outw, outh)) return "";
  /* report an absolute path so the user can actually find the file */
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd))) return std::string(cwd) + "/" + name;
  return std::string(name);
}

/* PHASE6-UI: top toolbar — view switch, run controls, key stats. Fixed,
 * full width, non-movable, drawn on top of the background plot. */
void draw_top_toolbar(AppState &app) {
  const float bar_h = 0.0f; (void)bar_h;
  ImGui::SetNextWindowPos(ImVec2(0, 0));
  ImGui::SetNextWindowSize(ImVec2(static_cast<float>(app.window_width), 0));
  ImGuiWindowFlags flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
                           ImGuiWindowFlags_NoTitleBar |
                           ImGuiWindowFlags_NoScrollbar |
                           ImGuiWindowFlags_NoCollapse |
                           ImGuiWindowFlags_NoBringToFrontOnFocus;
  ImGui::Begin("##toolbar", nullptr, flags);

  /* Application title. */
  ImGui::TextColored(ImVec4(0.55f, 0.95f, 0.65f, 1.0f), "Dynsys");
  ImGui::SameLine();
  if (ImGui::SmallButton(app.capture_request ? "Saving..." : "Save PNG")) {
    /* Defer the actual pixel capture to AFTER this frame is rendered (see the
     * main loop). Capturing here would read the previous/empty back buffer. */
    app.capture_request = true;
  }
  ImGui::SameLine();
  if (app.screenshot_msg_timer > 0) {
    --app.screenshot_msg_timer;
    ImGui::TextColored(ImVec4(0.6f, 0.95f, 0.7f, 1.0f), "%s", app.screenshot_msg.c_str());
    ImGui::SameLine();
  }
  ImGui::TextUnformatted("|"); ImGui::SameLine();

  if (ImGui::SmallButton("Reset view")) reset_current_view(app);
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("Reset zoom/pan (or camera) for the current view to a sane framing.");
  ImGui::SameLine();
  if (ImGui::SmallButton("Reset all")) {
    reset_current_view(app);
    reset_simulation(app);
    /* also drop every 2D view back to auto framing so nothing stays zoomed out */
    app.phase_auto_bounds = true; app.phase_bounds_valid = false;
    app.bif_view_valid = false; app.cont_view_valid = false;
    app.fractal_dirty = true; app.basin_dirty = true; app.scan_dirty = true;
    app.scan_view_init = false;
  }
  if (ImGui::IsItemHovered())
    ImGui::SetTooltip("Reset the simulation and all view framings.");
  ImGui::SameLine();
  ImGui::TextUnformatted("|"); ImGui::SameLine();

  if (ImGui::Button(app.show_side_panel ? "<< panel" : "panel >>"))
    app.show_side_panel = !app.show_side_panel;
  ImGui::SameLine();
  ImGui::TextUnformatted("|"); ImGui::SameLine();

  /* View switch — offer only the views valid for this system. */
  ensure_valid_view(app);
  auto view_radio = [&](const char *label, AppState::ActiveView v) {
    const bool ok = view_valid(app, v);
    if (!ok) ImGui::BeginDisabled();
    if (ImGui::RadioButton(label, app.active_view == v) && ok) app.active_view = v;
    if (!ok) ImGui::EndDisabled();
    /* Explain *why* a view is unavailable (or what it shows) on hover, so
     * the applicability rules are visible instead of a mystery grey-out.
     * IsItemHovered with AllowWhenDisabled still fires for disabled items. */
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
      if (ok) ImGui::SetTooltip("%s", view_requirement(v));
      else ImGui::SetTooltip("unavailable — needs %s", view_requirement(v));
    }
    /* Wrap to the next line instead of overflowing off the right edge of
     * the (non-scrolling) toolbar — otherwise later views (continuation,
     * IFS, ...) become invisible on narrower windows. */
    const float avail = ImGui::GetWindowContentRegionMax().x;
    const float next_w = ImGui::CalcTextSize(label).x + ImGui::GetStyle().FramePadding.x * 4.0f + 40.0f;
    if (ImGui::GetItemRectMax().x + next_w < avail) ImGui::SameLine();
  };
  view_radio("1D", AppState::ActiveView::Line1D);
  view_radio("2D phase", AppState::ActiveView::Phase2D);
  view_radio("3D", AppState::ActiveView::Scene3D);
  view_radio("bifurcation", AppState::ActiveView::Bifurcation);
  view_radio("fractal", AppState::ActiveView::Fractal);
  view_radio("3D bridge", AppState::ActiveView::Scene3DBridge);
  view_radio("basins", AppState::ActiveView::Basins);
  view_radio("param scan", AppState::ActiveView::ParamScan2D);
  view_radio("continuation", AppState::ActiveView::Continuation);
  view_radio("limit cycle", AppState::ActiveView::LimitCycle);
  ImGui::Separator();

  if (ImGui::Button(app.paused ? "Play" : "Pause")) app.paused = !app.paused;
  ImGui::SameLine();
  if (ImGui::Button("Reset")) reset_simulation(app);
  ImGui::SameLine();
  if (app.active_view == AppState::ActiveView::Phase2D) {
    if (ImGui::Button("Clear orbits")) { app.phase_trajectories.clear(); app.separatrix_curves.clear(); }
    ImGui::SameLine();
  }
  if (app.active_view == AppState::ActiveView::Bifurcation && !app.params.empty()) {
    ImGui::SetNextItemWidth(90);
    if (ImGui::BeginCombo("##bifp", app.bif_param)) {
      for (const auto &pr : app.params)
        if (ImGui::Selectable(pr.name.c_str(), pr.name == app.bif_param)) {
          std::snprintf(app.bif_param, sizeof(app.bif_param), "%s", pr.name.c_str());
          if (pr.has_range) { app.bif_start = pr.min_value; app.bif_end = pr.max_value; app.bif_view_valid = false; }
        }
      ImGui::EndCombo();
    }
    ImGui::SameLine(); ImGui::SetNextItemWidth(64);
    ImGui::InputDouble("##bifs", &app.bif_start, 0, 0, "%.3g"); ImGui::SameLine();
    ImGui::SetNextItemWidth(64);
    ImGui::InputDouble("##bife", &app.bif_end, 0, 0, "%.3g"); ImGui::SameLine();
    if (ImGui::Button("Run bifurcation")) run_bifurcation(app);
    ImGui::SameLine();
  }
  if (app.active_view == AppState::ActiveView::Continuation && !app.params.empty()) {
    if (app.cont_param[0] == '\0') std::snprintf(app.cont_param, sizeof(app.cont_param), "%s", app.params[0].name.c_str());
    ImGui::SetNextItemWidth(90);
    if (ImGui::BeginCombo("##contp", app.cont_param)) {
      for (const auto &pr : app.params)
        if (ImGui::Selectable(pr.name.c_str(), pr.name == app.cont_param)) {
          std::snprintf(app.cont_param, sizeof(app.cont_param), "%s", pr.name.c_str());
          if (pr.has_range) { app.cont_p_min = pr.min_value; app.cont_p_max = pr.max_value; }
        }
      ImGui::EndCombo();
    }
    ImGui::SameLine(); ImGui::SetNextItemWidth(56);
    ImGui::InputDouble("##contpmin", &app.cont_p_min, 0, 0, "%.3g"); ImGui::SameLine();
    ImGui::SetNextItemWidth(56);
    ImGui::InputDouble("##contpmax", &app.cont_p_max, 0, 0, "%.3g"); ImGui::SameLine();
    ImGui::SetNextItemWidth(70);
    if (ImGui::BeginCombo("##conty", app.state_names.empty() ? "y" : app.state_names[std::min((size_t)app.cont_y_index, app.state_names.size()-1)].c_str())) {
      for (int i = 0; i < (int)app.state_names.size(); ++i)
        if (ImGui::Selectable(app.state_names[i].c_str(), i == app.cont_y_index)) { app.cont_y_index = i; app.cont_view_valid = false; }
      ImGui::EndCombo();
    }
    ImGui::SameLine();
    if (ImGui::Button("Continue")) run_continuation(app);
    ImGui::SameLine();
  }
  if (app.active_view == AppState::ActiveView::Phase2D && app.mode == SystemMode::IFS) {
    ImGui::Text("IFS model: %zu maps", app.ifs_maps.size());
    ImGui::SameLine(); ImGui::SetNextItemWidth(120);
    int it = (int)app.ifs_iterations;
    if (ImGui::SliderInt("iters##ifs", &it, 50000, 2000000)) { app.ifs_iterations = it; app.ifs_dirty = true; }
    ImGui::SameLine();
  }
  if (app.active_view == AppState::ActiveView::LimitCycle && !app.params.empty()) {
    if (app.lcc_param[0] == '\0') std::snprintf(app.lcc_param, sizeof(app.lcc_param), "%s", app.params[0].name.c_str());
    ImGui::SetNextItemWidth(90);
    if (ImGui::BeginCombo("##lccp", app.lcc_param)) {
      for (const auto &pr : app.params)
        if (ImGui::Selectable(pr.name.c_str(), pr.name == app.lcc_param)) {
          std::snprintf(app.lcc_param, sizeof(app.lcc_param), "%s", pr.name.c_str());
          if (pr.has_range) { app.lcc_p_min = pr.min_value; app.lcc_p_max = pr.max_value; }
        }
      ImGui::EndCombo();
    }
    ImGui::SameLine(); ImGui::SetNextItemWidth(54);
    ImGui::InputDouble("##lccmin", &app.lcc_p_min, 0, 0, "%.3g"); ImGui::SameLine();
    ImGui::SetNextItemWidth(54);
    ImGui::InputDouble("##lccmax", &app.lcc_p_max, 0, 0, "%.3g"); ImGui::SameLine();
    ImGui::SetNextItemWidth(90);
    ImGui::SliderInt("slices##lcc", &app.lcc_slices, 16, 240); ImGui::SameLine();
    if (ImGui::Button("Sweep")) run_lc_continuation(app);
    ImGui::SameLine();
  }
  ImGui::TextDisabled("  %d FPS | %s | %dD%s", app.fps, mode_name(app.mode),
                      effective_dimension(app),
                      app.active_view == AppState::ActiveView::Phase2D
                          ? "  | click: add orbit, drag: pan, wheel/+/-: zoom"
                          : "  | drag: rotate, wheel/+/-: zoom");
  app.window_toolbar_h = ImGui::GetWindowHeight();
  ImGui::End();
}

void draw_gui(AppState &app) {
  /* PHASE-A: render the active full-area view behind the panels. (3D is
   * already drawn straight to the backbuffer in the main loop.) */
  switch (app.active_view) {
    case AppState::ActiveView::Line1D:      render_line1d_background(app); break;
    case AppState::ActiveView::Phase2D:
      /* An IFS system is a planar object; it displays in the phase view via
       * the chaos game (just as an ODE shows its trajectory here). */
      if (app.mode == SystemMode::IFS) render_ifs_background(app);
      else render_phase_background(app);
      break;
    case AppState::ActiveView::Bifurcation: render_bifurcation_background(app); break;
    case AppState::ActiveView::Fractal:     render_fractal_background(app); break;
    case AppState::ActiveView::Scene3DBridge: /* drawn to backbuffer in main loop */ break;
    case AppState::ActiveView::Basins:      render_basin_background(app); break;
    case AppState::ActiveView::ParamScan2D: render_scan_background(app); break;
    case AppState::ActiveView::Continuation: render_continuation_background(app); break;
    case AppState::ActiveView::IFS: render_ifs_background(app); break;
    case AppState::ActiveView::LimitCycle: render_lc_continuation_background(app); break;
    case AppState::ActiveView::Scene3D:     /* drawn in main loop */ break;
  }

  /* During a clean capture frame we suppress all chrome (top toolbar + side
   * controls panel) so the saved PNG shows the FULL plot at full size, with
   * nothing hidden behind the menus. */
  if (app.capturing_clean) return;

  draw_top_toolbar(app);

  if (!app.show_side_panel) {
    /* Side panel hidden: still need the analysis plots window? No — the
     * background is the plot. Show nothing else. */
    return;
  }

  const float top = app.window_toolbar_h > 0 ? app.window_toolbar_h : 32.0f;
  ImGui::SetNextWindowPos(ImVec2(0, top), ImGuiCond_Always);
  ImGui::SetNextWindowSize(ImVec2(440, static_cast<float>(app.window_height) - top), ImGuiCond_FirstUseEver);
  ImGuiWindowFlags side_flags = ImGuiWindowFlags_NoMove;
  ImGui::Begin("controls", nullptr, side_flags);
  /* Record the panel's on-screen rect so Save PNG can crop the controls out
   * and save only the plotted system. */
  {
    ImVec2 pp = ImGui::GetWindowPos();
    ImVec2 ps = ImGui::GetWindowSize();
    app.panel_x0 = pp.x; app.panel_y0 = pp.y;
    app.panel_x1 = pp.x + ps.x; app.panel_y1 = pp.y + ps.y;
  }

  ImGui::Text("FPS: %d | mode: %s | integrator: %s | dim: %zu", app.fps, mode_name(app.mode), integrator_name(app.integrator), app.state_names.size());
  ImGui::Checkbox("auto performance governor", &app.perf_governor_enabled);
  if (app.perf_governor_enabled) {
    ImGui::SameLine();
    if (app.perf_throttle_active)
      ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f),
                         "throttling %.0f%% (%.0f ms/frame) — easing load to stay responsive",
                         app.perf_throttle * 100.0, app.perf_frame_ms);
    else
      ImGui::TextColored(ImVec4(0.5f, 0.85f, 0.6f, 1.0f), "headroom OK (%.0f ms/frame)", app.perf_frame_ms);
  }
  std::string state_line = "t " + std::to_string(app.current.t);
  for (size_t i = 0; i < app.state_names.size(); ++i) {
    char buf[96];
    std::snprintf(buf, sizeof(buf), ", %s %.5g", app.state_names[i].c_str(), state_at(app.current, i));
    state_line += buf;
  }
  ImGui::TextWrapped("%s", state_line.c_str());
  if (app.lyapunov_enabled) {
    ImGui::Text("largest Lyapunov estimate: %.6g (%lld samples)", app.lyapunov_estimate, app.lyapunov_samples);
  }
  ImGui::Separator();

  /* PHASE-A: group the (previously 12) control sections into a few tabs
   * so nothing is buried in one long scroll. Within each tab the old
   * collapsing headers are kept. */
  if (ImGui::BeginTabBar("control_tabs")) {
  if (ImGui::BeginTabItem("Setup")) {

  /* PHASE6-UI: phase-plane controls (axes, layers, seeding) — the plot
   * itself is the window background; this drives it. Shown in 2D view. */
  if (app.active_view == AppState::ActiveView::Phase2D &&
      app.state_names.size() >= 2) {
    if (ImGui::CollapsingHeader("Phase plane", ImGuiTreeNodeFlags_DefaultOpen)) {
      const int ns = static_cast<int>(app.state_names.size());
      app.phase_x_index = std::max(0, std::min(app.phase_x_index, ns - 1));
      app.phase_y_index = std::max(0, std::min(app.phase_y_index, ns - 1));
      if (ImGui::BeginCombo("x axis", app.state_names[(size_t)app.phase_x_index].c_str())) {
        for (int i = 0; i < ns; ++i)
          if (ImGui::Selectable(app.state_names[(size_t)i].c_str(), i == app.phase_x_index)) app.phase_x_index = i;
        ImGui::EndCombo();
      }
      if (ImGui::BeginCombo("y axis", app.state_names[(size_t)app.phase_y_index].c_str())) {
        for (int i = 0; i < ns; ++i)
          if (ImGui::Selectable(app.state_names[(size_t)i].c_str(), i == app.phase_y_index)) app.phase_y_index = i;
        ImGui::EndCombo();
      }
      ImGui::Checkbox("vector field", &app.phase_show_vector_field); ImGui::SameLine();
      ImGui::Checkbox("nullclines", &app.phase_show_nullclines);
      ImGui::Checkbox("nullcline arrows", &app.phase_show_nullcline_arrows);
      ImGui::Checkbox("trajectories", &app.phase_show_trajectory); ImGui::SameLine();
      ImGui::Checkbox("auto fixed points", &app.phase_auto_fixed_points);
      ImGui::Checkbox("manifolds", &app.phase_show_manifolds); ImGui::SameLine();
      ImGui::Checkbox("separatrices", &app.phase_show_separatrices);
      ImGui::Checkbox("equal aspect", &app.phase_equal_aspect); ImGui::SameLine();
      ImGui::Checkbox("normalize arrows", &app.phase_normalize_vectors);
      ImGui::Checkbox("auto-fit bounds", &app.phase_auto_bounds);
      ImGui::SliderInt("grid density", &app.phase_grid, 8, 48);
      if (app.phase_auto_fixed_points && app.mode == SystemMode::ODE)
        ImGui::TextDisabled("equilibria in view: %zu (auto-classified)", app.phase_fixed_points.size());
      ImGui::Text("orbits: %zu", app.phase_trajectories.size());
      /* Seeding uses the current view bounds for the active axis pair. */
      PlotBounds vb = current_phase_bounds(app, (size_t)app.phase_x_index, (size_t)app.phase_y_index);
      if (ImGui::Button("Seed 5x5 grid")) seed_phase_grid(app, vb, (size_t)app.phase_x_index, (size_t)app.phase_y_index, 5, 5);
      ImGui::SameLine();
      if (ImGui::Button("Seed circle")) seed_phase_circle(app, vb, (size_t)app.phase_x_index, (size_t)app.phase_y_index, 16);
      ImGui::SameLine();
      if (ImGui::Button("Clear orbits")) { app.phase_trajectories.clear(); app.separatrix_curves.clear(); }

      /* PHASE6-UI: add an orbit from EXACT typed initial conditions. */
      if (ImGui::TreeNode("Add custom orbit")) {
        if (app.custom_orbit_ic.size() != app.state_names.size()) {
          app.custom_orbit_ic.assign(app.state_names.size(), 0.0f);
          for (size_t i = 0; i < app.state_names.size(); ++i)
            app.custom_orbit_ic[i] = static_cast<float>(state_at(app.current, i));
        }
        for (size_t i = 0; i < app.state_names.size(); ++i) {
          ImGui::PushID(static_cast<int>(i));
          ImGui::InputFloat(app.state_names[i].c_str(), &app.custom_orbit_ic[i], 0.0f, 0.0f, "%.5g");
          ImGui::PopID();
        }
        if (ImGui::Button("Add orbit at these values")) {
          State seed = app.current;
          resize_state(seed, app.state_names.size());
          for (size_t i = 0; i < app.state_names.size(); ++i)
            set_state_at(seed, i, app.custom_orbit_ic[i]);
          add_phase_trajectory(app, seed);
        }
        ImGui::SameLine();
        if (ImGui::Button("Use current state")) {
          for (size_t i = 0; i < app.state_names.size(); ++i)
            app.custom_orbit_ic[i] = static_cast<float>(state_at(app.current, i));
        }
        ImGui::TreePop();
      }
    }
  }

  if (ImGui::Button(app.paused ? "Resume" : "Pause")) app.paused = !app.paused;
  ImGui::SameLine();
  if (ImGui::Button("Reset")) reset_simulation(app);
  ImGui::SameLine();
  if (ImGui::Button("Quit")) {
    GLFWwindow *window = glfwGetCurrentContext();
    if (window != nullptr) glfwSetWindowShouldClose(window, GLFW_TRUE);
  }

  /* PHASE C: fractal controls (shown when the fractal view is active). */
  if (app.active_view == AppState::ActiveView::Fractal &&
      app.mode == SystemMode::Map &&
      (!app.params.empty() || app.state_names.size() >= 2)) {
    /* State-space (Julia) needs 2 state variables. For a 1-D map only
     * parameter-space is meaningful, so force it and disable the toggle. */
    const bool can_state_space = (app.state_names.size() >= 2);
    const bool can_buddha = (app.state_names.size() >= 2 && app.params.size() >= 2);
    if (!can_state_space) app.fractal_mode = AppState::FractalMode::ParameterSpace;
    if (ImGui::CollapsingHeader("Fractal", ImGuiTreeNodeFlags_DefaultOpen)) {
      int mode = (int)app.fractal_mode;
      const char *modes[] = {"parameter space (Mandelbrot-type)", "state space (Julia-type)",
                             "Buddhabrot (trajectory density)"};
      const int nmodes = can_buddha ? 3 : 2;
      if (!can_state_space) ImGui::BeginDisabled();
      if (ImGui::Combo("mode", &mode, modes, nmodes)) {
        app.fractal_mode = (AppState::FractalMode)mode;
        app.fractal_dirty = true;
      }
      if (!can_state_space) {
        ImGui::EndDisabled();
        if (ImGui::IsItemHovered())
          ImGui::SetTooltip("State-space (Julia) needs 2 state variables; this 1-D map only has parameter-space.");
      }
      if (app.fractal_mode == AppState::FractalMode::Buddhabrot)
        ImGui::TextWrapped("Buddhabrot: accumulates the orbits of ESCAPING points into a density "
                           "map. It builds up and sharpens the longer you leave it. Best on the "
                           "complex-quadratic (Mandelbrot) system.");
      else
      ImGui::TextWrapped("Mandelbrot = which parameter values keep the start orbit bounded. "
                         "Julia = which initial points stay bounded for fixed parameters. "
                         "The logistic bifurcation diagram is the real-axis slice of this.");
      if (app.fractal_mode == AppState::FractalMode::ParameterSpace && app.params.size() >= 1) {
        if (ImGui::BeginCombo("re axis = param", app.params[std::min((size_t)app.fractal_param_cx_index, app.params.size()-1)].name.c_str())) {
          for (int i = 0; i < (int)app.params.size(); ++i)
            if (ImGui::Selectable(app.params[i].name.c_str(), i == app.fractal_param_cx_index)) { app.fractal_param_cx_index = i; app.fractal_dirty = true; }
          ImGui::EndCombo();
        }
        if (app.params.size() >= 2 && ImGui::BeginCombo("im axis = param", app.params[std::min((size_t)app.fractal_param_cy_index, app.params.size()-1)].name.c_str())) {
          for (int i = 0; i < (int)app.params.size(); ++i)
            if (ImGui::Selectable(app.params[i].name.c_str(), i == app.fractal_param_cy_index)) { app.fractal_param_cy_index = i; app.fractal_dirty = true; }
          ImGui::EndCombo();
        }
        if (app.params.size() >= 2 && app.fractal_param_cx_index == app.fractal_param_cy_index)
          ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f),
                             "both axes are the same parameter — the image collapses to a diagonal; pick two different parameters.");
      }
      if (ImGui::SliderInt("max iterations", &app.fractal_max_iter, 20, 2000)) app.fractal_dirty = true;
      if (ImGui::InputDouble("escape radius", &app.fractal_escape_r, 0.5, 1.0, "%.3g")) app.fractal_dirty = true;
      if (ImGui::Checkbox("smooth coloring", &app.fractal_smooth)) app.fractal_dirty = true;
      if (ImGui::Button("Reset view")) {
        app.fractal_xmin = -2.5; app.fractal_xmax = 1.0;
        app.fractal_ymin = -1.5; app.fractal_ymax = 1.5;
        app.fractal_dirty = true;
      }
      ImGui::SameLine();
      if (ImGui::Button("Recompute")) app.fractal_dirty = true;
      ImGui::TextDisabled("Tip: load the 'Complex quadratic' preset for the classic Mandelbrot/Julia.");
    }
  }

  /* PHASE C+: controls for the 3D bridge scene. */
  if (app.active_view == AppState::ActiveView::Scene3DBridge) {
    if (ImGui::CollapsingHeader("3D bridge", ImGuiTreeNodeFlags_DefaultOpen)) {
      bool changed = false;
      /* The projection-solid bridge is only meaningful when the LOADED system
       * has a matching family (logistic/quadratic -> z^2+c, cubic map -> z^3+c,
       * sine map -> c*sin z). For systems with no such bridge (Henon, Lozi,
       * Ikeda, tent, Gauss, the ODE flows, ...) the button is disabled and we
       * say why. "current system" always works; "side-by-side" is the special
       * Mandelbrot/logistic illustration. */
      const int sys_fam = bridge_family_for_system(app);
      const bool proj_ok = (sys_fam >= 0);
      const char *fam_name = (sys_fam == 1) ? "cubic  z\u00b3+c"
                           : (sys_fam == 2) ? "sine  c\u00b7sin(z)"
                           : (sys_fam == 0) ? "quadratic  z\u00b2+c" : "(none)";
      /* keep the family in sync with the loaded system */
      if (sys_fam == 0) app.bridge_family = AppState::BridgeFamily::Quadratic;
      else if (sys_fam == 1) app.bridge_family = AppState::BridgeFamily::Cubic;
      else if (sys_fam == 2) app.bridge_family = AppState::BridgeFamily::Sine;
      if (!proj_ok && app.bridge_mode == AppState::BridgeMode::ProjectionSolid) {
        app.bridge_mode = AppState::BridgeMode::CurrentSystem;
        app.bridge_built = false; app.bridge_cam_init = false;
      }

      int mode = (app.bridge_mode == AppState::BridgeMode::ProjectionSolid) ? 2
               : (app.bridge_mode == AppState::BridgeMode::Classic)         ? 0 : 1;
      const int mode_was = mode;
      ImGui::TextUnformatted("mode:");
      ImGui::SameLine();
      if (!proj_ok) ImGui::BeginDisabled();
      if (ImGui::RadioButton("projection solid", mode == 2)) { mode = 2; }
      if (!proj_ok) ImGui::EndDisabled();
      if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
        if (proj_ok)
          ImGui::SetTooltip("ONE 3D object whose two shadows are this system's fractal set and its "
                            "bifurcation diagram. Family: %s (matched to the loaded system).", fam_name);
        else
          ImGui::SetTooltip("No unified bridge exists for this system. It only exists when the 1-D "
                            "map is the real slice of an iterated complex polynomial (logistic, "
                            "cubic, sine) or the complex-quadratic system. Load one of those, or "
                            "use 'current system'.");
      }
      ImGui::SameLine();
      if (ImGui::RadioButton("side-by-side", mode == 0)) { mode = 0; }
      if (ImGui::IsItemHovered())
        ImGui::SetTooltip("The Mandelbrot set (flat) + logistic bifurcation rising along its real "
                          "axis. A fixed Mandelbrot/logistic illustration, any system loaded.");
      ImGui::SameLine();
      if (ImGui::RadioButton("current system", mode == 1)) { mode = 1; }
      if (ImGui::IsItemHovered())
        ImGui::SetTooltip("Lifts the CURRENTLY LOADED system's bifurcation diagram into 3D.");
      if (mode != mode_was) { changed = true; app.bridge_built = false; app.bridge_cam_init = false; }
      app.bridge_mode = (mode == 2 && proj_ok) ? AppState::BridgeMode::ProjectionSolid
                      : (mode == 0) ? AppState::BridgeMode::Classic
                      : (mode == 1) ? AppState::BridgeMode::CurrentSystem
                                    : app.bridge_mode;

      if (app.bridge_mode == AppState::BridgeMode::ProjectionSolid) {
        ImGui::TextWrapped("One 3D object built by iterating this system's complex map over the "
                           "plane. Its footprint (look down the vertical axis) is the fractal "
                           "connectedness set; its real-axis silhouette (look from the front) is "
                           "the bifurcation diagram. Two projections of the SAME structure.");
        ImGui::Text("family (from loaded system): %s", fam_name);
        if (app.bridge_family == AppState::BridgeFamily::Sine)
          ImGui::TextDisabled("(sine is experimental: sparse and largely chaotic.)");
        changed |= ImGui::Checkbox("color by period", &app.bridge_color_by_period);
        changed |= ImGui::SliderInt("resolution", &app.bridge_mandel_res, 120, 700);
        ImGui::TextDisabled("Tip: look straight down -> the set; from the front -> the bifurcation "
                            "diagram. Left-drag rotate, right-drag move, wheel zoom.");
      } else if (app.bridge_mode == AppState::BridgeMode::Classic) {
        ImGui::TextWrapped("The Mandelbrot set (flat) with the logistic bifurcation diagram "
                           "rising along its real axis. They align because z^2+c and the "
                           "logistic map are conjugate. (This exact correspondence is special "
                           "to the quadratic family.)");
        changed |= ImGui::Checkbox("show Mandelbrot", &app.bridge_show_mandelbrot);
        ImGui::SameLine();
        changed |= ImGui::Checkbox("show bifurcation", &app.bridge_show_bifurcation);
        changed |= ImGui::Checkbox("color by period (show bulb \u2194 window correspondence)", &app.bridge_color_by_period);
        if (app.bridge_color_by_period)
          ImGui::TextDisabled("Each Mandelbrot bulb and the cascade window above it share a color "
                              "keyed to the attracting cycle's period (1 amber, 2 blue, 3 green, ...).");
        changed |= ImGui::SliderInt("Mandelbrot resolution", &app.bridge_mandel_res, 120, 700);
      } else {
        const int sf = bridge_family_for_system(app);
        if (sf < 0)
          ImGui::TextColored(ImVec4(1.0f, 0.78f, 0.35f, 1.0f),
                             "This system has no unified fractal<->bifurcation bridge, so this is "
                             "NOT the Mandelbrot-style 3D object -- it is just this system's "
                             "bifurcation diagram lifted into 3D.");
        ImGui::TextWrapped("Lifts THIS system's bifurcation diagram into 3D: the primary "
                           "parameter along x, the observable's value as height, and an "
                           "optional second parameter stacked along z. Works for any system.");
        /* primary parameter + range come from the bifurcation view's settings */
        ImGui::InputText("x param (primary)", app.bif_param, sizeof(app.bif_param));
        changed |= ImGui::InputDouble("x min", &app.bif_start, 0.0, 0.0, "%.4g");
        changed |= ImGui::InputDouble("x max", &app.bif_end, 0.0, 0.0, "%.4g");
        ImGui::InputText("observable", app.bif_observable, sizeof(app.bif_observable));
        ImGui::Separator();
        ImGui::InputText("z param (optional 2nd)", app.bridge_param2, sizeof(app.bridge_param2));
        if (app.bridge_param2[0] != '\0') {
          changed |= ImGui::InputDouble("z min", &app.bridge_p2_min, 0.0, 0.0, "%.4g");
          changed |= ImGui::InputDouble("z max", &app.bridge_p2_max, 0.0, 0.0, "%.4g");
          changed |= ImGui::SliderInt("z slices (stacked diagrams)", &app.bridge_p2_slices, 1, 80);
        } else {
          ImGui::TextDisabled("(leave z param empty for a single diagram lifted to 3D)");
        }
      }
      changed |= ImGui::SliderInt("bifurcation slices", &app.bridge_bif_slices, 60, 1500);
      changed |= ImGui::SliderInt("attractor points/slice", &app.bridge_bif_keep, 20, 400);
      changed |= ImGui::SliderFloat("bifurcation height", &app.bridge_height, 5.0f, 40.0f, "%.1f");
      if (ImGui::Button("Rebuild bridge") || changed) { app.bridge_built = false; app.bridge_cam_init = false; }
      ImGui::TextDisabled("%d points. drag to rotate, wheel to zoom.", app.bridge_point_count);
    }
  }

  /* PHASE B/C: basin-of-attraction controls. */
  if (app.active_view == AppState::ActiveView::Basins && app.state_names.size() >= 2) {
    if (ImGui::CollapsingHeader("Basins of attraction", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::TextWrapped("Each initial condition is integrated until it settles; pixels are "
                         "colored by which attractor they reach. The basin boundaries are "
                         "often fractal. Try a system with several attractors (e.g. Newton's "
                         "method, or the Duffing oscillator).");
      if (ImGui::SliderInt("steps per cell", &app.basin_steps, 100, 6000)) app.basin_dirty = true;
      if (ImGui::InputDouble("cluster tolerance", &app.basin_cluster_tol, 1e-3, 1e-2, "%.1e")) app.basin_dirty = true;
      if (ImGui::Checkbox("shade by convergence speed", &app.basin_shade_speed)) app.basin_dirty = true;
      if (ImGui::Button("Recompute basins")) app.basin_dirty = true;
      ImGui::TextDisabled("%d basins found. Pan/zoom shares the phase-plane view.", app.basin_attractor_count);
    }
  }

  /* PHASE B: 2-parameter scan controls. */
  if (app.active_view == AppState::ActiveView::ParamScan2D && app.params.size() >= 2) {
    if (ImGui::CollapsingHeader("2-parameter scan", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::TextWrapped("Largest Lyapunov exponent over a grid of TWO parameters: dark/blue = "
                         "periodic, bright/warm = chaotic.");
      ImGui::TextWrapped("When is this useful? To map a system's behaviour across its parameter "
                         "plane at a glance -- where it's periodic vs chaotic, where period-"
                         "doubling cascades and 'shrimp' windows sit, and how regimes border each "
                         "other. It answers 'for which parameter pairs is this system stable / "
                         "chaotic?' in one picture, instead of scanning one parameter at a time. "
                         "Best on maps (Henon, Tinkerbell, the standard map).");
      if (ImGui::BeginCombo("x param", app.params[std::min((size_t)app.scan_px_index, app.params.size()-1)].name.c_str())) {
        for (int i = 0; i < (int)app.params.size(); ++i)
          if (ImGui::Selectable(app.params[i].name.c_str(), i == app.scan_px_index)) { app.scan_px_index = i; app.scan_view_init = false; }
        ImGui::EndCombo();
      }
      if (ImGui::BeginCombo("y param", app.params[std::min((size_t)app.scan_py_index, app.params.size()-1)].name.c_str())) {
        for (int i = 0; i < (int)app.params.size(); ++i)
          if (ImGui::Selectable(app.params[i].name.c_str(), i == app.scan_py_index)) { app.scan_py_index = i; app.scan_view_init = false; }
        ImGui::EndCombo();
      }
      if (app.scan_px_index == app.scan_py_index)
        ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f),
                           "x and y are the same parameter — pick two different parameters for a 2-D scan.");
      ImGui::InputDouble("x min", &app.scan_xmin, 0.0, 0.0, "%.4g"); ImGui::SameLine();
      ImGui::InputDouble("x max", &app.scan_xmax, 0.0, 0.0, "%.4g");
      ImGui::InputDouble("y min", &app.scan_ymin, 0.0, 0.0, "%.4g"); ImGui::SameLine();
      ImGui::InputDouble("y max", &app.scan_ymax, 0.0, 0.0, "%.4g");
      ImGui::SliderInt("transient", &app.scan_transient, 50, 4000);
      ImGui::SliderInt("iterations", &app.scan_iterations, 100, 4000);
      if (ImGui::Button("Recompute scan")) app.scan_dirty = true;
      ImGui::TextDisabled("Note: scans are heavy; the grid is capped for responsiveness.");
    }
  }

  if (ImGui::CollapsingHeader("System editor")) {
    const char *current_preset_name = kPresets[std::max(0, std::min(app.selected_preset, kPresetCount - 1))].name;
    if (ImGui::BeginCombo("Preset", current_preset_name)) {
      for (int i = 0; i < kPresetCount; ++i) {
        char label[128];
        std::snprintf(label, sizeof(label), "%s / %s", kPresets[i].category, kPresets[i].name);
        const bool selected = app.selected_preset == i;
        if (ImGui::Selectable(label, selected)) {
          apply_preset(app, i);
          if (compile_system(app, app.system_input, &app.parse_error)) reset_simulation(app);
        }
        if (selected) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }
    /* Equation editor with a line-number gutter on the left. Systems are
     * short and typically don't scroll, so a simple aligned gutter is enough;
     * we don't reach into ImGui internals to chase the editor's scroll. */
    {
      const float ed_h = 280.0f;
      int nlines = 1;
      for (const char *p = app.system_input; *p; ++p) if (*p == '\n') ++nlines;
      ImGui::BeginChild("##gutter", ImVec2(34.0f, ed_h), false,
                        ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
      for (int i = 1; i <= nlines; ++i)
        ImGui::TextColored(ImVec4(0.45f, 0.50f, 0.58f, 1.0f), "%3d", i);
      ImGui::EndChild();
      ImGui::SameLine(0.0f, 4.0f);
      ImGui::InputTextMultiline("##system_input", app.system_input, sizeof(app.system_input),
                                ImVec2(-FLT_MIN, ed_h), ImGuiInputTextFlags_AllowTabInput);
    }
    if (ImGui::Button("Apply system")) {
      if (compile_system(app, app.system_input, &app.parse_error)) reset_simulation(app);
    }
    ImGui::SameLine();
    ImGui::TextDisabled("%zu defs, %zu params, %zu observables", app.definitions.size(), app.params.size(), app.observables.size());

    /* PHASE6-UI: dimensionality. Show what was auto-detected and let the
     * user override it. This controls whether the 3D view is offered and
     * how the system is treated, fixing "2D things treated as 3D" while
     * keeping the user in charge. */
    ImGui::Text("dimension: %d state vars, detected %d effective%s",
                static_cast<int>(app.state_names.size()), app.detected_dim,
                app.detected_dim < static_cast<int>(app.state_names.size())
                    ? " (has a trivial/dummy dimension)"
                    : "");
    int dim_mode = static_cast<int>(app.dim_override);
    const char *dim_modes[] = {"Auto (detected)", "Force 2D (phase only)",
                               "Force 3D"};
    if (ImGui::Combo("treat as", &dim_mode, dim_modes, 3))
      app.dim_override = static_cast<AppState::DimOverride>(dim_mode);
    ImGui::SameLine();
    ImGui::TextDisabled("now: %dD", effective_dimension(app));

    if (ImGui::TreeNode("Syntax")) {
      ImGui::TextWrapped("state x, y, z, w declares any dimension. ODEs use dx/dt or dx for every state variable; maps use x_next or x'.");
      ImGui::TextWrapped("plot3d = x, y, z chooses the 3D projection. param rho = 28 [0,100] creates a live slider.");
      ImGui::TextWrapped("observe r = sqrt(x*x+y*y+z*z) creates a plot variable. Poincare: section = z - 27; section_direction = positive; section_plot = x, y.");
      ImGui::TextWrapped("Builtins: sin cos tan asin acos atan exp log log10 sqrt abs floor ceil sign pow min max mod clamp select.");
      ImGui::TreePop();
    }
    if (!app.parse_error.empty()) ImGui::TextColored(ImVec4(1,0.35f,0.35f,1), "%s", app.parse_error.c_str());
    if (!app.runtime_error.empty()) ImGui::TextColored(ImVec4(1,0.35f,0.35f,1), "Runtime: %s", app.runtime_error.c_str());
  }

  if (app.mode == SystemMode::IFS &&
      ImGui::CollapsingHeader("IFS maps", ImGuiTreeNodeFlags_DefaultOpen)) {
    /* Surface the IFS's own contents: each affine map and its coefficients
     * x' = a x + b y + e,  y' = c x + d y + f,  with selection weight p.
     * Every coefficient gets a real SLIDER. When all coefficients are
     * constants the sliders edit them directly; when they're parameter
     * expressions the sliders are shown disabled (drive them via Parameters). */
    if (app.ifs_maps.empty()) evaluate_ifs_maps(app, false);
    ImGui::TextDisabled("x' = a x + b y + e      y' = c x + d y + f      (p = selection weight)");
    const bool has_param_coef = !app.ifs_maps_editable;
    if (has_param_coef)
      ImGui::TextColored(ImVec4(0.95f,0.8f,0.4f,1.0f),
                         "orange coefficients are parameter-driven; drag one to take manual control (it becomes a fixed value)");

    struct CoefSpec { const char *label; float lo, hi; };
    const CoefSpec specs[7] = {
      {"a", -1.0f, 1.0f}, {"b", -1.0f, 1.0f}, {"c", -1.0f, 1.0f}, {"d", -1.0f, 1.0f},
      {"e", -5.0f, 5.0f}, {"f", -5.0f, 5.0f}, {"p", 0.0f, 1.0f},
    };
    for (size_t mi = 0; mi < app.ifs_maps.size(); ++mi) {
      ImGui::PushID((int)mi);
      auto &m = app.ifs_maps[mi];
      double *coef[7] = {&m.a,&m.b,&m.c,&m.d,&m.e,&m.f,&m.p};
      ImGui::SeparatorText((std::string("map ") + std::to_string(mi + 1)).c_str());
      bool changed = false;
      for (int k = 0; k < 7; ++k) {
        /* every coefficient is draggable. If it's currently parameter-driven
         * (not literal), colour it orange; dragging it converts that slot to
         * a fixed literal so it stops being overwritten by the expression. */
        const bool param_driven = (mi < app.ifs_coef_literal.size())
                                   ? !app.ifs_coef_literal[mi][(size_t)k] : false;
        float lo = specs[k].lo, hi = specs[k].hi;
        if ((float)*coef[k] < lo) lo = (float)*coef[k];
        if ((float)*coef[k] > hi) hi = (float)*coef[k];
        float v = (float)*coef[k];
        ImGui::SetNextItemWidth(150);
        if (param_driven) ImGui::PushStyleColor(ImGuiCol_Text, IM_COL32(240,180,90,255));
        if (ImGui::SliderFloat(specs[k].label, &v, lo, hi, "%.4f")) {
          *coef[k] = (double)v;
          changed = true;
          if (param_driven && mi < app.ifs_coef_literal.size())
            app.ifs_coef_literal[mi][(size_t)k] = true; /* take manual control */
        }
        if (param_driven) ImGui::PopStyleColor();
        if (k % 2 == 0 && k < 6) ImGui::SameLine();
      }
      if (changed) app.ifs_dirty = true;
      ImGui::PopID();
    }

    ImGui::Separator();
    if (ImGui::SmallButton("+ add map")) {
      app.ifs_maps.push_back(dynsys::analysis::AffineMap{0.5,0,0,0.5,0,0,0.0});
      app.ifs_coef_literal.push_back({true,true,true,true,true,true,true});
      app.ifs_dirty = true;
    }
    ImGui::SameLine();
    if (ImGui::SmallButton("- remove last") && app.ifs_maps.size() > 1) {
      app.ifs_maps.pop_back();
      if (!app.ifs_coef_literal.empty()) app.ifs_coef_literal.pop_back();
      app.ifs_dirty = true;
    }
    ImGui::SameLine();
    ImGui::TextDisabled("(slider edits drive the attractor live)");
  }

  if (ImGui::CollapsingHeader("Parameters", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (app.params.empty()) ImGui::TextDisabled("No param declarations in this system.");
    for (auto &param : app.params) {
      double old = param.value;
      /* Every parameter gets a slider. If the model declared [min,max] we
       * use it; otherwise we synthesize an editable range from the value's
       * magnitude so there's always a draggable control (the user asked for
       * sliders on everything, not bare input boxes). */
      double lo, hi;
      if (param.has_range) { lo = param.min_value; hi = param.max_value; }
      else {
        if (!param.ui_range_init) { auto_slider_range(param.value, &param.ui_lo, &param.ui_hi); param.ui_range_init = true; }
        lo = param.ui_lo; hi = param.ui_hi;
      }
      float v = static_cast<float>(param.value);
      ImGui::SetNextItemWidth(-90);
      if (ImGui::SliderFloat(param.name.c_str(), &v, static_cast<float>(lo), static_cast<float>(hi), "%.6g"))
        param.value = static_cast<double>(v);
      ImGui::SameLine();
      if (ImGui::SmallButton((std::string(param.ui_show_exact ? "edit-" : "edit+") + "##e" + param.name).c_str()))
        param.ui_show_exact = !param.ui_show_exact;
      ImGui::SameLine();
      if (ImGui::SmallButton((std::string("reset##") + param.name).c_str())) {
        param.value = param.default_value;
        param.ui_range_init = false; /* re-fit the auto-range to the default */
        sync_param_values(app);
      }
      /* expander: type an exact value and/or adjust the (auto) range */
      if (param.ui_show_exact) {
        ImGui::Indent();
        ImGui::SetNextItemWidth(140);
        ImGui::InputDouble((std::string("value##v") + param.name).c_str(), &param.value, 0.01, 0.1, "%.8g");
        if (!param.has_range) {
          ImGui::SetNextItemWidth(90);
          ImGui::InputDouble((std::string("min##lo") + param.name).c_str(), &param.ui_lo, 0, 0, "%.4g");
          ImGui::SameLine(); ImGui::SetNextItemWidth(90);
          ImGui::InputDouble((std::string("max##hi") + param.name).c_str(), &param.ui_hi, 0, 0, "%.4g");
          ImGui::SameLine();
          if (ImGui::SmallButton((std::string("auto##ar") + param.name).c_str())) { auto_slider_range(param.value, &param.ui_lo, &param.ui_hi); }
        }
        ImGui::Unindent();
      }
      if (old != param.value) {
        app.poincare_points.clear();
        reset_lyapunov(app);
        sync_param_values(app);
        /* nudge cached views so the slider feels live there too */
        app.fractal_dirty = true;
        app.basin_dirty = true;
        app.ifs_dirty = true;
      }
    }
  }

  /* PHASE: initial-condition sliders. The starting point of the orbit is a
   * "variable value" the user can tweak; previously it had no control.
   * Each state variable gets a slider (auto-ranged); changing it re-seeds
   * the simulation so the orbit moves live. (Not for IFS, which has no
   * single trajectory.) */
  if (app.mode != SystemMode::IFS &&
      ImGui::CollapsingHeader("Initial conditions", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (app.state_names.empty()) ImGui::TextDisabled("No state variables.");
    static std::vector<double> ic_lo, ic_hi;
    static std::vector<unsigned char> ic_init;
    const size_t n = app.state_names.size();
    if (ic_lo.size() != n) { ic_lo.assign(n, -1.0); ic_hi.assign(n, 1.0); ic_init.assign(n, 0); }
    bool changed = false;
    for (size_t i = 0; i < n; ++i) {
      double cur = state_at(app.start, i);
      if (!ic_init[i]) { auto_slider_range(cur, &ic_lo[i], &ic_hi[i]); ic_init[i] = 1; }
      /* keep the range covering the current value */
      if (cur < ic_lo[i]) ic_lo[i] = cur;
      if (cur > ic_hi[i]) ic_hi[i] = cur;
      float v = static_cast<float>(cur);
      ImGui::SetNextItemWidth(-90);
      if (ImGui::SliderFloat((app.state_names[i] + "(0)").c_str(), &v,
                             static_cast<float>(ic_lo[i]), static_cast<float>(ic_hi[i]), "%.5g")) {
        set_state_at(app.start, i, static_cast<double>(v));
        changed = true;
      }
      ImGui::SameLine();
      if (ImGui::SmallButton((std::string("auto##ic") + std::to_string(i)).c_str()))
        ic_init[i] = 0;
    }
    if (changed) reset_simulation(app);
  }

  if (ImGui::CollapsingHeader("Simulation")) {
    /* Integrator and dt only apply to ODEs. Maps iterate discretely (no step
     * size), and IFS uses the chaos game — so those controls are hidden for
     * them and a short note explains why, rather than offering meaningless
     * choices. */
    if (app.mode == SystemMode::ODE) {
      int integrator_idx = static_cast<int>(app.integrator);
      const char *integrators[] = {"Euler", "RK2 midpoint", "Heun (RK2)", "RK4",
                                   "RK 3/8", "RKF45 (adaptive)", "Dormand-Prince (adaptive)"};
      if (ImGui::Combo("integrator", &integrator_idx, integrators, 7)) app.integrator = static_cast<Integrator>(integrator_idx);
      const bool adaptive = (app.integrator == Integrator::RKF45 || app.integrator == Integrator::DOPRI45);
      ImGui::InputDouble(adaptive ? "dt (output step)" : "dt", &app.dt, 0.001, 0.01, "%.8f");
      if (adaptive) {
        ImGui::InputDouble("tolerance", &app.adaptive_tol, 1e-7, 1e-6, "%.1e");
        ImGui::InputDouble("min substep", &app.adaptive_dt_min, 1e-7, 1e-6, "%.1e"); ImGui::SameLine();
        ImGui::InputDouble("max substep", &app.adaptive_dt_max, 1e-3, 1e-2, "%.1e");
        ImGui::TextDisabled("adaptive: subdivides each dt to meet the tolerance (last substep %.2e)", app.last_adaptive_dt);
      }
    } else if (app.mode == SystemMode::Map) {
      ImGui::TextDisabled("Maps iterate discretely — no integrator or step size.");
    } else {
      ImGui::TextDisabled("IFS uses the chaos game — no integrator or step size.");
    }
    ImGui::SliderInt("steps/frame", &app.steps_per_frame, 1, 2000);
    resize_state(app.start, app.state_names.size());
    for (size_t i = 0; i < app.state_names.size(); ++i) {
      std::string label = "start " + app.state_names[i];
      ImGui::InputDouble(label.c_str(), &app.start.v[i], 0.1, 1.0, "%.6f");
    }
    ImGui::InputDouble("start t", &app.start.t, 0.1, 1.0, "%.6f");
    int next_num_points = app.num_points;
    if (ImGui::SliderInt("3D point buffer", &next_num_points, 100, 300000)) {
      app.num_points = next_num_points;
      allocate_points(app);
      reset_simulation(app);
      upload_buffers(app);
    }
    ImGui::Text("history buffer: %d (tied to 3D point buffer)", app.history_limit);
  }

  ImGui::EndTabItem(); } /* Setup */
  if (ImGui::BeginTabItem("Analysis")) {

  if (app.mode == SystemMode::IFS) {
    /* The Lyapunov-spectrum / equilibria / limit-cycle / Poincare tools are
     * for flows and iterated maps; an IFS is a set of contractions explored by
     * the chaos game and characterized by its fractal dimension instead. Show
     * that rather than buttons that would produce meaningless results. */
    ImGui::TextWrapped("This is an IFS (iterated function system). The flow/map analysis tools "
                       "(Lyapunov spectrum, equilibria, limit cycles, Poincare sections) don't "
                       "apply to an IFS.");
    ImGui::TextWrapped("An IFS is characterized by its attractor's fractal dimension. Use the "
                       "Fractal/IFS view and the box-counting dimension there. Edit the maps and "
                       "weights in the Setup tab's \"IFS maps\" panel.");
  } else {
  if (ImGui::CollapsingHeader("Analysis", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("estimate largest Lyapunov exponent", &app.lyapunov_enabled);
    ImGui::InputDouble("Lyapunov epsilon", &app.lyapunov_epsilon, 1e-6, 1e-5, "%.1e");
    ImGui::SameLine();
    if (ImGui::Button("Reset Lyapunov")) reset_lyapunov(app);

    /* PHASE B: full Lyapunov spectrum + Kaplan-Yorke (fractal) dimension */
    ImGui::SeparatorText("Lyapunov spectrum (all exponents)");
    ImGui::SliderInt("transient##lyap", &app.lyap_spec_transient, 0, 20000);
    ImGui::SliderInt("average steps##lyap", &app.lyap_spec_steps, 1000, 200000);
    if (ImGui::Button("Compute Lyapunov spectrum")) run_lyapunov_spectrum(app);
    if (app.lyap_spectrum_ready && !app.lyap_spectrum.empty()) {
      std::string s = "exponents:";
      for (double l : app.lyap_spectrum) { char b[48]; std::snprintf(b, sizeof(b), " %.4f", l); s += b; }
      ImGui::TextWrapped("%s", s.c_str());
      ImGui::Text("sum = %.4f  (phase-space contraction rate)", app.lyap_spectrum_sum);
      ImGui::TextColored(ImVec4(0.55f, 0.85f, 1.0f, 1.0f),
                         "Kaplan-Yorke (fractal) dimension = %.4f", app.lyap_kaplan_yorke);
      /* quick interpretation */
      if (!app.lyap_spectrum.empty()) {
        const double l0 = app.lyap_spectrum[0];
        const char *verdict = l0 > 1e-3 ? "chaotic (largest exponent > 0)"
                             : (l0 < -1e-3 ? "stable (all exponents < 0)"
                                           : "marginal / periodic (largest ~ 0)");
        ImGui::TextDisabled("%s", verdict);
      }
    } else if (!app.lyap_spectrum_msg.empty()) {
      ImGui::TextDisabled("%s", app.lyap_spectrum_msg.c_str());
    }

    /* PHASE B/C: box-counting fractal dimension of the on-screen set */
    ImGui::SeparatorText("Box-counting (fractal) dimension");
    if (ImGui::Button("Measure box-counting dimension")) run_box_dimension(app);
    if (app.boxdim_ready) {
      ImGui::TextColored(ImVec4(0.55f, 0.85f, 1.0f, 1.0f),
                         "D_box = %.4f   (R^2 = %.3f, %ld points)",
                         app.boxdim_value, app.boxdim_r2, app.boxdim_n_points);
      ImGui::TextDisabled("%s",
        app.mode == SystemMode::Map ? "sampled the map attractor"
                                    : "measured from the current trajectory");
    } else if (!app.boxdim_msg.empty() && app.boxdim_msg != "ok") {
      ImGui::TextDisabled("%s", app.boxdim_msg.c_str());
    }

    /* PHASE D step 2 (foundation): limit-cycle period & amplitude */
    ImGui::SeparatorText("Limit cycle (period & amplitude)");
    if (ImGui::Button("Measure limit cycle")) run_limit_cycle(app);
    if (app.lc_ready) {
      ImGui::TextColored(ImVec4(0.55f, 0.85f, 1.0f, 1.0f),
                         "period T = %.5g   amplitude (peak-to-peak) = %.5g",
                         app.lc_period, app.lc_amplitude);
      ImGui::TextDisabled("measured on the '%s' component of the current orbit",
                          app.state_names[(size_t)std::max(0, std::min(app.phase_y_index, (int)app.state_names.size()-1))].c_str());
    } else if (!app.lc_msg.empty() && app.lc_msg != "ok") {
      ImGui::TextDisabled("%s", app.lc_msg.c_str());
    }

    ImGui::Text("Poincare points: %zu", app.poincare_points.size());
    if (ImGui::Button("Clear Poincare")) app.poincare_points.clear();
    ImGui::Separator();
    if (ImGui::Button("Find fixed point")) find_fixed_point(app);
    /* Exact equilibria via the CAS: finds ALL equilibria exactly (not just the
     * one Newton lands on), for polynomial fields with rational parameters. */
    if (dynsys::cas::is_available() && app.mode != SystemMode::IFS) {
      ImGui::SameLine();
      if (ImGui::Button("Find ALL equilibria exactly (CAS)")) run_exact_equilibria(app);
      if (app.cas_equi_ready) {
        if (!app.cas_equi_is_poly)
          ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.4f, 1.0f), "%s", app.cas_equi_msg.c_str());
        else {
          ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.6f, 1.0f), "exact equilibria:");
          for (const auto &s : app.cas_equi_lines) ImGui::BulletText("%s", s.c_str());
          if (app.state_names.size() > 1)
            ImGui::TextDisabled("(N-D: shown as the Groebner basis of the solution ideal — each generator = 0)");
        }
      } else if (!app.cas_equi_msg.empty()) {
        ImGui::TextDisabled("exact equilibria: %s", app.cas_equi_msg.c_str());
      }
    }
    if (app.fixed_ready) {
      std::string fp = "fixed point:";
      for (size_t i = 0; i < app.state_names.size(); ++i) {
        char buf[96];
        std::snprintf(buf, sizeof(buf), " %s=%.8g", app.state_names[i].c_str(), state_at(app.fixed_point, i));
        fp += buf;
      }
      ImGui::TextWrapped("%s | residual %.3g", fp.c_str(), app.fixed_residual);
      if (!app.fixed_classification.empty()) ImGui::Text("classification: %s", app.fixed_classification.c_str());
      if (!app.fixed_eigenvalues.empty()) {
        std::string ev = "eigenvalues:";
        for (const auto &z : app.fixed_eigenvalues) {
          char buf[96];
          if (std::fabs(z.second) < 1e-12) std::snprintf(buf, sizeof(buf), " %.6g", z.first);
          else std::snprintf(buf, sizeof(buf), " %.6g%+.6gi", z.first, z.second);
          ev += buf;
        }
        ImGui::TextWrapped("%s", ev.c_str());
      }

      /* PHASE1: general N-D classification + full complex spectrum. */
      if (app.fixed_general_ready) {
        ImGui::Separator();
        ImGui::TextWrapped("general classification: %s",
                           app.fixed_general.label.c_str());
        std::string spec = "spectrum:";
        for (const auto &z : app.fixed_general.eigenvalues) {
          char buf[96];
          if (std::fabs(z.imag()) < 1e-12)
            std::snprintf(buf, sizeof(buf), " %.6g", z.real());
          else
            std::snprintf(buf, sizeof(buf), " %.6g%+.6gi", z.real(), z.imag());
          spec += buf;
        }
        ImGui::TextWrapped("%s", spec.c_str());
        ImGui::Text("stable: %d   unstable: %d   marginal: %d",
                    app.fixed_general.n_stable, app.fixed_general.n_unstable,
                    app.fixed_general.n_center);
        if (app.fixed_general.hopf_candidate)
          ImGui::TextColored(ImVec4(0.5f, 0.8f, 1.0f, 1.0f),
                             "Hopf candidate (complex pair near imaginary axis)");
        if (app.fixed_general.fold_candidate)
          ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f),
                             "fold candidate (real eigenvalue near zero)");

        /* Normal-form classification of the Hopf point: the first Lyapunov
         * coefficient (the MatCont quantity). Offered whenever an equilibrium
         * is in hand — it self-checks for a pure-imaginary pair and reports if
         * this isn't actually a Hopf point. */
        if (app.fixed_ready) {
          if (ImGui::Button("Classify Hopf (first Lyapunov coeff)")) {
            AppState::Param *bp = app.params.empty() ? nullptr : find_param(app, app.cont_param);
            if (bp == nullptr && !app.params.empty()) bp = &app.params[0];
            if (bp != nullptr) {
              dynsys::analysis::Model model = build_model(app, bp);
              std::vector<double> xeq(app.state_names.size());
              for (size_t i = 0; i < xeq.size(); ++i) xeq[i] = state_at(app.fixed_point, i);
              std::string e;
              app.hopf_l1_ready = dynsys::analysis::hopf_first_lyapunov(
                  model, xeq, bp->value, &app.hopf_l1, &app.hopf_omega, &e);
              app.hopf_l1_msg = app.hopf_l1_ready ? "ok" : e;
            } else {
              app.hopf_l1_ready = false;
              app.hopf_l1_msg = "no parameter to evaluate the vector field";
            }
          }
          if (app.hopf_l1_ready) {
            const char *verdict = app.hopf_l1 < 0
                ? "SUPERcritical Hopf (l1<0): a stable limit cycle is born"
                : (app.hopf_l1 > 0
                     ? "SUBcritical Hopf (l1>0): an unstable cycle; hard loss of stability"
                     : "degenerate (l1~0): Bautin / generalized-Hopf codim-2 point");
            ImGui::TextColored(app.hopf_l1 < 0 ? ImVec4(0.5f,1.0f,0.6f,1.0f)
                                               : ImVec4(1.0f,0.7f,0.4f,1.0f),
                               "first Lyapunov coeff l1 = %.4g  (omega = %.4g)", app.hopf_l1, app.hopf_omega);
            ImGui::TextWrapped("%s", verdict);
          } else if (!app.hopf_l1_msg.empty() && app.hopf_l1_msg != "ok") {
            ImGui::TextDisabled("Hopf classification: %s", app.hopf_l1_msg.c_str());
          }

          /* Fold (limit point) normal-form coefficient — the LP counterpart of
           * the Hopf l1, also a MatCont quantity. a != 0 => genuine quadratic
           * fold; a ~ 0 => cusp (codim-2). Self-checks for a near-zero
           * eigenvalue and says when the point isn't actually a fold. */
          if (ImGui::Button("Classify fold (normal-form coeff a)")) {
            AppState::Param *bp = app.params.empty() ? nullptr : find_param(app, app.cont_param);
            if (bp == nullptr && !app.params.empty()) bp = &app.params[0];
            if (bp != nullptr) {
              dynsys::analysis::Model model = build_model(app, bp);
              std::vector<double> xeq(app.state_names.size());
              for (size_t i = 0; i < xeq.size(); ++i) xeq[i] = state_at(app.fixed_point, i);
              std::string e;
              app.fold_a_ready = dynsys::analysis::fold_normal_form(
                  model, xeq, bp->value, &app.fold_a, &app.fold_lambda0, &e);
              app.fold_a_msg = app.fold_a_ready ? "ok" : e;
            } else {
              app.fold_a_ready = false;
              app.fold_a_msg = "no parameter to evaluate the vector field";
            }
          }
          if (app.fold_a_ready) {
            const bool at_fold = app.fold_lambda0 >= 0 && app.fold_lambda0 < 1e-2;
            const char *verdict = (std::fabs(app.fold_a) < 1e-3)
                ? "a ~ 0: degenerate fold -> CUSP (codim-2)"
                : "a != 0: non-degenerate quadratic fold (saddle-node)";
            ImGui::TextColored(std::fabs(app.fold_a) < 1e-3 ? ImVec4(1.0f,0.7f,0.4f,1.0f)
                                                            : ImVec4(0.5f,0.85f,1.0f,1.0f),
                               "fold coefficient a = %.4g   (|smallest eigenvalue| = %.2e)",
                               app.fold_a, app.fold_lambda0);
            ImGui::TextWrapped("%s", verdict);
            if (!at_fold)
              ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.3f, 1.0f),
                                 "note: no eigenvalue is near zero here, so this equilibrium isn't at a fold — the coefficient is only meaningful AT a fold point.");
          } else if (!app.fold_a_msg.empty() && app.fold_a_msg != "ok") {
            ImGui::TextDisabled("fold classification: %s", app.fold_a_msg.c_str());
          }
        }
      }

      /* PHASE E: exact, certificate-backed eigenvalues from the Sangaku CAS.
       * Only offered when the CAS is reachable; degrades silently otherwise so
       * the numeric spectrum above always stands on its own. */
      if (dynsys::cas::is_available() && !app.fixed_jacobian.empty()) {
        ImGui::Separator();
        if (ImGui::Button("Certified eigenvalues (CAS)")) {
          const int n = (int)app.state_names.size();
          if ((int)app.fixed_jacobian.size() == n * n) {
            bool exact = false;
            app.cas_eig = dynsys::cas::eigen_report_from_doubles(
                app.fixed_jacobian.data(), n, &exact);
            app.cas_eig_exact = exact;
            app.cas_eig_ready = true;
          }
        }
        ImGui::SameLine();
        ImGui::TextDisabled("(exact spectrum + stability via Sangaku)");

        if (app.cas_eig_ready) {
          const auto &R = app.cas_eig;
          if (!R.ok) {
            ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.3f, 1.0f),
                               "CAS: %s", R.message.c_str());
          } else {
            if (!app.cas_eig_exact)
              ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.4f, 1.0f),
                                 "approximate: Jacobian not exactly rational "
                                 "(showing nearest rational system)");
            else
              ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.6f, 1.0f),
                                 "exact (rational Jacobian) -- kernel-checked");
            /* exact eigenvalue strings */
            for (const auto &s : R.eigen_strings)
              ImGui::TextWrapped("  lambda: %s", s.c_str());
            /* certified stability verdict from exact Routh-Hurwitz */
            const char *verdict =
                R.stable ? "STABLE (all eigenvalues in the left half-plane)"
                         : (R.n_imag > 0
                                ? "MARGINAL (eigenvalues on the imaginary axis)"
                                : "UNSTABLE (an eigenvalue in the right half-plane)");
            ImVec4 col = R.stable ? ImVec4(0.5f, 1.0f, 0.6f, 1.0f)
                       : (R.n_imag > 0 ? ImVec4(0.5f, 0.8f, 1.0f, 1.0f)
                                       : ImVec4(1.0f, 0.6f, 0.5f, 1.0f));
            ImGui::TextColored(col, "%s", verdict);
            ImGui::Text("Re<0: %d   Re>0: %d   Re=0: %d",
                        R.n_lhp, R.n_rhp, R.n_imag);
            if (R.n_imag == 2 && R.n_rhp == 0)
              ImGui::TextColored(ImVec4(0.5f, 0.8f, 1.0f, 1.0f),
                                 "exact Hopf point: a pair sits precisely on the "
                                 "imaginary axis (Re = 0 exactly)");
          }
        }
      }
      ImGui::InputDouble("separatrix epsilon", &app.separatrix_epsilon, 1e-5, 1e-4, "%.1e");
      ImGui::SliderInt("separatrix steps", &app.separatrix_steps, 100, 20000);
      if (ImGui::Button("Trace saddle separatrices")) trace_separatrices(app);
      ImGui::SameLine();
      if (ImGui::Button("Clear separatrices")) app.separatrix_curves.clear();
      ImGui::Text("Jacobian rows:");
      const size_t n = app.state_names.size();
      for (size_t r = 0; r < n && r < 12; ++r) {
        std::string row = "[";
        for (size_t c = 0; c < n; ++c) {
          char buf[64];
          const double v = (r * n + c) < app.fixed_jacobian.size() ? app.fixed_jacobian[r * n + c] : 0.0;
          std::snprintf(buf, sizeof(buf), "%s%.4g", c ? " " : "", v);
          row += buf;
        }
        row += "]";
        ImGui::TextUnformatted(row.c_str());
      }
    }
    if (!app.analysis_message.empty()) ImGui::TextWrapped("%s", app.analysis_message.c_str());
  }

  if (ImGui::CollapsingHeader("Equilibrium continuation", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::TextWrapped(
        "Pseudo-arclength continuation of the located equilibrium as one "
        "parameter varies, with fold and Hopf detection. Find a fixed point "
        "first, then continue.");
    ImGui::InputText("continuation parameter", app.cont_param,
                     sizeof(app.cont_param));
    ImGui::InputDouble("param min", &app.cont_p_min, 0.1, 1.0, "%.4g");
    ImGui::SameLine();
    ImGui::InputDouble("param max", &app.cont_p_max, 0.1, 1.0, "%.4g");
    ImGui::InputDouble("initial step h0", &app.cont_h0, 0.01, 0.1, "%.4g");
    ImGui::SliderInt("max points", &app.cont_max_points, 50, 5000);
    int dir_idx = app.cont_direction > 0 ? 0 : 1;
    const char *dirs[] = {"increasing first", "decreasing first"};
    if (ImGui::Combo("direction", &dir_idx, dirs, 2))
      app.cont_direction = dir_idx == 0 ? 1 : -1;
    ImGui::Checkbox("detect fold", &app.cont_detect_fold);
    ImGui::SameLine();
    ImGui::Checkbox("detect Hopf", &app.cont_detect_hopf);

    /* Two-parameter continuation: trace a fold or Hopf CURVE in the
     * (cont_param, 2nd param) plane. This is the codim-1-curve-in-2-params
     * feature — e.g. where does a fold persist as two parameters vary. */
    if (app.params.size() >= 2) {
      ImGui::SeparatorText("two-parameter curves (codim-1 in 2 params)");
      ImGui::TextDisabled("x axis = continuation param '%s'; choose a 2nd param for the y axis:", app.cont_param);
      if (ImGui::BeginCombo("2nd param (y)", app.twopar_p2[0] ? app.twopar_p2 : "(pick)")) {
        for (const auto &pr : app.params)
          if (ImGui::Selectable(pr.name.c_str(), pr.name == app.twopar_p2))
            std::snprintf(app.twopar_p2, sizeof(app.twopar_p2), "%s", pr.name.c_str());
        ImGui::EndCombo();
      }
      if (ImGui::Button("Trace fold curve")) run_two_param_curve(app, dynsys::analysis::TwoParamKind::Fold);
      ImGui::SameLine();
      if (ImGui::Button("Trace Hopf curve")) run_two_param_curve(app, dynsys::analysis::TwoParamKind::Hopf);
      if (!app.twopar_msg.empty())
        ImGui::TextWrapped("%s", app.twopar_msg.c_str());
      if (app.twopar_ready && app.twopar_curve.points.size() > 1) {
        /* small (p,q) plot of the curve */
        ImDrawList *dl = ImGui::GetWindowDrawList();
        ImVec2 o = ImGui::GetCursorScreenPos();
        const float W = ImGui::GetContentRegionAvail().x, H = 150.0f;
        ImVec2 e(o.x + W, o.y + H);
        dl->AddRectFilled(o, e, IM_COL32(16,16,20,220));
        dl->AddRect(o, e, IM_COL32(110,110,120,200));
        double pmn=1e300,pmx=-1e300,qmn=1e300,qmx=-1e300;
        for (auto &pt : app.twopar_curve.points){ pmn=std::min(pmn,pt.p);pmx=std::max(pmx,pt.p);qmn=std::min(qmn,pt.q);qmx=std::max(qmx,pt.q);}
        const double pr=(pmx-pmn)>1e-12?(pmx-pmn):1.0, qr=(qmx-qmn)>1e-12?(qmx-qmn):1.0;
        auto sc=[&](double p,double q){ return ImVec2((float)(o.x+6+(W-12)*(p-pmn)/pr),(float)(e.y-6-(H-12)*(q-qmn)/qr)); };
        const ImU32 col = app.twopar_curve.kind==dynsys::analysis::SpecialPointKind::Hopf
                          ? IM_COL32(120,200,255,235) : IM_COL32(255,210,80,235);
        for (size_t i=1;i<app.twopar_curve.points.size();++i){
          auto&a=app.twopar_curve.points[i-1]; auto&b=app.twopar_curve.points[i];
          /* break the line across the two-direction join (large jump) */
          if (std::fabs(a.p-b.p)+std::fabs(a.q-b.q) < 0.4*(pr+qr))
            dl->AddLine(sc(a.p,a.q), sc(b.p,b.q), col, 2.0f);
        }
        /* codim-2 markers detected ON the curve (cusp / BT / generalized Hopf) */
        for (size_t idx : app.twopar_curve.special_indices) {
          const auto &sp = app.twopar_curve.points[idx];
          ImVec2 c2 = sc(sp.p, sp.q); const float r = 5.0f;
          ImU32 mc; const char *ml;
          switch (sp.special) {
            case dynsys::analysis::SpecialPointKind::Cusp:            mc=IM_COL32(255,160,80,255);  ml="CP"; break;
            case dynsys::analysis::SpecialPointKind::BogdanovTakens:  mc=IM_COL32(255,120,120,255); ml="BT"; break;
            case dynsys::analysis::SpecialPointKind::GeneralizedHopf: mc=IM_COL32(180,140,255,255); ml="GH"; break;
            default: mc=IM_COL32(230,230,230,255); ml="*"; break;
          }
          dl->AddRectFilled(ImVec2(c2.x-r,c2.y-r), ImVec2(c2.x+r,c2.y+r), mc);
          dl->AddRect(ImVec2(c2.x-r,c2.y-r), ImVec2(c2.x+r,c2.y+r), IM_COL32(20,20,25,255), 0,0,1.5f);
          dl->AddText(ImVec2(c2.x+7,c2.y-6), mc, ml);
        }
        char lab[192]; std::snprintf(lab,sizeof(lab),"%s  (x: %s, y: %s)%s",
                       app.twopar_curve.kind==dynsys::analysis::SpecialPointKind::Hopf?"Hopf curve":"fold curve",
                       app.cont_param, app.twopar_p2,
                       app.twopar_curve.special_indices.empty()?"":"  CP=cusp BT=Bogdanov-Takens GH=gen.Hopf");
        dl->AddText(ImVec2(o.x+8,o.y+6), col, lab);
        ImGui::Dummy(ImVec2(W, H+4));

        /* text readout of each codim-2 point: refined location + normal form */
        for (size_t idx : app.twopar_curve.special_indices) {
          const auto &sp = app.twopar_curve.points[idx];
          const char *kn = sp.special==dynsys::analysis::SpecialPointKind::BogdanovTakens ? "Bogdanov-Takens"
                         : sp.special==dynsys::analysis::SpecialPointKind::Cusp ? "Cusp"
                         : sp.special==dynsys::analysis::SpecialPointKind::GeneralizedHopf ? "generalized Hopf" : "codim-2";
          if (sp.special==dynsys::analysis::SpecialPointKind::BogdanovTakens && sp.has_codim2_nf)
            ImGui::BulletText("%s at (%s=%.5g, %s=%.5g)  —  BT normal form: a=%.4g, b=%.4g  (a,b != 0 => non-degenerate)",
                              kn, app.cont_param, sp.p2, app.twopar_p2, sp.q2, sp.bt_a, sp.bt_b);
          else if (sp.special==dynsys::analysis::SpecialPointKind::Cusp && sp.has_codim2_nf)
            ImGui::BulletText("%s at (%s=%.5g, %s=%.5g)  —  cusp cubic coefficient c=%.4g  (c != 0 => non-degenerate)",
                              kn, app.cont_param, sp.p2, app.twopar_p2, sp.q2, sp.cusp_c);
          else if (sp.special==dynsys::analysis::SpecialPointKind::GeneralizedHopf && sp.has_codim2_nf)
            ImGui::BulletText("%s at (%s=%.5g, %s=%.5g)  —  2nd Lyapunov coeff l2=%.4g (sign %s => fold-of-cycles opens %s)",
                              kn, app.cont_param, sp.p2, app.twopar_p2, sp.q2, sp.gh_l2,
                              sp.gh_l2 < 0 ? "-" : "+", sp.gh_l2 < 0 ? "supercritically" : "subcritically");
          else
            ImGui::BulletText("%s at (%s=%.5g, %s=%.5g)  [refined]",
                              kn, app.cont_param, sp.p2, app.twopar_p2, sp.q2);
        }
      }
    }

    /* Periodic-orbit continuation by collocation: solves the cycle + its
     * period as a BVP and follows it (tracks unstable cycles and folds of
     * cycles), distinct from the simulate-and-measure sweep in the LimitCycle
     * view. Needs a 2+ D ODE with a stable cycle to seed from. */
    if (app.mode == SystemMode::ODE && app.state_names.size() >= 2 && !app.params.empty()) {
      ImGui::SeparatorText("periodic-orbit continuation (collocation)");
      ImGui::TextDisabled("continues in '%s' from a cycle seeded at the current value", app.cont_param);
      if (ImGui::Button("Continue limit cycle (collocation)")) run_cycle_collocation(app);
      if (!app.cyc_msg.empty()) ImGui::TextWrapped("%s", app.cyc_msg.c_str());

      /* Homoclinic orbit: connection from a saddle back to itself, solved as a
       * truncated boundary-value problem with projection boundary conditions. */
      ImGui::SeparatorText("homoclinic orbit (saddle connection)");
      ImGui::TextDisabled("finds a saddle, seeds an excursion off it, solves the truncated BVP");
      if (ImGui::Button("Find homoclinic orbit")) run_homoclinic(app);
      if (!app.homoclinic_msg.empty()) ImGui::TextWrapped("%s", app.homoclinic_msg.c_str());
      if (app.homoclinic_ready && app.homoclinic.orbit.size() > 2) {
        const auto &H = app.homoclinic;
        ImGui::Text("saddle at (%.4g, %.4g)  |  half-time T = %.3g  |  peak deviation = %.3g  |  %d unstable dir%s",
                    H.saddle.size()>0?H.saddle[0]:0.0, H.saddle.size()>1?H.saddle[1]:0.0,
                    H.T, H.amplitude, H.n_unstable, H.n_unstable==1?"":"s");
        /* phase-plane plot of the orbit (first two coordinates) + the saddle */
        ImDrawList *dl = ImGui::GetWindowDrawList();
        ImVec2 o = ImGui::GetCursorScreenPos();
        const float W = ImGui::GetContentRegionAvail().x, Hh = 200.0f;
        ImVec2 e(o.x + W, o.y + Hh);
        dl->AddRectFilled(o, e, IM_COL32(16,16,20,220));
        dl->AddRect(o, e, IM_COL32(110,110,120,200));
        double xmn=1e300,xmx=-1e300,ymn=1e300,ymx=-1e300;
        for (auto &pt : H.orbit) { xmn=std::min(xmn,pt[0]);xmx=std::max(xmx,pt[0]);
          if (pt.size()>1){ymn=std::min(ymn,pt[1]);ymx=std::max(ymx,pt[1]);} }
        if (H.orbit[0].size()<2) { ymn=-1; ymx=1; }
        const double xr=(xmx-xmn)>1e-9?(xmx-xmn):1.0, yr=(ymx-ymn)>1e-9?(ymx-ymn):1.0;
        auto sc=[&](double x,double y){ return ImVec2((float)(o.x+8+(W-16)*(x-xmn)/xr),(float)(e.y-8-(Hh-16)*(y-ymn)/yr)); };
        for (size_t i=1;i<H.orbit.size();++i){
          double y0 = H.orbit[i-1].size()>1?H.orbit[i-1][1]:0.0;
          double y1 = H.orbit[i].size()>1?H.orbit[i][1]:0.0;
          dl->AddLine(sc(H.orbit[i-1][0],y0), sc(H.orbit[i][0],y1), IM_COL32(140,230,255,255), 2.0f);
        }
        if (H.saddle.size()>=1) {
          ImVec2 sp = sc(H.saddle[0], H.saddle.size()>1?H.saddle[1]:0.0);
          dl->AddLine(ImVec2(sp.x-5,sp.y-5),ImVec2(sp.x+5,sp.y+5),IM_COL32(255,120,120,255),2.0f);
          dl->AddLine(ImVec2(sp.x-5,sp.y+5),ImVec2(sp.x+5,sp.y-5),IM_COL32(255,120,120,255),2.0f);
          dl->AddText(ImVec2(sp.x+7,sp.y-6),IM_COL32(255,120,120,255),"saddle");
        }
        dl->AddText(ImVec2(o.x+8,o.y+6), IM_COL32(200,200,210,220),
                    "homoclinic orbit (cyan) leaving & returning to the saddle (x-y projection)");
        ImGui::Dummy(ImVec2(W, Hh+4));

        /* two-parameter continuation of the homoclinic locus */
        ImGui::TextDisabled("continue the homoclinic in 2 params (%s vs %s):", app.cont_param, app.twopar_p2);
        if (ImGui::Button("Continue homoclinic curve")) run_homoclinic_curve(app);
        if (app.homoclinic_curve_ready && app.homoclinic_curve.points.size() > 1) {
          const auto &HC = app.homoclinic_curve;
          ImGui::Text("homoclinic locus: %zu points", HC.points.size());
          ImDrawList *dl2 = ImGui::GetWindowDrawList();
          ImVec2 o2 = ImGui::GetCursorScreenPos();
          const float W2 = ImGui::GetContentRegionAvail().x, H2 = 170.0f;
          ImVec2 e2(o2.x + W2, o2.y + H2);
          dl2->AddRectFilled(o2, e2, IM_COL32(16,16,20,220));
          dl2->AddRect(o2, e2, IM_COL32(110,110,120,200));
          double qmn=1e300,qmx=-1e300,pmn2=1e300,pmx2=-1e300;
          for (auto &pt : HC.points){ qmn=std::min(qmn,pt.q);qmx=std::max(qmx,pt.q);pmn2=std::min(pmn2,pt.p);pmx2=std::max(pmx2,pt.p);}
          const double qr=(qmx-qmn)>1e-12?(qmx-qmn):1.0, pr2=(pmx2-pmn2)>1e-12?(pmx2-pmn2):1.0;
          auto scC=[&](double q,double p){ return ImVec2((float)(o2.x+8+(W2-16)*(q-qmn)/qr),(float)(e2.y-8-(H2-16)*(p-pmn2)/pr2)); };
          /* sort points by q for a clean polyline */
          std::vector<const dynsys::analysis::HomoclinicCurvePoint*> sorted;
          for (auto &pt : HC.points) sorted.push_back(&pt);
          std::sort(sorted.begin(), sorted.end(), [](auto a, auto b){ return a->q < b->q; });
          for (size_t i=1;i<sorted.size();++i)
            dl2->AddLine(scC(sorted[i-1]->q,sorted[i-1]->p), scC(sorted[i]->q,sorted[i]->p), IM_COL32(255,180,120,255), 2.0f);
          for (auto *pt : sorted) dl2->AddCircleFilled(scC(pt->q,pt->p), 2.5f, IM_COL32(255,210,150,255));
          dl2->AddText(ImVec2(o2.x+8,o2.y+6), IM_COL32(200,200,210,220), "homoclinic locus (orange): primary param (y) vs secondary (x)");
          ImGui::Dummy(ImVec2(W2, H2+4));
        }
      }
      if (app.cyc_ready && app.cyc_branch.samples.size() > 1) {
        /* plot period and amplitude vs parameter */
        ImDrawList *dl = ImGui::GetWindowDrawList();
        ImVec2 o = ImGui::GetCursorScreenPos();
        const float W = ImGui::GetContentRegionAvail().x, H = 150.0f;
        ImVec2 e(o.x + W, o.y + H);
        dl->AddRectFilled(o, e, IM_COL32(16,16,20,220));
        dl->AddRect(o, e, IM_COL32(110,110,120,200));
        double pmn=1e300,pmx=-1e300,amn=1e300,amx=-1e300,Tmn=1e300,Tmx=-1e300;
        for (auto &smp : app.cyc_branch.samples){ pmn=std::min(pmn,smp.p);pmx=std::max(pmx,smp.p);
          amn=std::min(amn,smp.amplitude);amx=std::max(amx,smp.amplitude);
          Tmn=std::min(Tmn,smp.period);Tmx=std::max(Tmx,smp.period);}
        const double pr=(pmx-pmn)>1e-12?(pmx-pmn):1.0;
        const double ar=(amx-amn)>1e-12?(amx-amn):1.0, Tr=(Tmx-Tmn)>1e-12?(Tmx-Tmn):1.0;
        auto scA=[&](double p,double a){ return ImVec2((float)(o.x+6+(W-12)*(p-pmn)/pr),(float)(e.y-6-(H-12)*(a-amn)/ar)); };
        auto scT=[&](double p,double T){ return ImVec2((float)(o.x+6+(W-12)*(p-pmn)/pr),(float)(e.y-6-(H-12)*(T-Tmn)/Tr)); };
        for (size_t i=1;i<app.cyc_branch.samples.size();++i){
          auto&a=app.cyc_branch.samples[i-1]; auto&b=app.cyc_branch.samples[i];
          /* amplitude: green if stable, red if not; period: faint blue */
          dl->AddLine(scT(a.p,a.period), scT(b.p,b.period), IM_COL32(110,160,220,150), 1.5f);
          const ImU32 col = b.stable ? IM_COL32(120,220,135,235) : IM_COL32(255,130,95,235);
          dl->AddLine(scA(a.p,a.amplitude), scA(b.p,b.amplitude), col, 2.0f);
        }
        /* mark fold-of-cycles (LPC) points on the branch */
        for (const auto &smp : app.cyc_branch.samples)
          if (smp.is_fold) {
            ImVec2 cc = scA(smp.p, smp.amplitude);
            dl->AddCircleFilled(cc, 5.0f, IM_COL32(255,210,80,255));
            dl->AddText(ImVec2(cc.x+7,cc.y-6), IM_COL32(255,210,80,255), "LPC");
          }
        /* mark period-doubling (PD) and Neimark-Sacker (NS / torus) points */
        for (const auto &smp : app.cyc_branch.samples) {
          if (smp.is_pd) {
            ImVec2 cc = scA(smp.p, smp.amplitude);
            dl->AddCircleFilled(cc, 4.5f, IM_COL32(255,120,200,255));
            dl->AddText(ImVec2(cc.x+7,cc.y+4), IM_COL32(255,120,200,255), "PD");
          }
          if (smp.is_ns) {
            ImVec2 cc = scA(smp.p, smp.amplitude);
            dl->AddCircleFilled(cc, 4.5f, IM_COL32(150,200,255,255));
            dl->AddText(ImVec2(cc.x+7,cc.y+4), IM_COL32(150,200,255,255), "NS");
          }
        }
        dl->AddText(ImVec2(o.x+8,o.y+6), IM_COL32(200,200,210,220),
                    "amplitude (green=stable, red=unstable) + period (blue); yellow=LPC, pink=PD, cyan=NS");
        if (app.cyc_branch.turned)
          dl->AddText(ImVec2(o.x+8,o.y+22), IM_COL32(255,210,80,220),
                      "branch turned around a fold of cycles (arclength)");
        ImGui::Dummy(ImVec2(W, H+4));

        /* Floquet multipliers of the cycle at a chosen sample, drawn on the unit
         * circle: multipliers inside = stable directions, crossing the circle =
         * a bifurcation (+1 fold, -1 period-doubling, complex pair = torus). */
        {
          static int fl_idx = -1;
          const int nsamp = (int)app.cyc_branch.samples.size();
          if (fl_idx < 0 || fl_idx >= nsamp) fl_idx = nsamp / 2;
          ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x);
          ImGui::SliderInt("##flsample", &fl_idx, 0, nsamp - 1, "Floquet @ sample %d");
          const auto &fs = app.cyc_branch.samples[(size_t)fl_idx];
          ImGui::Text("parameter %s = %.5g   period = %.5g   |largest non-trivial multiplier| = %.4f  (%s)",
                      app.cont_param, fs.p, fs.period, fs.max_nontrivial_mult,
                      fs.stable ? "stable" : "unstable");
          ImDrawList *fd = ImGui::GetWindowDrawList();
          ImVec2 fo = ImGui::GetCursorScreenPos();
          const float FH = 150.0f, FW = ImGui::GetContentRegionAvail().x;
          ImVec2 fe(fo.x + FW, fo.y + FH);
          fd->AddRectFilled(fo, fe, IM_COL32(16,16,20,220));
          fd->AddRect(fo, fe, IM_COL32(110,110,120,200));
          const ImVec2 ctr(fo.x + FW*0.5f, fo.y + FH*0.5f);
          const float R = (FH*0.5f - 12.0f);
          /* unit circle + axes */
          fd->AddCircle(ctr, R, IM_COL32(120,120,140,220), 64, 1.5f);
          fd->AddLine(ImVec2(ctr.x - R - 8, ctr.y), ImVec2(ctr.x + R + 8, ctr.y), IM_COL32(80,80,90,200));
          fd->AddLine(ImVec2(ctr.x, ctr.y - R - 8), ImVec2(ctr.x, ctr.y + R + 8), IM_COL32(80,80,90,200));
          fd->AddText(ImVec2(ctr.x + R - 4, ctr.y + 4), IM_COL32(120,120,140,220), "+1");
          fd->AddText(ImVec2(ctr.x - R - 14, ctr.y + 4), IM_COL32(120,120,140,220), "-1");
          /* plot each multiplier */
          for (size_t i = 0; i < fs.floquet_re.size(); ++i) {
            const double re = fs.floquet_re[i], im = fs.floquet_im[i];
            const float X = ctr.x + (float)re * R, Y = ctr.y - (float)im * R;
            const double mag = std::hypot(re, im);
            const bool trivial = std::hypot(re - 1.0, im) < 0.06;
            ImU32 col = trivial ? IM_COL32(150,150,160,220)
                                : (mag < 1.0 ? IM_COL32(120,220,135,255) : IM_COL32(255,130,95,255));
            fd->AddCircleFilled(ImVec2(X, Y), 4.0f, col);
          }
          fd->AddText(ImVec2(fo.x+8, fo.y+6), IM_COL32(200,200,210,220),
                      "Floquet multipliers (grey=trivial ~1, green=|.|<1, red=|.|>1); circle = unit |.|=1");
          ImGui::Dummy(ImVec2(FW, FH+4));
        }
      }
      /* Two-parameter fold-of-cycles (LPC) curve. Needs a 2nd parameter. */
      if (app.params.size() >= 2) {
        if (ImGui::Button("Trace fold-of-cycles curve (2 params)")) run_lpc_curve(app);
        ImGui::SameLine();
        ImGui::TextDisabled("(x: %s, y: %s)", app.cont_param, app.twopar_p2[0]?app.twopar_p2:"pick 2nd param above");
        if (app.lpc_ready && !app.lpc_msg.empty()) ImGui::TextWrapped("%s", app.lpc_msg.c_str());
        else if (!app.lpc_ready && !app.lpc_msg.empty()) ImGui::TextDisabled("%s", app.lpc_msg.c_str());
        if (app.lpc_ready && app.lpc_curve_data.points.size() > 1) {
          ImDrawList *dl = ImGui::GetWindowDrawList();
          ImVec2 o = ImGui::GetCursorScreenPos();
          const float W = ImGui::GetContentRegionAvail().x, H = 150.0f;
          ImVec2 e(o.x+W,o.y+H);
          dl->AddRectFilled(o,e,IM_COL32(16,16,20,220)); dl->AddRect(o,e,IM_COL32(110,110,120,200));
          double pmn=1e300,pmx=-1e300,qmn=1e300,qmx=-1e300;
          for(auto&pt:app.lpc_curve_data.points){pmn=std::min(pmn,pt.p);pmx=std::max(pmx,pt.p);qmn=std::min(qmn,pt.q);qmx=std::max(qmx,pt.q);}
          const double pr=(pmx-pmn)>1e-12?(pmx-pmn):1.0, qr=(qmx-qmn)>1e-12?(qmx-qmn):1.0;
          auto sc=[&](double p,double q){ return ImVec2((float)(o.x+6+(W-12)*(p-pmn)/pr),(float)(e.y-6-(H-12)*(q-qmn)/qr)); };
          for(size_t i=1;i<app.lpc_curve_data.points.size();++i){
            auto&a=app.lpc_curve_data.points[i-1]; auto&b=app.lpc_curve_data.points[i];
            dl->AddLine(sc(a.p,a.q), sc(b.p,b.q), IM_COL32(255,180,90,235), 2.0f);
          }
          char lab[256]; std::snprintf(lab,sizeof(lab),"fold-of-cycles (LPC) curve  (x: %s, y: %s)", app.cont_param, app.twopar_p2);
          dl->AddText(ImVec2(o.x+8,o.y+6), IM_COL32(255,180,90,235), lab);
          ImGui::Dummy(ImVec2(W,H+4));
        }
      }
    }

    if (ImGui::Button("Continue equilibrium")) run_equilibrium_continuation(app);
    ImGui::SameLine();
    if (ImGui::Button("Clear branch")) {
      app.cont_branch = dynsys::analysis::Branch{};
      app.cont_ready = false;
      app.cont_switched_branch = dynsys::analysis::Branch{};
      app.cont_switched_ready = false;
    }

    if (app.cont_ready && !app.cont_branch.points.empty()) {
      ImGui::Text("branch: %zu points", app.cont_branch.points.size());
      /* Special points table with codim-2 classification + branch switching. */
      if (!app.cont_branch.special_indices.empty()) {
        ImGui::Text("special points:");
        for (size_t i : app.cont_branch.special_indices) {
          const auto &bp = app.cont_branch.points[i];
          using K = dynsys::analysis::SpecialPointKind;
          const char *kind = "end";
          switch (bp.special) {
            case K::Fold:            kind = "FOLD (LP)"; break;
            case K::Hopf:            kind = "HOPF"; break;
            case K::BogdanovTakens:  kind = "BOGDANOV-TAKENS (codim-2)"; break;
            case K::Cusp:            kind = "CUSP (codim-2)"; break;
            case K::GeneralizedHopf: kind = "GENERALIZED HOPF / Bautin (codim-2)"; break;
            case K::BranchPoint:     kind = "BRANCH POINT"; break;
            default: kind = "end"; break;
          }
          ImGui::BulletText("%s at %s = %.6g", kind, app.cont_param, bp.p);
          /* show the normal-form coefficient where we have one */
          if (bp.has_normal_form) {
            if (bp.special == K::Hopf || bp.special == K::GeneralizedHopf)
              ImGui::TextDisabled("      first Lyapunov coeff l1 = %.4g (%s)", bp.lyapunov1,
                                  bp.lyapunov1 < 0 ? "supercritical" : (bp.lyapunov1 > 0 ? "subcritical" : "degenerate"));
            else if (bp.special == K::Fold || bp.special == K::Cusp)
              ImGui::TextDisabled("      fold coeff a = %.4g (%s)", bp.fold_a,
                                  std::fabs(bp.fold_a) < 1e-2 ? "cusp/degenerate" : "non-degenerate");
          }
          /* offer branch switching at a branch point */
          if (bp.special == K::BranchPoint) {
            ImGui::SameLine();
            if (ImGui::SmallButton((std::string("switch branch##") + std::to_string(i)).c_str())) {
              AppState::Param *bpar = app.params.empty() ? nullptr : find_param(app, app.cont_param);
              if (bpar == nullptr && !app.params.empty()) bpar = &app.params[0];
              if (bpar != nullptr) {
                dynsys::analysis::Model model = build_model(app, bpar);
                dynsys::analysis::ContinuationSettings s;
                s.p_min = app.cont_p_min; s.p_max = app.cont_p_max;
                s.h0 = app.cont_h0; s.max_points = app.cont_max_points;
                s.detect_fold = app.cont_detect_fold; s.detect_hopf = app.cont_detect_hopf;
                app.cont_switched_branch = dynsys::analysis::switch_branch(model, bp, s);
                app.cont_switched_ready = app.cont_switched_branch.ok;
                app.analysis_message = app.cont_switched_branch.ok
                    ? std::string("branch switched: ") + app.cont_switched_branch.message
                    : std::string("branch switch: ") + app.cont_switched_branch.message;
                sync_param_values(app);
              }
            }
          }
        }
        if (app.cont_switched_ready)
          ImGui::TextColored(ImVec4(0.5f,1.0f,0.6f,1.0f),
                             "switched branch: %zu points (drawn in cyan on the diagram)",
                             app.cont_switched_branch.points.size());
      }
      /* Simple stability-vs-parameter strip: draw the branch as a
       * polyline of (p, first-coordinate), colored by stability. */
      ImDrawList *draw = ImGui::GetWindowDrawList();
      ImVec2 p0 = ImGui::GetCursorScreenPos();
      const float bw = ImGui::GetContentRegionAvail().x;
      const float bh = 160.0f;
      ImVec2 p1(p0.x + bw, p0.y + bh);
      draw->AddRectFilled(p0, p1, IM_COL32(18, 18, 22, 220));
      draw->AddRect(p0, p1, IM_COL32(110, 110, 120, 200));
      double pmn = 1e300, pmx = -1e300, xmn = 1e300, xmx = -1e300;
      for (const auto &bp : app.cont_branch.points) {
        pmn = std::min(pmn, bp.p);
        pmx = std::max(pmx, bp.p);
        const double xv = bp.x.empty() ? 0.0 : bp.x[0];
        xmn = std::min(xmn, xv);
        xmx = std::max(xmx, xv);
      }
      /* extend bounds over the switched branch too, so it's not clipped */
      if (app.cont_switched_ready)
        for (const auto &bp : app.cont_switched_branch.points) {
          pmn = std::min(pmn, bp.p); pmx = std::max(pmx, bp.p);
          const double xv = bp.x.empty() ? 0.0 : bp.x[0];
          xmn = std::min(xmn, xv); xmx = std::max(xmx, xv);
        }
      const double pr = (pmx - pmn) > 1e-12 ? (pmx - pmn) : 1.0;
      const double xr = (xmx - xmn) > 1e-12 ? (xmx - xmn) : 1.0;
      auto to_screen = [&](double p, double x) {
        return ImVec2(
            static_cast<float>(p0.x + 6 + (bw - 12) * (p - pmn) / pr),
            static_cast<float>(p1.y - 6 - (bh - 12) * (x - xmn) / xr));
      };
      for (size_t i = 1; i < app.cont_branch.points.size(); ++i) {
        const auto &a = app.cont_branch.points[i - 1];
        const auto &b = app.cont_branch.points[i];
        const ImU32 col = b.stable ? IM_COL32(120, 220, 135, 235)
                                   : IM_COL32(255, 130, 95, 235);
        draw->AddLine(to_screen(a.p, a.x.empty() ? 0 : a.x[0]),
                      to_screen(b.p, b.x.empty() ? 0 : b.x[0]), col, 2.0f);
      }
      for (size_t i : app.cont_branch.special_indices) {
        const auto &bp = app.cont_branch.points[i];
        const ImU32 col =
            bp.special == dynsys::analysis::SpecialPointKind::Hopf
                ? IM_COL32(120, 210, 255, 255)
                : IM_COL32(255, 210, 80, 255);
        draw->AddCircle(to_screen(bp.p, bp.x.empty() ? 0 : bp.x[0]), 5.0f, col,
                        16, 2.0f);
      }
      /* the switched (crossing) branch, drawn in cyan; stability still shown by
       * brightness so the user can read both branches at once. */
      if (app.cont_switched_ready) {
        for (size_t i = 1; i < app.cont_switched_branch.points.size(); ++i) {
          const auto &a = app.cont_switched_branch.points[i - 1];
          const auto &b = app.cont_switched_branch.points[i];
          const ImU32 col = b.stable ? IM_COL32(90, 230, 230, 235)
                                     : IM_COL32(120, 160, 200, 200);
          draw->AddLine(to_screen(a.p, a.x.empty() ? 0 : a.x[0]),
                        to_screen(b.p, b.x.empty() ? 0 : b.x[0]), col, 2.0f);
        }
      }
      draw->AddText(ImVec2(p0.x + 8, p0.y + 6), IM_COL32(200, 200, 210, 220),
                    "green: stable  red: unstable  cyan: switched branch  (vertical: first state var)");
      ImGui::Dummy(ImVec2(bw, bh + 4));
    }
  }

  if (ImGui::CollapsingHeader("Bifurcation diagram", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (app.params.empty()) {
      ImGui::TextDisabled("declare a param to sweep (e.g. param r = 3.9 [0,4]).");
    } else {
      /* One-click: set the sweep to this system's first parameter and its
       * declared range, observable = first state var. This is the fix for
       * "I can't set the right params" — it can't be wrong for the loaded
       * system. */
      if (ImGui::Button("Auto-configure for this system")) {
        const auto &p0 = app.params[0];
        std::snprintf(app.bif_param, sizeof(app.bif_param), "%s", p0.name.c_str());
        if (p0.has_range) { app.bif_start = p0.min_value; app.bif_end = p0.max_value; }
        else { app.bif_start = p0.value - 1.0; app.bif_end = p0.value + 1.0; }
        if (!app.state_names.empty())
          std::snprintf(app.bif_observable, sizeof(app.bif_observable), "%s", app.state_names[0].c_str());
        if (app.mode == SystemMode::Map) { app.bif_discard = 1000; app.bif_keep = 200; app.bif_slices = 800; }
        else { app.bif_discard = 4000; app.bif_keep = 2000; app.bif_slices = 400; }
        app.bif_view_valid = false;
      }
      /* parameter picker: when chosen, also snap the range to its [min,max] */
      if (ImGui::BeginCombo("sweep param", app.bif_param)) {
        for (const auto &p : app.params)
          if (ImGui::Selectable(p.name.c_str(), p.name == app.bif_param)) {
            std::snprintf(app.bif_param, sizeof(app.bif_param), "%s", p.name.c_str());
            if (p.has_range) { app.bif_start = p.min_value; app.bif_end = p.max_value; app.bif_view_valid = false; }
          }
        ImGui::EndCombo();
      }
      ImGui::InputText("observable", app.bif_observable, sizeof(app.bif_observable));
      ImGui::InputDouble("start", &app.bif_start, 0.0, 0.0, "%.4g"); ImGui::SameLine();
      ImGui::InputDouble("end", &app.bif_end, 0.0, 0.0, "%.4g");
      ImGui::SliderInt("param slices", &app.bif_slices, 50, 1200);
      ImGui::SliderInt("discard (transient)", &app.bif_discard, 0, 5000); ImGui::SameLine();
      ImGui::SliderInt("keep", &app.bif_keep, 10, 5000);
      if (app.mode != SystemMode::Map)
        ImGui::TextDisabled("ODE: records local maxima of the observable (a Poincare-style section)");
      ImGui::Checkbox("compute Lyapunov vs param", &app.bif_compute_lyapunov);
      if (ImGui::Button("Run bifurcation")) { run_bifurcation(app); app.active_view = AppState::ActiveView::Bifurcation; }
      ImGui::SameLine();
      if (ImGui::Button("View in full window")) app.active_view = AppState::ActiveView::Bifurcation;
      ImGui::SameLine();
      if (ImGui::Button("Clear")) { app.bifurcation_points.clear(); app.bifurcation_lyapunov.clear(); }
      ImGui::TextDisabled("sweeping %s in [%.4g, %.4g], observing %s",
                          app.bif_param, app.bif_start, app.bif_end, app.bif_observable);
    }
  }

  } /* else (non-IFS Analysis) */
  ImGui::EndTabItem(); } /* Analysis */
  if (ImGui::BeginTabItem("Data & view")) {

  if (ImGui::CollapsingHeader("Exports")) {
    ImGui::InputText("CSV path", app.export_path, sizeof(app.export_path));
    if (ImGui::Button("Export trajectory CSV")) export_trajectory_csv(app);
    ImGui::SameLine();
    if (ImGui::Button("Export Poincare CSV")) export_poincare_csv(app);
  }

  if (ImGui::CollapsingHeader("3D view controls")) {
    if (!system_is_3d(app)) {
      ImGui::TextColored(ImVec4(0.6f, 0.8f, 1.0f, 1.0f),
                         "This system is %dD — the 3D view tab is hidden.",
                         effective_dimension(app));
      ImGui::TextDisabled(
          "Use the 'phase plane' tab. To force a 3D view, set 'treat as' to "
          "Force 3D in the System editor.");
    } else {
      ImGui::TextDisabled(
          "The 3D scene is the '3D view' tab (drag to rotate, wheel to zoom). "
          "These are fine controls for the same camera.");
    }
    ImGui::Checkbox("axes", &app.show_axes); ImGui::SameLine(); ImGui::Checkbox("center cross", &app.show_center_cross);
    if (ImGui::Checkbox("orthographic 3D", &app.orthographic_3d)) update_projection(app);
    if (ImGui::SliderFloat("zoom", &app.zoom, 0.05f, 100.0f, "%.3g")) update_projection(app);
    ImGui::SliderFloat("scene scale", &app.scene_scale, 0.0005f, 2.0f, "%.4g");
    ImGui::SliderFloat("camera distance", &app.camera_distance, 5.0f, 500.0f, "%.3g");
    if (ImGui::Button("Reset view")) {
      app.angle_x = 0.0f; app.angle_y = 0.0f; app.zoom = 1.0f; app.scene_scale = 0.05f;
      app.camera_distance = 50.0f; app.center_x = app.center_y = app.center_z = 0.0f;
      update_projection(app);
    }
    ImGui::SliderFloat("angle x", &app.angle_x, -360.0f, 360.0f);
    ImGui::SliderFloat("angle y", &app.angle_y, -360.0f, 360.0f);
    ImGui::InputFloat("center x", &app.center_x, 0.1f, 1.0f, "%.3f");
    ImGui::InputFloat("center y", &app.center_y, 0.1f, 1.0f, "%.3f");
    ImGui::InputFloat("center z", &app.center_z, 0.1f, 1.0f, "%.3f");
  }

  if (ImGui::CollapsingHeader("Keyboard / mouse help")) {
    ImGui::TextWrapped("2D phase: left-click adds an orbit through the point, "
                       "shift+click restarts the live orbit, left-drag pans, "
                       "wheel or +/- zooms, double-click auto-fits. "
                       "3D: drag rotates, wheel or +/- zooms. "
                       "Space play/pause, C reset, Esc quit.");
  }

  /* PHASE6-UI: time series and Poincare live in the side panel now that
   * the phase/3D plots are the window background. */
  if (ImGui::CollapsingHeader("Time series")) {
    for (const auto &name : app.state_names) {
      std::string title = name + "(t)";
      draw_series_plot(app, title.c_str(), name.c_str(), ImVec2(-FLT_MIN, 90));
    }
    for (const auto &obs : app.observables)
      draw_series_plot(app, obs.name.c_str(), obs.name.c_str(), ImVec2(-FLT_MIN, 90));
  }
  if (ImGui::CollapsingHeader("Poincare section", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (app.poincare_points.empty()) {
      ImGui::TextWrapped("No section points yet. Enable a section in the system text, e.g. for a "
                         "3-D flow: 'section = z - 27; section_direction = positive; "
                         "section_plot = x, y'. Then run the simulation: each time the trajectory "
                         "crosses the plane z=27 going upward, the (x,y) point is recorded here.");
      ImGui::TextDisabled("For a flow with no explicit section, local maxima of the observable are "
                          "recorded as a Poincare-style section automatically.");
    } else {
      ImGui::Text("%zu section points", app.poincare_points.size());
      /* a large square plot so the section structure is actually legible */
      const float side = std::min(ImGui::GetContentRegionAvail().x, 560.0f);
      draw_scatter_plot("Poincare section", app.poincare_points, ImVec2(side, side));
      ImGui::TextDisabled("A Poincare section turns a continuous flow into a map: structure here "
                          "(curves, islands, scatter) reveals periodic / quasi-periodic / chaotic "
                          "motion. Export via the button above.");
    }
  }

  ImGui::EndTabItem(); } /* Data & view */
  ImGui::EndTabBar();
  } /* control_tabs */

  ImGui::End(); /* controls side panel */
}

bool init_glfw_window(AppState &app, GLFWwindow **out_window) {
  if (!glfwInit()) {
    std::fprintf(stderr, "Failed to initialize GLFW\n");
    return false;
  }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  GLFWwindow *window = glfwCreateWindow(app.window_width, app.window_height, "Dynsys", nullptr, nullptr);
  if (window == nullptr) {
    std::fprintf(stderr, "Failed to create GLFW window\n");
    glfwTerminate();
    return false;
  }
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);
  glfwSetWindowUserPointer(window, &app);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetScrollCallback(window, scroll_callback);
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK) {
    std::fprintf(stderr, "Failed to initialize GLEW\n");
    glfwDestroyWindow(window);
    glfwTerminate();
    return false;
  }
  glGetError();
  *out_window = window;
  return true;
}

void init_imgui(GLFWwindow *window) {
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
  ImGui::StyleColorsDark();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(kGlslVersion);
}

void shutdown_imgui() {
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
}

} // namespace

/* Headless driver: load a .dyn file (or use the default Lorenz),
 * step it N times, print the final state. Used for differential
 * testing (was the old simulation = new simulation?) and for
 * microbenchmarking the integrator hot loop in environments
 * without a display. */
int run_headless(int argc, char **argv) {
  AppState app{};
  const char *path = nullptr;
  long long steps = 1000;
  bool dump_each = false;
  bool use_ast = false;
  for (int i = 2; i < argc; ++i) {
    if (std::strcmp(argv[i], "--steps") == 0 && i + 1 < argc) {
      steps = std::strtoll(argv[++i], nullptr, 10);
    } else if (std::strcmp(argv[i], "--dump") == 0) {
      dump_each = true;
    } else if (std::strcmp(argv[i], "--use-ast") == 0) {
      use_ast = true;
    } else if (path == nullptr) {
      path = argv[i];
    } else {
      std::fprintf(stderr, "unexpected arg: %s\n", argv[i]);
      return EXIT_FAILURE;
    }
  }

  std::string source;
  if (path) {
    std::FILE *f = std::fopen(path, "rb");
    if (!f) { std::perror(path); return EXIT_FAILURE; }
    char buf[4096];
    while (size_t n = std::fread(buf, 1, sizeof buf, f)) source.append(buf, n);
    std::fclose(f);
  } else {
    set_lorenz(app);
    source.assign(app.system_input);
  }
  copy_system_input(app, source.c_str());

  std::string err;
  if (!compile_system(app, app.system_input, &err)) {
    std::fprintf(stderr, "compile failed: %s\n", err.c_str());
    return EXIT_FAILURE;
  }
  reset_simulation(app);
  app.use_ast_fallback = use_ast;

  std::printf("dim=%zu mode=%s integrator=%s dt=%.6f steps=%lld defs=%zu params=%zu path=%s\n",
              app.state_names.size(),
              mode_name(app.mode),
              integrator_name(app.integrator),
              app.dt, steps, app.definitions.size(), app.params.size(),
              use_ast ? "ast" : "ir");
  std::printf("initial: t=%.6f", app.current.t);
  for (size_t i = 0; i < app.state_names.size(); ++i) {
    std::printf(" %s=%.10f", app.state_names[i].c_str(), state_at(app.current, i));
  }
  std::printf("\n");

  char step_err[256] = {0};
  const auto t0 = std::chrono::steady_clock::now();
  for (long long s = 0; s < steps; ++s) {
    State next{};
    if (!step_state(app, app.current, &next, step_err, sizeof step_err)) {
      std::fprintf(stderr, "step %lld failed: %s\n", s, step_err);
      return EXIT_FAILURE;
    }
    app.current = next;
    if (dump_each) {
      std::printf("%lld t=%.6f", s, app.current.t);
      for (size_t i = 0; i < app.state_names.size(); ++i) {
        std::printf(" %s=%.10f", app.state_names[i].c_str(), state_at(app.current, i));
      }
      std::printf("\n");
    }
  }
  const auto t1 = std::chrono::steady_clock::now();
  const double elapsed_ns = std::chrono::duration<double, std::nano>(t1 - t0).count();

  std::printf("final:   t=%.6f", app.current.t);
  for (size_t i = 0; i < app.state_names.size(); ++i) {
    std::printf(" %s=%.10f", app.state_names[i].c_str(), state_at(app.current, i));
  }
  std::printf("\n");
  std::printf("elapsed: %.3f ms (%.1f ns/step)\n",
              elapsed_ns / 1e6, elapsed_ns / static_cast<double>(steps));

  if (app.arena_ready) arena_destroy(&app.system_arena);
  return EXIT_SUCCESS;
}

/* Adaptive performance governor. Called once per frame with the raw frame time
 * in milliseconds. Smooths it, then nudges app.perf_throttle toward 1 when
 * frames are too slow and toward 0 when there's headroom. The throttle is a
 * hysteretic ramp (not a hard switch) so the picture doesn't oscillate. The
 * heavy compute paths (fractal/IFS/integration substeps) read perf_throttle to
 * scale their budgets, and steps_per_frame is reduced directly here. This is
 * what keeps the program from spiraling into an unresponsive state on a very
 * heavy fractal or a huge step count. */
void update_perf_governor(AppState &app, double raw_frame_ms) {
  if (!app.perf_governor_enabled) {
    app.perf_throttle = 0.0;
    app.perf_throttle_active = false;
    return;
  }
  /* exponential smoothing so a single hitch doesn't jerk the throttle */
  const double alpha = 0.15;
  if (!std::isfinite(raw_frame_ms) || raw_frame_ms <= 0.0) raw_frame_ms = app.perf_frame_ms;
  raw_frame_ms = std::min(raw_frame_ms, 1000.0); /* clamp pathological spikes */
  app.perf_frame_ms = (1.0 - alpha) * app.perf_frame_ms + alpha * raw_frame_ms;

  /* budgets: comfortable below ~22 ms (~45 fps), distressed above ~40 ms
   * (~25 fps). Ramp the throttle within that band. */
  const double good_ms = 22.0, bad_ms = 40.0;
  if (app.perf_frame_ms > bad_ms) {
    app.perf_throttle = std::min(1.0, app.perf_throttle + 0.08);
  } else if (app.perf_frame_ms < good_ms) {
    app.perf_throttle = std::max(0.0, app.perf_throttle - 0.04);
  }
  app.perf_throttle_active = (app.perf_throttle > 0.02);
}

/* Headless single-frame render to a PNG. Used for visual verification / CI:
 *   dynsys --shot <preset-substring> <view> [--frames N] [--out FILE] [--mode M]
 * <view> is one of: line phase scene3d bifurcation fractal bridge basins
 *                   paramscan continuation ifs limitcycle
 * For the bridge, --mode is one of: projection | side | current
 * Requires a usable (possibly virtual, e.g. Xvfb) display. */
int run_shot(int argc, char **argv) {
  std::string preset_q, view_q = "phase", out_file, bridge_mode_q = "projection", bridge_family_q = "", fractalmode_q = "";
  int frames = 90;
  for (int i = 2; i < argc; ++i) {
    if (std::strcmp(argv[i], "--frames") == 0 && i + 1 < argc) frames = std::atoi(argv[++i]);
    else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) out_file = argv[++i];
    else if (std::strcmp(argv[i], "--mode") == 0 && i + 1 < argc) bridge_mode_q = argv[++i];
    else if (std::strcmp(argv[i], "--family") == 0 && i + 1 < argc) bridge_family_q = argv[++i];
    else if (std::strcmp(argv[i], "--fractalmode") == 0 && i + 1 < argc) fractalmode_q = argv[++i];
    else if (preset_q.empty()) preset_q = argv[i];
    else view_q = argv[i];
  }

  AppState app{};
  int idx = 0;
  if (!preset_q.empty()) {
    for (int i = 0; i < kPresetCount; ++i) {
      std::string n = kPresets[i].name;
      std::string lo; for (char c : n) lo += (char)std::tolower((unsigned char)c);
      std::string q; for (char c : preset_q) q += (char)std::tolower((unsigned char)c);
      if (lo.find(q) != std::string::npos) { idx = i; break; }
    }
  }
  apply_preset(app, idx);
  std::string err;
  if (!compile_system(app, app.system_input, &err)) {
    std::fprintf(stderr, "compile failed: %s\n", err.c_str()); return EXIT_FAILURE;
  }
  allocate_points(app);
  reset_simulation(app);

  auto pick = [&](const char *k){ return view_q == k; };
  if (pick("line")) app.active_view = AppState::ActiveView::Line1D;
  else if (pick("phase")) app.active_view = AppState::ActiveView::Phase2D;
  else if (pick("scene3d")) app.active_view = AppState::ActiveView::Scene3D;
  else if (pick("bifurcation")) app.active_view = AppState::ActiveView::Bifurcation;
  else if (pick("fractal")) app.active_view = AppState::ActiveView::Fractal;
  if (!fractalmode_q.empty()) {
    if (fractalmode_q == "buddhabrot" || fractalmode_q == "buddha") app.fractal_mode = AppState::FractalMode::Buddhabrot;
    else if (fractalmode_q == "state" || fractalmode_q == "julia") app.fractal_mode = AppState::FractalMode::StateSpace;
    else app.fractal_mode = AppState::FractalMode::ParameterSpace;
    /* frame on the preset window, then mark init done so entry won't reset mode */
    if (app.home_view_set) {
      app.fractal_xmin = app.home_x_min; app.fractal_xmax = app.home_x_max;
      app.fractal_ymin = app.home_y_min; app.fractal_ymax = app.home_y_max;
    }
    app.fractal_view_init = true;
    app.fractal_dirty = true;
  }
  else if (pick("bridge")) app.active_view = AppState::ActiveView::Scene3DBridge;
  else if (pick("basins")) app.active_view = AppState::ActiveView::Basins;
  else if (pick("paramscan")) app.active_view = AppState::ActiveView::ParamScan2D;
  else if (pick("continuation")) app.active_view = AppState::ActiveView::Continuation;
  else if (pick("ifs")) app.active_view = AppState::ActiveView::IFS;
  else if (pick("limitcycle")) app.active_view = AppState::ActiveView::LimitCycle;

  if (bridge_mode_q == "side") app.bridge_mode = AppState::BridgeMode::Classic;
  else if (bridge_mode_q == "current") app.bridge_mode = AppState::BridgeMode::CurrentSystem;
  else app.bridge_mode = AppState::BridgeMode::ProjectionSolid;
  /* enforce the same gating as the live UI: projection-solid only when the
   * loaded system actually has a unified bridge; otherwise lift the system. */
  if (app.bridge_mode == AppState::BridgeMode::ProjectionSolid &&
      bridge_family_for_system(app) < 0)
    app.bridge_mode = AppState::BridgeMode::CurrentSystem;
  if (bridge_family_q == "cubic") app.bridge_family = AppState::BridgeFamily::Cubic;
  else if (bridge_family_q == "sine") app.bridge_family = AppState::BridgeFamily::Sine;
  else if (bridge_family_q == "quadratic") app.bridge_family = AppState::BridgeFamily::Quadratic;
  /* else: leave whatever compile_system auto-detected for the loaded system */

  app.window_width = 1280; app.window_height = 1000; /* larger canvas for crisp shots */
  GLFWwindow *window = nullptr;
  if (!init_glfw_window(app, &window)) { std::fprintf(stderr, "no GL window (need a display, e.g. Xvfb)\n"); return EXIT_FAILURE; }
  init_imgui(window);
  app.shader_program = create_shader_program(kVertexShaderSrc, kFragmentShaderSrc);
  glGenVertexArrays(1, &app.vao); glGenBuffers(1, &app.vbo); glGenBuffers(1, &app.cbo);
  upload_buffers(app); update_projection(app);

  for (int f = 0; f < frames; ++f) {
    glfwPollEvents();
    if (!app.paused && app.mode != SystemMode::IFS) {
      char e[256] = {0};
      for (int s = 0; s < app.steps_per_frame; ++s) {
        State nx{}; if (!step_state(app, app.current, &nx, e, sizeof e)) break;
        app.current = nx; app.history.push_back(app.current);
      }
      upload_buffers(app);
    }
    ImGui_ImplOpenGL3_NewFrame(); ImGui_ImplGlfw_NewFrame(); ImGui::NewFrame();
    ensure_valid_view(app);
    if (app.active_view == AppState::ActiveView::Scene3D && view_valid(app, AppState::ActiveView::Scene3D))
      render_scene_background(app);
    else if (app.active_view == AppState::ActiveView::Scene3DBridge)
      render_bridge_scene(app);
    else
      draw_scene(app);
    if (f == frames - 1) app.capturing_clean = true;
    draw_gui(app);
    ImGui::Render();
    glViewport(0, 0, app.window_width, app.window_height);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    if (f == frames - 1) {
      std::string p = capture_screenshot_png(app);
      app.capturing_clean = false;
      if (!p.empty() && !out_file.empty()) { std::rename(p.c_str(), out_file.c_str()); p = out_file; }
      std::printf("shot: preset='%s' view=%s -> %s\n", kPresets[idx].name, view_q.c_str(), p.c_str());
    }
    glfwSwapBuffers(window);
  }
  ImGui_ImplOpenGL3_Shutdown(); ImGui_ImplGlfw_Shutdown(); ImGui::DestroyContext();
  glfwDestroyWindow(window); glfwTerminate();
  return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
  if (argc >= 2 && std::strcmp(argv[1], "--headless") == 0) {
    return run_headless(argc, argv);
  }
  if (argc >= 2 && std::strcmp(argv[1], "--shot") == 0) {
    return run_shot(argc, argv);
  }
  AppState app{};
  set_lorenz(app);
  if (!compile_system(app, app.system_input, &app.parse_error)) {
    std::fprintf(stderr, "%s\n", app.parse_error.c_str());
    return EXIT_FAILURE;
  }
  allocate_points(app);
  reset_simulation(app);
  GLFWwindow *window = nullptr;
  if (!init_glfw_window(app, &window)) {
    if (app.arena_ready) arena_destroy(&app.system_arena);
    return EXIT_FAILURE;
  }
  init_imgui(window);
  app.shader_program = create_shader_program(kVertexShaderSrc, kFragmentShaderSrc);
  glGenVertexArrays(1, &app.vao);
  glGenBuffers(1, &app.vbo);
  glGenBuffers(1, &app.cbo);
  upload_buffers(app);
  update_projection(app);
  double last_fps_time = glfwGetTime();
  int frame_count = 0;
  double prev_frame_time = glfwGetTime();
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();
    process_keyboard(app, window);
    /* measure the previous frame's wall time and update the governor BEFORE
     * doing this frame's heavy work, so the throttle reflects recent load. */
    const double frame_start = glfwGetTime();
    update_perf_governor(app, (frame_start - prev_frame_time) * 1000.0);
    prev_frame_time = frame_start;
    if (!app.paused && app.mode != SystemMode::IFS) {
      /* throttle integration substeps: at full throttle do ~15% of the
       * requested steps (min 1). Keeps trajectory cost bounded when the
       * frame budget is blown. */
      int eff_steps = app.steps_per_frame;
      if (app.perf_throttle > 0.02) {
        eff_steps = (int)std::ceil(app.steps_per_frame * (1.0 - 0.85 * app.perf_throttle));
        if (eff_steps < 1) eff_steps = 1;
      }
      for (int i = 0; i < eff_steps; ++i) {
        if (!calculate_next_point(app)) break;
      }
    }
    const double now = glfwGetTime();
    ++frame_count;
    if (now - last_fps_time >= 1.0) {
      app.fps = frame_count;
      frame_count = 0;
      last_fps_time = now;
    }
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    /* PHASE-A: render the active full-area view. 3D goes straight to the
     * backbuffer; the others paint via the ImGui background draw list
     * inside draw_gui. Snap to a valid view first (e.g. after switching
     * to a system of different dimension). */
    ensure_valid_view(app);
    if (app.active_view == AppState::ActiveView::Scene3D &&
        view_valid(app, AppState::ActiveView::Scene3D)) {
      render_scene_background(app);
    } else if (app.active_view == AppState::ActiveView::Scene3DBridge) {
      render_bridge_scene(app);
    } else {
      draw_scene(app);
    }
    draw_gui(app);
    ImGui::Render();
    glViewport(0, 0, app.window_width, app.window_height);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    /* Serviced here. We render ONE extra frame with the chrome suppressed
     * (capturing_clean) so the PNG is the full-size plot with nothing hidden
     * behind the toolbar/panel, then capture that clean backbuffer. */
    if (app.capture_request) {
      app.capture_request = false;
      app.capturing_clean = true;
      /* re-render this frame without the GUI chrome */
      ImGui_ImplOpenGL3_NewFrame();
      ImGui_ImplGlfw_NewFrame();
      ImGui::NewFrame();
      ensure_valid_view(app);
      if (app.active_view == AppState::ActiveView::Scene3D &&
          view_valid(app, AppState::ActiveView::Scene3D)) {
        render_scene_background(app);
      } else if (app.active_view == AppState::ActiveView::Scene3DBridge) {
        render_bridge_scene(app);
      } else {
        draw_scene(app);
      }
      draw_gui(app); /* returns early (chrome suppressed) but paints backgrounds */
      ImGui::Render();
      glViewport(0, 0, app.window_width, app.window_height);
      ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
      std::string p = capture_screenshot_png(app);
      app.capturing_clean = false;
      app.screenshot_msg = p.empty() ? "PNG save failed (check write permissions in the run directory)"
                                      : ("saved " + p);
      app.screenshot_msg_timer = 240; /* ~4s at 60fps */
    }
    glfwSwapBuffers(window);
  }
  glDeleteVertexArrays(1, &app.vao);
  glDeleteBuffers(1, &app.vbo);
  glDeleteBuffers(1, &app.cbo);
  glDeleteProgram(app.shader_program);
  if (app.scene_fbo) glDeleteFramebuffers(1, &app.scene_fbo);
  if (app.scene_tex) glDeleteTextures(1, &app.scene_tex);
  if (app.fractal_tex) glDeleteTextures(1, &app.fractal_tex);
  if (app.basin_tex) glDeleteTextures(1, &app.basin_tex);
  if (app.scan_tex) glDeleteTextures(1, &app.scan_tex);
  if (app.ifs_tex) glDeleteTextures(1, &app.ifs_tex);
  if (app.ifs_arena_init) { arena_destroy(&app.ifs_arena); app.ifs_arena_init = false; }
  if (app.bridge_vbo) glDeleteBuffers(1, &app.bridge_vbo);
  if (app.bridge_cbo) glDeleteBuffers(1, &app.bridge_cbo);
  if (app.bridge_vao) glDeleteVertexArrays(1, &app.bridge_vao);
  if (app.scene_depth) glDeleteRenderbuffers(1, &app.scene_depth);
  shutdown_imgui();
  glfwDestroyWindow(window);
  glfwTerminate();
  if (app.arena_ready) arena_destroy(&app.system_arena);
  return EXIT_SUCCESS;
}
