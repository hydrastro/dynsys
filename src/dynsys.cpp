#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <cglm/call.h>
#include <cglm/cglm.h>

#include <algorithm>
#include <array>
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
#include <deque>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
#include "arena.h"
#include "ast.h"
#include "pratt.h"
}

#include "expr_ir.h"

#include "analysis.h"
#include "expr_ir_ad.h"

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

const char *mode_name(SystemMode mode) {
  return mode == SystemMode::Map ? "map" : "ode";
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
  int window_width = 1100;
  int window_height = 820;
  std::string screenshot_msg;     /* transient "saved <path>" toast */
  int screenshot_msg_timer = 0;

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
  enum class ActiveView { Line1D, Phase2D, Scene3D, Bifurcation, Fractal, Scene3DBridge, Basins, ParamScan2D, Continuation };
  ActiveView active_view = ActiveView::Phase2D;
  bool show_side_panel = true;
  float window_toolbar_h = 32.0f;

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
  enum class FractalMode { ParameterSpace /*Mandelbrot-like*/, StateSpace /*Julia-like*/ };
  FractalMode fractal_mode = FractalMode::ParameterSpace;
  double fractal_xmin = -2.5, fractal_xmax = 1.0;   /* view window (re) */
  double fractal_ymin = -1.5, fractal_ymax = 1.5;   /* view window (im) */
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
  int bridge_mandel_res = 360;   /* c-plane sampling resolution */
  int bridge_bif_slices = 600;   /* bifurcation parameter slices */
  int bridge_bif_keep = 120;     /* attractor points kept per slice */
  float bridge_height = 18.0f;   /* vertical scale of the bifurcation */
  bool bridge_show_mandelbrot = true;
  bool bridge_show_bifurcation = true;

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
  int scan_px_index = 0;     /* parameter on the horizontal axis */
  int scan_py_index = 1;     /* parameter on the vertical axis */
  double scan_xmin = 0, scan_xmax = 4;   /* param-x range */
  double scan_ymin = 0, scan_ymax = 4;   /* param-y range */
  int scan_transient = 500;
  int scan_iterations = 500;
  bool scan_view_init = false;
  double scan_lyap_min = 0, scan_lyap_max = 0; /* observed range, for the legend */

  bool fixed_ready = false;
  State fixed_point{};
  double fixed_residual = 0.0;
  std::vector<double> fixed_jacobian;
  std::string fixed_classification;
  std::vector<std::pair<double, double>> fixed_eigenvalues;
  std::vector<PhaseTrajectory> separatrix_curves;
  int separatrix_steps = 2500;
  double separatrix_epsilon = 1e-4;

  /* PHASE1/2: general (N-D) equilibrium analysis + equilibrium
   * continuation. fixed_eigenvalues above is kept for the legacy 2D
   * panel; this richer classification works in any dimension. */
  dynsys::analysis::Classification fixed_general{};
  bool fixed_general_ready = false;
  char cont_param[64] = "rho";
  double cont_p_min = -10.0;
  double cont_p_max = 110.0;
  double cont_h0 = 0.05;
  int cont_max_points = 600;
  int cont_direction = 1;
  bool cont_detect_fold = true;
  bool cont_detect_hopf = true;
  dynsys::analysis::Branch cont_branch{};
  bool cont_ready = false;

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
    if (name == "t" || name == "pi" || name == "e" || is_builtin_function_name(name)) {
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
  return name == "t" || name == "pi" || name == "e";
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
    if (std::strcmp(name, "pi") == 0) { *out = M_PI; return true; }
    if (std::strcmp(name, "e") == 0) { *out = M_E; return true; }

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
      if (!is_identifier(name) || is_reserved_value_name(name) || is_state_value_name(next_state_names, name) || is_builtin_function_name(name)) {
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
      else {
        *error = "line " + std::to_string(line_no + 1) + ": mode must be 'ode' or 'map'";
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
  }

  /* PHASE0: equation-tied analysis is always stale after a recompile —
   * the located equilibrium, its Jacobian/spectrum, and the saddle
   * separatrices were computed for the previous vector field. Clear
   * them. */
  app.fractal_dirty = true; /* PHASE C: recompute the fractal image too */
  app.basin_dirty = true;   /* PHASE B/C: recompute basins for the new system */
  app.scan_dirty = true; app.scan_view_init = false; /* PHASE B: recompute 2-param scan */
  app.fixed_ready = false;
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
  }
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
    {"Three-Scroll Unified", "Continuous 3D", "# Three-Scroll Unified Chaotic System\nmode = ode\nintegrator = rk4\nparam a = 32.48 [0,60]\nparam b = 45.84 [0,70]\nparam c = 1.18 [0,5]\nparam d = 0.13 [0,1]\nparam e = 0.57 [0,2]\nparam f = 14.7 [0,30]\ndx = a * (y - x) + d * x * z\ndy = b * x - x * z + f * y\ndz = c * z + x * y - e * pow(x, 2)\n", {0.1,0.1,0.1,0}, 0.0005, 16, Integrator::RK4,1.0f,0,0,0},
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
param e = 0.5 [0,2]
param m = 0.3 [0,2]
view2d = 0, 12, 0, 4
dx = r * x * (1 - x / K) - a * x * y / (1 + a * h * x)
dy = e * a * x * y / (1 + a * h * x) - m * y
)dyn", {5.0,1.0,0.0,0}, 0.02, 4, Integrator::RK4,1.0f,0,0,0},
    {"Simple pendulum (undamped)", "Phase-plane examples", R"dyn(# Undamped pendulum: nested closed orbits + separatrix
state theta, omega
mode = ode
integrator = rk4
view2d = -6.5, 6.5, -3, 3
dtheta = omega
domega = 0 - sin(theta)
)dyn", {0.5,0.0,0.0,0}, 0.02, 4, Integrator::RK4,1.0f,0,0,0},
    {"Henon map", "Discrete maps", "# Hénon map\nmode = map\nparam a = 1.4 [0,2]\nparam b = 0.3 [0,1]\nx_next = 1 - a * x * x + y\ny_next = b * x\nz_next = z\nobserve radius = sqrt(x*x + y*y)\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Lozi map", "Discrete maps", "# Lozi map (piecewise-linear Henon cousin; chaotic)\nstate x, y\nmode = map\nparam a = 1.7 [0,2]\nparam b = 0.5 [0,1]\nview2d = -2, 2, -1, 1\nx_next = 1 - a * abs(x) + y\ny_next = b * x\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Ikeda map", "Discrete maps", "# Ikeda map (laser cavity; spiral chaotic attractor)\nstate x, y\nmode = map\nparam u = 0.918 [0,1]\nview2d = -0.5, 2, -2.5, 1\nx_next = 1 + u * (x * cos(0.4 - 6 / (1 + x*x + y*y)) - y * sin(0.4 - 6 / (1 + x*x + y*y)))\ny_next = u * (x * sin(0.4 - 6 / (1 + x*x + y*y)) + y * cos(0.4 - 6 / (1 + x*x + y*y)))\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Logistic map", "Discrete maps", "# Logistic map x_{n+1} = r x (1-x)  (genuine 1D map)\n# The classic period-doubling route to chaos.\nstate x\nmode = map\nparam r = 3.9 [0,4]\nview2d = 0, 1, 0, 1\nx_next = r * x * (1 - x)\nobserve x_value = x\n", {0.2,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Tent map", "Discrete maps", "# Tent map: x_next = mu * min(x, 1-x)\n# 1D map; chaotic for mu near 2, x stays in [0,1].\nstate x\nmode = map\nparam mu = 1.9 [0,2]\nview2d = 0, 1, 0, 1\nx_next = mu * min(x, 1 - x)\nobserve x_value = x\n", {0.35,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Sine map", "Discrete maps", "# Sine map x_{n+1} = r sin(pi x)  (1D, period-doubling like logistic)\nstate x\nmode = map\nparam r = 0.9 [0,1]\nview2d = 0, 1, 0, 1\nx_next = r * sin(3.141592653589793 * x)\nobserve x_value = x\n", {0.4,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Gauss / mouse map", "Discrete maps", "# Gauss (mouse) map x_{n+1} = exp(-alpha x^2) + beta  (1D)\nstate x\nmode = map\nparam alpha = 6.2 [1,10]\nparam beta = 0 - 0.5 [-1,1]\nview2d = -1, 1.5, -1, 1.5\nx_next = exp(0 - alpha * x * x) + beta\nobserve x_value = x\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Cubic map", "Discrete maps", "# Cubic map x_{n+1} = a x - x^3  (1D, symmetric period-doubling)\nstate x\nmode = map\nparam a = 2.8 [0,3]\nview2d = -2, 2, -2, 2\nx_next = a * x - x*x*x\nobserve x_value = x\n", {0.2,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Gingerbreadman map", "Discrete maps", "# Gingerbreadman map (piecewise-linear, chaotic)\nstate x, y\nmode = map\nview2d = -8, 8, -8, 8\nx_next = 1 - y + abs(x)\ny_next = x\n", {0.5,3.7,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Standard (Chirikov) map", "Discrete maps", "# Chirikov standard map on the cylinder\nstate theta, p\nmode = map\nparam K = 0.97 [0,5]\nview2d = 0, 6.2832, -3.1416, 3.1416\np_next = p + K * sin(theta)\ntheta_next = theta + p + K * sin(theta)\n", {3.0,0.1,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Ikeda map", "Discrete maps", "# Ikeda map\nmode = map\nparam u = 0.918 [0,1]\ntheta = 0.4 - 6 / (1 + x*x + y*y)\nx_next = 1 + u * (x * cos(theta) - y * sin(theta))\ny_next = u * (x * sin(theta) + y * cos(theta))\nz_next = z\n", {0.1,0.1,0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Tinkerbell map", "Discrete maps", "# Tinkerbell map\nmode = map\nparam a = 0.9 [0,2]\nparam b = 0 - 0.6013 [-2,2]\nparam c = 2.0 [0,4]\nparam d = 0.5 [0,2]\nx_next = x*x - y*y + a*x + b*y\ny_next = 2*x*y + c*x + d*y\nz_next = z\n", {0.1,0.1,0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Newton fractal (z^3 - 1)", "Fractals", "# Newton's method for z^3 = 1, as a 2D map.\n# Three roots -> three basins with fractal boundaries.\n# Switch to the 'basins' view after loading.\nstate x, y\nmode = map\nview2d = -2, 2, -2, 2\n# z^2 = a + i b\nx_next = x - ( ( (x*(x*x-y*y) - y*(2*x*y)) - 1 )*(3*(x*x-y*y)) + ( x*(2*x*y) + y*(x*x-y*y) )*(3*(2*x*y)) ) / ( 9*((x*x-y*y)*(x*x-y*y) + (2*x*y)*(2*x*y)) )\ny_next = y - ( ( x*(2*x*y) + y*(x*x-y*y) )*(3*(x*x-y*y)) - ( (x*(x*x-y*y) - y*(2*x*y)) - 1 )*(3*(2*x*y)) ) / ( 9*((x*x-y*y)*(x*x-y*y) + (2*x*y)*(2*x*y)) )\n", {0.5,0.5,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Complex quadratic (Mandelbrot/Julia)", "Fractals", "# z -> z^2 + c  (x+iy)^2 + (cx+i cy)\n# Fractal view: parameter space = Mandelbrot, state space = Julia.\nstate x, y\nmode = map\nparam cx = 0 [-2,2]\nparam cy = 0 [-2,2]\nview2d = -2, 2, -2, 2\nx_next = x*x - y*y + cx\ny_next = 2*x*y + cy\n", {0.0,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Henon map 2D", "N-dimensional examples", "# Hénon map as a true 2D state\nstate x, y\nmode = map\nview2d = -1.5, 1.5, -0.5, 0.5\nparam a = 1.4 [0,2]\nparam b = 0.3 [0,1]\nplot3d = x, y, 0\nx_next = 1 - a * x * x + y\ny_next = b * x\nobserve radius = sqrt(x*x + y*y)\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"4D coupled oscillator demo", "N-dimensional examples", "# Four-dimensional demo system\nstate x, y, z, w\nmode = ode\nintegrator = rk4\nparam a = 0.1 [0,1]\nplot3d = x, y, w\nobserve r4 = sqrt(x*x + y*y + z*z + w*w)\nsection = w\nsection_direction = positive\nsection_plot = x, z\ninitial x = 0.1\ninitial y = 0\ninitial z = 0\ninitial w = 0.2\ndx = y\ndy = 0 - x - a*y + z\ndz = w\ndw = 0 - z - a*w + x\n", {0.1,0.0,0.0,0}, 0.01, 3, Integrator::RK4,1.0f,0,0,0},
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
void load_logistic_bifurcation_demo(AppState &app) {
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
  app->zoom *= static_cast<float>(std::pow(1.1, yoffset));
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

/* PHASE-A: which full-area views make sense for the current system.
 *   1D line/cobweb : exactly 1 effective dimension
 *   2D phase       : >= 2 effective dimensions
 *   3D scene       : >= 3 effective dimensions
 *   bifurcation    : any system that has at least one parameter
 * (effective dimension respects the Force2D/Force3D override.) */
bool view_valid(const AppState &app, AppState::ActiveView v) {
  const int d = effective_dimension(app);
  switch (v) {
    case AppState::ActiveView::Line1D:      return d == 1;
    case AppState::ActiveView::Phase2D:     return d >= 2;
    case AppState::ActiveView::Scene3D:     return system_is_3d(app);
    case AppState::ActiveView::Bifurcation: return !app.params.empty();
    case AppState::ActiveView::Fractal:     return app.mode == SystemMode::Map && app.state_names.size() >= 2;
    case AppState::ActiveView::Scene3DBridge:
      /* The Mandelbrot<->bifurcation bridge is meaningful for 1D maps (the
       * logistic family), whose period-doubling cascade embeds in the
       * Mandelbrot needle. It isn't relevant to general ODEs/maps, so only
       * offer it there. */
      return app.mode == SystemMode::Map && effective_dimension(app) == 1;
    case AppState::ActiveView::Basins:      return app.state_names.size() >= 2;
    case AppState::ActiveView::ParamScan2D: return app.params.size() >= 2;
    case AppState::ActiveView::Continuation: return app.mode == SystemMode::ODE && !app.params.empty();
  }
  return false;
}

/* A sensible default view for a freshly compiled system. */
AppState::ActiveView default_view_for(const AppState &app) {
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
      app.zoom *= std::pow(1.15f, io.MouseWheel);
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

void zoom_phase_bounds(AppState &app, PlotBounds b, const Point2 &anchor_in, double factor) {
  /* Guard against a non-finite anchor (which could come from a prior bad
   * frame) — fall back to the view center so zoom can always recover. */
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

void pan_phase_bounds(AppState &app, PlotBounds b, double dx, double dy) {
  b.xmin += dx; b.xmax += dx;
  b.ymin += dy; b.ymax += dy;
  set_manual_phase_bounds(app, b);
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
      const double factor = std::pow(0.8, static_cast<double>(io.MouseWheel));
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
  const Point2 cp = plot_to_screen(bounds, p0, w, h, state_at(app.current, ix), state_at(app.current, iy));
  draw->AddCircleFilled(ImVec2(static_cast<float>(cp.x), static_cast<float>(cp.y)), 4.0f, IM_COL32(255, 255, 255, 255));

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
      zoom_phase_bounds(app, b, mouse_data, std::pow(0.8, (double)io.MouseWheel));
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
      zoom_phase_bounds(app, bounds, mouse_data, std::pow(0.8, static_cast<double>(io.MouseWheel)));
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

  const Point2 cp = plot_to_screen(bounds, p0, w, h, state_at(app.current, ix), state_at(app.current, iy));
  draw->AddCircleFilled(ImVec2((float)cp.x, (float)cp.y), 4.0f, IM_COL32(255, 255, 255, 255));

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

void compute_fractal_image(AppState &app, int W, int H, std::vector<uint32_t> &out, int step = 1) {
  if (W < 2) W = 2;
  if (H < 2) H = 2;
  if (step < 1) step = 1;
  out.assign((size_t)W * H, 0xff000000u);
  const size_t n = app.state_names.size();
  if (n < 2) return; /* needs a 2D plane */
  if (app.mode != SystemMode::Map) return; /* escape-time is for maps */

  const int maxit = std::max(10, app.fractal_max_iter);
  const double R2 = app.fractal_escape_r * app.fractal_escape_r;
  const double x0 = app.fractal_xmin, x1 = app.fractal_xmax;
  const double y0 = app.fractal_ymin, y1 = app.fractal_ymax;

  const bool param_mode = (app.fractal_mode == AppState::FractalMode::ParameterSpace);
  const int cxi = std::max(0, std::min(app.fractal_param_cx_index, (int)app.params.size() - 1));
  const int cyi = std::max(0, std::min(app.fractal_param_cy_index, (int)app.params.size() - 1));

  /* Save params we will temporarily override (param mode). */
  std::vector<double> saved = app.param_values;

  /* Reused scratch states (avoid per-pixel heap allocation). */
  State s = make_state_like(n, 0.0);
  State cur = make_state_like(n, 0.0);
  State nx = make_state_like(n, 0.0);

  char err[128] = {0};
  for (int py = 0; py < H; py += step) {
    const double b_im = y0 + (y1 - y0) * (double)py / (H - 1);
    for (int px = 0; px < W; px += step) {
      const double a_re = x0 + (x1 - x0) * (double)px / (W - 1);

      if (param_mode) {
        /* plane = (param cx, param cy); initial state = system start */
        for (size_t i = 0; i < n; ++i) set_state_at(s, i, state_at(app.start, i));
        if ((size_t)cxi < app.param_values.size()) app.param_values[(size_t)cxi] = a_re;
        if ((size_t)cyi < app.param_values.size()) app.param_values[(size_t)cyi] = b_im;
      } else {
        /* plane = initial (x, y); params fixed at current values */
        if (px == 0 && py == 0) app.param_values = saved;
        set_state_at(s, 0, a_re);
        set_state_at(s, 1, b_im);
        for (size_t i = 2; i < n; ++i) set_state_at(s, i, state_at(app.start, i));
      }

      int it = 0;
      double r2 = 0.0;
      for (size_t i = 0; i < n; ++i) set_state_at(cur, i, state_at(s, i));
      for (; it < maxit; ++it) {
        if (!step_map_state(app, cur, &nx, err, sizeof(err))) { it = maxit; break; }
        const double xx = state_at(nx, 0), yy = state_at(nx, 1);
        r2 = xx * xx + yy * yy;
        if (!std::isfinite(r2) || r2 > R2) { ++it; break; }
        std::swap(cur.v, nx.v); cur.t = nx.t;
      }

      uint32_t color = 0xff101014u; /* "in the set" */
      if (it < maxit && r2 > R2) {
        double mu = it;
        if (app.fractal_smooth && r2 > 1.0) {
          /* continuous (smooth) iteration count */
          mu = it + 1.0 - std::log(std::log(std::sqrt(r2))) / std::log(2.0);
        }
        const double t = std::fmod(mu * 0.04, 1.0);
        color = fractal_palette(t < 0 ? t + 1.0 : t);
      }
      /* fill the step x step block so a coarse pass covers the whole image */
      for (int by = py; by < py + step && by < H; ++by)
        for (int bx = px; bx < px + step && bx < W; ++bx)
          out[(size_t)by * W + bx] = color;
    }
  }
  app.param_values = saved; /* restore */
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
    app.fractal_dirty = true;
  }

  /* interaction: pan/zoom the complex window; mark dirty to recompute */
  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse) {
    auto invx = [&](float pxf){ return app.fractal_xmin + (double)pxf / w * (app.fractal_xmax - app.fractal_xmin); };
    auto invy = [&](float pyf){ return app.fractal_ymin + (double)pyf / h * (app.fractal_ymax - app.fractal_ymin); };
    if (io.MouseWheel != 0.0f) {
      const double mx = invx(io.MousePos.x), my = invy(io.MousePos.y);
      const double f = std::pow(0.8, (double)io.MouseWheel);
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
  draw->AddText(ImVec2(14, 12), IM_COL32(235, 235, 240, 235), hud);
  if (!io.WantCaptureMouse)
    draw->AddText(ImVec2(14, h - 24), IM_COL32(150, 150, 160, 200),
                  "drag: pan   wheel: zoom   (set mode/params/iterations in the Setup tab)");
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
  if (label == -1) return IM_COL32(8, 8, 12, 255);     /* diverged to infinity */
  if (label == -2) return IM_COL32(64, 64, 72, 255);   /* did not settle (chaotic attractor) */
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

void compute_basin_image(AppState &app, int W, int H, std::vector<uint32_t> &out) {
  if (W < 2) W = 2;
  if (H < 2) H = 2;
  out.assign((size_t)W * H, 0xff000000u);
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
  opt.width = W; opt.height = H;
  opt.max_steps = std::max(50, app.basin_steps);
  opt.cluster_tol = app.basin_cluster_tol;
  opt.settle_tol = (app.mode == SystemMode::Map) ? 1e-6 : 1e-5;

  dynsys::analysis::BasinResult R = dynsys::analysis::compute_basins(advance, opt);
  app.integrator = saved_integrator; /* restore the user's choice */
  app.basin_attractor_count = (int)R.attractors.size();
  app.basin_n_converged = R.n_converged;
  app.basin_n_diverged = R.n_diverged;
  app.basin_n_nonconvergent = R.n_nonconvergent;
  if (!R.ok) return;
  const int nl = (int)R.attractors.size();
  for (int j = 0; j < H; ++j)
    for (int i = 0; i < W; ++i) {
      const size_t idx = (size_t)j * W + i;
      out[idx] = basin_color(R.cell_attractor[idx], nl, R.cell_speed[idx], app.basin_shade_speed);
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
   * a coordinate frame; any view change marks the image dirty. */
  const size_t ix = (size_t)app.phase_x_index, iy = (size_t)app.phase_y_index;
  PlotBounds b = sanitize_bounds(current_phase_bounds(app, ix, iy));
  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse) {
    const ImVec2 p0(0, 0);
    const Point2 m = screen_to_plot(b, p0, w, h, io.MousePos);
    if (io.MouseWheel != 0.0f) { zoom_phase_bounds(app, b, m, std::pow(0.8, (double)io.MouseWheel)); app.basin_dirty = true; }
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) && !ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      const double ddx = (double)io.MouseDelta.x / w * (b.xmax - b.xmin);
      const double ddy = (double)io.MouseDelta.y / h * (b.ymax - b.ymin);
      pan_phase_bounds(app, b, -ddx, ddy);
      app.basin_dirty = true;
    }
    if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) { app.phase_bounds_valid = false; app.basin_dirty = true; }
  }

  const int CW = std::min(app.window_width, 700);
  const int CH = std::min(app.window_height, 560);
  const bool bforce = (app.basin_tex == 0 || app.basin_tex_w != CW || app.basin_tex_h != CH);
  if (app.basin_dirty) { app.basin_settle = 6; app.basin_dirty = false; }
  bool bcompute = bforce;
  if (app.basin_settle > 0) { if (--app.basin_settle == 0) bcompute = true; }
  if (bcompute) {
    std::vector<uint32_t> img;
    compute_basin_image(app, CW, CH, img);
    if (app.basin_tex == 0) glGenTextures(1, &app.basin_tex);
    glBindTexture(GL_TEXTURE_2D, app.basin_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, CW, CH, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    app.basin_tex_w = CW; app.basin_tex_h = CH;
    app.basin_dirty = false;
  }
  if (app.basin_tex != 0)
    draw->AddImage((ImTextureID)(uintptr_t)app.basin_tex, ImVec2(0, 0), ImVec2(w, h));

  char hud[320];
  std::snprintf(hud, sizeof(hud),
                "Basins of attraction — %d basin(s)  |  %s vs %s  |  converged %ld, escaped %ld, non-convergent %ld",
                app.basin_attractor_count, app.state_names[ix].c_str(), app.state_names[iy].c_str(),
                app.basin_n_converged, app.basin_n_diverged, app.basin_n_nonconvergent);
  draw->AddText(ImVec2(14, 12), IM_COL32(235, 235, 240, 235), hud);

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

void compute_scan_image(AppState &app, int W, int H, std::vector<uint32_t> &out) {
  if (W < 2) W = 2;
  if (H < 2) H = 2;
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

  for (int j = 0; j < H; ++j) {
    const double py = app.scan_ymin + (app.scan_ymax - app.scan_ymin) * (double)j / (H - 1);
    for (int i = 0; i < W; ++i) {
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
      ly[(size_t)j * W + i] = (float)val;
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
      const double f = std::pow(0.8, (double)io.MouseWheel);
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

  /* scans are expensive — keep the grid modest and debounce recompute */
  const int CW = std::min(app.window_width, 480);
  const int CH = std::min(app.window_height, 380);
  const bool sforce = (app.scan_tex == 0 || app.scan_tex_w != CW || app.scan_tex_h != CH);
  if (app.scan_dirty) { app.scan_settle = 8; app.scan_dirty = false; }
  bool scompute = sforce;
  if (app.scan_settle > 0) { if (--app.scan_settle == 0) scompute = true; }
  if (scompute) {
    std::vector<uint32_t> img;
    compute_scan_image(app, CW, CH, img);
    if (app.scan_tex == 0) glGenTextures(1, &app.scan_tex);
    glBindTexture(GL_TEXTURE_2D, app.scan_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, CW, CH, 0, GL_RGBA, GL_UNSIGNED_BYTE, img.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    app.scan_tex_w = CW; app.scan_tex_h = CH;
    app.scan_dirty = false;
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
  draw->AddText(ImVec2(14, 12), IM_COL32(235, 235, 240, 235), hud);
  draw->AddText(ImVec2(14, 30), IM_COL32(170, 170, 180, 220),
                "blue/dark = periodic (shrimps);  warm/bright = chaotic.  drag: pan  wheel: zoom");
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
}  // namespace bridge_detail

void build_bridge_geometry(AppState &app) {
  std::vector<Point> pts;
  std::vector<float> cols;
  pts.reserve(200000);
  cols.reserve(600000);

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
        /* draw interior solidly (deep blue) and a thin escaped halo */
        if (inside) {
          pts.push_back(Point{wx(re), 0.0f, wz(im)});
          cols.push_back(0.10f); cols.push_back(0.18f); cols.push_back(0.42f);
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
   * Map logistic r in [2.6,4.0] onto the Mandelbrot real-axis span where
   * the matching cascade lives: r=2.6..4.0 corresponds (via conjugacy) to
   * c roughly in [-0.75, -2.0]. We place the diagram along z=0 (im=0). */
  if (app.bridge_show_bifurcation) {
    const int slices = std::max(100, app.bridge_bif_slices);
    const int discard = 800;
    const int keep = std::max(20, app.bridge_bif_keep);
    /* r in [2.6,4]; c(r) = (2r - r^2)/4 is the logistic<->quadratic map's
     * parameter conjugacy (c = r/2 - r^2/4). This lands the cascade onto
     * the Mandelbrot needle. */
    auto c_of_r = [](double r) { return r / 2.0 - r * r / 4.0; };
    for (int s = 0; s < slices; ++s) {
      const double r = 2.6 + (4.0 - 2.6) * s / (slices - 1);
      double x = 0.5;
      for (int k = 0; k < discard; ++k) x = r * x * (1.0 - x);
      const double cre = c_of_r(r);
      for (int k = 0; k < keep; ++k) {
        x = r * x * (1.0 - x);
        if (!std::isfinite(x)) break;
        /* world position: x-axis from the conjugate c, z=0 (real axis),
         * y = attractor value lifted by bridge_height */
        const float wy = (float)x * app.bridge_height;
        pts.push_back(Point{wx(cre), wy, wz(0.0)});
        /* color by height: cyan low -> yellow high */
        const float t = (float)x;
        cols.push_back(0.2f + 0.8f * t);
        cols.push_back(0.7f + 0.3f * t);
        cols.push_back(0.9f - 0.6f * t);
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

  glEnable(GL_DEPTH_TEST);
  glClearColor(0.04f, 0.04f, 0.06f, 1.0f);
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

  /* camera interaction (same scheme as the attractor scene) */
  ImGuiIO &io = ImGui::GetIO();
  if (!io.WantCaptureMouse) {
    if (ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
      app.angle_y += io.MouseDelta.x * 0.4f;
      app.angle_x += io.MouseDelta.y * 0.4f;
    }
    if (io.MouseWheel != 0.0f) {
      app.zoom *= std::pow(1.15f, io.MouseWheel);
      app.zoom = std::max(0.05f, std::min(app.zoom, 100.0f));
    }
  }

  ImDrawList *draw = ImGui::GetBackgroundDrawList();
  draw->AddText(ImVec2(14, 12), IM_COL32(235, 235, 240, 235),
                "3D bridge — Mandelbrot set (flat) + logistic bifurcation (rising along the real axis)");
  draw->AddText(ImVec2(14, 30), IM_COL32(170, 170, 180, 220),
                "they line up because z^2+c and the logistic map are conjugate.  drag: rotate   wheel: zoom");
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

  /* fold (circle) and Hopf (diamond) markers */
  for (size_t i = 0; i < app.cont_special.size(); ++i) {
    if (app.cont_special[i] == 0) continue;
    const ImVec2 c(sx(app.cont_pp[i]), sy(app.cont_yy[i]));
    if (app.cont_special[i] == 1) {
      draw->AddCircle(c, 6.0f, IM_COL32(255, 220, 90, 255), 16, 2.0f);
      draw->AddText(ImVec2(c.x + 8, c.y - 6), IM_COL32(255, 220, 90, 255), "Fold (LP)");
    } else {
      const float r = 6.0f;
      draw->AddQuad(ImVec2(c.x, c.y-r), ImVec2(c.x+r, c.y), ImVec2(c.x, c.y+r), ImVec2(c.x-r, c.y),
                    IM_COL32(120, 200, 255, 255), 2.0f);
      draw->AddText(ImVec2(c.x + 8, c.y - 6), IM_COL32(120, 200, 255, 255), "Hopf (H)");
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
  draw->AddText(ImVec2(14, 12), IM_COL32(235, 235, 240, 235), hud);
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

  /* Restore the live parameter — the model closures already restore it
   * per call, but be explicit. */
  sync_param_values(app);

  size_t n_fold = 0, n_hopf = 0;
  for (size_t i : app.cont_branch.special_indices) {
    if (app.cont_branch.points[i].special ==
        dynsys::analysis::SpecialPointKind::Fold)
      ++n_fold;
    else if (app.cont_branch.points[i].special ==
             dynsys::analysis::SpecialPointKind::Hopf)
      ++n_hopf;
  }
  char msg[256];
  std::snprintf(msg, sizeof(msg),
                "continuation: %zu points, %zu fold(s), %zu Hopf point(s) — %s",
                app.cont_branch.points.size(), n_fold, n_hopf,
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
      if (pt.special == dynsys::analysis::SpecialPointKind::Fold) { sp = 1; ++app.cont_n_fold; }
      else if (pt.special == dynsys::analysis::SpecialPointKind::Hopf) { sp = 2; ++app.cont_n_hopf; }
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
      app.zoom *= std::pow(1.15f, io.MouseWheel);
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
  glReadPixels(0, 0, w, h, GL_RGBA, /*GL_UNSIGNED_BYTE*/ 0x1401, buf.data());
  /* flip vertically into an RGB buffer */
  std::vector<unsigned char> rgb((size_t)w * h * 3);
  for (int y = 0; y < h; ++y) {
    const unsigned char *src = &buf[(size_t)(h - 1 - y) * w * 4];
    unsigned char *dst = &rgb[(size_t)y * w * 3];
    for (int x = 0; x < w; ++x) {
      dst[x * 3 + 0] = src[x * 4 + 0];
      dst[x * 3 + 1] = src[x * 4 + 1];
      dst[x * 3 + 2] = src[x * 4 + 2];
    }
  }
  char path[256];
  std::time_t t = std::time(nullptr);
  std::tm *lt = std::localtime(&t);
  std::snprintf(path, sizeof(path), "dynsys_%04d%02d%02d_%02d%02d%02d.png",
                lt->tm_year + 1900, lt->tm_mon + 1, lt->tm_mday,
                lt->tm_hour, lt->tm_min, lt->tm_sec);
  if (png_write_rgb(path, rgb.data(), w, h)) return std::string(path);
  return "";
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

  /* Visible build stamp: if you don't see this line at the top of the
   * window, you are running an OLD binary — rebuild (make clean && make).*/
  ImGui::TextColored(ImVec4(0.45f, 0.95f, 0.55f, 1.0f), "dynsys NEW-UI " __DATE__);
  ImGui::SameLine();
  if (ImGui::SmallButton("Logistic demo")) load_logistic_bifurcation_demo(app);
  ImGui::SameLine();
  if (ImGui::SmallButton("Save PNG")) {
    std::string p = capture_screenshot_png(app);
    app.screenshot_msg = p.empty() ? "screenshot failed" : ("saved " + p);
    app.screenshot_msg_timer = 180; /* show ~3s at 60fps */
  }
  ImGui::SameLine();
  if (app.screenshot_msg_timer > 0) {
    --app.screenshot_msg_timer;
    ImGui::TextColored(ImVec4(0.6f, 0.95f, 0.7f, 1.0f), "%s", app.screenshot_msg.c_str());
    ImGui::SameLine();
  }
  ImGui::TextDisabled("(full-window plot)"); ImGui::SameLine();
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
    ImGui::SameLine();
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
  ImGui::TextUnformatted("|"); ImGui::SameLine();

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
    case AppState::ActiveView::Phase2D:     render_phase_background(app); break;
    case AppState::ActiveView::Bifurcation: render_bifurcation_background(app); break;
    case AppState::ActiveView::Fractal:     render_fractal_background(app); break;
    case AppState::ActiveView::Scene3DBridge: /* drawn to backbuffer in main loop */ break;
    case AppState::ActiveView::Basins:      render_basin_background(app); break;
    case AppState::ActiveView::ParamScan2D: render_scan_background(app); break;
    case AppState::ActiveView::Continuation: render_continuation_background(app); break;
    case AppState::ActiveView::Scene3D:     /* drawn in main loop */ break;
  }

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

  ImGui::Text("FPS: %d | mode: %s | integrator: %s | dim: %zu", app.fps, mode_name(app.mode), integrator_name(app.integrator), app.state_names.size());
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
      app.mode == SystemMode::Map && app.state_names.size() >= 2) {
    if (ImGui::CollapsingHeader("Fractal", ImGuiTreeNodeFlags_DefaultOpen)) {
      int mode = (int)app.fractal_mode;
      const char *modes[] = {"parameter space (Mandelbrot-type)", "state space (Julia-type)"};
      if (ImGui::Combo("mode", &mode, modes, 2)) {
        app.fractal_mode = (AppState::FractalMode)mode;
        app.fractal_dirty = true;
      }
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
      ImGui::TextWrapped("The Mandelbrot set (flat) with the logistic bifurcation diagram "
                         "rising along its real axis. They align because z^2+c and the "
                         "logistic map are conjugate.");
      bool changed = false;
      changed |= ImGui::Checkbox("show Mandelbrot", &app.bridge_show_mandelbrot);
      ImGui::SameLine();
      changed |= ImGui::Checkbox("show bifurcation", &app.bridge_show_bifurcation);
      changed |= ImGui::SliderInt("Mandelbrot resolution", &app.bridge_mandel_res, 120, 700);
      changed |= ImGui::SliderInt("bifurcation slices", &app.bridge_bif_slices, 100, 1500);
      changed |= ImGui::SliderInt("attractor points/slice", &app.bridge_bif_keep, 20, 400);
      changed |= ImGui::SliderFloat("bifurcation height", &app.bridge_height, 5.0f, 40.0f, "%.1f");
      if (ImGui::Button("Rebuild bridge") || changed) app.bridge_built = false;
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
      ImGui::TextWrapped("Largest Lyapunov exponent over a grid of two parameters. Periodic "
                         "windows (dark/blue) form the shrimp shapes inside the chaotic sea "
                         "(bright/warm). Best on maps (Henon, Tinkerbell).");
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

  if (ImGui::CollapsingHeader("System editor", ImGuiTreeNodeFlags_DefaultOpen)) {
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
    ImGui::InputTextMultiline("##system_input", app.system_input, sizeof(app.system_input), ImVec2(-FLT_MIN, 280.0f), ImGuiInputTextFlags_AllowTabInput);
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

  if (ImGui::CollapsingHeader("Parameters", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (app.params.empty()) ImGui::TextDisabled("No param declarations in this system.");
    for (auto &param : app.params) {
      double old = param.value;
      if (param.has_range) {
        float v = static_cast<float>(param.value);
        if (ImGui::SliderFloat(param.name.c_str(), &v, static_cast<float>(param.min_value), static_cast<float>(param.max_value), "%.6g")) {
          param.value = static_cast<double>(v);
        }
      } else {
        ImGui::InputDouble(param.name.c_str(), &param.value, 0.01, 0.1, "%.8g");
      }
      if (old != param.value) {
        app.poincare_points.clear();
        reset_lyapunov(app);
        sync_param_values(app);
      }
      ImGui::SameLine();
      if (ImGui::SmallButton((std::string("reset##") + param.name).c_str())) {
        param.value = param.default_value;
        sync_param_values(app);
      }
    }
  }

  if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
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

    ImGui::Text("Poincare points: %zu", app.poincare_points.size());
    if (ImGui::Button("Clear Poincare")) app.poincare_points.clear();
    ImGui::Separator();
    if (ImGui::Button("Find fixed point")) find_fixed_point(app);
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

    if (ImGui::Button("Continue equilibrium")) run_equilibrium_continuation(app);
    ImGui::SameLine();
    if (ImGui::Button("Clear branch")) {
      app.cont_branch = dynsys::analysis::Branch{};
      app.cont_ready = false;
    }

    if (app.cont_ready && !app.cont_branch.points.empty()) {
      ImGui::Text("branch: %zu points", app.cont_branch.points.size());
      /* Special points table. */
      if (!app.cont_branch.special_indices.empty()) {
        ImGui::Text("special points:");
        for (size_t i : app.cont_branch.special_indices) {
          const auto &bp = app.cont_branch.points[i];
          const char *kind =
              bp.special == dynsys::analysis::SpecialPointKind::Fold ? "FOLD"
              : bp.special == dynsys::analysis::SpecialPointKind::Hopf
                  ? "HOPF"
                  : "end";
          ImGui::BulletText("%s at %s = %.6g", kind, app.cont_param, bp.p);
        }
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
      draw->AddText(ImVec2(p0.x + 8, p0.y + 6), IM_COL32(200, 200, 210, 220),
                    "green: stable  red: unstable  (vertical: first state var)");
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
  if (ImGui::CollapsingHeader("Poincare section")) {
    draw_scatter_plot("Poincare section", app.poincare_points, ImVec2(-FLT_MIN, 300));
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
  GLFWwindow *window = glfwCreateWindow(app.window_width, app.window_height, "dynsys laboratory — TPCAS + Dear ImGui", nullptr, nullptr);
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
              app.mode == SystemMode::Map ? "map" : "ode",
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

int main(int argc, char **argv) {
  if (argc >= 2 && std::strcmp(argv[1], "--headless") == 0) {
    return run_headless(argc, argv);
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
  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();
    process_keyboard(app, window);
    if (!app.paused) {
      for (int i = 0; i < app.steps_per_frame; ++i) {
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
