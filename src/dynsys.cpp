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
  RK2,
  RK4,
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
  case Integrator::RK4: return "RK4";
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
  enum class ActiveView { Phase2D, Scene3D };
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

  char export_path[256] = "dynsys_export.csv";
  char bif_param[64] = "rho";
  char bif_observable[128] = "x";
  double bif_start = 0.0;
  double bif_end = 100.0;
  int bif_slices = 160;
  int bif_discard = 1000;
  int bif_keep = 80;
  std::vector<Point2> bifurcation_points;

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

bool step_ode_state(AppState &app, const State &in, State *out, char *err,
                    size_t err_cap) {
  const size_t dim = app.state_names.size();
  if (app.integrator == Integrator::Euler) {
    if (!eval_rhs(app, in, &app.scratch_k1, err, err_cap)) return false;
    add_scaled_into(in, app.scratch_k1, app.dt, dim, out);
    return true;
  }
  if (app.integrator == Integrator::RK2) {
    if (!eval_rhs(app, in, &app.scratch_k1, err, err_cap)) return false;
    add_scaled_into(in, app.scratch_k1, app.dt * 0.5, dim, &app.scratch_mid);
    if (!eval_rhs(app, app.scratch_mid, &app.scratch_k2, err, err_cap)) return false;
    add_scaled_into(in, app.scratch_k2, app.dt, dim, out);
    return true;
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
      else if (rhs == "rk4") app.integrator = Integrator::RK4;
      else {
        *error = "line " + std::to_string(line_no + 1) + ": integrator must be euler, rk2, or rk4";
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
    {"Henon map", "Discrete maps", "# Hénon map\nmode = map\nparam a = 1.4 [0,2]\nparam b = 0.3 [0,1]\nx_next = 1 - a * x * x + y\ny_next = b * x\nz_next = z\nobserve radius = sqrt(x*x + y*y)\n", {0.1,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Logistic map", "Discrete maps", "# Logistic map embedded in 3D\nmode = map\nparam r = 3.9 [0,4]\nx_next = r * x * (1 - x)\ny_next = x\nz_next = z\nobserve x_value = x\n", {0.2,0.0,0.0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Ikeda map", "Discrete maps", "# Ikeda map\nmode = map\nparam u = 0.918 [0,1]\ntheta = 0.4 - 6 / (1 + x*x + y*y)\nx_next = 1 + u * (x * cos(theta) - y * sin(theta))\ny_next = u * (x * sin(theta) + y * cos(theta))\nz_next = z\n", {0.1,0.1,0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
    {"Tinkerbell map", "Discrete maps", "# Tinkerbell map\nmode = map\nparam a = 0.9 [0,2]\nparam b = 0 - 0.6013 [-2,2]\nparam c = 2.0 [0,4]\nparam d = 0.5 [0,2]\nx_next = x*x - y*y + a*x + b*y\ny_next = 2*x*y + c*x + d*y\nz_next = z\n", {0.1,0.1,0,0}, 1.0, 1, Integrator::RK4,1.0f,0,0,0},
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
    /* backward */
    app.dt = -std::fabs(dt0);
    State s = start;
    for (int i = 0; i < half; ++i) {
      State nx{};
      if (!step_ode_state(app, s, &nx, err, sizeof(err))) break;
      bool finite = true;
      for (double v : nx.v) finite = finite && std::isfinite(v);
      if (!finite) break;
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
      if (!finite) break;
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
    b.xmin = DBL_MAX; b.xmax = -DBL_MAX;
    b.ymin = DBL_MAX; b.ymax = -DBL_MAX;
    auto include = [&](const State &s) {
      const double x = state_at(s, ix);
      const double y = state_at(s, iy);
      if (std::isfinite(x) && std::isfinite(y)) {
        b.xmin = std::min(b.xmin, x); b.xmax = std::max(b.xmax, x);
        b.ymin = std::min(b.ymin, y); b.ymax = std::max(b.ymax, y);
      }
    };
    for (const State &s : app.history) include(s);
    include(app.current);
    for (const auto &traj : app.phase_trajectories) {
      if (!traj.visible) continue;
      for (const State &s : traj.history) include(s);
    }
    for (const auto &curve : app.separatrix_curves) {
      if (!curve.visible) continue;
      for (const State &s : curve.history) include(s);
    }
    if (app.fixed_ready) include(app.fixed_point);
    if (b.xmin == DBL_MAX) {
      const double x = state_at(app.current, ix);
      const double y = state_at(app.current, iy);
      b.xmin = x - 10.0; b.xmax = x + 10.0;
      b.ymin = y - 10.0; b.ymax = y + 10.0;
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
        const double max_len_data = 0.34 * std::min(bx, by) / grid;
        double nvx = vx, nvy = vy;
        if (app.phase_normalize_vectors) { nvx = vx / len * max_len_data; nvy = vy / len * max_len_data; }
        else { const double s = max_len_data / (1.0 + len); nvx *= s; nvy *= s; }
        /* speed-colored arrows: slow=blue, fast=warm */
        const double tcol = std::min(1.0, len / (4.0 * (max_len_data / std::max(1e-9, max_len_data))) * 0.0 + std::min(1.0, len / 5.0));
        const int rr = static_cast<int>(90 + 150 * tcol);
        const int gg = static_cast<int>(140 - 40 * tcol);
        const int bb = static_cast<int>(190 - 120 * tcol);
        const Point2 a = plot_to_screen(bounds, p0, w, h, x - 0.5 * nvx, y - 0.5 * nvy);
        const Point2 b = plot_to_screen(bounds, p0, w, h, x + 0.5 * nvx, y + 0.5 * nvy);
        draw_arrow(draw, ImVec2((float)a.x, (float)a.y), ImVec2((float)b.x, (float)b.y),
                   IM_COL32(rr, gg, bb, 150));
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

[[maybe_unused]] void run_bifurcation(AppState &app) {
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
  const double old_value = param->value;
  char err[256] = {0};
  const int slices = std::max(2, app.bif_slices);
  for (int i = 0; i < slices; ++i) {
    const double u = static_cast<double>(i) / static_cast<double>(slices - 1);
    param->value = app.bif_start + u * (app.bif_end - app.bif_start);
    State s = app.start;
    for (int j = 0; j < app.bif_discard; ++j) {
      State next{};
      if (!step_state(app, s, &next, err, sizeof(err))) { app.analysis_message = err; param->value = old_value; arena_destroy(&tmp_arena); return; }
      s = next;
    }
    for (int j = 0; j < app.bif_keep; ++j) {
      State next{};
      if (!step_state(app, s, &next, err, sizeof(err))) { app.analysis_message = err; param->value = old_value; arena_destroy(&tmp_arena); return; }
      s = next;
      double val = 0.0;
      if (eval_expr_at(app, obs, s, &val, err, sizeof(err))) {
        app.bifurcation_points.push_back(Point2{param->value, val});
      }
    }
  }
  param->value = old_value;
  arena_destroy(&tmp_arena);
  app.analysis_message = "bifurcation scan completed";
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
  ImGui::TextColored(ImVec4(0.45f, 0.95f, 0.55f, 1.0f), "dynsys NEW-UI");
  ImGui::SameLine();
  ImGui::TextDisabled("(full-window plot)"); ImGui::SameLine();
  ImGui::TextUnformatted("|"); ImGui::SameLine();

  if (ImGui::Button(app.show_side_panel ? "<< panel" : "panel >>"))
    app.show_side_panel = !app.show_side_panel;
  ImGui::SameLine();
  ImGui::TextUnformatted("|"); ImGui::SameLine();

  /* View switch — 3D only offered when the system is 3D. */
  const bool can3d = system_is_3d(app);
  if (ImGui::RadioButton("2D phase", app.active_view == AppState::ActiveView::Phase2D))
    app.active_view = AppState::ActiveView::Phase2D;
  ImGui::SameLine();
  if (can3d) {
    if (ImGui::RadioButton("3D", app.active_view == AppState::ActiveView::Scene3D))
      app.active_view = AppState::ActiveView::Scene3D;
  } else {
    ImGui::BeginDisabled();
    bool dummy = false;
    ImGui::RadioButton("3D", dummy);
    ImGui::EndDisabled();
    app.active_view = AppState::ActiveView::Phase2D;
  }
  ImGui::SameLine(); ImGui::TextUnformatted("|"); ImGui::SameLine();

  if (ImGui::Button(app.paused ? "Play" : "Pause")) app.paused = !app.paused;
  ImGui::SameLine();
  if (ImGui::Button("Reset")) reset_simulation(app);
  ImGui::SameLine();
  if (app.active_view == AppState::ActiveView::Phase2D) {
    if (ImGui::Button("Clear orbits")) { app.phase_trajectories.clear(); app.separatrix_curves.clear(); }
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
  /* PHASE6-UI: render the active background plot first (behind panels). */
  if (!(system_is_3d(app) && app.active_view == AppState::ActiveView::Scene3D))
    render_phase_background(app);

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
    const char *integrators[] = {"Euler", "RK2 midpoint", "RK4"};
    if (ImGui::Combo("integrator", &integrator_idx, integrators, 3)) app.integrator = static_cast<Integrator>(integrator_idx);
    ImGui::InputDouble("dt", &app.dt, 0.001, 0.01, "%.8f");
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

  if (ImGui::CollapsingHeader("Analysis", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("estimate largest Lyapunov exponent", &app.lyapunov_enabled);
    ImGui::InputDouble("Lyapunov epsilon", &app.lyapunov_epsilon, 1e-6, 1e-5, "%.1e");
    ImGui::SameLine();
    if (ImGui::Button("Reset Lyapunov")) reset_lyapunov(app);
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
              app.integrator == Integrator::Euler ? "euler" :
              app.integrator == Integrator::RK2 ? "rk2" : "rk4",
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
    /* PHASE6-UI: the active plot is the full-window background. For 3D we
     * render geometry straight to the backbuffer now; for 2D we clear and
     * the phase plane paints itself via the ImGui background draw list
     * inside draw_gui. */
    const bool want_3d = system_is_3d(app) &&
                         app.active_view == AppState::ActiveView::Scene3D;
    if (want_3d) {
      render_scene_background(app);
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
  if (app.scene_depth) glDeleteRenderbuffers(1, &app.scene_depth);
  shutdown_imgui();
  glfwDestroyWindow(window);
  glfwTerminate();
  if (app.arena_ready) arena_destroy(&app.system_arena);
  return EXIT_SUCCESS;
}
