#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <cglm/call.h>
#include <cglm/cglm.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

extern "C" {
#include "arena.h"
#include "ast.h"
#include "pratt.h"
}

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

struct KeyState {
  bool pressed = false;
  double last_pressed_time = 0.0;
};

struct AppState {
  arena_t equation_arena{};
  bool arena_ready = false;
  node_t *equations[3] = {nullptr, nullptr, nullptr};

  char dx_input[256] = "10 * (y - x)";
  char dy_input[256] = "x * (28 - z) - y";
  char dz_input[256] = "x * y - (8 / 3) * z";
  std::string parse_error;
  std::string runtime_error;

  double current_x = 0.1;
  double current_y = 0.1;
  double current_z = 0.1;
  double start_x = 0.1;
  double start_y = 0.1;
  double start_z = 0.1;
  double dt = 0.01;

  int num_points = 10000;
  int current_point = 0;
  int steps_per_frame = 1;
  bool paused = false;

  float center_x = 0.0f;
  float center_y = 0.0f;
  float center_z = 0.0f;
  float angle_x = 0.0f;
  float angle_y = 0.0f;
  float zoom = 1.0f;
  bool show_axes = true;
  bool show_center_cross = true;

  int fps = 0;
  int window_width = 800;
  int window_height = 800;

  GLuint vbo = 0;
  GLuint vao = 0;
  GLuint cbo = 0;
  GLuint shader_program = 0;
  mat4 projection{};

  std::vector<Point> points;
  std::vector<float> colors;

  KeyState keys[GLFW_KEY_LAST + 1]{};
  bool mouse_left_pressed = false;
  double last_mouse_x = 0.0;
  double last_mouse_y = 0.0;
};

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

bool eval_ast(AppState &app, const node_t *ast, double *out, char *err,
              size_t err_cap);

bool eval_call(AppState &app, const char *name, size_t argc, node_t **args,
               double *out, char *err, size_t err_cap) {
  if (argc == 1) {
    double x = 0.0;
    if (!eval_ast(app, args[0], &x, err, err_cap)) {
      return false;
    }
    if (std::strcmp(name, "sin") == 0) {
      *out = std::sin(x);
      return true;
    }
    if (std::strcmp(name, "cos") == 0) {
      *out = std::cos(x);
      return true;
    }
    if (std::strcmp(name, "tan") == 0) {
      *out = std::tan(x);
      return true;
    }
    if (std::strcmp(name, "exp") == 0) {
      *out = std::exp(x);
      return true;
    }
    if (std::strcmp(name, "log") == 0) {
      *out = std::log(x);
      return true;
    }
    if (std::strcmp(name, "sqrt") == 0) {
      *out = std::sqrt(x);
      return true;
    }
    if (std::strcmp(name, "abs") == 0) {
      *out = std::fabs(x);
      return true;
    }
  }

  if (argc == 2) {
    double a = 0.0;
    double b = 0.0;
    if (!eval_ast(app, args[0], &a, err, err_cap) ||
        !eval_ast(app, args[1], &b, err, err_cap)) {
      return false;
    }
    if (std::strcmp(name, "pow") == 0) {
      *out = std::pow(a, b);
      return true;
    }
    if (std::strcmp(name, "min") == 0) {
      *out = std::fmin(a, b);
      return true;
    }
    if (std::strcmp(name, "max") == 0) {
      *out = std::fmax(a, b);
      return true;
    }
  }

  set_error(err, err_cap, "unsupported function call: %s/%zu", name, argc);
  return false;
}

bool eval_ast(AppState &app, const node_t *ast, double *out, char *err,
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

  case NODE_VAR:
    if (std::strcmp(ast->var.name, "x") == 0) {
      *out = app.current_x;
      return true;
    }
    if (std::strcmp(ast->var.name, "y") == 0) {
      *out = app.current_y;
      return true;
    }
    if (std::strcmp(ast->var.name, "z") == 0) {
      *out = app.current_z;
      return true;
    }
    if (std::strcmp(ast->var.name, "pi") == 0) {
      *out = M_PI;
      return true;
    }
    if (std::strcmp(ast->var.name, "e") == 0) {
      *out = M_E;
      return true;
    }
    set_error(err, err_cap, "unknown variable: %s", ast->var.name);
    return false;

  case NODE_APP:
    if (ast->app.head->kind == NODE_CONST && ast->app.head->cnst.op != nullptr) {
      const op_info_t *op = ast->app.head->cnst.op;
      if (ast->app.argc == 2) {
        double a = 0.0;
        double b = 0.0;
        if (!eval_ast(app, ast->app.args[0], &a, err, err_cap) ||
            !eval_ast(app, ast->app.args[1], &b, err, err_cap)) {
          return false;
        }
        if (op == &OP_ADD) {
          *out = a + b;
          return true;
        }
        if (op == &OP_SUB) {
          *out = a - b;
          return true;
        }
        if (op == &OP_MUL) {
          *out = a * b;
          return true;
        }
        if (op == &OP_DIV) {
          *out = a / b;
          return true;
        }
      }
      set_error(err, err_cap, "unsupported operator application: %s/%zu",
                op->name, ast->app.argc);
      return false;
    }

    if (ast->app.head->kind == NODE_VAR) {
      return eval_call(app, ast->app.head->var.name, ast->app.argc,
                       ast->app.args, out, err, err_cap);
    }

    set_error(err, err_cap, "unsupported application near byte range [%zu,%zu)",
              ast->span.start, ast->span.end);
    return false;

  case NODE_LAMBDA:
  case NODE_FORALL:
  case NODE_EXISTS:
    set_error(err, err_cap, "binders are not valid inside numeric ODE equations");
    return false;
  }

  set_error(err, err_cap, "unsupported AST node");
  return false;
}

std::string format_parse_error(const char *label, const char *input,
                               parse_result_t result) {
  char buf[768];
  std::snprintf(buf, sizeof(buf),
                "%s parse error at line %zu, column %zu: %s\n%s\n%*s^",
                label, result.err_span.line, result.err_span.col,
                result.err_msg ? result.err_msg : "unknown parse error", input,
                result.err_span.col > 0 ? static_cast<int>(result.err_span.col - 1)
                                        : 0,
                "");
  return std::string(buf);
}

bool compile_equations(AppState &app, const char *dx, const char *dy,
                       const char *dz, std::string *error) {
  arena_t next_arena{};
  arena_init(&next_arena, 65536);

  const char *inputs[3] = {dx, dy, dz};
  const char *labels[3] = {"dx/dt", "dy/dt", "dz/dt"};
  node_t *next_equations[3] = {nullptr, nullptr, nullptr};

  for (int i = 0; i < 3; ++i) {
    parse_result_t result = parse(inputs[i], &next_arena);
    if (!result.ok) {
      if (error != nullptr) {
        *error = format_parse_error(labels[i], inputs[i], result);
      }
      arena_destroy(&next_arena);
      return false;
    }
    next_equations[i] = result.ast;
  }

  if (app.arena_ready) {
    arena_destroy(&app.equation_arena);
  }
  app.equation_arena = next_arena;
  app.arena_ready = true;
  app.equations[0] = next_equations[0];
  app.equations[1] = next_equations[1];
  app.equations[2] = next_equations[2];
  app.parse_error.clear();
  app.runtime_error.clear();
  return true;
}

void set_lorenz(AppState &app) {
  std::snprintf(app.dx_input, sizeof(app.dx_input), "%s", "10 * (y - x)");
  std::snprintf(app.dy_input, sizeof(app.dy_input), "%s", "x * (28 - z) - y");
  std::snprintf(app.dz_input, sizeof(app.dz_input), "%s",
                "x * y - (8 / 3) * z");
  app.start_x = 0.1;
  app.start_y = 0.1;
  app.start_z = 0.1;
}

void set_rossler(AppState &app) {
  std::snprintf(app.dx_input, sizeof(app.dx_input), "%s", "0 - (y + z)");
  std::snprintf(app.dy_input, sizeof(app.dy_input), "%s", "x + (1 / 5) * y");
  std::snprintf(app.dz_input, sizeof(app.dz_input), "%s",
                "(1 / 5) + z * (x - (57 / 10))");
  app.start_x = 0.1;
  app.start_y = 0.1;
  app.start_z = 0.1;
}

void allocate_points(AppState &app) {
  if (app.num_points < 2) {
    app.num_points = 2;
  }
  app.points.assign(static_cast<size_t>(app.num_points), Point{0.0f, 0.0f, 0.0f});
  app.colors.assign(static_cast<size_t>(app.num_points) * 3U, 0.0f);
  app.current_point = 0;
}

void reset_simulation(AppState &app) {
  app.current_point = 0;
  app.current_x = app.start_x;
  app.current_y = app.start_y;
  app.current_z = app.start_z;
  for (int i = 0; i < app.num_points; ++i) {
    app.points[static_cast<size_t>(i)] = Point{0.0f, 0.0f, 0.0f};
    app.colors[static_cast<size_t>(i) * 3U + 0U] = 1.0f;
    app.colors[static_cast<size_t>(i) * 3U + 1U] = 0.5f;
    app.colors[static_cast<size_t>(i) * 3U + 2U] = 0.0f;
  }
  app.runtime_error.clear();
}

void upload_buffers(AppState &app) {
  glBindVertexArray(app.vao);

  glBindBuffer(GL_ARRAY_BUFFER, app.vbo);
  glBufferData(GL_ARRAY_BUFFER,
               static_cast<GLsizeiptr>(app.points.size() * sizeof(Point)),
               app.points.data(), GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Point),
                        reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, app.cbo);
  glBufferData(GL_ARRAY_BUFFER,
               static_cast<GLsizeiptr>(app.colors.size() * sizeof(float)),
               app.colors.data(), GL_DYNAMIC_DRAW);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float),
                        reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(1);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

bool calculate_next_point(AppState &app) {
  if (app.equations[0] == nullptr || app.equations[1] == nullptr ||
      app.equations[2] == nullptr) {
    app.runtime_error = "equations are not compiled";
    app.paused = true;
    return false;
  }

  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  char err[256] = {0};
  if (!eval_ast(app, app.equations[0], &dx, err, sizeof(err)) ||
      !eval_ast(app, app.equations[1], &dy, err, sizeof(err)) ||
      !eval_ast(app, app.equations[2], &dz, err, sizeof(err))) {
    app.runtime_error = err;
    app.paused = true;
    return false;
  }

  app.current_x += dx * app.dt;
  app.current_y += dy * app.dt;
  app.current_z += dz * app.dt;

  const auto idx = static_cast<size_t>(app.current_point);
  app.points[idx] = Point{static_cast<float>(app.current_x),
                          static_cast<float>(app.current_y),
                          static_cast<float>(app.current_z)};

  app.colors[idx * 3U + 0U] = static_cast<float>(app.current_point % 256) / 256.0f;
  app.colors[idx * 3U + 1U] = 0.5f;
  app.colors[idx * 3U + 2U] = static_cast<float>(255 - (app.current_point % 256)) / 256.0f;

  app.current_point = (app.current_point + 1) % app.num_points;
  return true;
}

bool key_pressed(AppState &app, GLFWwindow *window, int key, double debounce_time) {
  if (key < 0 || key > GLFW_KEY_LAST) {
    return false;
  }
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
  if (io.WantCaptureKeyboard) {
    return;
  }

  if (key_pressed(app, window, GLFW_KEY_W, 0.015)) {
    app.angle_x += 5.0f;
  }
  if (key_pressed(app, window, GLFW_KEY_S, 0.015)) {
    app.angle_x -= 5.0f;
  }
  if (key_pressed(app, window, GLFW_KEY_A, 0.015)) {
    app.angle_y -= 5.0f;
  }
  if (key_pressed(app, window, GLFW_KEY_D, 0.015)) {
    app.angle_y += 5.0f;
  }
  if (key_pressed(app, window, GLFW_KEY_EQUAL, 0.015)) {
    app.zoom *= 0.9f;
  }
  if (key_pressed(app, window, GLFW_KEY_MINUS, 0.015)) {
    app.zoom /= 0.9f;
  }
  if (key_pressed(app, window, GLFW_KEY_SPACE, 0.1)) {
    app.paused = !app.paused;
  }
  if (key_pressed(app, window, GLFW_KEY_X, 0.1)) {
    app.show_axes = !app.show_axes;
  }
  if (key_pressed(app, window, GLFW_KEY_Z, 0.1)) {
    app.show_center_cross = !app.show_center_cross;
  }
  if (key_pressed(app, window, GLFW_KEY_PERIOD, 0.05)) {
    app.steps_per_frame += 1;
  }
  if (key_pressed(app, window, GLFW_KEY_SLASH, 0.05)) {
    app.steps_per_frame *= 2;
  }
  if (key_pressed(app, window, GLFW_KEY_COMMA, 0.05)) {
    app.steps_per_frame = std::max(1, app.steps_per_frame - 1);
  }
  if (key_pressed(app, window, GLFW_KEY_M, 0.05)) {
    app.steps_per_frame = std::max(1, app.steps_per_frame / 2);
  }
  if (key_pressed(app, window, GLFW_KEY_V, 0.05)) {
    app.dt *= 0.5;
  }
  if (key_pressed(app, window, GLFW_KEY_B, 0.05)) {
    app.dt *= 2.0;
  }
  if (key_pressed(app, window, GLFW_KEY_C, 0.015)) {
    reset_simulation(app);
  }
  if (key_pressed(app, window, GLFW_KEY_J, 0.015)) {
    app.center_x += 0.5f;
  }
  if (key_pressed(app, window, GLFW_KEY_L, 0.015)) {
    app.center_x -= 0.5f;
  }
  if (key_pressed(app, window, GLFW_KEY_I, 0.015)) {
    app.center_y += 0.5f;
  }
  if (key_pressed(app, window, GLFW_KEY_K, 0.015)) {
    app.center_y -= 0.5f;
  }
  if (key_pressed(app, window, GLFW_KEY_P, 0.015)) {
    app.center_z += 0.5f;
  }
  if (key_pressed(app, window, GLFW_KEY_SEMICOLON, 0.015)) {
    app.center_z -= 0.5f;
  }
  if (key_pressed(app, window, GLFW_KEY_ESCAPE, 0.01)) {
    glfwSetWindowShouldClose(window, GLFW_TRUE);
  }
}

void check_compile_errors(GLuint shader, const char *type) {
  GLint success = 0;
  GLchar info_log[1024] = {0};
  if (std::strcmp(type, "PROGRAM") == 0) {
    glGetProgramiv(shader, GL_LINK_STATUS, &success);
    if (!success) {
      glGetProgramInfoLog(shader, 1024, nullptr, info_log);
      std::fprintf(stderr, "ERROR::PROGRAM_LINKING_ERROR of type: %s\n%s\n", type,
                   info_log);
    }
  } else {
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
      glGetShaderInfoLog(shader, 1024, nullptr, info_log);
      std::fprintf(stderr, "ERROR::SHADER_COMPILATION_ERROR of type: %s\n%s\n", type,
                   info_log);
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
  glm_perspective(glm_rad(45.0f), aspect, 1.0f, 500.0f, app.projection);
}

void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
  auto *app = static_cast<AppState *>(glfwGetWindowUserPointer(window));
  if (app == nullptr) {
    return;
  }
  app->window_width = width;
  app->window_height = height;
  glViewport(0, 0, width, height);
  update_projection(*app);
}

void mouse_button_callback(GLFWwindow *window, int button, int action, int mods) {
  (void)mods;
  auto *app = static_cast<AppState *>(glfwGetWindowUserPointer(window));
  if (app == nullptr || ImGui::GetIO().WantCaptureMouse) {
    return;
  }
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
  if (app == nullptr || ImGui::GetIO().WantCaptureMouse) {
    return;
  }
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
  if (app == nullptr || ImGui::GetIO().WantCaptureMouse) {
    return;
  }
  app->zoom *= static_cast<float>(1.0 - yoffset * 0.1);
  app->zoom = std::max(0.05f, app->zoom);
}

void draw_attractor(AppState &app) {
  glBindVertexArray(app.vao);
  if (app.current_point + 1 < app.num_points) {
    glDrawArrays(GL_LINE_STRIP, app.current_point + 1,
                 app.num_points - app.current_point - 1);
  }
  if (app.current_point > 1) {
    glDrawArrays(GL_LINE_STRIP, 0, app.current_point);
  }
  glBindVertexArray(0);
}

void draw_axes(AppState &app) {
  if (!app.show_axes) {
    return;
  }

  const float axes_vertices[] = {
      -20.0f, 0.0f,  0.0f, 1.0f, 0.0f, 0.0f,
       20.0f, 0.0f,  0.0f, 1.0f, 0.0f, 0.0f,
        0.0f,-20.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 20.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f,-20.0f, 0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 20.0f, 0.0f, 0.0f, 1.0f,
  };

  GLuint axes_vbo = 0;
  GLuint axes_vao = 0;
  glGenVertexArrays(1, &axes_vao);
  glGenBuffers(1, &axes_vbo);
  glBindVertexArray(axes_vao);
  glBindBuffer(GL_ARRAY_BUFFER, axes_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(axes_vertices), axes_vertices, GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float),
                        reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float),
                        reinterpret_cast<void *>(3 * sizeof(float)));
  glEnableVertexAttribArray(1);
  glDrawArrays(GL_LINES, 0, 6);
  glDeleteVertexArrays(1, &axes_vao);
  glDeleteBuffers(1, &axes_vbo);
}

void draw_scene(AppState &app) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  glUseProgram(app.shader_program);

  mat4 model;
  mat4 view;
  glm_mat4_identity(model);
  glm_mat4_identity(view);

  vec3 translate{};
  set_vec3(translate, app.center_x, app.center_y, app.center_z);
  glm_translate(model, translate);

  vec3 scale{};
  set_vec3(scale, 0.05f, 0.05f, 0.05f);
  glm_scale(model, scale);

  vec3 axis_x{};
  set_vec3(axis_x, 1.0f, 0.0f, 0.0f);
  glm_rotate(model, glm_rad(app.angle_x), axis_x);

  vec3 axis_y{};
  set_vec3(axis_y, 0.0f, 1.0f, 0.0f);
  glm_rotate(model, glm_rad(app.angle_y), axis_y);

  vec3 view_translation{};
  set_vec3(view_translation, 0.0f, 0.0f, -50.0f * app.zoom);
  glm_translate(view, view_translation);

  const GLint model_loc = glGetUniformLocation(app.shader_program, "model");
  const GLint view_loc = glGetUniformLocation(app.shader_program, "view");
  const GLint proj_loc = glGetUniformLocation(app.shader_program, "projection");
  glUniformMatrix4fv(model_loc, 1, GL_FALSE, reinterpret_cast<float *>(model));
  glUniformMatrix4fv(view_loc, 1, GL_FALSE, reinterpret_cast<float *>(view));
  glUniformMatrix4fv(proj_loc, 1, GL_FALSE, reinterpret_cast<float *>(app.projection));

  glBindBuffer(GL_ARRAY_BUFFER, app.vbo);
  glBufferSubData(GL_ARRAY_BUFFER, 0,
                  static_cast<GLsizeiptr>(app.points.size() * sizeof(Point)),
                  app.points.data());
  glBindBuffer(GL_ARRAY_BUFFER, app.cbo);
  glBufferSubData(GL_ARRAY_BUFFER, 0,
                  static_cast<GLsizeiptr>(app.colors.size() * sizeof(float)),
                  app.colors.data());

  draw_attractor(app);
  draw_axes(app);
}

void draw_center_cross(const AppState &app) {
  if (!app.show_center_cross) {
    return;
  }
  ImDrawList *draw_list = ImGui::GetForegroundDrawList();
  const ImVec2 center(static_cast<float>(app.window_width) * 0.5f,
                      static_cast<float>(app.window_height) * 0.5f);
  const ImU32 color = IM_COL32(180, 180, 180, 180);
  draw_list->AddLine(ImVec2(center.x - 8.0f, center.y),
                     ImVec2(center.x + 8.0f, center.y), color, 1.0f);
  draw_list->AddLine(ImVec2(center.x, center.y - 8.0f),
                     ImVec2(center.x, center.y + 8.0f), color, 1.0f);
}

void draw_gui(AppState &app) {
  ImGui::SetNextWindowPos(ImVec2(12, 12), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(430, 620), ImGuiCond_FirstUseEver);
  ImGui::Begin("dynsys control panel");

  ImGui::Text("FPS: %d", app.fps);
  ImGui::Text("Current: x %.5f, y %.5f, z %.5f", app.current_x, app.current_y,
              app.current_z);
  ImGui::Separator();

  if (ImGui::Button(app.paused ? "Resume" : "Pause")) {
    app.paused = !app.paused;
  }
  ImGui::SameLine();
  if (ImGui::Button("Reset simulation")) {
    reset_simulation(app);
  }
  ImGui::SameLine();
  if (ImGui::Button("Quit")) {
    GLFWwindow *window = glfwGetCurrentContext();
    if (window != nullptr) {
      glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
  }

  if (ImGui::CollapsingHeader("Equations", ImGuiTreeNodeFlags_DefaultOpen)) {
    if (ImGui::Button("Lorenz")) {
      set_lorenz(app);
      if (compile_equations(app, app.dx_input, app.dy_input, app.dz_input,
                            &app.parse_error)) {
        reset_simulation(app);
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Rössler")) {
      set_rossler(app);
      if (compile_equations(app, app.dx_input, app.dy_input, app.dz_input,
                            &app.parse_error)) {
        reset_simulation(app);
      }
    }

    ImGui::InputText("dx/dt", app.dx_input, sizeof(app.dx_input));
    ImGui::InputText("dy/dt", app.dy_input, sizeof(app.dy_input));
    ImGui::InputText("dz/dt", app.dz_input, sizeof(app.dz_input));
    if (ImGui::Button("Apply equations")) {
      if (compile_equations(app, app.dx_input, app.dy_input, app.dz_input,
                            &app.parse_error)) {
        reset_simulation(app);
      }
    }
    ImGui::TextWrapped("TPCAS infix mode: x, y, z, pi, e; + - * /; sin/cos/tan/exp/log/sqrt/abs/pow/min/max.");

    if (!app.parse_error.empty()) {
      ImGui::Separator();
      ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "%s",
                         app.parse_error.c_str());
    }
    if (!app.runtime_error.empty()) {
      ImGui::Separator();
      ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "Runtime: %s",
                         app.runtime_error.c_str());
    }
  }

  if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::InputDouble("dt", &app.dt, 0.001, 0.01, "%.8f");
    ImGui::SliderInt("steps/frame", &app.steps_per_frame, 1, 500);
    ImGui::InputDouble("start x", &app.start_x, 0.1, 1.0, "%.6f");
    ImGui::InputDouble("start y", &app.start_y, 0.1, 1.0, "%.6f");
    ImGui::InputDouble("start z", &app.start_z, 0.1, 1.0, "%.6f");

    int next_num_points = app.num_points;
    if (ImGui::SliderInt("point buffer", &next_num_points, 100, 200000)) {
      app.num_points = next_num_points;
      allocate_points(app);
      reset_simulation(app);
      upload_buffers(app);
    }
  }

  if (ImGui::CollapsingHeader("View", ImGuiTreeNodeFlags_DefaultOpen)) {
    ImGui::Checkbox("axes", &app.show_axes);
    ImGui::SameLine();
    ImGui::Checkbox("center cross", &app.show_center_cross);
    ImGui::SliderFloat("angle x", &app.angle_x, -360.0f, 360.0f);
    ImGui::SliderFloat("angle y", &app.angle_y, -360.0f, 360.0f);
    ImGui::SliderFloat("zoom", &app.zoom, 0.05f, 10.0f);
    ImGui::InputFloat("center x", &app.center_x, 0.1f, 1.0f, "%.3f");
    ImGui::InputFloat("center y", &app.center_y, 0.1f, 1.0f, "%.3f");
    ImGui::InputFloat("center z", &app.center_z, 0.1f, 1.0f, "%.3f");
  }

  if (ImGui::CollapsingHeader("Keyboard / mouse help")) {
    ImGui::TextWrapped("Mouse drag rotates, mouse wheel zooms. Keyboard: WASD rotate, +/- zoom, space pause, C reset, X axes, Z cross, comma/period speed, M/slash halve/double speed, V/B halve/double dt, IJKL/P; translate, Esc quit.");
  }

  ImGui::End();

  draw_center_cross(app);
}

bool init_glfw_window(AppState &app, GLFWwindow **out_window) {
  if (!glfwInit()) {
    std::fprintf(stderr, "Failed to initialize GLFW\n");
    return false;
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow *window = glfwCreateWindow(app.window_width, app.window_height,
                                        "dynsys — TPCAS + Dear ImGui", nullptr,
                                        nullptr);
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
  glGetError(); // GLEW may emit a harmless GL_INVALID_ENUM during init.

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

int main(int argc, char **argv) {
  (void)argc;
  (void)argv;

  AppState app{};
  set_lorenz(app);
  if (!compile_equations(app, app.dx_input, app.dy_input, app.dz_input,
                         &app.parse_error)) {
    std::fprintf(stderr, "%s\n", app.parse_error.c_str());
    return EXIT_FAILURE;
  }

  allocate_points(app);
  reset_simulation(app);

  GLFWwindow *window = nullptr;
  if (!init_glfw_window(app, &window)) {
    if (app.arena_ready) {
      arena_destroy(&app.equation_arena);
    }
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
        if (!calculate_next_point(app)) {
          break;
        }
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

    draw_scene(app);
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

  shutdown_imgui();
  glfwDestroyWindow(window);
  glfwTerminate();
  if (app.arena_ready) {
    arena_destroy(&app.equation_arena);
  }

  return EXIT_SUCCESS;
}
