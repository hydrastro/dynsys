#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <cglm/call.h>
#include <cglm/cglm.h>
#include <errno.h>
#include <ft2build.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include FT_FREETYPE_H

#include "arena.h"
#include "ast.h"
#include "pratt.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#ifndef DYNSYS_DEFAULT_FONT_PATH
#define DYNSYS_DEFAULT_FONT_PATH "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"
#endif

typedef struct {
  bool pressed;
  uint64_t last_pressed_time;
} key_state_t;

key_state_t key_states[GLFW_KEY_LAST];

uint64_t get_time_millis() { return (uint64_t)(glfwGetTime() * 1000); }

bool is_key_pressed(GLFWwindow *window, int key, uint64_t debounce_time) {
  uint64_t current_time = get_time_millis();
  if (glfwGetKey(window, key) == GLFW_PRESS) {
    if (!key_states[key].pressed ||
        (current_time - key_states[key].last_pressed_time) > debounce_time) {
      key_states[key].pressed = true;
      key_states[key].last_pressed_time = current_time;
      return true;
    }
  } else {
    key_states[key].pressed = false;
  }
  return false;
}

const char *vertex_shader_src =
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
    "}\0";

const char *fragment_shader_src = "#version 330 core\n"
                                  "in vec3 our_color;\n"
                                  "out vec4 frag_color;\n"
                                  "void main()\n"
                                  "{\n"
                                  "    frag_color = vec4(our_color, 1.0);\n"
                                  "}\n\0";

double sigma = 10.0, beta = 8.0 / 3.0, rho = 28.0;
double dt = 0.01;

float center_x = 0.0f, center_y = 0.0f, center_z = 0.0f;

float angle_x = 0.0f, angle_y = 0.0f;
float zoom = 1.0f;
int paused = 0;
bool show_axes = 1;
int animation_speed_ms = 1;
int render_speed = 1;
bool center_cross = true;

double current_x = 0.1, current_y = 0.0, current_z = 0.0;
double start_x = 0.1, start_y = 0.1, start_z = 0.1;

int num_points = 10000;
int current_point = 0;
typedef struct point {
  double x, y, z;
} point_t;
point_t *lorenz_points;
float *colors;

clock_t last_time = 0;
int fps = 0;

node_t *equations[3];
arena_t equation_arena;
int window_width = 400;
int window_height = 400;

double evaluate_ast(node_t *ast);

static bool parse_number_literal(const char *text, double *out) {
  if (text == NULL || *text == '\0')
    return false;
  errno = 0;
  char *end = NULL;
  double value = strtod(text, &end);
  if (errno != 0 || end == text || *end != '\0')
    return false;
  *out = value;
  return true;
}

static void die_bad_expression(const char *message, const node_t *ast) {
  fprintf(stderr, "%s", message);
  if (ast != NULL) {
    fprintf(stderr, " near byte range [%zu, %zu)", ast->span.start,
            ast->span.end);
  }
  fputc('\n', stderr);
  exit(EXIT_FAILURE);
}

static double eval_call(const char *name, size_t argc, node_t **args) {
  if (argc == 1) {
    double x = evaluate_ast(args[0]);
    if (strcmp(name, "sin") == 0)
      return sin(x);
    if (strcmp(name, "cos") == 0)
      return cos(x);
    if (strcmp(name, "tan") == 0)
      return tan(x);
    if (strcmp(name, "exp") == 0)
      return exp(x);
    if (strcmp(name, "log") == 0)
      return log(x);
    if (strcmp(name, "sqrt") == 0)
      return sqrt(x);
    if (strcmp(name, "abs") == 0)
      return fabs(x);
  }
  if (argc == 2) {
    double a = evaluate_ast(args[0]);
    double b = evaluate_ast(args[1]);
    if (strcmp(name, "pow") == 0)
      return pow(a, b);
    if (strcmp(name, "min") == 0)
      return fmin(a, b);
    if (strcmp(name, "max") == 0)
      return fmax(a, b);
  }

  fprintf(stderr, "Unsupported function call: %s/%zu\n", name, argc);
  exit(EXIT_FAILURE);
}

double evaluate_ast(node_t *ast) {
  switch (ast->kind) {
  case NODE_CONST: {
    double value = 0.0;
    if (parse_number_literal(ast->cnst.name, &value))
      return value;
    fprintf(stderr, "Unsupported constant in numeric equation: %s\n",
            ast->cnst.name);
    exit(EXIT_FAILURE);
  }
  case NODE_VAR: {
    if (strcmp(ast->var.name, "x") == 0)
      return current_x;
    if (strcmp(ast->var.name, "y") == 0)
      return current_y;
    if (strcmp(ast->var.name, "z") == 0)
      return current_z;
    if (strcmp(ast->var.name, "pi") == 0)
      return M_PI;
    if (strcmp(ast->var.name, "e") == 0)
      return M_E;
    fprintf(stderr, "Unknown variable: %s\n", ast->var.name);
    exit(EXIT_FAILURE);
  }
  case NODE_APP: {
    if (ast->app.head->kind == NODE_CONST && ast->app.head->cnst.op != NULL) {
      const op_info_t *op = ast->app.head->cnst.op;
      if (op == &OP_ADD && ast->app.argc == 2)
        return evaluate_ast(ast->app.args[0]) + evaluate_ast(ast->app.args[1]);
      if (op == &OP_SUB && ast->app.argc == 2)
        return evaluate_ast(ast->app.args[0]) - evaluate_ast(ast->app.args[1]);
      if (op == &OP_MUL && ast->app.argc == 2)
        return evaluate_ast(ast->app.args[0]) * evaluate_ast(ast->app.args[1]);
      if (op == &OP_DIV && ast->app.argc == 2)
        return evaluate_ast(ast->app.args[0]) / evaluate_ast(ast->app.args[1]);
    }

    if (ast->app.head->kind == NODE_VAR) {
      return eval_call(ast->app.head->var.name, ast->app.argc, ast->app.args);
    }

    die_bad_expression("Unsupported application in numeric equation", ast);
  }
  case NODE_LAMBDA:
  case NODE_FORALL:
  case NODE_EXISTS:
    die_bad_expression("Binders are not valid inside numeric ODE equations", ast);
  }

  die_bad_expression("Unsupported AST node in numeric equation", ast);
  return 0.0;
}

void reset_simulation() {
  current_point = 0;
  current_x = start_x;
  current_y = start_y;
  current_z = start_z;

  for (int i = 0; i < num_points; ++i) {
    lorenz_points[i].x = 0.0;
    lorenz_points[i].y = 0.0;
    lorenz_points[i].z = 0.0;
  }
}

void calculate_next_point() {
  double dx, dy, dz;
  dx = evaluate_ast(equations[0]);
  dy = evaluate_ast(equations[1]);
  dz = evaluate_ast(equations[2]);

  current_x += dx * dt;
  current_y += dy * dt;
  current_z += dz * dt;

  lorenz_points[current_point].x = current_x;
  lorenz_points[current_point].y = current_y;
  lorenz_points[current_point].z = current_z;

  colors[current_point * 3] = (current_point % 256) / 256.0f;
  colors[current_point * 3 + 1] = 0.5f;
  colors[current_point * 3 + 2] = (255 - (current_point % 256)) / 256.0f;

  current_point = (current_point + 1) % num_points;
}

GLuint vbo, vao, cbo, shader_program;

void draw_attractor() {

  glBindVertexArray(vao);

  glDrawArrays(GL_LINE_STRIP, current_point + 1,
               num_points - current_point - 1);
  glDrawArrays(GL_LINE_STRIP, 0, current_point);

  glBindVertexArray(0);
}

void draw_axes() {
  if (!show_axes)
    return;

  float axes_vertices[] = {
      -20.0f, 0.0f,  0.0f, 1.0f, 0.0f,   0.0f,  20.0f, 0.0f, 0.0f,
      1.0f,   0.0f,  0.0f, 0.0f, -20.0f, 0.0f,  0.0f,  1.0f, 0.0f,
      0.0f,   20.0f, 0.0f, 0.0f, 1.0f,   0.0f,  0.0f,  0.0f, -20.0f,
      0.0f,   0.0f,  1.0f, 0.0f, 0.0f,   20.0f, 0.0f,  0.0f, 1.0f};

  GLuint axes_vbo, axes_vao;
  glGenVertexArrays(1, &axes_vao);
  glGenBuffers(1, &axes_vbo);
  glBindVertexArray(axes_vao);
  glBindBuffer(GL_ARRAY_BUFFER, axes_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(axes_vertices), axes_vertices,
               GL_STATIC_DRAW);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float),
                        (void *)(3 * sizeof(float)));
  glEnableVertexAttribArray(1);

  glBindVertexArray(axes_vao);
  glDrawArrays(GL_LINES, 0, 6);

  glDeleteVertexArrays(1, &axes_vao);
  glDeleteBuffers(1, &axes_vbo);
}

void draw_center_cross() {
  if (!center_cross)
    return;

  float cross_vertices[] = {0.48f, 0.5f,  0.0f, 0.52f, 0.5f,  0.0f,
                            0.5f,  0.48f, 0.0f, 0.5f,  0.52f, 0.0f};

  float cross_color[] = {0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
                         0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f};

  GLuint cross_vao, cross_vbo, cross_cbo;
  glGenVertexArrays(1, &cross_vao);
  glGenBuffers(1, &cross_vbo);
  glGenBuffers(1, &cross_cbo);

  glBindVertexArray(cross_vao);

  glBindBuffer(GL_ARRAY_BUFFER, cross_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(cross_vertices), cross_vertices,
               GL_STATIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, cross_cbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(cross_color), cross_color,
               GL_STATIC_DRAW);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(1);

  mat4 model, view;
  glm_mat4_identity(model);
  glm_mat4_identity(view);

  mat4 projection;
  glm_ortho(0.0f, 1.0f, 0.0f, 1.0f, -1.0f, 1.0f, projection);

  glUseProgram(shader_program);
  GLint model_loc = glGetUniformLocation(shader_program, "model");
  GLint view_loc = glGetUniformLocation(shader_program, "view");
  GLint proj_loc = glGetUniformLocation(shader_program, "projection");
  glUniformMatrix4fv(model_loc, 1, GL_FALSE, (float *)model);
  glUniformMatrix4fv(view_loc, 1, GL_FALSE, (float *)view);
  glUniformMatrix4fv(proj_loc, 1, GL_FALSE, (float *)projection);

  glBindVertexArray(cross_vao);
  glDrawArrays(GL_LINES, 0, 4);

  glDeleteVertexArrays(1, &cross_vao);
  glDeleteBuffers(1, &cross_vbo);
  glDeleteBuffers(1, &cross_cbo);
}

mat4 projection;

typedef struct {
  GLuint texture_id;
  int size[2];
  int bearing[2];
  GLuint advance;
} character_t;

character_t characters[128];
GLuint text_vao, text_vbo, text_shader_program;

void load_font(const char *font_path) {
  FT_Library ft;
  if (FT_Init_FreeType(&ft)) {
    fprintf(stderr, "Could not init FreeType Library\n");
    exit(EXIT_FAILURE);
  }

  FT_Face face;
  if (FT_New_Face(ft, font_path, 0, &face)) {
    fprintf(stderr, "Failed to load font\n");
    exit(EXIT_FAILURE);
  }

  FT_Set_Pixel_Sizes(face, 0, 22);

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  for (unsigned char c = 0; c < 128; c++) {
    if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
      fprintf(stderr, "Failed to load Glyph\n");
      continue;
    }

    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, face->glyph->bitmap.width,
                 face->glyph->bitmap.rows, 0, GL_RED, GL_UNSIGNED_BYTE,
                 face->glyph->bitmap.buffer);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    character_t character = {
        texture,
        {face->glyph->bitmap.width, face->glyph->bitmap.rows},
        {face->glyph->bitmap_left, face->glyph->bitmap_top},
        (GLuint)face->glyph->advance.x};
    characters[c] = character;
  }

  FT_Done_Face(face);
  FT_Done_FreeType(ft);
}

void init_text_opengl() {
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glGenVertexArrays(1, &text_vao);
  glGenBuffers(1, &text_vbo);

  glBindVertexArray(text_vao);
  glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  const char *vertex_shader_source =
      "#version 330 core\n"
      "layout(location = 0) in vec4 vertex;\n"
      "out vec2 tex_coords;\n"
      "uniform mat4 projection;\n"
      "void main() {\n"
      "    gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);\n"
      "    tex_coords = vertex.zw;\n"
      "}\0";

  const char *fragment_shader_source =
      "#version 330 core\n"
      "in vec2 tex_coords;\n"
      "out vec4 color;\n"
      "uniform sampler2D text;\n"
      "uniform vec3 text_color;\n"
      "void main() {\n"
      "    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, tex_coords).r);\n"
      "    color = vec4(text_color, 1.0) * sampled;\n"
      "}\0";

  GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex_shader, 1, &vertex_shader_source, NULL);
  glCompileShader(vertex_shader);

  GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment_shader, 1, &fragment_shader_source, NULL);
  glCompileShader(fragment_shader);

  text_shader_program = glCreateProgram();
  glAttachShader(text_shader_program, vertex_shader);
  glAttachShader(text_shader_program, fragment_shader);
  glLinkProgram(text_shader_program);

  glDeleteShader(vertex_shader);
  glDeleteShader(fragment_shader);
}

void render_text(const char *text, float x, float y, float scale,
                 float color[3]) {
  glUseProgram(text_shader_program);
  glUniform3f(glGetUniformLocation(text_shader_program, "text_color"), color[0],
              color[1], color[2]);
  glActiveTexture(GL_TEXTURE0);
  glBindVertexArray(text_vao);

  mat4 projection;
  GLint proj_loc = glGetUniformLocation(text_shader_program, "projection");
  glm_ortho(0.0f, (float)window_width, 0.0f, (float)window_height, -1.0f, 1.0f,
            projection);
  glUniformMatrix4fv(proj_loc, 1, GL_FALSE, (float *)projection);

  const char *p;
  for (p = text; *p; p++) {
    character_t ch = characters[(unsigned char)*p];

    float xpos = x + ch.bearing[0] * scale;
    float ypos = y - (ch.size[1] - ch.bearing[1]) * scale;

    float w = ch.size[0] * scale;
    float h = ch.size[1] * scale;

    float vertices[6][4] = {
        {xpos, ypos + h, 0.0, 0.0}, {xpos, ypos, 0.0, 1.0},
        {xpos + w, ypos, 1.0, 1.0}, {xpos, ypos + h, 0.0, 0.0},
        {xpos + w, ypos, 1.0, 1.0}, {xpos + w, ypos + h, 1.0, 0.0}};

    glBindTexture(GL_TEXTURE_2D, ch.texture_id);
    glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glDrawArrays(GL_TRIANGLES, 0, 6);

    x += (ch.advance >> 6) * scale;
  }

  glBindVertexArray(0);
  glBindTexture(GL_TEXTURE_2D, 0);
}

void draw_hud() {
  char info1[50], info2[50], info3[50], info4[50], info5[50], info6[50],
      info7[50];
  snprintf(info1, sizeof(info1), "FPS: %d", fps);
  snprintf(info2, sizeof(info2), "Speed: %d", render_speed);
  snprintf(info3, sizeof(info3), "Points: %d", num_points);
  snprintf(info4, sizeof(info4), "Angles: %.2f %.2f", angle_x, angle_y);
  snprintf(info5, sizeof(info5), "Point: (%.2f, %.2f, %.2f)", center_x,
           center_y, center_z);
  snprintf(info6, sizeof(info6), "Zoom: %.2f", zoom);
  snprintf(info7, sizeof(info7), "dt: %f", dt);

  float color[3] = {1.0f, 1.0f, 1.0f};

  render_text(info1, 2.0f, (float)(window_height - 18 * 1), 1.0f, color);
  render_text(info2, 2.0f, (float)(window_height - 18 * 2), 1.0f, color);
  render_text(info3, 2.0f, (float)(window_height - 18 * 3), 1.0f, color);
  render_text(info4, 2.0f, (float)(window_height - 18 * 4), 1.0f, color);
  render_text(info5, 2.0f, (float)(window_height - 18 * 5), 1.0f, color);
  render_text(info6, 2.0f, (float)(window_height - 18 * 6), 1.0f, color);
  render_text(info7, 2.0f, (float)(window_height - 18 * 7), 1.0f, color);
}

void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
  (void)window;
  glViewport(0, 0, width, height);
  window_width = width;
  window_height = height;

  float aspect_ratio = (float)width / (float)height;
  glm_perspective(glm_rad(45.0f), aspect_ratio, 1.0f, 500.0f, projection);

  GLint proj_loc = glGetUniformLocation(shader_program, "projection");
  glUniformMatrix4fv(proj_loc, 1, GL_FALSE, (float *)projection);
}

void process_input(GLFWwindow *window) {
  if (is_key_pressed(window, GLFW_KEY_W, 15))
    angle_x += 5.0f;
  if (is_key_pressed(window, GLFW_KEY_S, 15))
    angle_x -= 5.0f;
  if (is_key_pressed(window, GLFW_KEY_A, 15))
    angle_y -= 5.0f;
  if (is_key_pressed(window, GLFW_KEY_D, 15))
    angle_y += 5.0f;
  if (is_key_pressed(window, GLFW_KEY_EQUAL, 15))
    zoom *= 0.9f;
  if (is_key_pressed(window, GLFW_KEY_MINUS, 15))
    zoom /= 0.9f;
  if (is_key_pressed(window, GLFW_KEY_SPACE, 100))
    paused = !paused;
  if (is_key_pressed(window, GLFW_KEY_X, 100))
    show_axes = !show_axes;
  if (is_key_pressed(window, GLFW_KEY_PERIOD, 50))
    render_speed += 1;
  if (is_key_pressed(window, GLFW_KEY_SLASH, 50)) {
    render_speed *= 2;
  }
  if (is_key_pressed(window, GLFW_KEY_COMMA, 50)) {
    render_speed -= 1;
    if (render_speed == 0)
      render_speed = 1;
  }
  if (is_key_pressed(window, GLFW_KEY_M, 50)) {
    render_speed /= 2;
    if (render_speed == 0)
      render_speed = 1;
  }
  if (is_key_pressed(window, GLFW_KEY_V, 50)) {
    dt /= 2;
    if (dt == 0)
      dt = 0.01;
  }
  if (is_key_pressed(window, GLFW_KEY_B, 50)) {
    dt *= 2;
  }
  if (is_key_pressed(window, GLFW_KEY_C, 15))
    reset_simulation();
  if (is_key_pressed(window, GLFW_KEY_J, 15))
    center_x += 0.5f;
  if (is_key_pressed(window, GLFW_KEY_L, 15))
    center_x -= 0.5f;
  if (is_key_pressed(window, GLFW_KEY_I, 15))
    center_y += 0.5f;
  if (is_key_pressed(window, GLFW_KEY_K, 15))
    center_y -= 0.5f;
  if (is_key_pressed(window, GLFW_KEY_P, 15))
    center_z += 0.5f;
  if (is_key_pressed(window, GLFW_KEY_SEMICOLON, 15))
    center_z -= 0.5f;
  if (is_key_pressed(window, GLFW_KEY_Z, 100))
    center_cross = !center_cross;
  if (is_key_pressed(window, GLFW_KEY_ESCAPE, 10)) {
    glfwSetWindowShouldClose(window, true);
  }
}

void display(GLFWwindow *window) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glUseProgram(shader_program);

  mat4 model, view;
  glm_mat4_identity(model);
  glm_translate(model, (vec3){center_x, center_y, center_z});
  glm_scale(model, (vec3){0.05f, 0.05f, 0.05f});
  glm_rotate(model, glm_rad(angle_x), (vec3){1.0f, 0.0f, 0.0f});
  glm_rotate(model, glm_rad(angle_y), (vec3){0.0f, 1.0f, 0.0f});

  glm_mat4_identity(view);
  glm_translate(view, (vec3){0.0f, 0.0f, -50.0f * zoom});

  GLint model_loc = glGetUniformLocation(shader_program, "model");
  GLint view_loc = glGetUniformLocation(shader_program, "view");
  GLint proj_loc = glGetUniformLocation(shader_program, "projection");
  glUniformMatrix4fv(model_loc, 1, GL_FALSE, (float *)model);
  glUniformMatrix4fv(view_loc, 1, GL_FALSE, (float *)view);
  glUniformMatrix4fv(proj_loc, 1, GL_FALSE, (float *)projection);

  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferSubData(GL_ARRAY_BUFFER, 0, num_points * sizeof(point_t),
                  lorenz_points);

  glBindBuffer(GL_ARRAY_BUFFER, cbo);
  glBufferSubData(GL_ARRAY_BUFFER, 0, num_points * 3 * sizeof(float), colors);

  draw_attractor();
  draw_axes();
  draw_center_cross();
  draw_hud();

  glfwSwapBuffers(window);
  glfwPollEvents();
}

void check_compile_errors(GLuint shader, const char *type) {
  GLint success;
  GLchar info_log[1024];
  if (strcmp(type, "PROGRAM") == 0) {
    glGetProgramiv(shader, GL_LINK_STATUS, &success);
    if (!success) {
      glGetProgramInfoLog(shader, 1024, NULL, info_log);
      fprintf(stderr, "ERROR::PROGRAM_LINKING_ERROR of type: %s\n%s\n", type,
              info_log);
    }
  } else {
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
      glGetShaderInfoLog(shader, 1024, NULL, info_log);
      fprintf(stderr, "ERROR::SHADER_COMPILATION_ERROR of type: %s\n%s\n", type,
              info_log);
    }
  }
}

GLuint create_shader_program(const char *vertex_src, const char *fragment_src) {
  GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex_shader, 1, &vertex_src, NULL);
  glCompileShader(vertex_shader);
  check_compile_errors(vertex_shader, "VERTEX");

  GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment_shader, 1, &fragment_src, NULL);
  glCompileShader(fragment_shader);
  check_compile_errors(fragment_shader, "FRAGMENT");

  GLuint shader_program = glCreateProgram();
  glAttachShader(shader_program, vertex_shader);
  glAttachShader(shader_program, fragment_shader);
  glLinkProgram(shader_program);
  check_compile_errors(shader_program, "PROGRAM");

  glDeleteShader(vertex_shader);
  glDeleteShader(fragment_shader);

  return shader_program;
}

bool mouse_left_pressed = false;
double last_mouse_x, last_mouse_y;

void mouse_button_callback(GLFWwindow *window, int button, int action,
                           int mods) {
  (void)mods;
  if (button == GLFW_MOUSE_BUTTON_LEFT) {
    if (action == GLFW_PRESS) {
      mouse_left_pressed = true;
      glfwGetCursorPos(window, &last_mouse_x, &last_mouse_y);
    } else if (action == GLFW_RELEASE) {
      mouse_left_pressed = false;
    }
  }
}

void cursor_position_callback(GLFWwindow *window, double xpos, double ypos) {
  (void)window;
  if (mouse_left_pressed) {
    double dx = xpos - last_mouse_x;
    double dy = ypos - last_mouse_y;
    angle_y += (float)dx * 0.1f;
    angle_x += (float)dy * 0.1f;
    last_mouse_x = xpos;
    last_mouse_y = ypos;
  }
}

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
  (void)window;
  (void)xoffset;
  zoom *= (1.0f - yoffset * 0.1f);
}


static void read_line_or_die(const char *prompt, char *buf, size_t buf_size) {
  fputs(prompt, stdout);
  fflush(stdout);
  if (fgets(buf, buf_size, stdin) == NULL) {
    fprintf(stderr, "Failed to read input.\n");
    exit(EXIT_FAILURE);
  }
  size_t n = strlen(buf);
  if (n > 0 && buf[n - 1] == '\n')
    buf[n - 1] = '\0';
}

static void print_parse_error(const char *label, const char *input,
                              parse_result_t result) {
  fprintf(stderr, "%s parse error at line %zu, column %zu: %s\n", label,
          result.err_span.line, result.err_span.col,
          result.err_msg ? result.err_msg : "unknown parse error");
  fprintf(stderr, "%s\n", input);
  for (size_t i = 1; i < result.err_span.col; ++i)
    fputc(' ', stderr);
  fputs("^\n", stderr);
}

static node_t *parse_equation_or_die(const char *label, const char *input) {
  parse_result_t result = parse(input, &equation_arena);
  if (!result.ok) {
    print_parse_error(label, input, result);
    exit(EXIT_FAILURE);
  }
  return result.ast;
}

int main(int argc, char **argv) {
  (void)argc;
  (void)argv;

  arena_init(&equation_arena, 65536);

  lorenz_points = (point_t *)malloc(num_points * sizeof(point_t));
  if (lorenz_points == NULL) {
    fprintf(stderr, "Failed to allocate memory for Lorenz points\n");
    arena_destroy(&equation_arena);
    return -1;
  }

  colors = (float *)malloc(num_points * 3 * sizeof(float));
  if (colors == NULL) {
    fprintf(stderr, "Failed to allocate memory for colors\n");
    arena_destroy(&equation_arena);
    free(lorenz_points);
    return -1;
  }

  int choice = 0;
  printf("Choose an option:\n");
  printf("1. Display Lorenz Attractor\n");
  printf("2. Display Rössler Attractor\n");
  printf("3. Enter system manually\n");
  printf("Enter your choice: ");
  if (scanf("%d", &choice) != 1) {
    fprintf(stderr, "Invalid choice. Exiting.\n");
    arena_destroy(&equation_arena);
    free(lorenz_points);
    free(colors);
    return -1;
  }

  int ch;
  while ((ch = getchar()) != '\n' && ch != EOF) {
  }

  if (choice == 1) {
    equations[0] = parse_equation_or_die("dx", "10 * (y - x)");
    equations[1] = parse_equation_or_die("dy", "x * (28 - z) - y");
    equations[2] = parse_equation_or_die("dz", "x * y - (8 / 3) * z");

    start_x = 0.1;
    start_y = 0.1;
    start_z = 0.1;
  } else if (choice == 2) {
    equations[0] = parse_equation_or_die("dx", "0 - (y + z)");
    equations[1] = parse_equation_or_die("dy", "x + (1 / 5) * y");
    equations[2] = parse_equation_or_die("dz", "(1 / 5) + z * (x - (57 / 10))");

    start_x = 0.1;
    start_y = 0.1;
    start_z = 0.1;
  } else if (choice == 3) {
    char input1[256], input2[256], input3[256];
    read_line_or_die("Enter dx/dt in TPCAS infix mode: ", input1,
                     sizeof(input1));
    read_line_or_die("Enter dy/dt in TPCAS infix mode: ", input2,
                     sizeof(input2));
    read_line_or_die("Enter dz/dt in TPCAS infix mode: ", input3,
                     sizeof(input3));

    equations[0] = parse_equation_or_die("dx", input1);
    equations[1] = parse_equation_or_die("dy", input2);
    equations[2] = parse_equation_or_die("dz", input3);

    printf("Enter the starting point (x y z): ");
    if (scanf("%lf %lf %lf", &start_x, &start_y, &start_z) != 3) {
      fprintf(stderr, "Invalid starting point. Exiting.\n");
      arena_destroy(&equation_arena);
      free(lorenz_points);
      free(colors);
      return -1;
    }
  } else {
    fprintf(stderr, "Invalid choice. Exiting.\n");
    arena_destroy(&equation_arena);
    free(lorenz_points);
    free(colors);
    return -1;
  }

  current_point = 0;
  for (int i = 0; i < num_points; ++i) {
    lorenz_points[i].x = 0.1;
    lorenz_points[i].y = 0.1;
    lorenz_points[i].z = 0.1;
    colors[i * 3] = 1.0f;
    colors[i * 3 + 1] = 0.5f;
    colors[i * 3 + 2] = 0.0f;
  }

  if (!glfwInit()) {
    fprintf(stderr, "Failed to initialize GLFW\n");
    arena_destroy(&equation_arena);
    free(lorenz_points);
    free(colors);
    return -1;
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow *window = glfwCreateWindow(
      800, 800, "Lorenz Attractor with Advanced Controls", NULL, NULL);
  if (!window) {
    fprintf(stderr, "Failed to create GLFW window\n");
    glfwTerminate();
    arena_destroy(&equation_arena);
    free(lorenz_points);
    free(colors);
    return -1;
  }

  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetScrollCallback(window, scroll_callback);

  if (glewInit() != GLEW_OK) {
    fprintf(stderr, "Failed to initialize GLEW\n");
    glfwDestroyWindow(window);
    glfwTerminate();
    free(lorenz_points);
    free(colors);
    return -1;
  }

  const char *font_path = getenv("DYNSYS_FONT_PATH");
  if (font_path == NULL || *font_path == '\0')
    font_path = DYNSYS_DEFAULT_FONT_PATH;
  load_font(font_path);

  init_text_opengl();

  shader_program =
      create_shader_program(vertex_shader_src, fragment_shader_src);

  glGenVertexArrays(1, &vao);
  glGenBuffers(1, &vbo);
  glGenBuffers(1, &cbo);

  glBindVertexArray(vao);

  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, num_points * sizeof(point_t), lorenz_points,
               GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, sizeof(point_t), (void *)0);
  glEnableVertexAttribArray(0);

  glBindBuffer(GL_ARRAY_BUFFER, cbo);
  glBufferData(GL_ARRAY_BUFFER, num_points * 3 * sizeof(float), colors,
               GL_DYNAMIC_DRAW);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
  glEnableVertexAttribArray(1);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  glm_perspective(glm_rad(45.0f), 1.0f, 1.0f, 500.0f, projection);

  double last_time = 0.0;
  int nb_frames = 0;

  while (!glfwWindowShouldClose(window)) {
    double current_time = glfwGetTime();
    nb_frames++;
    if (current_time - last_time >= 1.0) {
      fps = nb_frames;
      nb_frames = 0;
      last_time += 1.0;
    }

    if (!paused) {
      for (int i = 0; i < render_speed; i++) {
        calculate_next_point();
      }
    }

    process_input(window);
    display(window);
  }

  glDeleteVertexArrays(1, &vao);
  glDeleteBuffers(1, &vbo);
  glDeleteBuffers(1, &cbo);
  glDeleteProgram(shader_program);

  glfwDestroyWindow(window);
  glfwTerminate();
  arena_destroy(&equation_arena);
  free(lorenz_points);
  free(colors);

  return 0;
}
