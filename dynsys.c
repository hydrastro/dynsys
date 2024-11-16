#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <ds.h>
#include <lizard.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600
#define SCALE 10.0
#define MAX_POINTS 10000

double time_step = 0.01;

typedef struct {
  double x, y, z;
} vector3;

vector3 attractor_points[MAX_POINTS];
int point_count = 0;

double camera_theta = 0.0;
double camera_phi = 0.0;

void draw_text(Display *display, Window window, GC gc, int x, int y,
               const char *text) {
  XDrawString(display, window, gc, x, y, text, strlen(text));
}

vector3 apply_camera_rotation(vector3 p) {
  double sin_theta = sin(camera_theta), cos_theta = cos(camera_theta);
  double sin_phi = sin(camera_phi), cos_phi = cos(camera_phi);

  double x = p.x * cos_theta - p.z * sin_theta;
  double z = p.x * sin_theta + p.z * cos_theta;
  p.x = x;
  p.z = z;

  double y = p.y * cos_phi - p.z * sin_phi;
  z = p.y * sin_phi + p.z * cos_phi;
  p.y = y;
  p.z = z;

  return p;
}

void to_screen_coordinates(vector3 p, int *screen_x, int *screen_y) {
  *screen_x = (int)(p.x * SCALE) + WINDOW_WIDTH / 2;
  *screen_y = (int)(p.y * SCALE) + WINDOW_HEIGHT / 2;
}

double evaluate_ast(lizard_ast_node_t *ast, vector3 *state) {
  switch (ast->type) {
  case AST_NUMBER:
    return mpz_get_d(ast->number);
  case AST_SYMBOL: {
    if (strcmp(ast->variable, "x") == 0) {
      return state->x;
    }
    if (strcmp(ast->variable, "y") == 0) {
      return state->y;
    }
    if (strcmp(ast->variable, "z") == 0) {
      return state->z;
    }
    fprintf(stderr, "Unknown variable: %s\n", ast->variable);
    exit(1);
  }
  case AST_APPLICATION: {
    lizard_ast_node_t *op_node =
        &CAST(ast->application_arguments->head, lizard_ast_list_node_t)->ast;
    lizard_ast_node_t *arg1 =
        &CAST(ast->application_arguments->head->next, lizard_ast_list_node_t)
             ->ast;
    lizard_ast_node_t *arg2 =
        &CAST(ast->application_arguments->head->next->next,
              lizard_ast_list_node_t)
             ->ast;
    char *op = op_node->variable;
    double val1 = evaluate_ast(arg1, state);
    double val2 = evaluate_ast(arg2, state);
    if (strcmp(op, "+") == 0)
      return val1 + val2;
    if (strcmp(op, "*") == 0)
      return val1 * val2;
    if (strcmp(op, "-") == 0)
      return val1 - val2;
    if (strcmp(op, "/") == 0)
      return val1 / val2;
  } break;
  default:
    fprintf(stderr, "Unsupported AST node type: %d \n", ast->type);
    exit(1);
  }
}

void update_lorenz_attractor(vector3 *p, lizard_ast_node_t **equations) {
  double dx, dy, dz;

  dx = evaluate_ast(equations[0], p);
  dy = evaluate_ast(equations[1], p);
  dz = evaluate_ast(equations[2], p);

  p->x += dx * time_step;
  p->y += dy * time_step;
  p->z += dz * time_step;
}

void draw_axes(Display *display, Window window, GC gc) {
  vector3 origin = {0, 0, 0};
  vector3 x_axis = {5, 0, 0};
  vector3 y_axis = {0, 5, 0};
  vector3 z_axis = {0, 0, 5};

  origin = apply_camera_rotation(origin);
  x_axis = apply_camera_rotation(x_axis);
  y_axis = apply_camera_rotation(y_axis);
  z_axis = apply_camera_rotation(z_axis);

  int ox, oy, x1, y1, x2, y2, x3, y3;
  to_screen_coordinates(origin, &ox, &oy);
  to_screen_coordinates(x_axis, &x1, &y1);
  to_screen_coordinates(y_axis, &x2, &y2);
  to_screen_coordinates(z_axis, &x3, &y3);

  XDrawLine(display, window, gc, ox, oy, x1, y1);
  draw_text(display, window, gc, x1 + 5, y1, "X");
  XDrawLine(display, window, gc, ox, oy, x2, y2);
  draw_text(display, window, gc, x2 + 5, y2, "Y");
  XDrawLine(display, window, gc, ox, oy, x3, y3);
  draw_text(display, window, gc, x3 + 5, y3, "Z");
}

int main() {
  Display *display = XOpenDisplay(NULL);
  if (display == NULL) {
    fprintf(stderr, "Cannot open display\n");
    exit(1);
  }

  int screen = DefaultScreen(display);
  Window window = XCreateSimpleWindow(
      display, RootWindow(display, screen), 10, 10, WINDOW_WIDTH, WINDOW_HEIGHT,
      1, BlackPixel(display, screen), WhitePixel(display, screen));
  XSelectInput(display, window, ExposureMask | KeyPressMask);
  XMapWindow(display, window);

  GC gc = XCreateGC(display, window, 0, NULL);
  XSetForeground(display, gc, BlackPixel(display, screen));

  char *input1 = "( * 10 (- y x))";
  char *input2 = "( - (* x (- 28 z)) y)";
  char *input3 = "( - (* x y) (* (/ 8 3) z))";

  list_t *ast_list1 = lizard_parse(lizard_tokenize(input1));
  list_t *ast_list2 = lizard_parse(lizard_tokenize(input2));
  list_t *ast_list3 = lizard_parse(lizard_tokenize(input3));
  lizard_ast_node_t *equations[3] = {
      &CAST(ast_list1->head, lizard_ast_list_node_t)->ast,
      &CAST(ast_list2->head, lizard_ast_list_node_t)->ast,
      &CAST(ast_list3->head, lizard_ast_list_node_t)->ast};

  vector3 point = {0.1, 0.1, 0.1};
  point_count = 0;

  clock_t last_time = clock();
  int frame_count = 0;
  double fps = 0.0;

  while (1) {
    XEvent event;
    while (XPending(display)) {
      XNextEvent(display, &event);

      if (event.type == Expose) {
        XClearWindow(display, window);
      }

      if (event.type == KeyPress) {
        KeySym key = XLookupKeysym(&event.xkey, 0);
        if (key == XK_Escape) {
          XCloseDisplay(display);
          return 0;
        } else if (key == XK_w) {
          camera_phi -= 0.05;
        } else if (key == XK_s) {
          camera_phi += 0.05;
        } else if (key == XK_a) {
          camera_theta -= 0.05;
        } else if (key == XK_d) {
          camera_theta += 0.05;
        } else if (key == XK_plus || key == XK_equal) {
          time_step *= 1.1;
        } else if (key == XK_minus || key == XK_underscore) {
          time_step *= 0.9;
        } else if (key == XK_c) {
          point_count = 0;
          point.x = 0.1;
          point.y = 0.0;
          point.z = 0.0;
        }
      }
    }

    if (point_count < MAX_POINTS) {
      update_lorenz_attractor(&point, equations);
      attractor_points[point_count++] = point;
    }

    XClearWindow(display, window);
    draw_axes(display, window, gc);

    for (int i = 1; i < point_count; i++) {
      vector3 rotated_point1 = apply_camera_rotation(attractor_points[i - 1]);
      vector3 rotated_point2 = apply_camera_rotation(attractor_points[i]);

      int x1, y1, x2, y2;
      to_screen_coordinates(rotated_point1, &x1, &y1);
      to_screen_coordinates(rotated_point2, &x2, &y2);

      XDrawLine(display, window, gc, x1, y1, x2, y2);
    }

    frame_count++;
    clock_t current_time = clock();
    double elapsed_time = (double)(current_time - last_time) / CLOCKS_PER_SEC;
    if (elapsed_time >= 1.0) {
      fps = frame_count / elapsed_time;
      frame_count = 0;
      last_time = current_time;
    }

    char buffer[256];
    sprintf(buffer, "FPS: %.2f | Speed: %.4f | Points: %d", fps, time_step,
            point_count);
    draw_text(display, window, gc, 10, 20, buffer);

    usleep(1000);
  }

  XCloseDisplay(display);
  return 0;
}
