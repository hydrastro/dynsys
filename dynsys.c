#include <GL/glut.h>
#include <ds.h>
#include <lizard.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

clock_t last_time = 0;
int fps = 0;

lizard_ast_node_t *equations[3];

double evaluate_ast(lizard_ast_node_t *ast) {
  switch (ast->type) {
  case AST_NUMBER:
    return mpz_get_d(ast->number);
  case AST_SYMBOL: {
    if (strcmp(ast->variable, "x") == 0) {
      return current_x;
    }
    if (strcmp(ast->variable, "y") == 0) {
      return current_y;
    }
    if (strcmp(ast->variable, "z") == 0) {
      return current_z;
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
    double val1 = evaluate_ast(arg1);
    double val2 = evaluate_ast(arg2);
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

  current_point = (current_point + 1) % num_points;
}

void draw_attractor() {
  glLineWidth(1.0f);
  glBegin(GL_LINE_STRIP);
  for (int i = 0; i < current_point; ++i) {
    glColor3f((i % 256) / 256.0, 0.5, (255 - (i % 256)) / 256.0);
    glVertex3d(lorenz_points[i].x, lorenz_points[i].y, lorenz_points[i].z);
  }
  glEnd();
}

void drawAxes() {
  if (!show_axes)
    return;

  glLineWidth(1.0f);
  glBegin(GL_LINES);
  glColor3f(1.0f, 0.0f, 0.0f);
  glVertex3f(-20.0f, 0.0f, 0.0f);
  glVertex3f(20.0f, 0.0f, 0.0f);

  glColor3f(0.0f, 1.0f, 0.0f);
  glVertex3f(0.0f, -20.0f, 0.0f);
  glVertex3f(0.0f, 20.0f, 0.0f);

  glColor3f(0.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, -20.0f);
  glVertex3f(0.0f, 0.0f, 20.0f);
  glEnd();
}
void drawCenterCross() {
  if (!center_cross) {
    return;
  }
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, 1, 0, 1);

  glMatrixMode(GL_MODELVIEW);

  glColor3f(0.5f, 0.5f, 0.5f);

  glBegin(GL_LINES);
  glVertex2f(0.5f - 0.02f, 0.5f);
  glVertex2f(0.5f + 0.02f, 0.5f);
  glVertex2f(0.5f, 0.5f - 0.02f);
  glVertex2f(0.5f, 0.5f + 0.02f);

  glEnd();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
}

void draw_hud() {
  char info[128];
  snprintf(info, sizeof(info), "FPS: %d  Speed: %d  Points: %d", fps,
           render_speed, current_point);

  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, 800, 0, 800);
  glMatrixMode(GL_MODELVIEW);

  glColor3f(1.0f, 1.0f, 1.0f);
  glRasterPos2i(10, 780);
  for (char *c = info; *c != '\0'; c++) {
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *c);
  }

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glTranslatef(0.0f, 0.0f, -50.0f * zoom);
  glRotatef(angle_x, 1.0f, 0.0f, 0.0f);
  glRotatef(angle_y, 0.0f, 1.0f, 0.0f);

  glTranslatef(center_x, center_y, center_z);
  glScalef(0.05, 0.05, 0.05);

  drawAxes();
  draw_attractor();
  draw_hud();
  drawCenterCross();
  glutSwapBuffers();
}

void reshape(int w, int h) {
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, (float)w / (float)h, 1.0, 500.0);
  glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y) {
  switch (key) {
  case 'w':
    angle_x += 5.0f;
    break;
  case 's':
    angle_x -= 5.0f;
    break;
  case 'a':
    angle_y -= 5.0f;
    break;
  case 'd':
    angle_y += 5.0f;
    break;
  case '+':
    zoom *= 0.9f;
    break;
  case '-':
    zoom /= 0.9f;
    break;
  case ' ':
    paused = !paused;
    break;
  case 'x':
    show_axes = !show_axes;
    break;
  case '>':
    render_speed *= 2;
    break;
  case '<':
    render_speed /= 2;
    if (render_speed == 0) {
      render_speed = 1;
    }
    break;
  case 'c':
    reset_simulation();
    break;
  case 'j':
    center_x += 0.5f;
    break;
  case 'l':
    center_x -= 0.5f;
    break;
  case 'i':
    center_y += 0.5f;
    break;
  case 'k':
    center_y -= 0.5f;
    break;
  case 'p':
    center_z += 0.5f;
    break;
  case ';':
    center_z -= 0.5f;
    break;
  case 'z':
    center_cross = !center_cross;
    break;
  case 27:
    free(lorenz_points);
    exit(0);
    break;
  }
  glutPostRedisplay();
}

void timer(int value) {
  if (!paused) {
    for (int i = 0; i < render_speed; i++) {
      calculate_next_point();
    }
  }

  clock_t now = clock();
  fps = CLOCKS_PER_SEC / (now - last_time + 1);
  last_time = now;

  glutPostRedisplay();
  glutTimerFunc(animation_speed_ms, timer, 0);
}

void init_opengl() {
  glEnable(GL_DEPTH_TEST);
  glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0, 1.0, 1.0, 500.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

int main(int argc, char **argv) {
  lorenz_points = (point_t *)malloc(num_points * sizeof(point_t));
  char *input1 = "( * 10 (- y x))";
  char *input2 = "( - (* x (- 28 z)) y)";
  char *input3 = "( - (* x y) (* (/ 8 3) z))";
  list_t *ast_list1 = lizard_parse(lizard_tokenize(input1));
  list_t *ast_list2 = lizard_parse(lizard_tokenize(input2));
  list_t *ast_list3 = lizard_parse(lizard_tokenize(input3));
  equations[0] = &CAST(ast_list1->head, lizard_ast_list_node_t)->ast;
  equations[1] = &CAST(ast_list2->head, lizard_ast_list_node_t)->ast;
  equations[2] = &CAST(ast_list3->head, lizard_ast_list_node_t)->ast;
  start_x = 0.1;
  start_y = 0.1;
  start_z = 0.1;
  current_point = 0;
  for (int i = 0; i < num_points; ++i) {
    lorenz_points[i].x = 0.0;
    lorenz_points[i].y = 0.0;
    lorenz_points[i].z = 0.0;
  }

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(800, 800);
  glutCreateWindow("Lorenz Attractor with Advanced Controls");
  init_opengl();

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutTimerFunc(animation_speed_ms, timer, 0);

  glutMainLoop();
  return 0;
}
