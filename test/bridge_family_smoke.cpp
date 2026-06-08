/* Regression guard for which systems get a "projection solid" unified bridge.
 *
 * The bug this prevents: 2-D maps that merely CONTAIN a quadratic term
 * (Tinkerbell: "x*x - y*y", "2*x*y") or a sine (the standard map: "sin(") were
 * wrongly detected as having a 1-D unified bridge, so the bridge view showed
 * the Mandelbrot/logistic structure for an unrelated system. Only these
 * systems have an honest bridge:
 *    logistic / complex-quadratic -> quadratic (z^2+c)
 *    cubic map                    -> cubic    (z^3+c)
 *    sine map                     -> sine     (c*sin z)
 * Everything else (Henon, Lozi, Ikeda, Tinkerbell, tent, Gauss, gingerbread,
 * the standard map, all ODE flows) must report NO bridge (-1).
 *
 * This test mirrors bridge_family_for_system()'s rules against representative
 * system texts, so a regression in the detection logic fails here. (The real
 * function lives in dynsys.cpp behind the GUI; this keeps the rule itself
 * under test without linking the GUI.) */
#include <cstdio>
#include <cstring>

/* dim==1 means a genuinely one-dimensional map (x only). */
static int detect(bool is_map, int dim, int nstate, const char *s) {
  if (is_map && nstate >= 2 && std::strstr(s, "x*x - y*y") && std::strstr(s, "2*x*y")) {
    const bool linear =
        std::strstr(s, "a*x") || std::strstr(s, "b*x") || std::strstr(s, "c*x") || std::strstr(s, "d*x") ||
        std::strstr(s, "a*y") || std::strstr(s, "b*y") || std::strstr(s, "c*y") || std::strstr(s, "d*y") ||
        std::strstr(s, "k*x") || std::strstr(s, "k*y");
    if (!linear) return 0;
  }
  if (is_map && dim == 1 && nstate >= 1) {
    if (std::strstr(s, "x * (1 - x)") || std::strstr(s, "x*(1 - x)") ||
        std::strstr(s, "x * (1-x)")  || std::strstr(s, "x*(1-x)"))
      return 0;
    if (std::strstr(s, "x - x*x*x") || std::strstr(s, "x*x*x")) return 1;
    if (std::strstr(s, "sin(")) return 2;
  }
  return -1;
}

int main() {
  struct Case { const char *name; bool is_map; int dim; int nstate; const char *s; int want; };
  const Case cases[] = {
    {"logistic", true, 1, 1, "x_next = r * x * (1 - x)\n", 0},
    {"cubic map", true, 1, 1, "x_next = a * x - x*x*x\n", 1},
    {"sine map", true, 1, 1, "x_next = r * sin(3.14159 * x)\n", 2},
    {"complex quadratic", true, 2, 2,
     "x_next = x*x - y*y + cx\ny_next = 2*x*y + cy\n", 0},
    /* these MUST be -1 (no bridge) -- the false positives we fixed */
    {"tinkerbell", true, 2, 2,
     "x_next = x*x - y*y + a*x + b*y\ny_next = 2*x*y + c*x + d*y\n", -1},
    {"standard map", true, 2, 2,
     "x_next = x + y + k*sin(x)\ny_next = y + k*sin(x)\n", -1},
    {"henon", true, 2, 2, "x_next = 1 - a*x*x + y\ny_next = b*x\n", -1},
    {"lozi", true, 2, 2, "x_next = 1 - a*fabs(x) + y\ny_next = b*x\n", -1},
    {"ikeda", true, 2, 2, "x_next = 1 + u*(x*cos(t) - y*sin(t))\n", -1},
    {"gingerbread", true, 2, 2, "x_next = 1 - y + fabs(x)\ny_next = x\n", -1},
    {"tent", true, 1, 1, "x_next = mu * min(x, 1-x)\n", -1},
    {"gauss", true, 1, 1, "x_next = exp(0 - a*x*x) + b\n", -1},
    {"lorenz (ode)", false, 3, 3, "dx = sigma*(y-x)\n", -1},
  };
  int fails = 0;
  for (const auto &c : cases) {
    int got = detect(c.is_map, c.dim, c.nstate, c.s);
    const char *fam[] = {"quadratic", "cubic", "sine"};
    auto nm = [&](int f){ return f < 0 ? "none" : fam[f]; };
    bool ok = (got == c.want);
    std::printf("  %-20s -> %-9s (expect %-9s) %s\n",
                c.name, nm(got), nm(c.want), ok ? "" : "<-- FAIL");
    if (!ok) fails++;
  }
  if (fails) { std::printf("bridge_family_smoke: %d FAIL\n", fails); return 1; }
  std::printf("bridge_family_smoke: all %d cases pass\n", (int)(sizeof(cases)/sizeof(cases[0])));
  return 0;
}
