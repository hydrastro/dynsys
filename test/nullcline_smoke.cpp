/* Verifies the marching-squares contour logic used by the phase-plane
 * nullcline renderer (src/dynsys.cpp draw_nullcline_for_component).
 *
 * The renderer can't be unit-tested directly because it draws into an
 * ImGui draw list, but the geometric core — the 16-case contour of a
 * scalar field's zero set over a grid — is reproduced here verbatim
 * and checked against an analytic circle f(x,y)=x^2+y^2-4. Every
 * emitted vertex must lie on the circle to within the grid resolution.
 *
 *   make test-nullcline
 */

#include <cmath>
#include <cstdio>
#include <vector>

namespace {
struct P {
  double x, y;
};
double field(double x, double y) { return x * x + y * y - 4.0; }
}  // namespace

int main() {
  std::printf("=== nullcline marching-squares test ===\n");
  const int grid = 80;
  const double lo = -3.0, hi = 3.0;
  const double d = (hi - lo) / grid;
  const int nx = grid + 1;

  std::vector<double> f(static_cast<size_t>(nx) * nx);
  for (int j = 0; j <= grid; ++j)
    for (int i = 0; i <= grid; ++i)
      f[static_cast<size_t>(j) * nx + i] = field(lo + i * d, lo + j * d);

  auto node = [&](int i, int j) { return f[static_cast<size_t>(j) * nx + i]; };
  auto interp = [](double a, double fa, double b, double fb) {
    const double dn = fb - fa;
    return std::fabs(dn) < 1e-300 ? 0.5 * (a + b) : a + (b - a) * (-fa / dn);
  };

  std::vector<P> verts;
  int segs = 0;
  for (int j = 0; j < grid; ++j)
    for (int i = 0; i < grid; ++i) {
      const double x0 = lo + i * d, x1 = x0 + d, y0 = lo + j * d, y1 = y0 + d;
      const double f00 = node(i, j), f10 = node(i + 1, j),
                   f11 = node(i + 1, j + 1), f01 = node(i, j + 1);
      int c = 0;
      if (f00 > 0) c |= 1;
      if (f10 > 0) c |= 2;
      if (f11 > 0) c |= 4;
      if (f01 > 0) c |= 8;
      if (c == 0 || c == 15) continue;
      const P eb{interp(x0, f00, x1, f10), y0};
      const P er{x1, interp(y0, f10, y1, f11)};
      const P et{interp(x0, f01, x1, f11), y1};
      const P el{x0, interp(y0, f00, y1, f01)};
      auto add = [&](P a, P b) {
        verts.push_back(a);
        verts.push_back(b);
        ++segs;
      };
      switch (c) {
        case 1: case 14: add(el, eb); break;
        case 2: case 13: add(eb, er); break;
        case 3: case 12: add(el, er); break;
        case 4: case 11: add(er, et); break;
        case 6: case 9:  add(eb, et); break;
        case 7: case 8:  add(el, et); break;
        case 5:  add(el, eb); add(er, et); break;
        case 10: add(eb, er); add(el, et); break;
        default: break;
      }
    }

  double maxerr = 0.0;
  for (const P &p : verts)
    maxerr = std::fmax(maxerr, std::fabs(std::sqrt(p.x * p.x + p.y * p.y) - 2.0));

  const double cell_diag = std::sqrt(2.0) * d;
  std::printf("segments=%d vertices=%zu max|r-2|=%.4g (cell diag=%.4g)\n", segs,
              verts.size(), maxerr, cell_diag);
  const bool ok = segs > 100 && maxerr < cell_diag;
  std::printf("=== %s ===\n", ok ? "PASS" : "FAIL");
  return ok ? 0 : 1;
}
