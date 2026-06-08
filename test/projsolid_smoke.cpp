#include <cstdio>
#include <cmath>
#include <vector>
#include <set>
// Replicate the ProjectionSolid slice logic and verify both shadows.
static int mandel_escape(double c,double im,int maxit){
  double zr=0,zi=0; for(int i=0;i<maxit;i++){double r2=zr*zr-zi*zi+c,i2=2*zr*zi+im;zr=r2;zi=i2;if(zr*zr+zi*zi>16)return i;} return maxit; }
static double r_of_c(double c){double d=1-4*c;return d>0?1+std::sqrt(d):1;}
struct P{double x,y,z;};
static int fails=0; static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
int main(){
  const double c_lo=-2.05,c_hi=0.30,im_hi=1.20;
  const int NC=400,NIM=120,settle=1500,keep=80,maxit=200;
  std::vector<P> pts;
  for(int ic=0;ic<NC;ic++){
    double c=c_lo+(c_hi-c_lo)*ic/(NC-1);
    std::vector<double> imhits;
    for(int k=0;k<NIM;k++){double im=(double)k/(NIM-1)*im_hi; if(mandel_escape(c,im,maxit)>=maxit) imhits.push_back(im);}
    double r=r_of_c(c); std::vector<double> yf;
    double x=0.3137; for(int k=0;k<settle;k++)x=r*x*(1-x);
    for(int k=0;k<keep;k++){x=r*x*(1-x); if(std::isfinite(x))yf.push_back(x);}
    bool hs=!imhits.empty(),hf=!yf.empty();
    if(hs&&hf){ for(double im:imhits) for(double y:yf){ pts.push_back({c,y,im}); if(im>1e-6)pts.push_back({c,y,-im}); } }
    else if(hf){ for(double y:yf) pts.push_back({c,y,0}); }
    else if(hs){ for(double im:imhits){ pts.push_back({c,0,im}); if(im>1e-6)pts.push_back({c,0,-im}); } }
  }
  printf("built %zu points\n",pts.size());

  // SHADOW 1 (X-Z): for a grid of (c,im), a point with that (c,im) should exist IFF (c,im) in Mandelbrot.
  // Check agreement on a coarse grid.
  int agree=0,total=0;
  auto has_xz=[&](double c,double im,double tol)->bool{ for(auto&p:pts) if(std::fabs(p.x-c)<tol&&std::fabs(p.z-im)<tol) return true; return false; };
  for(double c=-2.0;c<=0.25;c+=0.05) for(double im=0;im<=1.0;im+=0.05){
    bool inset=mandel_escape(c,im,maxit)>=maxit;
    bool inproj=has_xz(c,im,0.03);
    // projection should contain set points; allow proj to be superset only on axis fallback rows
    if(inset){ total++; if(inproj)agree++; }
  }
  printf("X-Z shadow: %d/%d Mandelbrot grid points present in projection\n",agree,total);
  chk("X-Z shadow covers the Mandelbrot set", total>0 && agree>=0.9*total);

  // SHADOW 2 (X-Y): the set of y-values at each c equals the logistic fiber.
  // Check: at r=3.2 (period 2) there are ~2 distinct y bands; at r=3.5 (period 4) ~4.
  auto ybands=[&](double rtarget){ double c=rtarget/2-rtarget*rtarget/4; std::set<long> bands; for(auto&p:pts) if(std::fabs(p.x-c)<0.01){ bands.insert((long)std::llround(p.y*200)); } return bands.size(); };
  size_t b32=ybands(3.2), b35=ybands(3.5), b29=ybands(2.9);
  printf("X-Y shadow fiber widths: r=2.9 -> %zu band(s), r=3.2 -> %zu, r=3.5 -> %zu\n",b29,b32,b35);
  chk("r=2.9 is period-1 (1 narrow band)", b29<=3);
  chk("r=3.2 splits (period-2, more bands than r=2.9)", b32>b29);
  chk("r=3.5 splits further (period-4, more than r=3.2)", b35>=b32);

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
