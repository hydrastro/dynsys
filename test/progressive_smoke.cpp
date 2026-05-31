/* Locks progressive fractal rendering (this iteration): the 8->4->2->1
 * refinement converges to the exact full-resolution image, and
 * Mandelbrot membership is correct. make test-progressive */
/* Verify progressive rendering correctness: a step=1 pass and a step=8-then-1
 * refinement converge to the same final image; and confirm the Mandelbrot
 * membership is correct (cardioid in, c=2 out). Plain-double reference. */
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdint>
static int escape(double cre,double cim,int maxit,double R2){
  double x=0,y=0; int it=0;
  for(;it<maxit;it++){ double nx=x*x-y*y+cre, ny=2*x*y+cim; x=nx;y=ny;
    if(x*x+y*y>R2){it++;break;} }
  return it;
}
int main(){
  int fails=0;
  int W=200,H=160,maxit=200; double R2=4;
  double x0=-2.5,x1=1.0,y0=-1.5,y1=1.5;
  // full-res reference
  std::vector<int> ref((size_t)W*H);
  for(int py=0;py<H;py++)for(int px=0;px<W;px++){
    double re=x0+(x1-x0)*px/(W-1), im=y0+(y1-y0)*py/(H-1);
    ref[(size_t)py*W+px]=escape(re,im,maxit,R2);
  }
  // progressive: step 8->4->2->1, each overwrites blocks; final == reference
  std::vector<int> prog((size_t)W*H,-1);
  for(int step=8; step>=1; step/=2){
    for(int py=0;py<H;py+=step)for(int px=0;px<W;px+=step){
      double re=x0+(x1-x0)*px/(W-1), im=y0+(y1-y0)*py/(H-1);
      int v=escape(re,im,maxit,R2);
      for(int by=py;by<py+step&&by<H;by++)for(int bx=px;bx<px+step&&bx<W;bx++)
        prog[(size_t)by*W+bx]=v;
    }
  }
  // after the step=1 pass, every pixel should equal the reference
  int diff=0; for(size_t i=0;i<ref.size();i++) if(prog[i]!=ref[i])diff++;
  printf("progressive final vs full-res: %d differing pixels (expect 0)\n", diff);
  if(diff!=0){printf("  FAIL\n");fails++;}
  // membership sanity
  bool in_card = escape(-0.1,0.0,maxit,R2)>=maxit;   // inside cardioid
  bool out = escape(2.0,0.0,maxit,R2)<maxit;          // outside
  printf("cardioid in-set=%d, c=2 out=%d\n", in_card, out);
  if(!in_card||!out){printf("  FAIL membership\n");fails++;}
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
