/* Locks the dependency-free PNG writer (this iteration): writes RGB+RGBA
 * and checks the PNG signature. (Pixel-level decode is verified
 * separately with PIL during development.) make test-png */
#define PNG_WRITER_IMPLEMENTATION
#include "../src/png_writer.h"
#include <vector>
#include <cstdio>
int main(){
  int W=64,H=48;
  std::vector<unsigned char> rgb((size_t)W*H*3);
  for(int y=0;y<H;y++)for(int x=0;x<W;x++){
    size_t i=((size_t)y*W+x)*3;
    rgb[i+0]=(unsigned char)(x*4);     // R ramp in x
    rgb[i+1]=(unsigned char)(y*5);     // G ramp in y
    rgb[i+2]=128;                      // B const
  }
  bool ok=png_write_rgb("/tmp/pngtest.png", rgb.data(), W, H);
  printf("write rgb: %s\n", ok?"ok":"FAIL");
  // RGBA too
  std::vector<unsigned char> rgba((size_t)W*H*4,255);
  for(size_t p=0;p<(size_t)W*H;p++){rgba[p*4+0]=255;rgba[p*4+1]=0;rgba[p*4+2]=0;}
  printf("write rgba: %s\n", png_write_rgba("/tmp/pngtest_rgba.png", rgba.data(), W, H)?"ok":"FAIL");
    // verify the file starts with the PNG signature and is non-trivial
  FILE* fp=fopen("/tmp/pngtest.png","rb"); if(!fp){printf("no file\n");return 1;}
  unsigned char sig[8]; size_t n=fread(sig,1,8,fp); fclose(fp);
  const unsigned char exp[8]={137,80,78,71,13,10,26,10};
  bool sigok = (n==8); for(int i=0;i<8;i++) sigok = sigok && (sig[i]==exp[i]);
  printf("PNG signature: %s\n", sigok?"valid":"INVALID");
  printf("=== %s ===\n", (ok && sigok)?"PASS":"FAIL");
  return (ok && sigok)?0:1;
}
