/* Locks the dependency-free PNG writer (this iteration): writes RGB+RGBA
 * and checks the PNG signature. (Pixel-level decode is verified
 * separately with PIL during development.) make test-png */
#define PNG_WRITER_IMPLEMENTATION
#include "../src/png_writer.h"
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#ifdef _WIN32
#include <process.h>
#define test_getpid _getpid
#else
#include <unistd.h>
#define test_getpid getpid
#endif
int main(){
  const char *tmp = std::getenv("TMPDIR");
  std::string prefix = std::string(tmp ? tmp : ".") + "/dynsys-png-" +
                       std::to_string(test_getpid());
  std::string rgb_path = prefix + ".png";
  std::string rgba_path = prefix + "-rgba.png";
  int W=64,H=48;
  std::vector<unsigned char> rgb((size_t)W*H*3);
  for(int y=0;y<H;y++)for(int x=0;x<W;x++){
    size_t i=((size_t)y*W+x)*3;
    rgb[i+0]=(unsigned char)(x*4);     // R ramp in x
    rgb[i+1]=(unsigned char)(y*5);     // G ramp in y
    rgb[i+2]=128;                      // B const
  }
  bool ok=png_write_rgb(rgb_path.c_str(), rgb.data(), W, H);
  printf("write rgb: %s\n", ok?"ok":"FAIL");
  // RGBA too
  std::vector<unsigned char> rgba((size_t)W*H*4,255);
  for(size_t p=0;p<(size_t)W*H;p++){rgba[p*4+0]=255;rgba[p*4+1]=0;rgba[p*4+2]=0;}
  bool rgba_ok = png_write_rgba(rgba_path.c_str(), rgba.data(), W, H);
  printf("write rgba: %s\n", rgba_ok?"ok":"FAIL");
    // verify the file starts with the PNG signature and is non-trivial
  FILE* fp=fopen(rgb_path.c_str(),"rb"); if(!fp){printf("no file\n");std::remove(rgba_path.c_str());return 1;}
  unsigned char sig[8]; size_t n=fread(sig,1,8,fp); fclose(fp);
  const unsigned char exp[8]={137,80,78,71,13,10,26,10};
  bool sigok = (n==8); for(int i=0;i<8;i++) sigok = sigok && (sig[i]==exp[i]);
  printf("PNG signature: %s\n", sigok?"valid":"INVALID");
  std::remove(rgb_path.c_str());
  std::remove(rgba_path.c_str());
  printf("=== %s ===\n", (ok && rgba_ok && sigok)?"PASS":"FAIL");
  return (ok && rgba_ok && sigok)?0:1;
}
