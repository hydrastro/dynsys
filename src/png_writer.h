/* Minimal dependency-free PNG writer.
 *
 * Writes 8-bit RGB or RGBA images. Uses uncompressed (stored) DEFLATE
 * blocks wrapped in a zlib stream, so there is no compression-library
 * dependency — the file is valid PNG, just not size-optimized. Good enough
 * for screenshots / figure export.
 *
 * Public API:
 *   bool png_write_rgb(const char* path, const unsigned char* rgb,
 *                      int w, int h);   // 3 bytes/pixel, top row first
 *   bool png_write_rgba(const char* path, const unsigned char* rgba,
 *                       int w, int h);  // 4 bytes/pixel, top row first
 *
 * Header-only; define PNG_WRITER_IMPLEMENTATION in exactly one TU.
 */
#ifndef DYNSYS_PNG_WRITER_H
#define DYNSYS_PNG_WRITER_H

bool png_write_rgb(const char *path, const unsigned char *rgb, int w, int h);
bool png_write_rgba(const char *path, const unsigned char *rgba, int w, int h);

#ifdef PNG_WRITER_IMPLEMENTATION
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <vector>

namespace dynsys_png_detail {

inline uint32_t crc32_update(uint32_t crc, const unsigned char *buf, size_t len) {
  static uint32_t table[256];
  static bool init = false;
  if (!init) {
    for (uint32_t n = 0; n < 256; ++n) {
      uint32_t c = n;
      for (int k = 0; k < 8; ++k) c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
      table[n] = c;
    }
    init = true;
  }
  crc ^= 0xFFFFFFFFu;
  for (size_t i = 0; i < len; ++i) crc = table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
  return crc ^ 0xFFFFFFFFu;
}

inline uint32_t adler32(const unsigned char *data, size_t len) {
  uint32_t a = 1, b = 0;
  const uint32_t MOD = 65521;
  for (size_t i = 0; i < len; ++i) { a = (a + data[i]) % MOD; b = (b + a) % MOD; }
  return (b << 16) | a;
}

inline void put_u32(std::vector<unsigned char> &v, uint32_t x) {
  v.push_back((x >> 24) & 0xFF); v.push_back((x >> 16) & 0xFF);
  v.push_back((x >> 8) & 0xFF);  v.push_back(x & 0xFF);
}

inline void write_chunk(FILE *f, const char *type, const std::vector<unsigned char> &data) {
  unsigned char len[4] = { (unsigned char)((data.size() >> 24) & 0xFF),
                           (unsigned char)((data.size() >> 16) & 0xFF),
                           (unsigned char)((data.size() >> 8) & 0xFF),
                           (unsigned char)(data.size() & 0xFF) };
  fwrite(len, 1, 4, f);
  fwrite(type, 1, 4, f);
  if (!data.empty()) fwrite(data.data(), 1, data.size(), f);
  uint32_t crc = crc32_update(0, (const unsigned char *)type, 4);
  if (!data.empty()) crc = crc32_update(crc ^ 0xFFFFFFFFu, data.data(), data.size()) ; // continue
  /* simpler: recompute over type+data */
  std::vector<unsigned char> tmp;
  tmp.insert(tmp.end(), type, type + 4);
  tmp.insert(tmp.end(), data.begin(), data.end());
  uint32_t c = crc32_update(0, tmp.data(), tmp.size());
  unsigned char cb[4] = { (unsigned char)((c >> 24) & 0xFF), (unsigned char)((c >> 16) & 0xFF),
                          (unsigned char)((c >> 8) & 0xFF), (unsigned char)(c & 0xFF) };
  fwrite(cb, 1, 4, f);
}

inline bool write_png(const char *path, const unsigned char *pix, int w, int h, int channels) {
  if (w <= 0 || h <= 0 || (channels != 3 && channels != 4)) return false;
  FILE *f = std::fopen(path, "wb");
  if (!f) return false;
  const unsigned char sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
  fwrite(sig, 1, 8, f);

  /* IHDR */
  std::vector<unsigned char> ihdr;
  put_u32(ihdr, (uint32_t)w); put_u32(ihdr, (uint32_t)h);
  ihdr.push_back(8); /* bit depth */
  ihdr.push_back(channels == 4 ? 6 : 2); /* color type: 6=RGBA, 2=RGB */
  ihdr.push_back(0); ihdr.push_back(0); ihdr.push_back(0);
  write_chunk(f, "IHDR", ihdr);

  /* raw image data with filter byte (0) per scanline */
  const size_t stride = (size_t)w * channels;
  std::vector<unsigned char> raw;
  raw.reserve((stride + 1) * h);
  for (int y = 0; y < h; ++y) {
    raw.push_back(0);
    raw.insert(raw.end(), pix + (size_t)y * stride, pix + (size_t)y * stride + stride);
  }

  /* zlib stream: 2-byte header + stored DEFLATE blocks + adler32 */
  std::vector<unsigned char> z;
  z.push_back(0x78); z.push_back(0x01); /* CMF, FLG (no compression) */
  size_t pos = 0;
  const size_t MAX = 65535;
  while (pos < raw.size()) {
    size_t block = raw.size() - pos;
    if (block > MAX) block = MAX;
    const int final_block = (pos + block >= raw.size()) ? 1 : 0;
    z.push_back((unsigned char)final_block); /* BFINAL=final, BTYPE=00 stored */
    z.push_back((unsigned char)(block & 0xFF));
    z.push_back((unsigned char)((block >> 8) & 0xFF));
    z.push_back((unsigned char)(~block & 0xFF));
    z.push_back((unsigned char)((~block >> 8) & 0xFF));
    z.insert(z.end(), raw.begin() + pos, raw.begin() + pos + block);
    pos += block;
  }
  uint32_t ad = adler32(raw.data(), raw.size());
  z.push_back((ad >> 24) & 0xFF); z.push_back((ad >> 16) & 0xFF);
  z.push_back((ad >> 8) & 0xFF);  z.push_back(ad & 0xFF);

  write_chunk(f, "IDAT", z);
  std::vector<unsigned char> end;
  write_chunk(f, "IEND", end);
  std::fclose(f);
  return true;
}

}  // namespace dynsys_png_detail

bool png_write_rgb(const char *path, const unsigned char *rgb, int w, int h) {
  return dynsys_png_detail::write_png(path, rgb, w, h, 3);
}
bool png_write_rgba(const char *path, const unsigned char *rgba, int w, int h) {
  return dynsys_png_detail::write_png(path, rgba, w, h, 4);
}

#endif  /* PNG_WRITER_IMPLEMENTATION */
#endif  /* DYNSYS_PNG_WRITER_H */
