#include "types.h"

bool createBitmap(uint32_t width, uint32_t height, Bitmap* out);
void colorPixel(const Bitmap& bitmap, uint32_t x, uint32_t y, const Vec3& color);
bool saveAsPNG(const char* filename, Bitmap* bitmap);
bool saveAsPNGAndClose(const char* filename, Bitmap* bitmap);