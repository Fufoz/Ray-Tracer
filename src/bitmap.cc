#include "bitmap.h"
#include "logging.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

bool createBitmap(uint32_t width, uint32_t height, Bitmap* out)
{
	assert(width > 0);
	assert(height > 0);

	void* blob = malloc(width * height * 3);
	if(!blob)
	{
		hoth::log::error("Failed to allocate bitmap memory!");
		return false;
	}

	//set pixels start from bottom left corner
	int FLIP_IMAGE = 1;
	stbi_flip_vertically_on_write(FLIP_IMAGE);

	out->width = width;
	out->height = height;
	out->pixels = blob;
	return true;
}

bool saveAsPNG(const char* filename, Bitmap* bitmap)
{
	return stbi_write_png(filename, bitmap->width, bitmap->height, 3, bitmap->pixels, bitmap->width * 3);
}

bool saveAsPNGAndClose(const char* filename, Bitmap* bitmap)
{
	int status = stbi_write_png(filename, bitmap->width, bitmap->height, 3, bitmap->pixels, bitmap->width * 3);
	if(status != 1)
	{
		hoth::log::error("Failed to write bitmap as png file!");
		return false;
	}

	assert(bitmap->pixels);
	free(bitmap->pixels);
	bitmap->width = 0;
	bitmap->height = 0;
	bitmap->pixels = nullptr;
	return true;
}

void colorPixel(const Bitmap& bitmap, uint32_t x, uint32_t y, const Vec3& color)
{
	assert(x < bitmap.width);
	assert(y < bitmap.height);
	uint8_t* pixelPtr = (uint8_t*)bitmap.pixels;
	constexpr uint8_t numc = 3;
	pixelPtr[(x + y * bitmap.width) * numc + 0] = (uint8_t)color.R;
	pixelPtr[(x + y * bitmap.width) * numc + 1] = (uint8_t)color.G;
	pixelPtr[(x + y * bitmap.width) * numc + 2] = (uint8_t)color.B;
}