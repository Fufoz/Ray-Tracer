#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <random>
#include <vector>

#include "maths.h"
#include "timer.h"

struct Bitmap
{
	uint32_t width;
	uint32_t height;
	void* pixels;
};

struct Film
{
	float aspect;
	float pixelSize;
	float worldWidth;
	float worldHeight;
	uint32_t sampleCount;
};

struct Camera
{
	Vec3 u, v, w;
	Vec3 origin;
	float vfov;
	float rlens;
	float znear;
	float zfocal;
};

struct Ray
{
	Vec3 origin;
	Vec3 direction;
};

struct PointLight
{
	float intensity;
	Vec3 color;//specular && diffuse cols
	Vec3 position;
};

struct Plane
{
	Vec3 normal;
	Vec3 point;
	uint32_t matId;
};

struct Sphere
{
	Vec3 position;
	float radius;
	uint32_t matId;
};

struct Triangle
{
	Vec3 v0;
	Vec3 v1;
	Vec3 v2;
	uint32_t matId;
};

enum MaterialType
{
	MTYPE_MIRROR,
	MTYPE_DIFFUSE,
	MTYPE_DIELECTRIC
};

struct Material
{
	Vec3 diffuseColor;
	Vec3 specularColor;
	MaterialType mtype;
	int specularGlossiness;
	float kr;
	float ior;//intex of refraction
	bool pattern;
};

struct Surfel
{
	Vec3 normal;
	Vec3 position;
	Material material;
};

struct Scene
{
	std::vector<Plane> planes;
	std::vector<Sphere> spheres;
	std::vector<Triangle> triangles;
	std::vector<Material> materials;
	std::vector<PointLight> pointLights;
};

struct RandomCtx
{
	std::default_random_engine rengine;
	std::uniform_real_distribution<float> distr;
};

static const float EPSILON = 0.0001f;

#endif