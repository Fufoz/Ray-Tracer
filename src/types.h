#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
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

struct PinholeCamera
{
	Vec3 u, v, w;
	Vec3 origin;
	float znear;
	float vfov;
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

struct SurfaceRecord
{
	Vec3 normal;
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
	MTYPE_DIFFUSE,
	MTYPE_REFLECTIVE,
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
};

struct Surfel
{
	Vec3 normal;
	Vec3 position;
	Material material;
};

struct Scene
{
	Sphere* spheres;
	uint32_t sphereCount;

	Plane* planes;
	uint32_t planesCount;

	Triangle* triangles;
	uint32_t triangleCount;

	Material* materials;
	uint32_t materialCount;

	PointLight* pointLights;
	uint32_t pointLightCount;
};


#endif