#include <random>
#include <cassert>
#include "types.h"
#include "bitmap.h"
#include "logging.h"

static Ray getPrimaryRay(float x, float y, uint32_t width, uint32_t height, const Film& film, const PinholeCamera& camera)
{
	//pixel to world image plane
	Vec2 filmPoint = 
	{
		x * film.pixelSize - film.worldWidth / 2.f,
		y * film.pixelSize - film.worldHeight / 2.f,
	};

	Ray ray = {};
	ray.origin = camera.origin;
	ray.direction = normaliseVec3(filmPoint.x * camera.u + filmPoint.y * camera.v - camera.znear * camera.w);
	
	return ray;
}

static inline Vec3 reflectRay(const Vec3& incident, const Vec3& normal)
{
	return incident - 2.f * dotVec3(normal, incident) * normal;
}

static inline void fresnel(const Vec3& incident, const Vec3& normal, float ior, float* kr)
{
	float n1 = 1;//assume incoming medium is air
	float n2 = ior;

	float cosI = clamp(dotVec3(incident, normal), -1.f, 1.f);
	if(cosI > 0)
	{
		std::swap(n1, n2);
	}
	else
	{
		cosI = -cosI;//we need positive cosine
	}

	float sinT = n1 / n2 * sqrtf(std::max(0.f, 1.f - cosI * cosI));
	if(sinT >= 1)
	{
		*kr = 1;//total internal reflection
	}
	else
	{//Schlick's approximation
		float R0 = ((n1 - n2) * (n1 - n2)) / ((n1 + n2) * (n1 + n2));
		*kr = R0 + (1 - R0) * std::powf((1 - cosI), 5.f);
		assert((*kr) < 1.f);
	}
}

static inline Vec3 refractRay(const Vec3& incident, const Vec3& normal, float ior, bool* refracted)
{
	*refracted = false;
	Vec3 refractedRay = {};
	Vec3 norm = normal;

	float n1 = 1;//assume incoming medium is air
	float n2 = ior;

	float cosI = clamp(dotVec3(incident, norm), -1.f, 1.f);

	if(cosI < 0.f)
	{
		cosI = -cosI;//we need positive cosine
	}
	else //if normal and incident ray are collinear (i.e exiting the surface to air)
	{
		norm = -1.f * norm; //reverse normal
		std::swap(n1, n2); // swap refractive indices
	}

	float n = n1 / n2;
	
	float sinT2 = n * n * (1.f - cosI * cosI);
	if(sinT2 > 1.f)
	{
		return refractedRay;
	}

	//probably need to swap signs here : DONE
	refractedRay = n * incident + (n * cosI - sqrtf( 1 - sinT2)) * norm;
	*refracted = true;
	return refractedRay;
}

const float EPSILON = 0.0001f;

static bool rayPlaneIntersect(const Ray& ray, const Plane& plane, float* closestDistance)
{
	float numenator = dotVec3({plane.point - ray.origin}, plane.normal);
	
	float denom = dotVec3(ray.direction, plane.normal);

	//ray and plane normal has the same direction
	if(denom > 0)
	{
		return false;
	}

	//ray and plane are perpendicular
	if(std::abs(denom) < EPSILON)
	{
		return false;
	}

	float tmin = numenator / denom;
	//plane behind ray origin
	if(tmin < 0)
	{
		return false;
	}

	*closestDistance = tmin;

	return true;	
}

static bool raySphereIntersect(const Ray& ray, const Sphere& sphere, float* closestDistance)
{
	float a = dotVec3(ray.direction, ray.direction);
	//v - vector from sphere center to ray origin
	Vec3 v = ray.origin - sphere.position;
	float b = 2 * dotVec3(v, ray.direction);
	float c = dotVec3(v, v) - sphere.radius * sphere.radius;
	float disqSquared = b * b - 4 * a * c;
	if(disqSquared < 0)
	{
		return false;
	}

	float t1 = (-b - sqrt(disqSquared)) / (2.f * a);

	if(t1 >= EPSILON)
	{
		*closestDistance = t1;
		return true;
	}

	float t2 = (-b + sqrt(disqSquared)) / (2.f * a);
	if(t2 >= EPSILON)
	{
		*closestDistance = t2;
		return true;
	}

	return false;
}

static bool rayTriangleIntersect(const Ray& ray, const Triangle& triangle, float* closestDistance)
{
	Vec3 E1 = triangle.v1 - triangle.v0;
	Vec3 E2 = triangle.v2 - triangle.v0;

	Vec3 WxE2 = cross(ray.direction, E2);
	float det = dotVec3(WxE2, E1);
	if(det < EPSILON)
	{
		return false;
	}

	Vec3 tvec = ray.origin - triangle.v0;

	float u = dotVec3(WxE2, tvec);
	if(u < 0 || u > det)
	{
		return false;
	}

	Vec3 TxE1 = cross(tvec, E1);

	float v = dotVec3(TxE1, ray.direction);
	if(v < 0 || u + v > det)
	{
		return false;
	}
	
	float invDet = 1.f / det;
	
	*closestDistance = dotVec3(TxE1, E2) * invDet;
	u *= invDet;
	v *= invDet;
	return true;
}

static bool rayTriangleIntersectNaive(const Ray& ray, const Triangle& triangle, float* closestDistance)
{
	Vec3 v1v0 = triangle.v1 - triangle.v0;
	Vec3 v2v0 = triangle.v2 - triangle.v0;
	Vec3 normal = normaliseVec3(cross(v1v0, v2v0));

	Plane trianglePlane = {};
	trianglePlane.normal = normal;
	trianglePlane.point = triangle.v0;

	float t = {};
	if(!rayPlaneIntersect(ray, trianglePlane, &t))
	{
		return false;
	}

	Vec3 p = ray.origin + t * ray.direction;
	
	Vec3 pv0 = p - triangle.v0;
	Vec3 pv1 = p - triangle.v1;
	Vec3 pv2 = p - triangle.v2;

	Vec3 v2v1 = triangle.v2 - triangle.v1;
	Vec3 v0v2 = triangle.v0 - triangle.v2;

	if(dotVec3(normal, cross(v1v0, pv0)) < 0)
	{
		return false;
	}

	if(dotVec3(normal, cross(v2v1, pv1)) < 0)
	{
		return false;
	}

	if(dotVec3(normal, cross(v0v2, pv2)) < 0)
	{
		return false;
	}

	return true;
}

static inline Vec3 phongShade(const Ray& ray, const Surfel& surfel, const PointLight& pointLight)
{
//	Vec3 ambientTerm = hadamard(surfel.material.diffuseColor, Vec3{0.5f, 0.5f, 0.5f});
	Vec3 lightVector = normaliseVec3(surfel.position - pointLight.position);
	Vec3 diffuseTerm = surfel.material.diffuseColor * std::max<float>(std::abs(dotVec3(lightVector, surfel.normal)), 0.f);
	Vec3 reflectedRay = reflectRay(lightVector, surfel.normal);
	Vec3 viewVector = -1.f * ray.direction;
	Vec3 specularTerm = surfel.material.specularColor * std::powf(std::max(dotVec3(viewVector, reflectedRay), 0.f), surfel.material.specularGlossiness);

	diffuseTerm = pointLight.color ^ diffuseTerm;
	specularTerm =  pointLight.color ^ specularTerm;

	Vec3 color = diffuseTerm + specularTerm;

	return color;
}

static bool visible(const Scene& scene, const Vec3& from, const Vec3& to)
{
	const float epsilonOffset = 0.0001f;
	Ray shadowRay = {};
	shadowRay.direction = normaliseVec3(to - from);
	shadowRay.origin = from + epsilonOffset * shadowRay.direction;
	const float maxDistance = lengthVec3(to - from);
	
	for(uint32_t i = 0; i < scene.sphereCount; i++)
	{
		float currentDistance = {};
		if(raySphereIntersect(shadowRay, scene.spheres[i], &currentDistance) && currentDistance < maxDistance)
		{
			return false;
		}
	}//spheres count

	for(uint32_t i = 0; i < scene.planesCount; i++)
	{
		float currentDistance = {};
		if(rayPlaneIntersect(shadowRay, scene.planes[i], &currentDistance) && currentDistance < maxDistance)
		{
			return false;
		}
	}//planes count

	for(uint32_t i = 0; i < scene.triangleCount; i++)
	{
		float currentDistance = {};
		if(rayTriangleIntersect(shadowRay, scene.triangles[i], &currentDistance) && currentDistance < maxDistance)
		{
			return false;
		}
	}//triangle count

	return true;
}

static const uint32_t MAX_DEPTH = 4;

static Vec3 rayCast(const Scene& scene, const Ray& ray, int depth = 0)
{
	Vec3 color = {};
	const Vec3 background = {0.447f, 0.941f, 0.972f};
//	color = background;
	if(depth > MAX_DEPTH) return color;

	Surfel surfel = {};
	float minDistance = std::numeric_limits<float>::max();
	bool intersectionFound = false;
	int matId = -1;
	int sphereId = -1;
	int planeId = -1;
	int triangleId = -1;
	bool isShadowRay = false;
	Ray currentRay = ray;
	const float bounceCoeff = 0.001f;
	
	for(uint32_t i = 0; i < scene.sphereCount; i++)
	{
		float currentDistance = {};
		if(raySphereIntersect(currentRay, scene.spheres[i], &currentDistance))
		{
			if(currentDistance < minDistance)
			{
				minDistance = currentDistance;
				matId = scene.spheres[i].matId;
				intersectionFound = true;
				sphereId = i;
			}
		}
	}//spheres count

	for(uint32_t i = 0; i < scene.planesCount; i++)
	{
		float currentDistance = {};
		if(rayPlaneIntersect(currentRay, scene.planes[i], &currentDistance))
		{
			if(currentDistance < minDistance)
			{
				minDistance = currentDistance;
				matId = scene.planes[i].matId;
				intersectionFound = true;
				planeId = i;
			}
		}
	}//planes count

	for(uint32_t i = 0; i < scene.triangleCount; i++)
	{
		float currentDistance = {};
		if(rayTriangleIntersect(currentRay, scene.triangles[i], &currentDistance))
		{
			if(currentDistance < minDistance)
			{
				minDistance = currentDistance;
				matId = scene.triangles[i].matId;
				intersectionFound = true;
				triangleId = i;
			}
		}
	}//triangle

	if(intersectionFound)
	{
		Vec3 X = currentRay.origin + minDistance * currentRay.direction;

		if(triangleId >= 0)
		{
			Triangle closestTriangle = scene.triangles[triangleId];
			surfel.position = X;
			Vec3 edge1 = closestTriangle.v1 - closestTriangle.v0;
			Vec3 edge2 = closestTriangle.v2 - closestTriangle.v0;
			surfel.normal =	normaliseVec3(cross(edge1, edge2));
			surfel.material = scene.materials[closestTriangle.matId];
		}
		else if(planeId >= 0)
		{
			Plane closestPlane = scene.planes[planeId];
			surfel.position = X;
			surfel.normal = closestPlane.normal;
			surfel.material = scene.materials[closestPlane.matId];
		}
		else if(sphereId >= 0)
		{
			Sphere closestSphere = scene.spheres[sphereId];
			surfel.position = X;
			surfel.normal = normaliseVec3(X - closestSphere.position);
			surfel.material = scene.materials[closestSphere.matId];
		}

		switch(surfel.material.mtype)
		{
			case MTYPE_REFLECTIVE : {
				Vec3 reflectRayDir = reflectRay(currentRay.direction, surfel.normal);
				Vec3 rayOrigin = X + bounceCoeff * reflectRayDir;
				Ray newRay = {rayOrigin, reflectRayDir};
				color += 0.8f * rayCast(scene, newRay, depth + 1);
				break;
			}
			case MTYPE_DIELECTRIC: {

				bool isOutside = dotVec3(currentRay.direction, surfel.normal) < 0;
				float kr = {};
				fresnel(currentRay.direction, surfel.normal, surfel.material.ior, &kr);

				//cast reflection ray
				Vec3 reflectRayDir = reflectRay(currentRay.direction, surfel.normal);
				Vec3 rayOrigin = isOutside ? X + bounceCoeff * surfel.normal : X - bounceCoeff * surfel.normal;
				Ray newRay = {rayOrigin, reflectRayDir};
				color += kr * rayCast(scene, newRay, depth + 1);
				assert(kr > 0.f && kr < 1.f);
				//cast refracted ray
				bool isRefracted = false;
				Vec3 refractRayDir = refractRay(currentRay.direction, surfel.normal, surfel.material.ior, &isRefracted);
				if(isRefracted)
				{
					Vec3 refractRayOrigin = isOutside ? X - bounceCoeff * surfel.normal : X + bounceCoeff * surfel.normal;
					Ray newRefractRay = {refractRayOrigin, refractRayDir};
					color += (1.f - kr) * rayCast(scene, newRefractRay, depth + 1);
				}
				break;
			}
			case MTYPE_DIFFUSE : {
				
				//ambient constant
				Vec3 ambientTerm = surfel.material.diffuseColor ^ Vec3{0.5f, 0.5f, 0.5f};
				color = ambientTerm;

				for(uint32_t i = 0; i < scene.pointLightCount; i++)
				{
					// bool vis = visible(scene, surfel.position, scene.pointLights[i].position);
					// assert(vis && !"WTF IS THIS");
					if(visible(scene, surfel.position, scene.pointLights[i].position))
					{
						color += phongShade(currentRay, surfel, scene.pointLights[i]);
					}
				}

				if(surfel.material.pattern)
				{
					float x = X.x;
					float checkerSize = 0.3f;
					if((int)std::floor(x / checkerSize) % 2 == 0)
					{
						color = RGB_BLACK;
					}
				}

				break;
			}
			default : {
				assert(!"WTF");
				return color;
			}
		};
	}
	else//if no intersection were found 
	{
		color = background;
	}

	return clamp(color, RGB_BLACK, RGB_WHITE);
}

int main(int arc, char** argv)
{
	hoth::log::setSeverityMask(hoth::log::MASK_ALL);
	hoth::log::info("Hello from hoth!");

	Bitmap bitmap = {};
	if(!createBitmap(1920u, 1080u, &bitmap))
	{
		return -1;
	}

	PinholeCamera camera = {};
	//absolute distance to image plane
	camera.znear = 1.f;
	camera.vfov = toRad(90.f);
	camera.origin = {0.f, 0.5f, 1.f};
	Vec3 UpDir = {0.f, 0.5f, 0.f};
	Vec3 lookAt = {0.f, 0.f, -1.f};
	Vec3 zAxis  = normaliseVec3(camera.origin - lookAt);
	Vec3 xAxis = normaliseVec3(cross(UpDir, zAxis));
	Vec3 yAxis = cross(zAxis, xAxis);
	camera.u = xAxis;
	camera.v = yAxis;
	camera.w = zAxis;

	Film film = {};
	film.aspect = bitmap.width / (float)bitmap.height;
	film.worldHeight = 2.f * camera.znear * tanf(camera.vfov / 2.f);
	film.worldWidth = film.worldHeight * film.aspect;
	film.pixelSize = film.worldWidth / (float)bitmap.width;
	film.sampleCount = 16;

	PointLight lights[2] = {};
	lights[0].color = {0.96f, 1.f, 0.46f};
	lights[0].position = {1.f, 1.f, 1.f};

	lights[1].color = {0.6f, 0.16f, 0.16f};
	lights[1].position = {-1.f, 1.f, 1.f};
	
	Material materials[10] = {};
	materials[0].diffuseColor = {0.46f, 0.69f, 0.77f};
	materials[0].specularColor = {0.5f, 0.5f, 0.5f};
	materials[0].specularGlossiness = 256;
	materials[0].mtype = MTYPE_DIELECTRIC;
	materials[0].ior = 1.7f;

	materials[1].diffuseColor = {0.98f, 0.73f, 0.14f};
	materials[1].specularColor = {0.5f, 0.5f, 0.5f};
	materials[1].specularGlossiness = 256;
	materials[1].mtype = MTYPE_DIFFUSE;

	materials[2].diffuseColor = {0.59f, 0.47f, 0.69f};
	materials[2].specularColor = {0.1f, 0.1f, 0.1f};
	materials[2].specularGlossiness = 4;
	materials[2].mtype = MTYPE_DIFFUSE;
	materials[2].pattern = true;

	Sphere spheres[10] = {};
	spheres[0].position = {0.f, 0.5f, -2.f};
	spheres[0].radius = 1.f;
	spheres[0].matId = 0;

	spheres[1].position = {1.0f, 0.5f, -0.5f};
	spheres[1].radius = 0.6f;
	spheres[1].matId = 1;

	Plane planes[10] = {};
	planes[0].normal = {0.f, 1.f, 0.f};
	planes[0].point = {0.f, -0.8f, 0.f};
	planes[0].matId = 2;

	planes[1].normal = {0.f, 0.f, 1.f};
	planes[1].point = {0.f, 0.f, -4.f};
	planes[1].matId = 2;


	Scene scene = {};
	scene.spheres = spheres;
	scene.sphereCount = 2;
	scene.planes = planes;
	scene.planesCount = 1;
	//scene.triangles = &triangle;
	//scene.triangleCount = 1;
	scene.materials = materials;
	scene.materialCount = 3;
	scene.pointLights = lights;
	scene.pointLightCount = 1;

	uint8_t* pixels = (uint8_t*)bitmap.pixels;

	Timer timer;
	timer.start();
	
	uint32_t sampleCount = sqrt(film.sampleCount);
	
	std::random_device device;
	std::default_random_engine engine(device());
	std::uniform_real_distribution<float> distr(0.f, 1.f);
	auto randomVal = [&]() {return distr(engine);};

	for(uint32_t i = 0; i < bitmap.height; i++)
	{
		for(uint32_t j = 0; j < bitmap.width; j++)
		{
			Vec3 color = {};

			for(uint32_t v = 0; v < sampleCount; v++)
			{
				for(uint32_t u = 0; u < sampleCount; u++)
				{
					const float samplePosU = j + (u + randomVal())/(float)sampleCount;
					const float samplePosV = i + (v + randomVal())/(float)sampleCount;
					Ray ray = getPrimaryRay(samplePosU, samplePosV, bitmap.width, bitmap.height, film, camera);
					color += rayCast(scene, ray);
				}
			}

			color = color / (float)film.sampleCount;

			color.R *= 255.f;
			color.G *= 255.f;
			color.B *= 255.f;
			colorPixel(bitmap, j, i, color);
		}
	}
	hoth::log::info("Ray tracing took {:.2f}ms", timer.stopMs());

	saveAsPNGAndClose("test.png", &bitmap);
	return 0;
}