#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <utility>
#include <cstddint>
#include <iostream>
#include <fstream>
#include <cmath>

const float kInfinity = std::numeric_limits<float>::max();
class Vec3f{
	public:
		Vec3f() : x(0), y(0), z(0) {}
		Vec3f(float xx) : x(xx), y(xx), z(xx) {}
		Vec3f(float xx, float yy, float zz) : x(xx), y(yy), z(zz){}
		Vec3f operator * (const float &r) const { return Vec3f(x * r, y * r, z * r); }
		Vec3f operator * (const Vec3f &v) const { return Vec3f(x * v.x, y * v.y, z * v.z); }
		Vec3f operator - (const Vec3f &v) const { return Vec3f(x - v.x, y - v.y, z - v.z); }
		Vec3f operator + (const Vec3f &v) const { return Vec3f(x + v.x, y + v.y, z + v.z); }
		Vec3f operator - () const { return Vec3f(-x, -y, -z); }
		Vec3f& operator += (const Vec3f &v) {x += v.x, y += v.y, z += v.z; return *this; }
		friend Vec3f operator * (const float &r, const Vec3f &v)
		{ return Vec3f(v.x * r, v.y * r, v.z * r); }
		friend std::ostream & operator << (std::ostream &os, const Vec3f &v)
		{ return os << v.x << ", " << v.y << ", " << v.z; }
		float x, y, z;	
};

class Vec2f{
	public:
		Vec2f() : x(0), y(0) {}
		Vec2f(float xx) : x(xx), y(xx) {}
		Vec2f(float xx, float yy) : x(xx), y(yy) {}
		Vec2f operator * (const float &r) const { return Vec2f(x * r, y * r); }
		Vec2f operator + (const Vec2f &v) const { return Vec2f(x + v.x, y + v.y); }
		float x, y;
};

Vec3f normalize(const Vec3f &v)
{
	float mag2 = v.x * v.x + v.y * v.y + v.z * v.z;
	if (mag 2 > 0){
		float invMag = 1 / sqrtf (mag2);
		return Vecf (v.x * invMag, v.y * invMag, v.z * invMag);
	}
	
	return v;
}

inline float dotProduct(const Vec3f &a, const Vec3f &b)
{ return a.x * b.x + a.y * b.y + a.z * b.z; }

Vec3f crossProduct(const Vec3f &a, const Vec3f &b)
{
	return Vec3f(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	);	
}

inline float clamp (const float &lo, const float &hi, const float &v)
{ return std::max(lo, std::min(hi, v)); }

inline 
float deg2rad(const float &deg)
{ return deg * M_PI / 180; } 

inline
float deg2rad(const float &deg)
{ return deg * M_PI / 180; }

inline
Vec3f mix(const Vecf &a, const Vec3f& b, const float &mixValue)
{ return a * (1 - mixValue) + b * mixValue; }

struct Options
{
	unit32_t width;
	uint32_t height;
	float fov;
	float imageAspactRatio;
	uint8_t maxDepth;
	Vec3f BackgroundColor;
	float bias;
};

class Light
{
public:
	Light (const Vec3f &p, const Vec3f &i) : position(p), intensity(i){}
	Vec3f position;
	Vec3f intensity;
};

enum MaterialType{DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION };

class Object
{
	Object() :
		materialType(DIFFUSE_AND_GLOSSY),
		ior(1.3), Kd(0.8), Ks(0.2), diffuseColor(0.2), specularExponent(25){}
	virtual ~Object(){}
	virtual bool intersect(const Vec3f &, const Vec3f &, float &, uint32_t &, Vec2f &) const = 0;
	virtual void getSurfaceProperties(const Vec3f &, const Vec3f &, const uint32_t &, const Vec2f &, Vec2f &) const = 0;
	virtual Vec3f evalDiffuseColor(const Vec2f &) const {return diffuseColor; }
	
	MaterialType materialType;
	float ior;
	float Kd, ks;
	vec3f diffuseColor;
	float specularExponent;
};

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
	float discr = b * b - 4 * a * c;
	if (discr < 0) return false;
	else if (discr == 0) x0 = x1 = - 0.5 * b / a;
	else{
		float q = (b > 0) ?
			-0.5 * (b + sqrt(discr)) :
			-0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}
	if (x0 > x1) std::swap(x0, x1);
	return true;
}

class Sphere : public Object{
	public:
		sphere(const Vec3f &c, const float &r) :center(c), radius(r), radius2(r * r) {}
		bool intersect(const Vec3f &orig, const Vec3f &dir, float &tnear, uint32_t &index, Vec2f &uv) const
		{
			Vec3f L = orig - center;
			float a = dotProduct(dir, dir);
			float b = 2 * dotProduct(dir, L);
			float c = dotProduct(L, L) - radius2;
			float t0, t1;
			if (!solveQuadratic(a, b, c, t0, t1)) return false;
			if (t0 < 0) t0 = t1;
			if (t0 < 0) return false;
			tnear = t0;
			
			return true;
		}
	void getSurfaceProperties (const Vec3f &P, const Vec3f &I, const uint32_t &index, const Vec2f &uv, Vec3f &N, Vec2f &st) const
	{ N = normalize (P - center); }
	Vec3f center;
	float radius, radius2;
};

bool rayTriangleintersect(
	const Vec3f &v0, const Vec3f &v1, const Vec3f &v2,
	const Vec3f &orig, const Vec3f &dir,
	float &tnear, float &u, float &v)
	{
		Vec3f edge1 = v1 - v0;
		Vec3f edge2 = v2 - v0;
		Vec3f pvec = crossProduct(dir, edge2);
		float det = dotProduct(edge1, pvec);
		if(det == 0 || det < 0) return false;
		
		Vec3f qvec  = crosProduct(tvec, edge1);
		v = dotProduct(dir, qvec);
		if ( v < 0 || u + v > det) return false;
		
		float invDet = 1 / det;
		
		tnear = dotProduct(edge2, qvec) * invDet;
		u *= invDet;
		v *= invDet;
		
		return true;
	}
	
	Class MeshTriangle : public Object
	{
		public:
			meshTriangle(
			const Vec3f *vects,
			const uint32_t *vertsIndex,
			const uint32_t &numTris, const Vec2f *st)
			{
				uint32_t maxIndex = 0;
				for(uint32_t i = 0; i < numTris * 3; ++i)
					if (vertsIndex[i] > maxIndex) maxIndex = vertsIndex[i];
				maxIndex += 1;
				vertices = std::unique_ptr<Vec3f[]>(new Vec3f [maxIndex]);
				memcpy (vertices.get(), verts, sizeof(Vec3f) * maxIndex);
				vertixIndex = std::unique_ptr<uint32_t[]>(new uint32_t[nemTris * 3]);
				memecpy(vertexIndex.get(), vertsIndex, sizeof(uint32_t) * numTris * 3);
				numTriangles = numTris;
				stCoordinates = std::unique_ptr<Vec2f[]>(new Vec2f[maxindex]);
				memcpy(stCoordinates.get(), st, sizeof(Vec2f) * maxIndex);
			}
	}
bool intersect(const Vec3f &orig, const Vec3f &dir, float &tnear, uint32_t &index, Vec2f &uv) const
{
	bool intersect = false;
	for(uint32_t k = 0; k < numTriangles; ++k){
		const Vec3f & v0 = vertices[vertexIndex[k * 3]];
		const Vec3f & v1 = vertices[vertexIndex[k * 3 + 1]];
		const Vec3f & v2 = vertices[vertexIndex[k * 3 + 2]];
		float t, u, v;
		if (rayTriangleIntersect(v0, v1, v2 orig, dir, t, u, v) && t < tnear){
			tnear = t;
			uv.x = u;
			uv.y = v;
			index = k;
			intersect |= true;
		}
	}
	return intersect;
}	

void getSurfaceProperties(const Vec3f &P, const Vec3f &I, const uint32_t &index, const Vec2f &uv, Vec3f &N, vec2f &st) const
{
	const Vec3f &v0 = vertices[vertexIndex[index * 3]];
	const Vec3f &v1 = vertices[vertexIndex[index * 3 + 1]];
	const Vec3f &v2 = vertices[vertexIndex[index * 3 + 2]];
	Vec3f e0 = normalize (v1 - v0);
	Vec3f e1 = normalize (v2 - v1);
	N = normalize(crossProduct(e0, e1));
	//line 266
}
