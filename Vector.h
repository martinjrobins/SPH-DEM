#include <iostream>
#include <math.h>

class Vector2D {
   public:
	real x,y;				// vector components
	
	Vector2D() { x=y=0; }			// constructor
	Vector2D(const real pX, const real pY)
	: x(pX), y(pY) {}

	// vector arithmatic

	Vector2D operator+(const Vector2D pV) const
	{ Vector2D out( x+pV.x, y+pV.y); return out; }

	Vector2D operator-(const Vector2D pV) const
	{ Vector2D out( x-pV.x, y-pV.y); return out; }

	Vector2D& operator+=(const Vector2D pV)
	{ x+=pV.x; y+=pV.y; return *this; }

	Vector2D& operator-=(const Vector2D pV)
	{ x-=pV.x; y-=pV.y; return *this; }

	real dot(const Vector2D pV) const			// dot product
	{ return x*pV.x + y*pV.y; }

	Vector2D operator*(const real pR) const		// * a scalar
	{ Vector2D out( x*pR, y*pR); return out; }

	friend Vector2D operator*(const real pR, const Vector2D pV)
	{ Vector2D out( pV.x*pR, pV.y*pR); return out; }

	Vector2D& operator*=(const real pR)
	{ x*=pR; y*=pR; return *this; }

	Vector2D& operator=(const real pR)
	{ x=pR; y=pR; return *this; }
	// magnitude

	friend real len2(const Vector2D pV)
	{ return pV.x*pV.x+pV.y*pV.y;}

	friend real len(const Vector2D pV) 
	{ return sqrt(len2(pV)); }

	// comparison
	int operator==( const Vector2D pV) const
	{ return x==pV.x && y==pV.y; }

	int operator!=( const Vector2D pV) const
	{ return x!=pV.x || y!=pV.y; }

	int operator<=( const Vector2D pV) const
	{ return x<=pV.x && y<=pV.y; }
	int operator>=( const Vector2D pV) const
	{ return x>=pV.x && y>=pV.y; }
	int operator<( const Vector2D pV) const
	{ return x<pV.x && y<pV.y; }
	int operator>( const Vector2D pV) const
	{ return x>pV.x && y>pV.y; }

	// input/output

	//friend ostream& operator<<(ostream& pStr, const Vector2D& pV)
	//{ return (pStr << '(' << pV.x << ',' << pV.y << ')'); }

	//friend istream& operator>>(istream pStr, Vector2D& pV)
	//{ return (pStr >> pV.x >> pV.y); }
};


class Vector3D {
   public:
	real x,y,z;				// vector components
	
	Vector3D() { x=y=z=0; }			// constructor
	Vector3D(const real pX, const real pY, const real pZ)
	: x(pX), y(pY), z(pZ) {}

	// vector arithmatic

	Vector3D operator+(const Vector3D pV) const
	{ Vector3D out( x+pV.x, y+pV.y, z+pV.z ); return out; }

	Vector3D operator-(const Vector3D pV) const
	{ Vector3D out( x-pV.x, y-pV.y, z-pV.z ); return out; }

	Vector3D& operator+=(const Vector3D pV)
	{ x+=pV.x; y+=pV.y; z+=pV.z; return *this; }

	Vector3D& operator-=(const Vector3D pV)
	{ x-=pV.x; y-=pV.y; z-=pV.z; return *this; }

	real dot(const Vector3D pV) const			// dot product
	{ return x*pV.x + y*pV.y + z*pV.z; }

	Vector3D operator*(const real pR) const		// * a scalar
	{ Vector3D out( x*pR, y*pR, z*pR ); return out; }

	friend Vector3D operator*(const real pR, const Vector3D pV)
	{ Vector3D out( pV.x*pR, pV.y*pR, pV.z*pR ); return out; }

	Vector3D& operator*=(const real pR)
	{ x*=pR; y*=pR; z*=pR; return *this; }

	Vector3D& operator=(const real pR)
	{ x=pR; y=pR; z=pR; return *this; }

	// magnitude

	friend real len2(const Vector3D pV)
	{ return pV.x*pV.x+pV.y*pV.y+pV.z*pV.z; }

	friend real len(const Vector3D pV) 
	{ return sqrt(len2(pV)); }

	// comparison
	int operator==( const Vector3D pV) const
	{ return x==pV.x && y==pV.y && z==pV.z; }

	int operator!=( const Vector3D pV) const
	{ return x!=pV.x || y!=pV.y || z!=pV.z; }

   int operator<=( const Vector3D pV) const
	{ return x<=pV.x && y<=pV.y && z<=pV.z; }
	int operator>=( const Vector3D pV) const
	{ return x>=pV.x && y>=pV.y && z>=pV.z; }
	int operator<( const Vector3D pV) const
	{ return x<pV.x && y<pV.y && z<pV.z; }
	int operator>( const Vector3D pV) const
	{ return x>pV.x && y>pV.y && z>pV.z; }


	// input/output

	//friend ostream& operator<<(ostream& pStr, const Vector3D& pV)
	//{ return (pStr << '(' << pV.x << ',' << pV.y << ',' << pV.z << ')'); }

	//friend istream& operator>>(istream pStr, Vector3D& pV)
	//{ return (pStr >> pV.x >> pV.y >> pV.z); }
};


