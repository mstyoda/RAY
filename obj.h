#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <set>
#include <vector>
#include <opencv2/opencv.hpp>

#define rep(i,l,r) for (i = l; i <= r; i++)
#define drep(i,l,r) for (i = l; i >= r; i--)
using namespace std;
typedef double db;

class UV
{
public:
	db u,v;
	UV(){}
	UV(db U,db V){u = U; v = V;}
};

class Point
{
public:
	db x,y,z;
	Point(){x = y = z = 0;}
	Point(db X, db Y, db Z){x = X; y = Y; z = Z;}
	db len(){return sqrt(x * x + y * y + z * z);}
	Point operator * (const db & k) const
	{
		return Point(x * k,y * k,z * k);
	}
	db operator * (const Point &B) const
	{
		return x * B.x + y * B.y + z * B.z;
	}
	Point operator + (const Point & B) const
	{
		return Point(x + B.x,y + B.y,z + B.z);
	}
	Point operator - (const Point &B) const
	{
		return Point(x - B.x,y - B.y,z - B.z);
	}
};

class Line
{
public:
	Point P0,Pd;
	Line(){}
	Line(Point p0,Point pd){P0 = p0; Pd = pd;}
	db operator * (const Line & B) const
	{
		return Pd * B.Pd;
	}
};

class Color
{
public:
	db r,g,b;
	Color(){}
	Color(db R,db G,db B){r = R; g = G; b = B;}
	Color operator * (const db & k)const
	{
		return Color(k * r,k * g,k * b);
	}
	Color operator * (const Color & B)const
	{
		return Color(r * B.r,g * B.g,b * B.b);
	}
	Color operator + (const Color & B)const
	{
		return Color(r + B.r,g + B.g,b + B.b);
	}
	db len(){return r * r + g * g + b * b;}
};

class Cp
{
public:
	Point P; int i,j,k; Line N,V; Color w,I;
};
class Mtr
{
public:
	int ar,at,ad;
	Color wr,wt,wm,Ka,Kd,Ks;
	db n1,n2;
};
class obj
{
public:
	Mtr mtr;db u1,u2,v1,v2;
	inline virtual Point get(UV uv)
	{

	}

	inline virtual Line getReflect(Line L,UV uv)
	{

	}

	inline virtual Line getN(Line L,UV uv)
	{

	}

	inline virtual Line getTrans(Line L,UV uv)
	{

	}

	inline virtual UV getCross(Line L)
	{
		
	}
};

class Plane : public obj
{
public:
	Line N; Point A,B;
	inline Plane(){}
	inline Plane(Line n,Point a,Point b)
	{
		N = n; A = a; B = b; mtr.n1 = 1.0; mtr.n2 = 1.6; 
		mtr.ar = 20; mtr.at = 0; mtr.ad = 80; 
		mtr.Kd = Color(0.75,0,25.25); mtr.Ka = Color(0.0,0.0,0.0); mtr.Ks = Color(0.0,0.0,0.0);
		mtr.wr = Color(0.0,0.0,0.0); mtr.wt = Color(0,0,0); mtr.wm = Color(0.75,0.25,0.25);
		u1 = -2; u2 = 2; v1 = -1.5; v2 = 2;
	}
	inline virtual Point get(UV uv)
	{
		return N.P0 + A * uv.u + B * uv.v;
	}

	inline virtual Line getReflect(Line L,UV uv)
	{
		return Line(get(uv),L.Pd - N.Pd * (2.0 * (L * N))); 
	}

	inline virtual Line getN(Line L,UV uv)
	{
		Point P = get(uv);
		return Line(P,N.Pd); 
	}

	inline virtual Line getTrans(Line L,UV uv)
	{
		db n; if (N * L < 0) n = mtr.n2 / mtr.n1; else n = mtr.n1 / mtr.n2;
		db cosI = -(N * L),cosT2 = 1.0 - (n * n) * (1 - cosI * cosI);
		if (cosT2 > 1e-9) return Line(get(uv),L.Pd * n + N.Pd * (n * cosI - sqrt(cosT2)));
		return getReflect(L,uv);  
	}

	inline virtual UV getCross(Line L)
	{
		db Z0 = (N.P0 - L.P0) * N.Pd;
		db k = -1,u,v;
		u = 1.2345; v = 5.4321;
		if (fabs(Z0) > 1e-9) k = Z0 / (L.Pd * N.Pd);
		if (k > 1e-9)
		{
			Point C = L.P0 + (L.Pd * k) - (N.P0),T1 = Point(11,3,7),T2 = Point(123,2,13);
			db a,b,c,d,e,f;
			a = A * T1; b = B * T1; c = C * T1;
			d = A * T2; e = B * T2; f = C * T2; 
			u = (b * f - c * e) / (b * d - a * e); v = (c * d - a * f) / (b * d - a * e);
		}
		if ((u < u1) || (u > u2) || (v < v1) || (v > v2)) {u = 1.2345; v = 5.4321;}																																																																																																																																																																																																																																																																																																																																																																																																																			
		//printf("Line (%lf,%lf,%lf) [%lf %lf %lf]\n",L.P0.x,L.P0.y,L.P0.z,L.Pd.x,L.Pd.y,L.Pd.z);
		return UV(u,v);
	}

};
class Bobj : public obj
{
public:
	Point P[5];
	int n;
	inline Point get(db u, db v)
	{
		int i,j; db pu[10],tu[10],F[10]; F[0] = pu[0] = tu[0] = 1.0;
		rep(i,1,n){pu[i] = pu[i - 1] * u; tu[i] = tu[i - 1] * (1.0 - u); F[i] = F[i - 1] * (db)i;}
		Point B; 
		rep(i,0,n)
		{
			B = B + P[i] * (F[n] / (F[i] * F[n - i]) * pu[i] * tu[n - i]);
		}
		v = v * 4.0 * atan2(1,0);
		return Point(B.y * cos(v),B.y * sin(v),B.z);
	}
};

class bowlout : public Bobj
{
public:
	bowlout()
	{
		P[0] = Point(0.0,0.0,0.0); P[1] = Point(0.0,1.1,0.0); 
		P[2] = Point(0.0,0.65,0.78); P[3] = Point(0.0,1.0,0.75);
		n = 3;
	}
};


class bowlin : public Bobj
{
public:
	bowlin()
	{
		P[3] = Point(0.0,0.0,0.0); P[2] = Point(0.0,1.1,0.3); 
		P[1] = Point(0.0,0.45,0.78); P[0] = Point(0.0,1.0,0.75);
		n = 3;
	}
};

class bowlbt2 : public Bobj
{
public:
	bowlbt2()
	{
		P[3] = Point(0.0,376.0,-221.0); P[1] = Point(0.0,499,-220); 
		P[2] = Point(0.0,377,-294); P[0] = Point(0.0,475,-319);
		int i; rep(i,0,3) P[i] = P[i] * (1.0/1150.0) + Point(0,0,0.09 + 0.03);
		n = 3;
	}
};

class bowlbt1 : public Bobj
{
public:
	bowlbt1()
	{
		P[3] = Point(0.0,361.0,-203.0); P[1] = Point(0.0,454,-279); 
		P[2] = Point(0.0,457,-202); P[0] = Point(0.0,364,-283);
		int i; rep(i,0,3) P[i] = P[i] * (1.0/1300.0) + Point(0,0,0.1 + 0.03);
		n = 3;
	}
};
class bowlbt0 : public Bobj
{
public:
	bowlbt0()
	{
		P[3] = Point(0.0,353.0,-240.0); P[1] = Point(0.0,353,-357); 
		P[2] = Point(0.0,334,-281); P[0] = Point(0.0,442,-399);
		int i; rep(i,0,3) P[i] = P[i] * (1.0/1300.0) + Point(0,0,0.175 + 0.03);
		n = 3;
	}
};
