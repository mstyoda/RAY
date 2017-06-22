#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <set>
#include <vector>
#include <opencv2/opencv.hpp>
#include "Eigen/Core"
#include "Eigen/Dense"

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
	Point P; int i,j,k; Line N,V; Color w,I; UV uv;
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
	Mtr mtr;db u1,u2,v1,v2; int row,col; db rc,cc;
	Color F[600][600];
	inline void Load(char *fname,db Rc,db Cc)
	{
		cv :: Mat img = cv :: imread(fname);
		db r,g,b; rc = Rc; Cc = cc;
		int i,j; uchar *o; rep(i,0,img.rows)
		{
			o = img.ptr<uchar>(i);
			rep(j,0,img.cols)
			{
				b = (db)(*o++); g = (db)(*o++); r = (db)(*o++);
				F[i][j] = Color(r,g,b) * (1.0 / 255.);
			}
		}
		mtr.Kd = Color(0,0,0); row = img.rows; col = img.cols; rc = -350.0; cc = 500.0;
	}
	inline Color getKd(UV uv)
	{
		Color Ans; int u,v;
		u = (int)((uv.u - u1) * (db)rc); v = (int)((uv.v - v1) * (db)cc);
		if (row){ u %= row; v %= col; }
		if (u < 0) u += row; if (v < 0) v += col;
		return mtr.Kd + F[u][v];
	}
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
		int i,j; rep(i,0,row) rep(j,0,col)F[i][j] = Color(0.,0.,0.);
	}
	inline virtual Point get(UV uv)
	{
		return N.P0 + A * uv.u + B * uv.v;
	}

	inline virtual Line getN(Line L,UV uv)
	{
		Point P = get(uv);
		return Line(P,N.Pd); 
	}
	
	inline virtual Line getReflect(Line L,UV uv)
	{
		return Line(get(uv),L.Pd - N.Pd * (2.0 * (L * N))); 
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
		if ((u < u1) || (u > u2) || (v < v1) || (v > v2)) {u = 1.2345; v=5.4321;}	
		//printf("Line (%lf,%lf,%lf) [%lf %lf %lf]\n",L.P0.x,L.P0.y,L.P0.z,L.Pd.x,L.Pd.y,L.Pd.z);
		return UV(u,v);
	}

};
class Bobj : public obj
{
public:
	Point P[5];
	int n; db ay0,ay1,ay2,ay3,az0,az1,az2,az3;
	inline virtual Point get(UV uv)
	{
		db u = uv.u,v = uv.v;
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
	
	inline db f(db u){return ay0 + u*(ay1 + u*(ay2 + u*(ay3)));}
	inline db g(db u){return az0 + u*(az1 + u*(az2 + u*(az3)));}
	inline db fd(db u){return ay1 + 2.*ay2*u + 3.*ay3*u*u;}
	inline db gd(db u){return az1 + 2.*az2*u + 3.*az3*u*u;}
	inline virtual Line getN(Line L,UV uv)
	{
		
	}
	inline virtual UV getCross(Line L)
	{
		//printf("L = (%lf,%lf,%lf) [%lf %lf %lf]\n",L.P0.x,L.P0.y,L.P0.z,L.Pd.x,L.Pd.y,L.Pd.z);
		//求解射线L和包围圆柱的交
		int i,j; db z1,z2,r2,k,bestk; r2 = 0.; bestk = 1e+9;
		rep(i,0,n) r2 = max(r2,(P[i].x * P[i].x + P[i].y * P[i].y));
		z2 = -1e+9; z1 = -z1; rep(i,0,n){z2 = max(z2,P[i].z); z1 = min(z1,P[i].z);}
		Point P0 = L.P0, Pd = L.Pd,C;
		std::vector<db> K; K.clear();
		if (Pd.z)
		{
			k = (z1 - P0.z) / Pd.z; if (k > 1e-6) K.push_back(k);
			k = (z2 - P0.z) / Pd.z; if (k > 1e-6) K.push_back(k);
		}
		//和柱面求交点，(P0.x + k * Pd.x)^2 + (P0.y + k * Pd.y)^2 = r^2
		//令a = P0.x,b = Pd.x,c = P0.y,d = Pd.y,变为(a+b*k)^2+(c+d*k)^2 = r2
		db a = P0.x,b = Pd.x,c = P0.y,d = Pd.y;
		if (b * b + d * d > 1e-6)
		{
			k = (-sqrt(-a*a*d*d + 2.*a*b*c*d + b*b*(-c*c) + b*b*r2 + d*d*r2) - a*b - c*d)/(b*b + d*d);
			if (k > 1e-6) K.push_back(k);
			//bestk = k > 1e-6 ? min(k,bestk) : bestk;
		}

		//牛顿迭代
		for (int l = 0; l < (int)K.size(); l++)
		{
			k = K[l]; C = P0 + Pd * k;
			printf("z1 = %lf z2 = %lf C.z = %lf\n",z1,z2,C.z);
			if ((C.z > z1 - 0.00001) && (C.z < z2 + 0.00001))
			{
				db u,v,t; int td = 0; t = k; u = 0.5 * (C.z - P[0].z) / (P[3].z - P[0].z); v = acos(C.x/sqrt(C.x*C.x+C.y*C.y));
				v = v / (4.0 * atan2(1,0));
				
				printf("u = %lf v = %lf t = %lf\n",u,v,t); v = 0.0;
				Eigen :: Matrix3f ma; td = 0;
				while (td <= 100)
				{
						C = (P0 + Pd * t) - get(UV(u,v));
						Point Ut = get(UV(u,v));
						printf("len = %.2lf C = (%.2lf,%.2lf,%.2lf)\n",sqrt(C*C),Ut.x,Ut.y,Ut.z);
						if (C*C < 1e-5) return UV(u,v);
						ma << (Pd.x),(-cos(v)*fd(u)),(f(u)*sin(v)),
							  (Pd.y),(-sin(v)*fd(u)),(-f(u)*cos(v)),
							  (Pd.z),(-gd(u)),(0.0);
						ma = ma.inverse().eval();
						C = (P0 + Pd * t) - get(UV(u,v));
						Eigen :: Vector3f V1(C.x,C.y,C.z);
						V1 = ma * V1;
						t = -V1[0] + t; u = -V1[1] + u; v = -V1[2] + v;
						td++;
				}
			}
		}
		printf("NO Solution\n");
		return UV(1.2345,5.4321);
	}
};

class bowlout : public Bobj
{
public:
	bowlout()
	{
		//y : 1 * 0 *  (1-u)^3 + 3 * 1.1 * u * (1-u)^2 + 3 *  0.65 * u^2 * (1-u) + 1 * 1 * u^3
		//z : 1 * 0 *  (1-u)^3 + 3 * 0 * u * (1-u)^2 + 3 *  0.78 * u^2 * (1-u) + 1 * 0.75 * u^3
		P[0] = Point(0.0,0.0,0.0); P[1] = Point(0.0,1.1,0.0); 
		P[2] = Point(0.0,0.65,0.78); P[3] = Point(0.0,1.0,0.75);
		ay0 = 0.; ay1 = 3.3; ay2 = 4.65; ay3 = 2.35;
		az0 = 0.; az1 = 0.0; az2 = 2.34; az3 = 1.59;
		n = 3;
		mtr.wm = Color(1.,1.,1.); mtr.wr = Color(0.,0.,0.); mtr.wt = Color(0.,0.,0.);
		mtr.Kd = Color(0.25,0.25,0.75);
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
