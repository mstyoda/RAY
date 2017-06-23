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

const db pi = 2. * atan2(1,0);

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
	Point det(const Point &B)const
	{
		return Point(y * B.z - B.y * z,x * B.z - B.x * z,x * B.y - B.x * y);
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
		row = img.rows; col = img.cols; rc = -350.0; cc = 500.0;
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
		mtr.Kd = Color(0.0,0.0,0.0);
		mtr.wr = Color(0.0,0.0,0.0); mtr.wt = Color(0,0,0); 
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
		Line N = getN(L,uv);
		Point Pd = L.Pd - (N.Pd * (2.0 * (L * N)));
		Pd = Pd * (1.0 /sqrt(Pd*Pd));
		return Line(get(uv),Pd);
	}
	
	inline virtual Line getTrans(Line L,UV uv)
	{
		db n; if (N * L < 0) n = mtr.n2 / mtr.n1; else n = mtr.n1 / mtr.n2;
		//L.Pd = L.Pd * (-1.0);
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
class Ball : public obj
{
public:
	Point P0; db r;
	Ball(){}
	Ball(Point p0,db R,Color WR,Color WT,Color KD)
	{
		P0 = p0; r = R; mtr.Kd = KD; 
		//mtr.wr = Color(0.999,0.999,0.999); mtr.wt = Color(0.8,0.8,0.8)*0;
		mtr.wr = WR; mtr.wt = WT;
	}
	inline virtual Point get(UV uv)
	{
		db u = uv.u,v = uv.v;
		u *= pi; v *= 2. * pi;
		return P0 + Point(r * sin(u) * cos(v),r * sin(u) * sin(v),r * cos(u));
	}
	inline virtual Line getN(Line L, UV uv)
	{
		Point P = get(uv),N; N = P - P0; N = N * (1.0 / sqrt(N * N));
		return Line(P,N);
	}
	inline virtual Line getReflect(Line L,UV uv)
	{
		Line N = getN(L,uv);
		Point Pd = L.Pd - (N.Pd * (2.0 * (L * N)));
		Pd = Pd * (1.0 /sqrt(Pd*Pd));
		return Line(get(uv),Pd);
	}
	inline virtual Line getTrans(Line L,UV uv)
	{
		Line N = getN(L,uv);
		db n; if (N * L < 0) n = 1./1.6; else {N.Pd = N.Pd * (-1); n = 1.6;}
		db cosI = -(N * L),cosT2 = 1.0 - (n * n) * (1 - cosI * cosI);
		if (cosT2 > 1e-9)
		{ 
			Point pd = (L.Pd * n) + (N.Pd * (n * cosI - sqrt(cosT2)));
			//printf("%lf %lf %lf cos = %lf\n",(pd-L.Pd)*(pd-L.Pd),L*L,N*N,cosI);
			return Line(get(uv),pd);
		}
		return getReflect(L,uv);
	}
	inline db Sqr(db x){return x * x;}

	inline virtual UV getCross(Line L)
	{
		db a,b,c,d,f,g,R; R = r * r;
		a = L.P0.x - P0.x; b = L.P0.y - P0.y; c = L.P0.z - P0.z;
		d = L.Pd.x; f = L.Pd.y; g = L.Pd.z;
		db t,delta,dfg;
		delta = Sqr(2*a*d + 2*b*f + 2*c*g) - 4*(d*d+f*f+g*g)*(a*a+b*b+c*c-R);
		dfg = (d*d + f*f + g*g);
		if ((delta > 1e-6) && (dfg > 1e-6))
		{
			db best = 1e+9;
			t = (-sqrt(delta) - 2*a*d - 2*b*f - 2*c*g) / (2 * dfg);
			if (t > 1e-6) best = min(best,t);
			t = (sqrt(delta) - 2*a*d - 2*b*f - 2*c*g) / (2 * dfg);
			if (t > 1e-6) best = min(best,t);
			if (best > 1e+9) return UV(1.2345,5.4321);
			Point C = L.P0 + L.Pd * best - P0;
			db u,v;
			u = acos(C.z); v = atan2(C.y,C.x);
			return UV(u/pi,v/(2.*pi));
		}
		return UV(1.2345,5.4321);
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
	
	inline db f(db u){return get(UV(u,0)).y;}
	inline db g(db u){return get(UV(u,0)).z;}
	inline db fd(db u){return (ay1 + 2.*ay2*u + 3.*ay3*u*u);}
	inline db gd(db u){return (az1 + 2.*az2*u + 3.*az3*u*u);}
	inline virtual Line getN(Line L,UV uv)
	{
		db u = uv.u,v = uv.v;
		db A,B,C;
			  /*d(x,y,z)/d(t,u,v)
			  (Pd.x),(-cos(v)*fd(u)),(f(u)*sin(v)),
			  (Pd.y),(-sin(v)*fd(u)),(-f(u)*cos(v)),
			  (Pd.z),(-gd(u)),(0.0);
			  */
		A = (-sin(v)*fd(u)) * ((0.0)) - (-f(u)*cos(v)) * (-gd(u)); 
		B = (-gd(u)) * (f(u)*sin(v)) - (0.0) * (-cos(v)*fd(u));
		C = (-cos(v)*fd(u)) * (-f(u)*cos(v)) - (f(u)*sin(v)) * (-sin(v)*fd(u));
		Line N; N.P0 = get(uv); N.Pd = Point(A,B,C) * (1.0 / sqrt(A*A + B*B + C*C)) * (1.);
		return N;
	}
	inline virtual Line getReflect(Line L,UV uv)
	{
		Line N = getN(L,uv);
		return Line(N.P0,L.Pd - N.Pd * (2.0 * (L * N))); 
	}
	inline UV ND(db u,db v,db t,Line L)
	{
		Eigen :: Matrix3f ma; int td = 0;
		Point P0 = L.P0, Pd = L.Pd,C;
		while (td <= 20)
		{
				Point C = (P0 + Pd * t) - get(UV(u,v));
				Point Ut = get(UV(u,v));
				if ((C*C < 1e-4) && (t > 1e-6)) return UV(u,v);
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
		return UV(1.2345,5.4321);
	}
	inline db Sqr(db x){return x * x;}
	inline UV Solve(db u1, db u2, Line L)
	{
		db z1,z2,r1,r2,k,bestk;
		Point P0 = L.P0, Pd = L.Pd,P1,P2;
		
		P1 = get(UV(u1,0)); P2 = get(UV(u2,0));
		z1 = P1.z; z2 = P2.z;
		r1 = P1.y*P1.y; r2 = P2.y*P2.y; 
		std::vector<db> K; K.clear();
		if (Pd.z)
		{
			k = (z1 - P0.z) / Pd.z; 
			if ((k > 1e-6) ) K.push_back(k);
			k = (z2 - P0.z) / Pd.z;
			if ((k > 1e-6)) K.push_back(k);
		}
		//和柱面求交点，(P0.x + k * Pd.x)^2 + (P0.y + k * Pd.y)^2 = r^2
		//令a = P0.x,b = Pd.x,c = P0.y,d = Pd.y,变为(a+b*k)^2+(c+d*k)^2 = r2
		db a = P0.x,b = Pd.x,c = P0.y,d = Pd.y;
		if (b * b + d * d > 1e-6)
		{
			k = (-sqrt(-a*a*d*d + 2.*a*b*c*d + b*b*(-c*c) + b*b*r1 + d*d*r1) - a*b - c*d)/(b*b + d*d);
			if (k > 1e-6) K.push_back(k);
			k = (sqrt(-a*a*d*d + 2.*a*b*c*d + b*b*(-c*c) + b*b*r1 + d*d*r1) - a*b - c*d)/(b*b + d*d);
			if (k > 1e-6) K.push_back(k);
		}
		if (b * b + d * d > 1e-6)
		{
			k = (-sqrt(-a*a*d*d + 2.*a*b*c*d + b*b*(-c*c) + b*b*r2 + d*d*r2) - a*b - c*d)/(b*b + d*d);
			if (k > 1e-6) K.push_back(k);
			k = (sqrt(-a*a*d*d + 2.*a*b*c*d + b*b*(-c*c) + b*b*r2 + d*d*r2) - a*b - c*d)/(b*b + d*d);
			if (k > 1e-6) K.push_back(k);
		}

		UV INF = UV(1.2345,5.4321),ans,cur; ans = INF;
		if (!K.size()) return INF;
		if (fabs(u2 - u1) < 1e-2)
		{
			int i; rep(i,0,K.size()-1)
			{
				k = K[i]; Point C = P0 + Pd * k; db phi;
				phi = atan2(C.y,C.x);
				if (phi < 0) phi += 2. * pi;
				ans = ND(u1,phi/(4.0 * atan2(1,0)),k,L);
				if (fabs(ans.u - INF.u) + fabs(ans.v - INF.v) > 1e-5) return ans;
			}
			return INF;
		}
		else
		{
			UV ANS; ANS = INF; db best = 1e+9;
			ans = Solve(u1,(u1 + u2) * 0.5,L);
			if (fabs(ans.u - INF.u) + fabs(ans.v - INF.v) > 1e-5)
			{
				Point cur = (get(ans) - L.P0); 
				if (cur * cur < best){ ANS = ans; best = cur * cur;}
			}
			ans = Solve((u1 + u2) * 0.5,u2,L);
			if (fabs(ans.u - INF.u) + fabs(ans.v - INF.v) > 1e-5)
			{
				Point cur = (get(ans) - L.P0); 
				if (cur * cur < best){ ANS = ans; best = cur * cur;}
			}
			return ANS;
		}
	}
	inline virtual UV getCross(Line L)
	{
		return Solve(0.,1.,L);
	}
};

class B1 : public Bobj
{
public:
	B1()
	{
		P[3] = Point(0.0,0.0,0.0); P[2] = Point(0.0,0.84,0.39); 
		P[1] = Point(0.0,0.76,0.67); P[0] = Point(0.0,1.0,1.0);
		n = 3;
		mtr.wr = Color(0.0,0.0,0.0); mtr.wt = Color(0.,0.,0.);
		mtr.Kd = Color(0.75,0.25,0.25) * 0.0;
		int i; rep(i,0,3) P[i] = P[i] * 0.35;
		ay3 = 0.266; ay2 = -0.714; ay1 = 0.798; ay0 = 0;
		az3 = 0.644; az2 = -0.9975; az1 = 0.7035; az0 = 0;

	}
};


class B2 : public Bobj
{
public:
	B2()
	{
		P[0] = Point(0.0,0,0); P[1] = Point(0,.99,.01); 
		P[2] = Point(0,.93,.59); P[3] = Point(0,1.,1.);
		n = 3;
		mtr.wm = Color(1.,1.,1.); mtr.wr = Color(0.,0.,0.); mtr.wt = Color(0.,0.,0.);
		mtr.Kd = Color(0.25,0.25,0.75); mtr.ad = 100;
		int i; rep(i,0,3) P[i] = P[i] * 0.35;
	}
};

