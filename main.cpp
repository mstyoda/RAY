#include "obj.h"

#define pb push_back

const int N = 3000000,W = 800,H = 600;
const db pi = 2.0 * atan2(1,0),rds = 0.17;

vector <obj * > OBJ;
vector <Point> Light;
int tCP,tot;
Cp CP[N];

Color IMG[H + 10][W + 10];
db cnt[H + 10][W + 10];
Point Eye,InfP;

struct arr
{
	int k,l,r,L,R; Cp now; db rd; Color flux;
}Tree[N];

inline int sn(db x){return x < -1e-9 ? -1 : x > 1e-9 ? 1 : 0;}
inline bool cmpx(const Cp & A,const Cp & B){ return sn(A.P.x - B.P.x) < 0;}
inline bool cmpy(const Cp & A,const Cp & B){ return sn(A.P.y - B.P.y) < 0;}
inline bool cmpz(const Cp & A,const Cp & B){ return sn(A.P.z - B.P.z) < 0;}
inline void Init()
{
	srand(time(0)); Light.clear(); 
	Light.pb(Point(1,1,1)); OBJ.clear(); 
	OBJ.pb(new Plane(Line(Point(0,0,0),Point(0,0,1)),Point(1,0,0),Point(0,1,0)));
	OBJ.pb(new Plane(Line(Point(0,0,2),Point(0,0,-1)),Point(1,0,0),Point(0,1,0)));
	OBJ.pb(new Plane(Line(Point(0,-1.5,0),Point(0,1,0)),Point(1,0,0),Point(0,0,1)));
	OBJ.pb(new Plane(Line(Point(0,1.5,0),Point(0,-1,0)),Point(1,0,0),Point(0,0,1)));
	OBJ.pb(new Plane(Line(Point(2,0,0),Point(-1,0,0)),Point(0,1,0),Point(0,0,1)));
	
	Mtr mtr;
	mtr.ar = 80; mtr.at = 0; mtr.ad = 10; 
	mtr.Kd = Color(0.0,0,0.0); mtr.Ka = Color(0.1,0.1,0.1); mtr.Ks = Color(0.5,0.5,0.5);
	mtr.wr = Color(0.9,0.9,0.9); mtr.wt = Color(0,0,0); mtr.wm = Color(0.1,0.1,0.1);
	OBJ[1]->mtr = mtr;
	OBJ[0]->v1 = OBJ[0]->u1 = -2; OBJ[0]->v2 = OBJ[0]->u2 = 2;
	OBJ[2]->mtr.Kd = Color(0.25,0.25,0.75);
	OBJ[3]->mtr.Kd = Color(0.25,0.75,0.25);
	OBJ[4]->mtr.Kd = Color(0.5,0.1,0.5);
	tCP = 0; tot = 0;
	int i,j; rep(i,0,H - 1) rep(j,0,W - 1){IMG[i][j] = Color(0,0,0); cnt[i][j] = 0.;}
}
inline db len(const Point &A,const Point &B) {return sqrt((A - B) * (A - B));}
inline db len(const UV &A,const UV &B) {return fabs(A.u - B.u) + fabs(A.v - B.v);}
inline int toInt(double x)
{
	return int(pow(1-exp(-x),1/2.2)*255+.5);
} 
inline void getCp(Line L,int i,int j,Color w,int deep)
{
	if (sn(w.len()) && (deep < 20))
	{
		int k,ck; db dr = 1e+10,cur; UV bestuv; ck = -1;
		rep(k,0,(int)OBJ.size() - 1)
		{
			UV uv = OBJ[k]->getCross(L),Inf = UV(1.2345,5.4321);
			Point C = OBJ[k]->get(uv);
			cur = len(C,L.P0);
			if ((sn(cur -  dr) < 0) && (sn(len(uv,Inf))))
			{
				dr = cur; bestuv = uv;  ck = k;
			}
		}
		if (ck != -1)//has cross
		{
			Color wr = OBJ[ck]->mtr.wr;
			if (sn(wr.len())) getCp(OBJ[ck]->getReflect(L,bestuv),i,j,w * wr,deep + 1);
			//Color wt = OBJ[ck]->mtr.wt;
			//if (sn(wt.len())) getCp(OBJ[ck]->getTrans(L,bestuv),i,j,w * wt,deep + 1);
			
			Color wm = OBJ[ck]->mtr.wm;
			Cp cp; cp.P = OBJ[ck]->get(bestuv); cp.i = i; cp.j = j; 
			cp.k = ck; cp.V = L; cp.N = OBJ[ck]->getN(L,bestuv);
			cp.w = w;
			CP[++tCP] = cp; 
		}
		//if (deep >= 1) printf("out\n");
	}
}

inline void Build(int x,int k,int L,int R)
{
	Tree[x].L = L; Tree[x].R = R; Tree[x].l = Tree[x].r = 0; 
	Tree[x].rd = rds; Tree[x].flux = Color(0,0,0);
	if (k == 0) sort(CP + L,CP + R + 1,cmpx);
	if (k == 1) sort(CP + L,CP + R + 1,cmpy);
	if (k == 2) sort(CP + L,CP + R + 1,cmpz);
	int mid = L + R >> 1; Tree[x].now = CP[mid];
	if (L <= mid - 1) Build(Tree[x].l = ++tot,(k + 1) % 3,L,mid - 1);
	if (mid + 1 <= R) Build(Tree[x].r = ++tot,(k + 1) % 3,mid + 1,R);
}
inline void Prepare()
{
	int i,j,k; Point P; Line L;
	Eye = Point(-1.5,W * 0.5 / 900.,H * 0.5 / 900.); 
	rep(i,0,W - 1) rep(j,0,H - 1)
	{
		P = Point(0,i,H - 1 - j) + (Point(0,1,1) * ((db)(rand() % 100) / 100.0));//random select
		P = P * (1.0 / 900.0);
		L = Line(Eye,P + Point(-1.2,0,0.0) - Eye); 
		//printf("Line (%lf,%lf,%lf) [%lf %lf %lf]\n",L.P0.x,L.P0.y,L.P0.z,L.Pd.x,L.Pd.y,L.Pd.z);
		getCp(L,j,i,Color(1.0,1.0,1.0),0);
	}
	printf("fuck\n");
	Build(++tot,0,1,tCP);
	printf("fuck\n");
}
inline db mi(db x,int m)
{
	db ans = 1.; while (m) {if (m & 1) ans = ans * x; m >>= 1; x = x * x;}
	return ans;
}
inline void Add(int x,Line L,Color w,UV uv,Cp cp)
{
	int i = cp.i,j = cp.j,k = cp.k;
	db g = (cnt[i][j] * 0.7 + 0.7) / (cnt[i][j] + 1.0);
	Tree[x].rd *= g; cnt[i][j] += 1.0;
	Tree[x].flux = (Tree[x].flux + (w * OBJ[k]->mtr.Kd) * (1./pi)) * g;
}
inline void search(int i,int k,Line L,Color w,UV uv,int ck)
{
	int cutr,cutl; Cp cp; Point C = OBJ[ck]->get(uv);
	cutl = cutr = 0; cp = Tree[i].now;
	if (len(cp.P,C) <= Tree[i].rd)
	{ 
		Add(i,L,w,uv,cp);
	}
	if (k == 0){ cutl = sn(C.x - cp.P.x - rds) > 0; cutr = sn(cp.P.x - C.x - rds) > 0; }
	if (k == 1){ cutl = sn(C.y - cp.P.y - rds) > 0; cutr = sn(cp.P.y - C.y - rds) > 0; }
	if (k == 2){ cutl = sn(C.z - cp.P.z - rds) > 0; cutr = sn(cp.P.z - C.z - rds) > 0; }	
	if ((Tree[i].l) && (!cutl)) {search(Tree[i].l,(k + 1) % 3,L,w,uv,ck);}
	if ((Tree[i].r) && (!cutr)) {search(Tree[i].r,(k + 1) % 3,L,w,uv,ck);}
}
inline void getph(Line L,Color w,int deep)
{
	int k,ck = -1; UV bestuv; db dr = 1e+10,cur;
	//printf("Line (%lf,%lf,%lf) [%lf %lf %lf]\n",L.P0.x,L.P0.y,L.P0.z,L.Pd.x,L.Pd.y,L.Pd.z);
	rep(k,0,(int)OBJ.size() - 1)
	{
		UV uv = OBJ[k]->getCross(L),Inf = UV(1.2345,5.4321);
		Point C = OBJ[k]->get(uv);
		cur = len(C,L.P0);
		if ((sn(cur -  dr) < 0) && (sn(len(uv,Inf))))
		{
			dr = cur; bestuv = uv;  ck = k;
		}
	}
	if (ck != -1)
	{
		Point C = OBJ[ck]->get(bestuv);
		if (sn(w.len()) && (deep <= 20))
		{
			Mtr mtr = OBJ[ck]->mtr;
			db ar = mtr.ar,at = mtr.at,ad = mtr.ad;
			if (sn(mtr.wr.len())) getph(OBJ[ck]->getReflect(L,bestuv),w * mtr.wr * ar * 0.01,deep + 1);
		//	if (sn(wt.len())) getph(OBJ[ck]->getTrans(L,bestuv),w * mtr.wt * ar * 0.01,deep + 1);
			if (ad) search(1,0,L,w * ad * 0.01,bestuv,ck);
		}
		else search(1,0,L,w * 0.01,bestuv,ck);
	}
}
Color sIMG[H + 10][W + 10];
inline void PUT()
{
	cv :: Mat img = cv :: Mat(cv :: Size(W,H),CV_8UC3,cv :: Scalar(0));
	uchar *o,*p; int i,j,k;
	rep(i,0,H - 1) rep(j,0,W - 1) sIMG[i][j] = Color(0,0,0);
	rep(k,1,tot)
	{
		arr hp = Tree[k];
		sIMG[hp.now.i][hp.now.j] = sIMG[hp.now.i][hp.now.j] + hp.flux *  (1.0 / (pi * hp.rd * 1000. * 1000.));
	}
	rep(i,0,H - 1) 
	{
		o = img.ptr<uchar>(i);
		rep(j,0,W - 1)
		{
			*o++ = (uchar)toInt(sIMG[i][j].b);
			*o++ = (uchar)toInt(sIMG[i][j].g); 
			*o++ = (uchar)toInt(sIMG[i][j].r);
		}
	}
	cv :: imwrite("1.jpg",img);
}
inline void Work()
{
	int i; Line L; Point pd,C0; db maxr = 0.5,r,phi,theta;
	for (db u = 0.0; u <= 1.0; u += 0.001)
	{
		printf("u = %lf\n",u);
		for (db v = 0.0; v <= 1.0; v += 0.001)
		{
			db du,dv; du = (rand() % 100)/100000.0; dv = (rand() % 100)/100000.0;
			r = maxr * (u + du); theta = 2.0 * pi * (v + dv);
			phi = pi * (u + du);
			getph(Line(Light[0],Point(sin(phi) * cos(theta),sin(phi) * sin(theta),cos(phi))),Color(2500,2500,2500) * 4 * pi,0);
		}
		PUT();
	}
}

int main()
{
	Init();
	Prepare();
	Work();
	PUT();
	return 0;
}
