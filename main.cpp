#include "obj.h"

#define pb push_back

const int N = 3000000,W = 800,H = 600;
const db pi = 2.0 * atan2(1,0),rds = 1e-1;

vector <obj * > OBJ;
vector <Point> Light;
int tCP,tot;
Cp CP[N];

Color IMG[H + 10][W + 10];
db cnt[H + 10][W + 10];
Point Eye,InfP;

struct arr
{
	int k,l,r,L,R; Cp now;
}Tree[N];

inline int sn(db x){return x < -1e-9 ? -1 : x > 1e-9 ? 1 : 0;}
inline bool cmpx(const Cp & A,const Cp & B){ return sn(A.P.x - B.P.x) < 0;}
inline bool cmpy(const Cp & A,const Cp & B){ return sn(A.P.y - B.P.y) < 0;}
inline bool cmpz(const Cp & A,const Cp & B){ return sn(A.P.z - B.P.z) < 0;}
inline void Init()
{
	srand(time(0)); Eye = Point(-1,0,1); Light.clear(); 
	Light.pb(Point(-1,-1,1));OBJ.clear(); 
	OBJ.pb(new Plane(Line(Point(0,0,0),Point(0,0,1)),Point(1,0,0),Point(0,1,0))); 
	tCP = 0; tot = 0;
	int i,j; rep(i,0,H - 1) rep(j,0,W - 1){IMG[i][j] = Color(0,0,0); cnt[i][j] = 0.;}
}
inline db len(const Point &A,const Point &B) {return sqrt((A - B) * (A - B));}
inline db len(const UV &A,const UV &B) {return fabs(A.u - B.u) + fabs(A.v - B.v);}
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
			Color wt = OBJ[ck]->mtr.wt;
			if (sn(wt.len())) getCp(OBJ[ck]->getTrans(L,bestuv),i,j,w * wt,deep + 1);
			
			Color wm = OBJ[ck]->mtr.wm;
			Cp cp; cp.P = OBJ[ck]->get(bestuv); cp.i = i; cp.j = j; 
			cp.k = ck; cp.V = L; cp.N = OBJ[ck]->getN(L,bestuv);
			cp.w = w * wm;
			CP[++tCP] = cp; 
			if (tCP % 10000 == 0)
			{
				printf("%d\n",tCP);
				printf("Line(%lf,%lf,%lf),[%lf %lf %lf]) cp = (%lf,%lf,%lf)\n",L.P0.x,L.P0.y,L.P0.z,L.Pd.x,L.Pd.y,L.Pd.z,cp.P.x,cp.P.y,cp.P.z);
				//system("pause");
			}
		}
	}
}

inline void Build(int x,int k,int L,int R)
{
	Tree[x].L = L; Tree[x].R = R; Tree[x].l = Tree[x].r = 0;
	if (k == 0) sort(CP + L,CP + R + 1,cmpx);
	if (k == 1) sort(CP + L,CP + R + 1,cmpy);
	if (k == 2) sort(CP + L,CP + R + 1,cmpz);
	int mid = L + R >> 1;
	if (L <= mid - 1) Build(Tree[x].l = ++tot,(k + 1) % 3,L,mid - 1);
	if (mid + 1 <= R) Build(Tree[x].r = ++tot,(k + 1) % 3,mid + 1,R);
}
inline void Prepare()
{
	int i,j,k; Point P; Line L;
	rep(i,0,H - 1) rep(j,0,W - 1)
	{
		P = Point(0,i,j) + (Point(0,1,1) * ((db)(rand() % 100) / 100.0));//random select
		P = P * (1.0 / 800.0);
		L = Line(Eye,P - Eye); 
		//printf("Line (%lf,%lf,%lf) [%lf %lf %lf]\n",L.P0.x,L.P0.y,L.P0.z,L.Pd.x,L.Pd.y,L.Pd.z);
		getCp(L,i,j,Color(1.0,1.0,1.0),0);
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
inline void Add(Line L,Color w,UV uv,Cp cp)
{
	int i = cp.i,j = cp.j,k = cp.k,n = 5;
	Color Ans = Color(0,0,0),Ia = Color(30,30,30),Ip = Color(255,255,255) * w;
	Color Ka = OBJ[k]->mtr.Ka, Kd = OBJ[k]->mtr.Kd, Ks = OBJ[k]->mtr.Ks;
	Line R = OBJ[k]->getReflect(L,uv);
	IMG[i][j] = (IMG[i][j] * cnt[i][j]) + (Ia * Ka) + (Ip * Kd * (L * cp.N)) + Ip * Ks * mi((R * cp.V),n);
	IMG[i][j] = IMG[i][j] * (1.0 / (1.0 + cnt[i][j]));
	cnt[i][j] += 1.0;

}
inline void search(int i,int k,Line L,Color w,UV uv,int ck)
{
	int cutr,cutl; Cp cp; Point C = OBJ[ck]->get(uv);
	cutl = cutr = 0; cp = Tree[i].now;
	if (len(cp.P,C) <= rds) Add(L,w,uv,cp);
	if (k == 0){ cutl = sn(C.x - cp.P.x - rds) > 0; cutr = sn(cp.P.x - C.x - rds) > 0; }
	if (k == 1){ cutl = sn(C.y - cp.P.y - rds) > 0; cutr = sn(cp.P.y - C.y - rds) > 0; }
	if (k == 2){ cutl = sn(C.z - cp.P.z - rds) > 0; cutr = sn(cp.P.z - C.z - rds) > 0; }	
	if ((Tree[i].l) && (!cutl)) {search(Tree[i].l,(k + 1) % 3,L,w,uv,ck);}
	if ((Tree[i].r) && (!cutr)) {search(Tree[i].r,(k + 1) % 3,L,w,uv,ck);}
}
inline void getph(Line L,Color w)
{
	int k,ck = -1; UV bestuv; db dr = 1e+10,cur;
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
		Point C = OBJ[k]->get(bestuv);
		printf("C %lf %lf %lf\n",C.x,C.y,C.z);
		if (sn(w.len()))
		{
			Mtr mtr = OBJ[ck]->mtr;
			int ar = mtr.ar,at = mtr.at,ad = mtr.ad;
			int Rdd = rand() % 100;
			if (Rdd <= ar)
			{
				Color wr = mtr.wr;
				if (sn(wr.len())) getph(OBJ[ck]->getReflect(L,bestuv),w * wr);
			}else if (Rdd <= ar + at)
			{
				Color wt = mtr.wt;
				if (sn(wt.len())) getph(OBJ[ck]->getTrans(L,bestuv),w * wt);
			}else if (Rdd <= ar + at + ad) search(1,0,L,w,bestuv,ck);
		}
		else search(1,0,L,w,bestuv,ck);
	}
}
inline void Work()
{
	int i; Line L; Point pd,C0; db maxr = 1.0,r,theta;
	for (db u = 0.0; u <= 1.0; u += 0.01)
		for (db v = 0.0; v <= 1.0; v += 0.01)
		{
			r = maxr * u + (rand() % 100)/10000.0; theta = 2.0 * pi * v * (rand() % 100)/10000.0;
			getph(Line(Light[0],Point(r * cos(theta),r * sin(theta),0) - Light[0]),Color(1.0,1.0,1.0));
		}
}
inline void PUT()
{
	cv :: Mat img = cv :: Mat(cv :: Size(W,H),CV_8UC3,cv :: Scalar(0));
	uchar *o,*p; int i,j;
	rep(i,0,H - 1) 
	{
		o = img.ptr<uchar>(i);
		rep(j,0,W - 1)
		{
			*o++ = (uchar)(int)(IMG[i][j].b / cnt[i][j]); 
			*o++ = (uchar)(int)(IMG[i][j].g / cnt[i][j]); 
			*o++ = (uchar)(int)(IMG[i][j].r / cnt[i][j]);
		}
	}
	cv :: imwrite("1.jpg",img);
}
int main()
{
	Init();
	Prepare();
	Work();
	PUT();
	return 0;
}
