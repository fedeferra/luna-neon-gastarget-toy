// Monte Carlo of photons from point sources on the beam axis: peak and total efficiency
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#define NT 1700
#define PI 3.14159265358979
static double EG[NT], MU[4][3][NT], ZEFF[4];   // [det][ph,co,pr][E]
typedef struct{uint64_t s[2];}rng_t;
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t nxt(rng_t*r){uint64_t s0=r->s[0],s1=r->s[1],res=s0+s1;s1^=s0;r->s[0]=rotl(s0,55)^s1^(s1<<14);r->s[1]=rotl(s1,36);return res;}
static inline double U(rng_t*r){return ((nxt(r)>>11)+0.5)*(1.0/9007199254740992.0);}
static void seed(rng_t*r,uint64_t a){uint64_t z=a+0x9E3779B97F4A7C15ULL;for(int i=0;i<2;i++){z+=0x9E3779B97F4A7C15ULL;uint64_t x=z;x=(x^(x>>30))*0xBF58476D1CE4E5B9ULL;x=(x^(x>>27))*0x94D049BB133111EBULL;r->s[i]=x^(x>>31);}for(int i=0;i<20;i++)nxt(r);}
typedef struct{double c[3],a[3],R,h,Rin;int hole;}geo_t;
static geo_t G[4];
static int BOXN=0;static double BX[8][6];   // box: x0,x1,y0,y1,z0,z1 (union of disjoint boxes, g = 0)
static double LOGE0,DLOGE;
static inline int eidx(double E,double*f){double x=(log(E)-LOGE0)/DLOGE;if(x<0){*f=0;return 0;}int i=(int)x;if(i>=NT-1){*f=1;return NT-2;}*f=x-i;return i;}
static inline double mu(int d,int k,double E){double f;int i=eidx(E,&f);return MU[d][k][i]*(1-f)+MU[d][k][i+1]*f;}
// interval [t0,t1] of a ray inside a finite cylinder (axis a, centre c, radius R, half length h)
static int cyl(const double*p,const double*d,const double*c,const double*a,double R,double h,double*t0,double*t1){
  double q[3]={p[0]-c[0],p[1]-c[1],p[2]-c[2]};double qa=q[0]*a[0]+q[1]*a[1]+q[2]*a[2],da=d[0]*a[0]+d[1]*a[1]+d[2]*a[2];
  double qp[3]={q[0]-qa*a[0],q[1]-qa*a[1],q[2]-qa*a[2]},dp[3]={d[0]-da*a[0],d[1]-da*a[1],d[2]-da*a[2]};
  double A=dp[0]*dp[0]+dp[1]*dp[1]+dp[2]*dp[2],B=qp[0]*dp[0]+qp[1]*dp[1]+qp[2]*dp[2],C=qp[0]*qp[0]+qp[1]*qp[1]+qp[2]*qp[2]-R*R;
  double lo=-1e30,hi=1e30;
  if(A<1e-14){if(C>=0)return 0;}else{double disc=B*B-A*C;if(disc<=0)return 0;double s=sqrt(disc);lo=(-B-s)/A;hi=(-B+s)/A;}
  if(fabs(da)<1e-14){if(fabs(qa)>h)return 0;}else{double s0=(-h-qa)/da,s1=(h-qa)/da;if(s0>s1){double t=s0;s0=s1;s1=t;}if(s0>lo)lo=s0;if(s1<hi)hi=s1;}
  if(hi<=lo)return 0;*t0=lo;*t1=hi;return 1;}
// segments of the solid (cylinder, possibly with a bore) along the ray, with t>0
static int boxint(const double*p,const double*d,const double*b,double*t0,double*t1){double lo=-1e30,hi=1e30;
  for(int k=0;k<3;k++){double a0=b[2*k],a1=b[2*k+1];if(fabs(d[k])<1e-14){if(p[k]<a0||p[k]>a1)return 0;}else{double s0=(a0-p[k])/d[k],s1=(a1-p[k])/d[k];if(s0>s1){double t=s0;s0=s1;s1=t;}if(s0>lo)lo=s0;if(s1<hi)hi=s1;}}
  if(hi<=lo||hi<=0)return 0;*t0=lo<0?0:lo;*t1=hi;return 1;}
// prisma esagonale (apotema a, semilunghezza hz lungo z): normali a 90 + 60k gradi, cioe facce piane in alto e in basso
static int HEX=0;
static int hexint(const double*p,const double*d,double a,double hz,double*t0,double*t1){double lo=-1e30,hi=1e30;
  for(int k=0;k<3;k++){double th=(90+60*k)*PI/180,nx=cos(th),ny=sin(th),pn=p[0]*nx+p[1]*ny,dn=d[0]*nx+d[1]*ny;
    if(fabs(dn)<1e-14){if(fabs(pn)>a)return 0;}else{double s0=(-a-pn)/dn,s1=(a-pn)/dn;if(s0>s1){double t=s0;s0=s1;s1=t;}if(s0>lo)lo=s0;if(s1<hi)hi=s1;}}
  if(fabs(d[2])<1e-14){if(fabs(p[2])>hz)return 0;}else{double s0=(-hz-p[2])/d[2],s1=(hz-p[2])/d[2];if(s0>s1){double t=s0;s0=s1;s1=t;}if(s0>lo)lo=s0;if(s1<hi)hi=s1;}
  if(hi<=lo)return 0;*t0=lo;*t1=hi;return 1;}
static int segs(int g,const double*p,const double*d,double*s,double*e){
  if(g==1&&HEX){double o0,o1;if(!hexint(p,d,10,14,&o0,&o1))return 0;if(o1<=0)return 0;if(o0<0)o0=0;int n=0;double i0,i1;
    if(hexint(p,d,3,28,&i0,&i1)&&i1>o0&&i0<o1){if(i0>o0){s[n]=o0;e[n]=i0;n++;}if(i1<o1){s[n]=i1;e[n]=o1;n++;}return n;}
    s[0]=o0;e[0]=o1;return 1;}
  if(g==0&&BOXN){int n=0;for(int b=0;b<BOXN;b++){double t0,t1;if(boxint(p,d,BX[b],&t0,&t1)){s[n]=t0;e[n]=t1;n++;}}return n;}
  double o0,o1;if(!cyl(p,d,G[g].c,G[g].a,G[g].R,G[g].h,&o0,&o1))return 0;if(o1<=0)return 0;if(o0<0)o0=0;
  int n=0;
  if(G[g].hole){double i0,i1;if(cyl(p,d,G[g].c,G[g].a,G[g].Rin,G[g].h*2,&i0,&i1)&&i1>o0&&i0<o1){if(i0>o0){s[n]=o0;e[n]=i0;n++;}if(i1<o1){s[n]=i1;e[n]=o1;n++;}return n;}}
  s[0]=o0;e[0]=o1;return 1;}
static void isodir(rng_t*r,double*d){double c=2*U(r)-1,ph=2*PI*U(r),s=sqrt(1-c*c);d[0]=s*cos(ph);d[1]=s*sin(ph);d[2]=c;}
static void rot(rng_t*r,double*d,double ct){double st=sqrt(fmax(0,1-ct*ct)),ph=2*PI*U(r),cp=cos(ph),sp=sin(ph);double ux=d[0],uy=d[1],uz=d[2];
  if(fabs(uz)>0.99999){d[0]=st*cp;d[1]=st*sp;d[2]=(uz>0?ct:-ct);return;}
  double sq=sqrt(1-uz*uz);double nx=st*cp*ux*uz/sq-st*sp*uy/sq+ct*ux,ny=st*cp*uy*uz/sq+st*sp*ux/sq+ct*uy,nz=-st*cp*sq+ct*uz;d[0]=nx;d[1]=ny;d[2]=nz;}
static double fwhm(int g,double E){double R=(g==0)?0.07:(g==1)?0.10:0.0015/0.0;if(g==2){double f=0.002*1.332*sqrt(E/1.332);return fmax(f,0.0015);} return R*sqrt(0.662*E);}
static inline double Yrad(int g,double T){double x=6e-4*ZEFF[g]*T;return x/(1+x);}
typedef struct{double p[3],d[3],E;}ph_t;
// one history: returns the total deposited energy
static double history(int g,double E0,double zs,rng_t*r){
  ph_t st[64];int n=0;st[0].p[0]=st[0].p[1]=0;st[0].p[2]=zs;isodir(r,st[0].d);st[0].E=E0;n=1;double dep=0;
  if(g==2){} 
  while(n>0){ph_t ph=st[--n];
    for(;;){
      if(ph.E<0.010){dep+=ph.E;break;}
      double s[16],e[16];int mt_[16];int ns=0;{double s1[8],e1[8];int n1=segs(g,ph.p,ph.d,s1,e1);for(int k=0;k<n1;k++){s[ns]=s1[k];e[ns]=e1[k];mt_[ns]=g;ns++;}
        int n2=segs(3,ph.p,ph.d,s1,e1);for(int k=0;k<n2;k++){s[ns]=s1[k];e[ns]=e1[k];mt_[ns]=3;ns++;}}
      if(ns==0)break;
      for(int a=1;a<ns;a++){double ss=s[a],ee=e[a];int mm=mt_[a];int b=a-1;while(b>=0&&s[b]>ss){s[b+1]=s[b];e[b+1]=e[b];mt_[b+1]=mt_[b];b--;}s[b+1]=ss;e[b+1]=ee;mt_[b+1]=mm;}
      int hit=0,gm=g;double tint=0;double mph=0,mco=0,mpr=0,mt=0;
      for(int k=0;k<ns;k++){double a1=mu(mt_[k],0,ph.E),a2=mu(mt_[k],1,ph.E),a3=mu(mt_[k],2,ph.E),at=a1+a2+a3;double l=-log(U(r))/at;if(l<e[k]-s[k]){tint=s[k]+l;hit=1;gm=mt_[k];mph=a1;mco=a2;mpr=a3;mt=at;break;}}
      if(!hit)break;
      for(int j=0;j<3;j++)ph.p[j]+=tint*ph.d[j];
      double u=U(r)*mt;
      if(u<mph){if(gm==g)dep+=ph.E;break;}
      else if(u<mph+mco){
        double k=ph.E/0.510999,eps0=1/(1+2*k),eps,ct;
        for(;;){double a1=log(1/eps0),a2=(1-eps0*eps0)*0.5,a=a1+a2;double ee;
          if(U(r)<a1/a)ee=exp(-U(r)*a1);else ee=sqrt(eps0*eps0+(1-eps0*eps0)*U(r));
          double t=(1-ee)/(k*ee);double sin2=t*(2-t);if(sin2<0)continue;
          if(U(r)<=1-ee*sin2/(1+ee*ee)){eps=ee;ct=1-t;break;}}
        double Ee=ph.E*eps,T=ph.E-Ee,y=(gm==g)?Yrad(g,T):0;
        if(y*T>0.01&&n<63){if(gm==g)dep+=T*(1-y);ph_t b;for(int j=0;j<3;j++)b.p[j]=ph.p[j];isodir(r,b.d);b.E=y*T;st[n++]=b;}else if(gm==g)dep+=T;
        rot(r,ph.d,ct);ph.E=Ee;
      }else{
        double T=ph.E-1.021998,y=(gm==g)?Yrad(g,0.5*T):0;
        if(y*T>0.01&&n<62){if(gm==g)dep+=T*(1-y);ph_t b;for(int j=0;j<3;j++)b.p[j]=ph.p[j];isodir(r,b.d);b.E=y*T;st[n++]=b;}else if(gm==g)dep+=T;
        if(n<62){ph_t a1;for(int j=0;j<3;j++)a1.p[j]=ph.p[j];isodir(r,a1.d);a1.E=0.510999;ph_t a2=a1;for(int j=0;j<3;j++)a2.d[j]=-a1.d[j];a2.E=0.510999;st[n++]=a1;st[n++]=a2;}else if(gm==g)dep+=1.021998;
        break;}
    }
  }
  return dep;}
int main(int argc,char**argv){
  long NH=argc>1?atol(argv[1]):100000;
  FILE*f=fopen("tables.bin","rb");if(!f)return 1;
  if(fread(EG,8,NT,f)!=NT)return 1;
  for(int d=0;d<4;d++){for(int k=0;k<3;k++)if(fread(MU[d][k],8,NT,f)!=NT)return 1;if(fread(&ZEFF[d],8,1,f)!=1)return 1;}
  fclose(f);LOGE0=log(EG[0]);DLOGE=(log(EG[NT-1])-LOGE0)/(NT-1);
  // geometries (cm): 0 = segmented NaI, 20 cm long; 1 = BGO 28 cm long, hexagonal bore (apothem 3 cm), 7 cm thick; 2 = HPGe 8x8 cm, transverse axis
  for(int g=0;g<4;g++){G[g].hole=0;G[g].c[0]=G[g].c[1]=G[g].c[2]=0;G[g].a[0]=G[g].a[1]=0;G[g].a[2]=1;}
  G[0].R=5+10;G[0].Rin=5;G[0].h=10;G[0].hole=1;
  { // segmented NaI: 6 crystals 10x10x20 cm in a brick arrangement around a 10x10 cm square bore, centred on z = 0 (the deposits in all crystals are summed)
    double B[6][4]={{-10,0,5,15},{0,10,5,15},{5,15,-5,5},{0,10,-15,-5},{-10,0,-15,-5},{-15,-5,-5,5}};BOXN=6;
    for(int b=0;b<6;b++){BX[b][0]=B[b][0];BX[b][1]=B[b][1];BX[b][2]=B[b][2];BX[b][3]=B[b][3];BX[b][4]=-10;BX[b][5]=10;}}
  G[1].R=3+7;G[1].Rin=3;G[1].h=14;G[1].hole=1;HEX=1;   // BGO: 6 prismi trapezoidali (altezza 7 cm) in contatto = guscio esagonale, apotema 3 cm (foro) e 10 cm, lungo 28 cm
  G[3].R=1.2;G[3].Rin=1.0;G[3].h=40;G[3].hole=1;
  G[2].R=4;G[2].h=4;G[2].c[0]=4+4;G[2].a[0]=1;G[2].a[2]=0;
  const int NE=28,NZ=21;double Ev[28]={0.02,0.05,0.1,0.15,0.2,0.3,0.4,0.511,0.662,0.8,1.0,1.332,1.7,2,2.6,3,4,5,6,7,8,9,10,12,14,16,18,20};
  double zv[21];for(int i=0;i<NZ;i++)zv[i]=1.5*i;
  static double pk[3][28][21],tot[3][28][21];
  int gmin=0,gmax=3;if(argc>2){gmin=argv[2][0]-'0';gmax=gmin+1;}
  for(int g=gmin;g<gmax;g++){
    #pragma omp parallel for schedule(dynamic)
    for(int idx=0;idx<NE*NZ;idx++){int ie=idx/NZ,iz=idx%NZ;rng_t r;seed(&r,1000003ULL*(g+1)+idx);
      long cp=0,ct=0;double E=Ev[ie],w=fwhm(g,E);
      for(long h=0;h<NH;h++){double d=history(g,E,zv[iz],&r);if(d>0.0005)ct++;
        // resolution: Gaussian with sigma = FWHM(d)/2.3548 on the deposited energy; peak = symmetric window E +/- FWHM(E)
        if(d>0.0005){double u1=U(&r),u2=U(&r),gs=sqrt(-2*log(u1))*cos(2*PI*u2);double ds=d+gs*fwhm(g,d)/2.3548;if(ds>=E-w&&ds<=E+w)cp++;}}
      pk[g][ie][iz]=(double)cp/NH;tot[g][ie][iz]=(double)ct/NH;}
    fprintf(stderr,"det %d done\n",g);}
  printf("{\"E\":[");for(int i=0;i<NE;i++)printf("%s%g",i?",":"",Ev[i]*1000);printf("],\"z\":[");for(int i=0;i<NZ;i++)printf("%s%g",i?",":"",zv[i]);printf("],\"N\":%ld",NH);
  const char*nm[3]={"nai","bgo","ge"};
  for(int g=gmin;g<gmax;g++){printf(",\"%s\":{",nm[g]);for(int t=0;t<2;t++){printf("%s\"%s\":[",t?",":"",t?"tot":"pk");for(int i=0;i<NE;i++){printf("%s[",i?",":"");for(int j=0;j<NZ;j++)printf("%s%.5f",j?",":"",t?tot[g][i][j]:pk[g][i][j]);printf("]");}printf("]");}printf("}");}
  printf("}\n");return 0;}
