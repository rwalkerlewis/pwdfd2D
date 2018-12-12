/* This file is automatically generated. DO NOT EDIT! */

#ifndef _ksutil_h
#define _ksutil_h


/*------------------------------------------------------------*/
/* CENTERED derivatives                                       */
/*------------------------------------------------------------*/
#define NOP 4 /* derivative operator half-size */
#define C1 4.0f/5.0f
#define C2 -1.0f/5.0f
#define C3 4.0f/105.0f
#define C4 -1.0f/280.0f
#define Dx(a,ix,iy,iz,s) (C4*(a[iy][ix+4][iz] - a[iy][ix-4][iz]) +	\
			  C3*(a[iy][ix+3][iz] - a[iy][ix-3][iz]) +	\
			  C2*(a[iy][ix+2][iz] - a[iy][ix-2][iz]) +	\
			  C1*(a[iy][ix+1][iz] - a[iy][ix-1][iz])  )*s
#define Dy(a,ix,iy,iz,s) (C4*(a[iy+4][ix][iz] - a[iy-4][ix][iz]) +	\
			  C3*(a[iy+3][ix][iz] - a[iy-3][ix][iz]) +	\
			  C2*(a[iy+2][ix][iz] - a[iy-2][ix][iz]) +	\
			  C1*(a[iy+1][ix][iz] - a[iy-1][ix][iz])  )*s
#define Dz(a,ix,iy,iz,s) (C4*(a[iy][ix][iz+4] - a[iy][ix][iz-4]) +	\
			  C3*(a[iy][ix][iz+3] - a[iy][ix][iz-3]) +	\
			  C2*(a[iy][ix][iz+2] - a[iy][ix][iz-2]) +	\
			  C1*(a[iy][ix][iz+1] - a[iy][ix][iz-1])  )*s


#define XX 0
#define YY 1
#define ZZ 2
#define XY 3
#define XZ 4
#define YZ 5


typedef struct fdm2 *fdm2d;


typedef struct fdm3 *fdm3d;


typedef struct lcoef2 *lint2d;


typedef struct lcoef3 *lint3d;


typedef struct abc2 *abcone2d;


typedef struct abc3 *abcone3d;


typedef struct spon *sponge;


typedef struct ofg *ofg2d;


typedef struct dft3 *dft3d;


typedef struct lps3 *lps3d;


typedef struct ksp3 *ksp3d;


typedef struct vksp3 *vksp3d;


typedef struct clr3 *clr3d;


typedef struct dat3 *dat3d;


typedef struct mut3 *mut3d;


struct fdm2{
    int nb;
    int   nz,nzpad;
    int   nx,nxpad;
    float oz,ozpad;
    float ox,oxpad;
    float dz;
    float dx;
    bool verb;
    bool free;
    int ompchunk;
};


struct fdm3{
    int nb;
    int   nz,nzpad;
    int   nx,nxpad;
    int   ny,nypad;
    float oz,ozpad;
    float ox,oxpad;
    float oy,oypad;
    float dz;
    float dx;
    float dy;
    bool verb;
    bool free;
    int ompchunk;
};


struct lcoef2{
    int n;
    float *w00;
    float *w01;
    float *w10;
    float *w11;
    int *jz;
    int *jx;
};


struct lcoef3{
    int n;
    float *w000;
    float *w001;
    float *w010;
    float *w011;
    float *w100;
    float *w101;
    float *w110;
    float *w111;
    int *jz;
    int *jx;
    int *jy;
};


struct abc2{
    bool free;
    float *bzl;
    float *bzh;
    float *bxl;
    float *bxh;
};


struct abc3{
    bool free;
    float**bzl;
    float**bzh;
    float**bxl;
    float**bxh;
    float**byl;
    float**byh;
};


struct spon{
    float *w;
};


struct ofg{
    float **tt;
};


struct dft3{
    int nkz;
    int nkx;
    int nky;
    float okz;
    float okx;
    float oky;
    float dkz;
    float dkx;
    float dky;
};


struct lps3{
    float ***flt;
};


struct ksp3{
    float***kz;
    float***kx;
    float***ky;
};


struct vksp3{
    float***kpa;
    float***ksa;
    float***kpb;
    float***ksb;
};


struct clr3{
    int *n2_start;
    int *n2_local;
    int **map;
    int n2_sum;
    int n2_max;
};


struct dat3{
    int   ns,nr,nc;
    int   nt,nx,ny;
    float ot,ox,oy;
    float dt,dx,dy;
};


struct mut3{
    float **wt;
};


/*------------------------------------------------------------*/
fdm2d fdutil_init(bool verb_,
		  bool free_,
		  sf_axis az_,
		  sf_axis ax_,
		  int     nb_,
		  int ompchunk_);
/*< init fdm utilities >*/


/*------------------------------------------------------------*/
fdm3d fdutil3d_init(bool verb_, 
		    bool free_,
		    sf_axis az_, 
		    sf_axis ax_, 
		    sf_axis ay_, 
		    int     nb_,
		    int ompchunk_);
/*< init fdm utilities >*/


/*------------------------------------------------------------*/
ofg2d offgrid_init(fdm2d fdm);
/*< init off-grid interpolation >*/


/*------------------------------------------------------------*/
void offgridfor(float **ti,
		ofg2d  ofg,
		fdm2d  fdm);
/*< forward off-grid interpolation (in place) >*/


/*------------------------------------------------------------*/
void offgridadj(float **ti,
		ofg2d  ofg,
		fdm2d  fdm);
/*< adjoint off-grid interpolation (in place) >*/


/*------------------------------------------------------------*/
void expand(float** a, 
	    float** b, 
	    fdm2d fdm);
/*< expand domain >*/


/*------------------------------------------------------------*/
void expand3d(float ***a, 
	      float ***b, 
	      fdm3d  fdm);
/*< expand domain >*/


/*------------------------------------------------------------*/
void expand2d(float** a, 
              float** b, 
	      fdm3d fdm);
/*< expand domain >*/


/*------------------------------------------------------------*/
void cut2d(float**  a,
	   float**  b,
	   fdm2d  fdm,
	   sf_axis cz, 
	   sf_axis cx);
/*< cut a rectangular wavefield subset >*/


/*------------------------------------------------------------*/
void cut3d(float*** a,
	   float*** b,
	   fdm3d  fdm,
	   sf_axis cz, 
	   sf_axis cx,
	   sf_axis cy);
/*< cut a rectangular wavefield subset >*/


/*------------------------------------------------------------*/
void cut3d_complex(sf_complex*** a,
	           sf_complex*** b,
                   fdm3d  fdm,
                   sf_axis cz, 
                   sf_axis cx,
                   sf_axis cy);
/*< cut a rectangular wavefield subset and convert to float >*/


/*------------------------------------------------------------*/
void cut3d_slice(float** a,
                 float** b,
                 fdm3d  fdm,
                 sf_axis cz,
                 sf_axis cx);
/*< cut a rectangular wavefield subset >*/


/*------------------------------------------------------------*/
void bfill(float** b, 
	   fdm2d fdm);
/*< fill boundaries >*/


/*------------------------------------------------------------*/
lint2d lint2d_make(int    na, 
		   pt2d*  aa, 
		   fdm2d fdm);
/*< init 2D linear interpolation >*/


/*------------------------------------------------------------*/
lint3d lint3d_make(int    na, 
		   pt3d*  aa, 
		   fdm3d fdm);
/*< init 3D linear interpolation >*/


/*------------------------------------------------------------*/
void lint2d_hold(float**uu,
		 float *ww,
		 lint2d ca);
/*< hold fixed value in field >*/


/*------------------------------------------------------------*/
void lint2d_inject(float**uu,
		   float *dd,
		   lint2d ca);
/*< inject into wavefield >*/


void lint2d_extract(float**uu,
		    float* dd,
		    lint2d ca);
/*< extract from wavefield >*/


/*------------------------------------------------------------*/
void lint3d_inject(float***uu,
		   float  *dd,
		   lint3d  ca);
/*< inject into wavefield >*/


void lint3d_extract(float***uu,
		    float  *dd,
		    lint3d  ca);
/*< extract from wavefield >*/


/*------------------------------------------------------------*/
void lint3d_inject_complex(sf_complex***uu,
                           sf_complex  *dd,
                           lint3d  ca);
/*< inject into wavefield >*/


void lint3d_extract_complex(sf_complex***uu,
                            sf_complex  *dd,
                            lint3d  ca);
/*< extract from wavefield >*/


/*------------------------------------------------------------*/
void lint3d_expl_complex(sf_complex***uz,
                         sf_complex***ux,
                         sf_complex***uy,
                         sf_complex **dd,
                         lint3d  ca);
/*< inject into wavefield >*/


/*------------------------------------------------------------*/
void lint2d_inject1(float**uu,
		    float  dd,
		    lint2d ca);
/*< inject into wavefield >*/


/*------------------------------------------------------------*/
void lint3d_inject1(float***uu,
		    float   dd,
		    lint3d  ca);
/*< inject into wavefield >*/


/*------------------------------------------------------------*/
void fdbell_init(int n);
/*< init bell taper >*/


/*------------------------------------------------------------*/
void fdbell3d_init(int n);
/*< init bell taper >*/


/*------------------------------------------------------------*/
void lint2d_bell(float**uu,
		 float *ww,
		 lint2d ca);
/*< apply bell taper >*/


/*------------------------------------------------------------*/
void lint3d_bell(float***uu,
		 float  *ww,
		 lint3d  ca);
/*< apply bell taper >*/


/*------------------------------------------------------------*/
void lint3d_bell_complex(sf_complex***uu,
                         sf_complex  *ww,
                         lint3d  ca);
/*< apply bell taper >*/


/*------------------------------------------------------------*/
void lint2d_bell1(float**uu,
		  float ww,
		  lint2d ca);
/*< apply bell taper >*/


/*------------------------------------------------------------*/
void lint3d_bell1(float***uu,
		  float  ww,
		  lint3d  ca);
/*< apply bell taper >*/


/*------------------------------------------------------------*/
abcone2d abcone2d_make(int     nop,
		       float    dt,
		       float**  vv,
		       bool   free, 
		       fdm2d   fdm);
/*< init 2D ABC >*/


/*------------------------------------------------------------*/
abcone3d abcone3d_make(int     nop,
		       float    dt,
		       float ***vv,
		       bool   free, 
		       fdm3d   fdm);
/*< init 3D ABC >*/


/*------------------------------------------------------------*/
void abcone2d_apply(float**   uo,
		    float**   um,
		    int      nop,
		    abcone2d abc,
		    fdm2d    fdm);
/*< apply 2D ABC >*/


/*------------------------------------------------------------*/
void abcone3d_apply(float  ***uo,
		    float  ***um,
		    int      nop,
		    abcone3d abc,
		    fdm3d    fdm);
/*< apply 3D ABC >*/


/*------------------------------------------------------------*/
sponge sponge_make(int nb);
/*< init boundary sponge >*/


sponge sponge_make2(int nb, float cb);
/*< init boundary sponge >*/


/*------------------------------------------------------------*/
void sponge2d_apply(float**   uu,
		    sponge   spo,
		    fdm2d    fdm);
/*< apply boundary sponge >*/


/*------------------------------------------------------------*/
void sponge3d_apply(float  ***uu,
		    sponge   spo,
		    fdm3d    fdm);
/*< apply boundary sponge >*/


void sponge3d_apply_complex(sf_complex  ***uu,
		            sponge   spo,
		            fdm3d    fdm);
/*< apply boundary sponge for complex-valued wavefield >*/


bool cfl_generic(
    float vpmin, float vpmax,
    float dx, float dy, float dz,
    float dt, float fmax, float safety,
    int intervals, char *wave);
/*< cfl check for both 2d and 3d acoustic fdcode >*/


bool cfl_elastic(
    float vpmin, float vpmax,
    float vsmin, float vsmax,
    float dx, float dy, float dz,
    float dt, float fmax, float safety,
    int intervals);
/*< cfl check for both 2d and 3d elastic fdcode >*/


bool cfl_acoustic(
    float vpmin, float vpmax,
    float dx, float dy, float dz,
    float dt, float fmax, float safety,
    int intervals);
/*< cfl check for acoustic wave equation >*/


/*------------------------------------------------------------*/
dft3d dft3d_init(int pad1,
                 bool rtoc,
                 bool padio,
                 fdm3d fdm);
/*< init 3d fft >*/


/*------------------------------------------------------------*/
void dft3d_finalize();
/*< finalize 3d fft >*/


/*------------------------------------------------------------*/
lps3d lps3d_init(float pcut /* pcut/2 is the width of tapered region w.r.t. 1 */,
                 float fcut /* cutoff frequency */,
                 float vmax /* maximum p-wave velocity */,
                 dft3d dft);
/*< initialize low-pass filter coefficients >*/


ksp3d ksp3d_make(dft3d dft);
/*< init k-space arrays for pseudo-spectral method >*/


/*------------------------------------------------------------*/
void ksp3d_apply(float *wavedx,
                 float *wavedy,
                 float *wavedz,
                 float *wave,
                 dft3d dft,
                 ksp3d ksp);
/*< apply k-space to calculate spatial derivatives (can be in-place) >*/


/*------------------------------------------------------------*/
void ksp3d_apply1(float *wavedx,
                  float *wavedy,
                  float *wavedz,
                  float *wavex,
                  float *wavey,
                  float *wavez,
                  dft3d dft,
                  ksp3d ksp);
/*< apply k-space to calculate spatial derivatives (can be in-place) >*/


/*------------------------------------------------------------*/
void ksp3d_apply2(float *wavedxx,
                  float *wavedyy,
                  float *wavedzz,
                  float *wavedxy,
                  float *wavedyz,
                  float *wavedzx,
                  float *wave,
                  dft3d dft,
                  ksp3d ksp);
/*< apply k-space to calculate spatial derivatives (can be in-place) >*/


/*------------------------------------------------------------*/
void ksp3d_finalize();
/*< free static memory allocated for ksp >*/


vksp3d vksp3d_make(float gpavg,
                   float gsavg,
                   dft3d dft,
                   lps3d lps);
/*< init k-space arrays for fractional-Laplacian pseudo-spectral method >*/


/*------------------------------------------------------------*/
void vksp3d_apply(float *wavea,
                  float *waveb,
                  int mode /* 1 -> p mode; 2 -> s mode */,
                  dft3d dft,
                  vksp3d vksp);
/*< apply k-space to calculate spatial derivatives (can be in-place) >*/


/*------------------------------------------------------------*/
void vksp3d_finalize();
/*< free static memory allocated for fft and ksp >*/


clr3d clr3d_make(int *n2s,
                 fdm3d fdm);
/*< prepare lowrank arrays for rite method >*/


clr3d clr3d_make2(int *n2s,
                  fdm3d fdm);
/*< prepare lowrank arrays for rite method >*/


void clr3d_init(fdm3d fdm,
                dft3d dft,
                clr3d clr);
/*< initialize computation array >*/


/*------------------------------------------------------------*/
void clr3d_apply(sf_complex **uo,
                 sf_complex **ui,
                 sf_complex **lt,
                 sf_complex **rt,
                 fdm3d fdm,
                 dft3d dft,
                 clr3d clr);
/*< apply lowrank matrices for time stepping (can be in-place) >*/


/*------------------------------------------------------------*/
void clr3d_apply2(sf_complex **u2,
                  sf_complex **u1,
                  sf_complex **u0,
                  sf_complex **lt,
                  sf_complex **rt,
                  fdm3d fdm,
                  dft3d dft,
                  clr3d clr);
/*< apply lowrank matrices for time stepping (can be in-place) -- two-step version >*/


/*------------------------------------------------------------*/
void clr3d_apply_dbg(sf_complex **uo,
                     sf_complex **ui,
                     sf_complex **lt,
                     sf_complex **rt,
                     fdm3d fdm,
                     dft3d dft,
                     clr3d clr);
/*< apply lowrank matrices for time stepping (can be in-place) >*/


/*------------------------------------------------------------*/
void clr3d_finalize();
/*< free static memory allocated for clr >*/


/*------------------------------------------------------------*/
double walltime( double *t0 );
/*< report runtime >*/


/*------------------------------------------------------------*/
dat3d dat3d_init(int ns_,
                 int nr_,
                 int nc_,
		 sf_axis at_,
		 sf_axis ax_,
                 sf_axis ay_);
/*< init 3D data >*/


/*------------------------------------------------------------*/
mut3d mut3d_make(float t0  /* source delay */,
                 float velw /* water velocity */,
                 float eps /* decay parameter */,
                 pt3d* ss,
                 pt3d* rr,
                 dat3d dat);
/*< init 3D mutting, currently only handles single source >*/


/*------------------------------------------------------------*/
void mut3d_apply(float ***dd,
                 dat3d dat,
                 mut3d mut);
/*< apply in-place muting >*/

#endif
