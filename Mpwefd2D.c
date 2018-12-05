/* Poroelastic time-domain FD modeling, 2D, supports different types of elastic.

   Porolastic wave equation finite difference modeling in both 2D and 3D, using an explicit time-domain solver.

   *** Please see the SConstruct in book/tutorial/ewe for a SConstruct that demonstrates how to use
   predefined functions for using this program. ***

   This program is designed to be as generic as possible, and allows you to use files
   with arbitrary models, and arbitrary source and receiver geometries.  Source types are
   as generic as possible.  Supports arbitrary types of anisotropy as well.

   The downside to the generality, is that the program is not as performant as dedicated solvers
   that are less flexible.  The program is parallelized using OpenMP, so be sure to use a compatible compiler to take
   advantage of the performance boost.
   =========== Rotated Staggered grid ==========================
   Ux,Uz=====================Ux,Uz
   ||            |             ||
   ||                          ||
   ||            |             ||
   ||             tij          ||
   ||- - - - - - |- - - - - - -||
   ||             Cij          ||
   ||            |             ||
   ||                          ||
   ||            |             ||
   Ux,Uz=====================Ux,Uz
   ===========  OPTIONS  =======================================
   ani - The type of anisotropy for this simulation.  Valid options:
   For 2D:
   ISO/HTI/VTI = 0
   TTI = 1

   For 3D:
   ISO/HTI/VTI = 0
   TTI    = 1

   VTI, HTI, and Isotropic media are special cases of ISO/HTI/VTI media.
   TTI media can be represented using TTI media.

   cfl   - Execute the CFL check.  If the CFL check fails, then it will cause the program to fail.
   The CFL check will check both the stability and accuracy conditions for both p-waves and
   s-waves. Depending on the type of anisotropy that you specify, the CFL condition will
   use a safety factor (that you can override if necessary).

   NOTE: the CFL condition will return both minimum and maximum
   constraints on the grid given your velocity model, desired frequency content, and other
   parameters.  IT IS POSSIBLE TO HAVE NO STABLE, AND ACCURATE SOLUTIONS FOR A GIVEN
   MODEL WITH GIVEN PARAMETERS. THE CFL CONDITION WILL WARN YOU IF THIS IS THE CASE.

   YOU MUST SPECIFY fmax Parameter as well!

   ----- STABILITY ------
   The stability condition is related to the maximum wave speed and minimum grid sampling
   as follows:

   dt < min(dx,dy,dz) / (sqrt(2)*vmax)

   Given a time sampling dt, it is possible to determine the minimum dx,dy,dz for stability.
   vmax is the MAXIMUM velocity for all waves in the model (usually P-wave velocity).

   For elastic FD, the P-wave most greatly influences the stability, as it moves fastest
   on the grid.

   The stability condition gives us a LOWER bound on the grid sampling for a given dt.

   ------ ACCURACY -------
   The accuracy condition is related to the number of gridpoints per wavelength.  Thus,

   safety*vmin / fmax > N * sqrt(dx^2+dy^2+dz^2)

   where vmin is the minimum wave velocity in the model (usually S-wave), fmax is some
   relative measure of the maximum frequency of your wavelet (usually 1.5*peak for Ricker),
   N is the number of points desired per wavelength (5), and safety is a safety factor that
   is dependent on the type of anisotropy.

   For elastic FD, the S-wave most greatly impacts the accuracy of the solution, as the S-wave
   is typically much higher frequency and travels at slower wave speeds, meaning shorter
   wavelengths.

   The accuracy condition places an UPPER bound on the grid sampling.

   ---- SAFETY FACTOR -----
   The safety factor depends on the type of anisotropy specified, and attempts to place a lower
   bound on the slowest S-wave velocity (guess):

   ISO/HTI/VTI - (3/4)
   TTI    - (1/2)

   You can also override the safety factor using the safety parameter.

   safety- Override the safety factor for the CFL condition.  This should be a floating point (0-1.0).

   fmax  - An estimate of the highest frequency content in your wavelet (for Ricker use 1.5*peak)

   fsrf  - Use a free surface at the top boundary (z=0).
   WARNING: The free surface condition appears to introduce numerical artifacts into the simulation.
   USE AT YOUR OWN RISK.

   snap  - Save snapshots of the wavefield every n iterations of the modeling program.

   jsnap - Number of iterations between snapshots of the wavefield.
   i.e. jsnap=10, means save a snapshot every 10 iterations.
   If you had 1000 total iterations, then you would have 100 snapshots total.
   The default, will output no snapshots.

   jdata - Number of time imterations between exporting the data at the receivers.
   i.e. jdata=1, means save a snapshot every iteration, which should be the default.
   This can be used to change the sampling of the data to something different from
   the wavelet/wavefield.

   verb  - Print useful information

   debug - Print debugging information.  This is more detailed than verbose.

   srctype - An integer which determines where the source wavelet is injected
   in the simulation.  Valid options are:
   0 - Acceleration source
   1 - Displacement source
   2 - Stress source
   3 - Tensor source
   The default option is 2: Acceleration source.
   For Stress, Displacement and Acceleration sources, your wavelet
   needs to have only 3 components (z,x,y).
   For a Tensor source, you must specify wavelet components for
   all 3 (2D) or 6 (3D) tensor components in the following order:
   2D: tzz, txx, tzx
   3D: tzz, txx, tyy, tyz, tzx, txy

   Hint:  To inject an acoustic source, use a stress source,
   with equal components on all three components.

   dabc  - Use a sponge layer to attenuate waves off the edge of the grid.  Use this in
   combination with the nb parameter.

   abcone- In addition to the sponge layer, using a severe ramp at the very edge of the expanded
   sponge layer to severely attenuate zero-incidence waves at the boundaries.
   It's not clear if this condition actually affects most computations.

   opot  - True: output is second spatial derivative of potentials; False: output wavefield.

   nbell - Size of gaussian used to linearly interpolate curves.  A value of 5 seems to work well.

   nb    - Not listed, but is an important parameter.  Allows you to control the size of the sponge
   layer for the absorbing boundary condition.  If you are getting reflections off the sides,
   with dabc=y, then make this number larger (int).  This pads the grid by this amount on all sides.
   For example:

   |--------------------------|
   |            ramp layer    |
   |r |--------------------|  |
   |a |        nb          |r |
   |m |      |~~~~~~~~|    |a |
   |p |      |  MODEL |    |m |
   |  |  nb  |  SPACE | nb |p |
   |  |      |~~~~~~~~|    |  |
   |  |         nb         |  |
   |  |--------------------|  |
   |         ramp layer       |
   |--------------------------|

   nqz, nqx, oqz, oqx, nqy, oqy, - Allows you to set the parameters for the axes.  Leave as defaults.

   =============BOUNDARY CONDITIONS ========================

   This code enforces a fixed reflecting boundary condition at the
   edge of the computational domain.  The absorbing sponge is used
   IN ADDITION to this condition.

   =============FILE DESCRIPTIONS   ========================

   Fdat.rsf - An RSF file containing your data in the following format:
   axis 1 - source location
   axis 2 - wavefield component (z,x,y) order
   axis 3 - Time

   Fwav.rsf - An RSF file containing your wavelet information.  For elastic modeling, the wavelet needs
   to have 3 samples on N1 one for each component Z-X-Y (or just Z-X for 2D).  The second
   axis describes the component as a function of time.  The sampling interval, origin time,
   and number of time samples will be used as the defaults for the modeling code.
   i.e. your wavelet needs to have the same length and parameters that you want to model with!
   Ex:
   1st axis    index
   Z component  0     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
   X component  1     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
   Y component  2     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
   2nd axis
   NOTE: For tensor sources, you must have an appropriate number of components.  See srctype for more information.

   cccc  - An N+1 dimensional RSF file that contains the values for the stiffness coefficients to be used
   as the model for the modeling code.  So, for 2D, this would be a 3D array of values concatenated
   together in the order as described in the anisotropy section.  Each coefficient file contains
   the value of that coefficient for every point in space.
   The axes for this file are: Axis 1: Z; Axis 2: X; Axis 3: Y;

   The stiffness tensor coefficients are defined uniformly as follows, where
   --x---y---z--(y)-----(y) describes how the coefficients depend on space.
   |C11 C12 C13 C14 C15 C16|
   |    C22 C23 C24 C25 C26|
   |        C33 C34 C35 C36|
   |            C44 C45 C46|
   |                C55 C56|
   |                    C66|

   The tensor is assumed to be symmetric.

   Order of the coefficients in the N+1 dimensional file...
   (First coefficient is the first 2D array in the 3D array).
   2D Anisotropy Modes:

   ISO/HTI/VTI: C11, C33, C55, C13
   "TTI:" C11, C13, C15, C33, C35, C55
   ***TTI basically allows access to all coefs in 2D, but is not really triclinic media
   ------------------------------------------------------------
   (First coefficient is the first 3D array in the 4D array).
   3D Anisotropy Modes:

   ISO/HTI/VTI: C11, C22, C33, C44, C55, C66, C12, C13, C23
   TTI: C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34,
   C35, C36, C44, C45, C46, C55, C56, C66


   sro      - An N dimensional RSF file that contains the values for the solid grain density to be used for the model.
   For 2D, this would be a 2D array.

   fro      - An N dimensional RSF file that contains the values for the fluid density to be used for the model.
   For 2D, this would be a 2D array.

   phi      - An N dimensional RSF file that contains the values for the porosity to be used for the model.
   For 2D, this would be a 2D array.

   kdr      - An N dimensional RSF file that contains the values for the drained bulk modulus to be used for the model.
   For 2D, this would be a 2D array.

   kfl      - An N dimensional RSF file that contains the values for the fluid modulus to be used for the model.
   For 2D, this would be a 2D array.

   ksg      - An N dimensional RSF file that contains the values for the solid grain modulus to be used for the model.
   For 2D, this would be a 2D array.

   shm      - An N dimensional RSF file that contains the values for the shear modulus to be used for the model.
   For 2D, this would be a 2D array.

   tor      - An N dimensional RSF file that contains the values for the tortuosity to be used for the model.
   For 2D, this would be a 2D array.

   sou, rec -The source and receiver RSF files respectively.
   The 1st axis contains the locations for the points like so:
   [x,y,z]
   The second axis is a concatenated list of all points in the list.
   So, for an array of receivers, it would look like:
   [x1,y1,z1]
   [x2,y2,z2]
   [x3,y3,z3]
   [x4,y4,z4]

   wfl     - The name of the file to save the wavefield snapshots to.  This will be an N+2
   dimensional file.  The file will be organized as follows:
   1-2(3) axes, spatial coordinates
   3(4) axis, wavefield components, in the Z,X,(Y) order
   4(5) axis, time, sequential snapshots
   ***The parentheses indicate what the axes will be for 3D models.

   dat     - The name of the file to save the receiver data to.  The data has the format of:
   spatial coordinates, then the data components of the elastic wavefield in the
   same order as the wavefield.  Lastly, time.

   ========== USEFUL COMMANDS  =============================

   To view the wavefield snapshots (2D case):
   sfwindow < Fwfl.rsf n3=1 f3=0 | sfgrey gainpanel=a pclip=100 | sfpen

   To view the data (2D case):
   sfwindow < Fdat.rsf n3=1 f3=0 | sfgrey gainpanel=a pclip=100 | sfpen

   ========== TROUBLESHOOTING ===============================

   If you aren't getting output, or your output is full of Nans, make sure
   that you have the proper dimensions for your wavelet files, and that
   your input files make sense.

   Make sure your source and receiver points are located inside the
   model space, otherwise you will get all NaNs and the simulation will
   run forever.

   ======= TIPS ========

   If the simulation seems to slow down as it's running, its a pretty
   good indication that the simulation has become unstable and is overflowing
   with NaNs.


   Modified by Robert Walker, University of Southern California, 2018-11-29.

*/

/*
  Copyright (C) 2008 Colorado School of Mines

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#include "omputil.h"
#endif

#include "fdutil.h"

#include <math.h>

/* ========================================================================== */
/* FD Coefficients                                                            */
/* ========================================================================== */
#define NOP 4 /* derivative operator half-size */

/* Coefficients are:  (-1)^(k+1)/k * (n!)^2 / (n+k)!(n-k)!
   for n = order of approximation i.e. eighth order
   and k = coefficient number */
/* CENTERED derivatives 8th order */
#define C1 +0.598144530000  /*1.196289060000 */
#define C2 -0.039876302100  /*0.079752604200 */
#define C3 +0.004785156250  /*0.009570312500 */
#define C4 -0.0003487723215 /*0.000697544643*/
#define C5 +0.299072265000  /*1.196289060000/4 */
#define C6 -0.019938151050  /*0.079752604200/4 */
#define C7 +0.002392578125  /*0.009570312500/4 */
#define C8 -0.00017438616075/*0.000697544643/4 */

#define Dx(a,ix,iz,s) (C4*(a[ix+4][iz+4] + a[ix+4][iz-3] -	\
                           a[ix-3][iz+4] - a[ix-3][iz-3]) +	\
		       C3*(a[ix+3][iz+3] + a[ix+3][iz-2] -	\
                           a[ix-2][iz+3] - a[ix-2][iz-2]) +	\
		       C2*(a[ix+2][iz+2] + a[ix+2][iz-1] -	\
                           a[ix-1][iz+2] - a[ix-1][iz-1]) +	\
		       C1*(a[ix+1][iz+1] + a[ix+1][iz  ] -	\
                           a[ix  ][iz+1] - a[ix  ][iz  ])  )*s
#define Dz(a,ix,iz,s) (C4*(a[ix+4][iz+4] - a[ix+4][iz-3] +	\
                           a[ix-3][iz+4] - a[ix-3][iz-3]) +	\
		       C3*(a[ix+3][iz+3] - a[ix+3][iz-2] +	\
                           a[ix-2][iz+3] - a[ix-2][iz-2]) +	\
		       C2*(a[ix+2][iz+2] - a[ix+2][iz-1] +	\
                           a[ix-1][iz+2] - a[ix-1][iz-1]) +	\
		       C1*(a[ix+1][iz+1] - a[ix+1][iz  ] +	\
                           a[ix  ][iz+1] - a[ix  ][iz  ])  )*s

#define Dx3(a,ix,iy,iz,s) (C8*(a[iy+4][ix+4][iz+4] + a[iy-3][ix+4][iz+4] + \
			       a[iy+4][ix+4][iz-3] + a[iy-3][ix+4][iz-3] - \
			       a[iy-3][ix-3][iz-3] - a[iy+4][ix-3][iz-3] - \
			       a[iy-3][ix-3][iz+4] - a[iy+4][ix-3][iz+4]) + \
			   C7*(a[iy+3][ix+3][iz+3] + a[iy-2][ix+3][iz+3] + \
			       a[iy+3][ix+3][iz-2] + a[iy-2][ix+3][iz-2] - \
			       a[iy-2][ix-2][iz-2] - a[iy+3][ix-2][iz-2] - \
			       a[iy-2][ix-2][iz+3] - a[iy+3][ix-2][iz+3]) + \
			   C6*(a[iy+2][ix+2][iz+2] + a[iy-1][ix+2][iz+2] + \
			       a[iy+2][ix+2][iz-1] + a[iy-1][ix+2][iz-1] - \
			       a[iy-1][ix-1][iz-1] - a[iy+2][ix-1][iz-1] - \
			       a[iy-1][ix-1][iz+2] - a[iy+2][ix-1][iz+2]) + \
			   C5*(a[iy+1][ix+1][iz+1] + a[iy  ][ix+1][iz+1] + \
			       a[iy+1][ix+1][iz  ] + a[iy  ][ix+1][iz  ] - \
			       a[iy  ][ix  ][iz  ] - a[iy+1][ix  ][iz  ] - \
			       a[iy  ][ix  ][iz+1] - a[iy+1][ix  ][iz+1])  )*s
#define Dy3(a,ix,iy,iz,s) (C8*(a[iy+4][ix+4][iz+4] + a[iy+4][ix-3][iz-3] + \
			       a[iy+4][ix+4][iz-3] + a[iy+4][ix-3][iz+4] - \
			       a[iy-3][ix-3][iz-3] - a[iy-3][ix+4][iz+4] - \
			       a[iy-3][ix-3][iz+4] - a[iy-3][ix+4][iz-3]) + \
			   C7*(a[iy+3][ix+3][iz+3] + a[iy+3][ix-2][iz-2] + \
			       a[iy+3][ix+3][iz-2] + a[iy+3][ix-2][iz+3] - \
			       a[iy-2][ix-2][iz-2] - a[iy-2][ix+3][iz+3] - \
			       a[iy-2][ix-2][iz+3] - a[iy-2][ix+3][iz-2]) + \
			   C6*(a[iy+2][ix+2][iz+2] + a[iy+2][ix-1][iz-1] + \
			       a[iy+2][ix+2][iz-1] + a[iy+2][ix-1][iz+2] - \
			       a[iy-1][ix-1][iz-1] - a[iy-1][ix+2][iz+2] - \
			       a[iy-1][ix-1][iz+2] - a[iy-1][ix+2][iz-1]) + \
			   C5*(a[iy+1][ix+1][iz+1] + a[iy+1][ix  ][iz  ] + \
			       a[iy+1][ix+1][iz  ] + a[iy+1][ix  ][iz+1] - \
			       a[iy  ][ix  ][iz  ] - a[iy  ][ix+1][iz+1] - \
			       a[iy  ][ix  ][iz+1] - a[iy  ][ix+1][iz  ])  )*s
#define Dz3(a,ix,iy,iz,s) (C8*(a[iy+4][ix+4][iz+4] + a[iy-3][ix+4][iz+4] + \
			       a[iy-3][ix-3][iz+4] + a[iy+4][ix-3][iz+4] - \
			       a[iy-3][ix-3][iz-3] - a[iy+4][ix-3][iz-3] - \
			       a[iy+4][ix+4][iz-3] - a[iy-3][ix+4][iz-3]) + \
			   C7*(a[iy+3][ix+3][iz+3] + a[iy-2][ix+3][iz+3] + \
			       a[iy-2][ix-2][iz+3] + a[iy+3][ix-2][iz+3] - \
			       a[iy-2][ix-2][iz-2] - a[iy+3][ix-2][iz-2] - \
			       a[iy+3][ix+3][iz-2] - a[iy-2][ix+3][iz-2]) + \
			   C6*(a[iy+2][ix+2][iz+2] + a[iy-1][ix+2][iz+2] + \
			       a[iy-1][ix-1][iz+2] + a[iy+2][ix-1][iz+2] - \
			       a[iy-1][ix-1][iz-1] - a[iy+2][ix-1][iz-1] - \
			       a[iy+2][ix+2][iz-1] - a[iy-1][ix+2][iz-1]) + \
			   C5*(a[iy+1][ix+1][iz+1] + a[iy  ][ix+1][iz+1] + \
			       a[iy  ][ix  ][iz+1] + a[iy+1][ix  ][iz+1] - \
			       a[iy  ][ix  ][iz  ] - a[iy+1][ix  ][iz  ] - \
			       a[iy+1][ix+1][iz  ] - a[iy  ][ix+1][iz  ])  )*s


enum AniType {ORTHORHOMBIC=0,TRICLINIC=1};
enum SourceType {STRESS=2, DISPLACEMENT=1, ACCELERATION=0, TENSOR=3};
/* ========================================================================== */
/* Aid Functions                                                              */
/* ========================================================================== */
void zero2d(int n1, int n2, float **array)
{
  int i1, i2;
  for (i2 = 0; i2 < n2; i2++){
    for (i1 = 0; i1 < n1; i1++){
      array[i2][i1] = 0.0;
    }
  }
  return;
}



// /* ========================================================================== */
// /* Implicit Finite Difference                                                 */
// /* ========================================================================== */
// void implicit_spatial_FD2D(float **p, float **q, float h, int L)
// {
//   float *r, *a;
//
//   a = sf_floatalloc1(L+1);
//
// }


/* ========================================================================== */
/* Main Loop                                                                  */
/* ========================================================================== */

int main(int argc, char* argv[])
{
  /* Declare RSF params */
  bool verb,fsrf,snap,dabc,debug,abcone,cfl,is2d,opot;
  int  jsnap,ntsnap,jdata,nbell;
  enum AniType type;
  enum SourceType srctype;
  float fmax,safety;

  /* I/O files */
  sf_file Fwav=NULL; /* wavelet   */
  sf_file Fsou=NULL; /* sources   */
  sf_file Frec=NULL; /* receivers */
  sf_file Fccc=NULL; /* velocity  */
  sf_file Ffro=NULL; /* fluid density   */
  sf_file Fsro=NULL; /* solid grain density   */
  sf_file Fdat=NULL; /* data      */
  sf_file Fwfl=NULL; /* wavefield */
  sf_file Fphi=NULL; /* porosity  */
  sf_file Fkdr=NULL; /* Drained Bulk Modulus */
  sf_file Fkfl=NULL; /* Fluid Bulk Modulus */
  sf_file Fksg=NULL; /* solid grain modulus */
  sf_file Fprm=NULL; /* permeability (tensor, or should be) */
  sf_file Ffvs=NULL; /* fluid viscosity */
  sf_file Fshm=NULL; /* shear modulus */
  sf_file Ftor=NULL; /* tortuosity */

  int     nt,nz,nx,ny=0,ns,nr,nc,nb;
  int     it,iz,ix,iy;
  float   dt,dz,dx,dy=0,idz,idx,idy;

  sf_axis at,ax,ay=NULL,az;
  sf_axis as,ar,ac,acq;

  /* init RSF */
  sf_init(argc,argv);


  int tsrc;
  if (! sf_getint("srctype",&tsrc)) tsrc=0; /* source type, see comments */

  srctype = tsrc;

  /* MAKE SURE USER HAS DECLARED ANISOTROPY TYPE*/
  int ttype;
  if (! sf_getint("ani", &ttype)) ttype=-1;/* Anisotropy type, see comments */
  if (ttype == -1 || ttype >2 || ttype < 0) sf_error("Invalid Anisotropy type");
  type = ttype;

  /* I/O arrays */
  float***ww=NULL;           /* wavelet   */
  float **dd=NULL;           /* data      */

  /*------------------------------------------------------------*/
  /* Temp values for solid strain */
  float    szz,   sxx,   syy,   sxy,   syz,   szx;
  /* Temp values for relative strain */
  float    ezz,   exx,   eyy,   exy,   eyz,   ezx;
  /*------------------------------------------------------------*/
  /* execution flags */
  if (! sf_getbool("verb",&verb)) verb=true; /* verbosity flag */
  if (! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
  if (! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
  if (! sf_getbool("dabc",&dabc)) dabc=true; /* use sponge absorbing BC */
  if (! sf_getbool("abcone",&abcone)) abcone=false; /* use sharp brake at end of boundary layer */
  if (! sf_getbool("debug",&debug)) debug=true; /* print debugging info */
  if (! sf_getbool("cfl",&cfl)) cfl=true; /* use CFL check, will cause program to fail if not satisfied */
  if (! sf_getbool("opot",&opot)) opot=false; /* output potentials */

  /*------------------------------------------------------------*/
  /* I/O files */
  Fwav = sf_input ("in" ); /* wavelet   */
  Fccc = sf_input ("ccc"); /* stiffness */
  Fsro = sf_input ("sro"); /* solid grain density   */
  Ffro = sf_input ("fro"); /* fluid density   */
  Fphi = sf_input ("phi"); /* porosity  */
  Fkdr = sf_input ("kdr"); /* drained bulk modulus */
  Fkfl = sf_input ("kfl"); /* fluid modulus */
  Fksg = sf_input ("ksg"); /* solid grain modulus */
  Ffvs = sf_input ("fvs"); /* fluid viscosity */
  Fshm = sf_input ("shm"); /* shear modulus */
  Ftor = sf_input ("tor"); /* tortuosity */
  Fprm = sf_input ("prm"); /* permeability (tensor?) */
  Fsou = sf_input ("sou"); /* sources   */
  Frec = sf_input ("rec"); /* receivers */
  Fwfl = sf_output("wfl"); /* wavefield */
  Fdat = sf_output("out"); /* data      */

  /* print debugging information */

  /*------------------------------------------------------------*/
  /* Determine if 3D or 2D, test for existence of 4th axis, n=1 if none */
  sf_axis test = sf_iaxa(Fccc,4);
  if (sf_n(test) == 1) is2d = true;
  else is2d = false;

  /* ---------------------------------------------------------------------- */
  /* warn if not using absorbing boundary */
  if (!dabc){
    fprintf(stderr,"*********************\n");
    fprintf(stderr,"NOT USING ABSORBING BOUNDARY CONDITIONS!!!\n");
    fprintf(stderr,"THIS WILL CAUSE LARGE REFLECTIONS AT EDGES\n");
    fprintf(stderr,"TO TURN ON ABSORBING, SET DABC=y\n");
    fprintf(stderr,"*********************\n");
  }

  /* ---------------------------------------------------------------------- */
  /* Make sure that we have the proper number of coefs in our stiffness file */
  sf_axis temp; int tempn;

  if (is2d) {
    temp = sf_iaxa(Fccc,3);
  } else {
    temp = sf_iaxa(Fccc,4);
  }
  tempn = sf_n(temp);
  sf_warning("Stiffness coefficient file has %d coefficients\n",tempn);

  switch (type) {
  case ORTHORHOMBIC:
    if (is2d && tempn != 4)
      sf_error("ISO/HTI/VTI requires 4 coefficients in 2D");
    else if ( !is2d && tempn != 9)
      sf_error("ISO/HTI/VTI requires 9 coefficients in 3D");
    break;
  case TRICLINIC:
    if (is2d && tempn != 6)
      sf_error("TTI requires 6 coefficients in 2D");
    else if ( !is2d && tempn != 21)
      sf_error("TTI requires 21 coefficients in 3D");
    break;
  }
  /*-------------------------------------------------------------*/
  /* Wavefield cut parameters */
  sf_axis   acz=NULL,acx=NULL,acy=NULL;
  int       nqz,nqx,nqy;
  float     oqz,oqx,oqy;
  float     dqz,dqx,dqy;

  /* axes */
  at = sf_iaxa(Fwav,3); sf_setlabel(at,"time");
  sf_setunit(at,"s"); if (verb) sf_raxa(at); /* time */

  az = sf_iaxa(Fccc,1); sf_setlabel(az,"space z");
  sf_setunit(az,"km");if (verb) sf_raxa(az); /* depth */

  ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"space x");
  sf_setunit(ax,"km");if (verb) sf_raxa(ax); /* space x */

  as = sf_iaxa(Fsou,2); sf_setlabel(as,"sources");
  sf_setunit(as,"km");if (verb) sf_raxa(as); /* sources */

  ar = sf_iaxa(Frec,2); sf_setlabel(ar,"receivers");
  sf_setunit(ar,"km");if (verb) sf_raxa(ar); /* receivers */

  /* number and increments */
  nt = sf_n(at); dt = sf_d(at);
  nz = sf_n(az); dz = sf_d(az);
  nx = sf_n(ax); dx = sf_d(ax);

  /* check if model is 2D */
  if (!is2d){
    ay = sf_iaxa(Fccc,3); sf_setlabel(ay,"space y");
    sf_setunit(ay,"km");if (verb) sf_raxa(ay); /* space y */
    ny = sf_n(ay); dy = sf_d(ay);
  }

  /* number of sources and receivers */
  ns = sf_n(as);
  nr = sf_n(ar);

  /* if wavefield snapshots are requested */
  if (snap){
    if (!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
    if (!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);

    if (!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
    if (!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);

    if (!is2d){
      if (!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay);
      if (!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay);
    }
  }
  /*------------------------------------------------------------*/
  /* other execution parameters */
  if (! sf_getint("nbell",&nbell)) nbell=1;  /* bell size */
  if (verb) sf_warning("nbell=%d",nbell);
  if (! sf_getint("jdata",&jdata)) jdata=1;

  if ( !sf_getint("nb",&nb)) nb=100; /* padding size for absorbing boundary */
  if (nb < NOP) nb = NOP;

  if (snap) {  /* save wavefield every *jsnap* time steps */
    if (! sf_getint("jsnap",&jsnap)) jsnap=nt;
  }
  /* CFL check parameters */
  if (cfl){
    if (! sf_getfloat("fmax",&fmax)) sf_error("CFL: you did not specify fmax");

    if (! sf_getfloat("safety",&safety) || safety < 0) {
      switch (type){
      case ORTHORHOMBIC:
	safety = 0.75f;
	break;
      case TRICLINIC:
	safety = 0.50f;
	break;
      }
    }
    sf_warning("CFL: safety margin %f",safety);
  }
  /* ====================================================================== */
  /* --------------2d code---------------- */
  /* IO Arrays */
  pt2d   *ss=NULL;           /* sources   */
  pt2d   *rr=NULL;           /* receivers */
  /* FDM structure */
  fdm2d    fdm=NULL;
  abcone2d abcp=NULL,abcs=NULL;
  sponge   spo=NULL;

  /*-------------- material parameters -------------------------------------*/
  float **tt=NULL;            /* work array */
  float **ro=NULL;            /* ro term for ABCs */
  float **fro=NULL;           /* fluid density   */
  float **sro=NULL;           /* solid grain density */
  float **bro=NULL;           /* bulk density   */
  float **mro=NULL;           /* momentum added density */
  float **kdr=NULL;           /* drained bulk modulus */
  float **kud=NULL;           /* undrained bulk modulus */
  float **kfl=NULL;           /* fluid modulus */
  float **ksg=NULL;           /* solid grain modulus */
  float **fud=NULL;           /* undrained bulk modulus */
  float **fvs=NULL;           /* fluid viscosity */
  float **alpha=NULL;         /* biot coefficient */
  float **biotmod=NULL;       /* biod modulus */
  float **skem=NULL;          /* skempton coefficient */
  float **lambdau=NULL;       /* undrained lame #1 */
  float **shm=NULL;           /* shear modulus */
  float **fmo=NULL;           /* fluid mobility */
  float **phi=NULL;           /* porosity */
  float **tor=NULL;           /* tortuosity */
  float **prm=NULL;           /* permeability (tensor?) */
  float **c11=NULL;           /* stiffness components */
  float **c33=NULL;
  float **c55=NULL;
  float **c13=NULL;
  float **c15=NULL;
  float **c35=NULL;
  float **r_bar=NULL;
  float **vp,**vs;
  float **qp=NULL,**qs=NULL;
  float **coeff_A=NULL;       /* FD Coefficients */
  float **coeff_B=NULL;
  float **coeff_C=NULL;
  float **coeff_D=NULL;
  float **coeff_E=NULL;
  float **coeff_F=NULL;
  float **coeff_G=NULL;
  float **coeff_H=NULL;
  float **coeff_I=NULL;
  float **coeff_J=NULL;
  float **coeff_K=NULL;

  /*------------------------------------------------------------*/
  /* particle displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
  float **umz,**uoz,**upz,**uaz,**utz;
  float **umx,**uox,**upx,**uax,**utx;
  /* particle velocity: vm = V @ t-1; vo = V @ t; vp = V @ t+1 */
  float **vmz,**voz,**vpz,**vaz,**vtz;
  float **vmx,**vox,**vpx,**vax,**vtx;
  /* filtration displacement: wm = W @ t-1; wo = W @ t; wp = W @ t+1 */
  float **wmz,**woz,**wpz,**waz,**wtz;
  float **wmx,**wox,**wpx,**wax,**wtx;
  /* filtration velocity: qm = Q @ t-1; uo = Q @ t; up = Q @ t+1 */
  float **qmz,**qoz,**qpz,**qaz,**qtz;
  float **qmx,**qox,**qpx,**qax,**qtx;
  /* */
  /* stress/strain tensor */
  float **tzz,**tzx,**txx;
  /* pore pressure */
  float **p;
  /*------------------------------------------------------------*/
  /* linear interpolation weights/indices */
  lint2d cs,cr;
  /* wavefield cut params */
  float     **uc=NULL;
  /*------------------------------------------------------------*/
  /* expand domain for FD operators and ABC */
  fdm=fdutil_init(verb,fsrf,az,ax,nb,1);
  fdbell_init(nbell);

  sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if (verb) sf_raxa(az);
  sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if (verb) sf_raxa(ax);
  /*------------------------------------------------------------*/
  /* 2D vector components */
  nc=2;
  ac=sf_maxa(nc,0,1);
  /*------------------------------------------------------------*/
  /* setup output data header */
  sf_oaxa(Fdat,ar,1);
  sf_oaxa(Fdat,ac,2);
  sf_setn(at,nt/jdata);
  sf_setd(at,dt*jdata);
  sf_oaxa(Fdat,at,3);
  /* setup output wavefield header */
  if (snap) {
    dqz=sf_d(az);
    dqx=sf_d(ax);

    acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
    acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
    /* TODO: check if the imaging window fits in the wavefield domain */

    uc=sf_floatalloc2(sf_n(acz),sf_n(acx));

    ntsnap=0;
    for (it=0; it<nt; it++) {
      if (it%jsnap==0) ntsnap++;
    }
    sf_setn(at,  ntsnap);
    sf_setd(at,dt*jsnap);
    if (verb) sf_raxa(at);

    sf_oaxa(Fwfl,acz,1);
    sf_oaxa(Fwfl,acx,2);
    sf_oaxa(Fwfl,ac, 3);
    sf_oaxa(Fwfl,at, 4);
  }

  /*------------------------------------------------------------*/
  /* source array */
  if (srctype == TENSOR) {
    ww=sf_floatalloc3(ns,3,nt);
    sf_floatread(ww[0][0],nt*3*ns,Fwav);
  } else {
    ww=sf_floatalloc3(ns,nc,nt);
    sf_floatread(ww[0][0],nt*nc*ns,Fwav);
  }

  /* data array */
  dd = sf_floatalloc2(nr,nc);
  /*------------------------------------------------------------*/
  /* setup source/receiver coordinates */
  ss = (pt2d*) sf_alloc(ns,sizeof(*ss));
  rr = (pt2d*) sf_alloc(nr,sizeof(*rr));

  /* read in source/receiver data */
  pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
  pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

  cs = lint2d_make(ns,ss,fdm); /* setup linear interpolation */
  cr = lint2d_make(nr,rr,fdm);

  /*------------------------------------------------------------*/
  /* setup FD coefficients */

  idz = 1/dz;
  idx = 1/dx;

  /* ---------------allocate space for material parameters ---------------- */
  tt = sf_floatalloc2(nz,nx);

  fro =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input fluid density */
  sro =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input solid grain density */
  kdr =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input drained bulk modulus */
  kfl =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input fluid bulk modulus */
  ksg =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input solid grain bulk modulus */
  fvs =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input fluid viscosity */
  shm =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input shear modulus */
  tor =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input tortuosity */
  prm =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input permeability (tensor?) */
  phi =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input porosity */

  /* --------------read in material parameters ---------------------------- */
  sf_floatread(tt[0],nz*nx,Ffro);     expand(tt,fro ,fdm);
  sf_floatread(tt[0],nz*nx,Fsro);     expand(tt,sro ,fdm);
  sf_floatread(tt[0],nz*nx,Fkdr);     expand(tt,kdr ,fdm);
  sf_floatread(tt[0],nz*nx,Fkfl);     expand(tt,kfl ,fdm);
  sf_floatread(tt[0],nz*nx,Fksg);     expand(tt,ksg ,fdm);
  sf_floatread(tt[0],nz*nx,Ffvs);     expand(tt,fvs ,fdm);
  sf_floatread(tt[0],nz*nx,Fshm);     expand(tt,shm ,fdm);
  sf_floatread(tt[0],nz*nx,Ftor);     expand(tt,tor ,fdm);
  sf_floatread(tt[0],nz*nx,Fprm);     expand(tt,tor ,fdm);
  sf_floatread(tt[0],nz*nx,Fphi);     expand(tt,tor ,fdm);

  /* ---------------allocate space for derived parameters------------------ */
  kud     =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* undrained bulk modulus */
  alpha   =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* biot coefficient */
  biotmod =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* biot modulus */
  bro     =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* bulk density */
  lambdau =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* undrained lame #1 */
  ro      =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* ro term for ABCs */
  mro     =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  r_bar   =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  fmo     =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* fluid mobility */

  /* --------------allocate space for FD Coefficients---------------------- */
  coeff_A =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_B =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_C =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_D =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_E =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_F =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_G =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_H =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_I =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_J =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  coeff_K =sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  /* ---------------load stiffness tensor---------------------------------- */
  switch (type){
  case ORTHORHOMBIC:
    c11=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    c33=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    c55=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    c13=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    /* input stiffness */
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c11,fdm);
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c33,fdm);
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c55,fdm);
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c13,fdm);
    break;
  case TRICLINIC:
    c11=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    c33=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    c55=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    c13=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    c15=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    c35=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    /* input stiffness */
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c11,fdm);
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c13,fdm);
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c15,fdm);
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c33,fdm);
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c35,fdm);
    sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c55,fdm);
    break;
  }
  /* free work vector */
  free(*tt); free(tt);

  /*------------------------------------------------------------*/
  /* calculate local parameters */
  for    (ix = NOP; ix < fdm->nxpad - NOP; ix++) {
    for    (iz = NOP; iz < fdm->nzpad - NOP; iz++) {
      /* bulk density */
      bro[ix][iz] = phi[ix][iz]*fro[ix][iz] + (1.0 - phi[ix][iz])*sro[ix][iz];
      /* biot coefficient */
      alpha[ix][iz] = 1.0 - kdr[ix][iz]/ksg[ix][iz];
      /* skempton coefficient */
      skem[ix][iz] = (1.0/kdr[ix][iz] - 1.0/ksg[ix][iz])/(1.0/kdr[ix][iz] - 1.0/ksg[ix][iz] + phi[ix][iz]*(1.0/kfl[ix][iz] - 1.0/ksg[ix][iz]));
      /* undrained bulk modulus */
      kud[ix][iz] = kdr[ix][iz] / (1.0 - skem[ix][iz]*(1.0 - kdr[ix][iz]/ksg[ix][iz]));
      /* biot modulus */
      biotmod[ix][iz] = skem[ix][iz] * kud[ix][iz]/alpha[ix][iz];
      /* undrained lame # 1 */
      lambdau[ix][iz] = kud[ix][iz] - 2.0*shm[ix][iz]/3.0;
      /* Poroelastodynamic Values */
      mro[ix][iz] = tor[ix][iz]*fro[ix][iz] / phi[ix][iz];
      r_bar[ix][iz] = mro[ix][iz]*bro[ix][iz] - fro[ix][iz]*fro[ix][iz];
      fmo[ix][iz] = fvs[ix][iz]/prm[ix][iz];
      /* ro value for ABCs */
      ro[ix][iz] = dt*dt / r_bar[ix][iz];
      //ro[ix][iz] = dt*dt/bro[ix][iz];
      /* Poroelastic FD Coefficients */
      coeff_A[ix][iz] = (mro[ix][iz]*dt)/r_bar[ix][iz];
      coeff_B[ix][iz] = (fro[ix][iz]*dt)/r_bar[ix][iz];
      coeff_C[ix][iz] = (fro[ix][iz]*fmo[ix][iz]*dt)/r_bar[ix][iz];
      coeff_D[ix][iz] = (2.0*r_bar[ix][iz] - bro[ix][iz]*fmo[ix][iz]*dt) / (2.0*r_bar[ix][iz] + bro[ix][iz]*fmo[ix][iz]*dt);
      coeff_E[ix][iz] = -1.0*(bro[ix][iz]*dt) / (r_bar[ix][iz] + bro[ix][iz]*dt);
      coeff_F[ix][iz] = -1.0*(bro[ix][iz]*dt) / (r_bar[ix][iz] + bro[ix][iz]*dt);
      coeff_G[ix][iz] = shm[ix][iz]*dt;
      coeff_H[ix][iz] = lambdau[ix][iz]*dt;
      coeff_I[ix][iz] = alpha[ix][iz]*biotmod[ix][iz]*dt;
      coeff_J[ix][iz] = -1.0*biotmod[ix][iz]*dt;
    }
  }

  /* ------------------------------------------------------------------------ */
  /* one-way abc setup   */
  vp = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  vs = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  float vpmax = 0.0f;
  float vpmin = 100000.0f;
  float vsmax = 0.0f;
  float vsmin = 100000.0f;
  float c11f;
  float c55f;
  for    (ix=0; ix<fdm->nxpad; ix++) {
    for  (iz=0; iz<fdm->nzpad; iz++) {
      c11f = c11[ix][iz];
      c55f = c55[ix][iz];
      if (c11f < 0 ) {
	vp[ix][iz] = sqrt( -c11f/ro[ix][iz] );
      } else {
	vp[ix][iz] = sqrt( c11f/ro[ix][iz] );
      }
      if (vp[ix][iz] > vpmax) vpmax = vp[ix][iz];
      if (vp[ix][iz] < vpmin) vpmin = vp[ix][iz];
      if (c55f < 0 ) {
	vs[ix][iz] = sqrt(-c55f/ro[ix][iz] );
      } else {
	vs[ix][iz] = sqrt( c55f/ro[ix][iz] );
      }
      if ( vs[ix][iz] > vsmax) vsmax = vs[ix][iz];
      if ( vs[ix][iz] < vsmin) vsmin = vs[ix][iz];
    }
  }
  /* CFL condition */
  if (cfl) {
    cfl_elastic( vpmin,vpmax,vsmin,vsmax,
		 dx,-1.0,dz,dt,fmax,safety,4);
  }
  /* Absorbing boundary conditions setup */
  if (abcone) {
    abcp = abcone2d_make(NOP,dt,vp,fsrf,fdm);
    abcs = abcone2d_make(NOP,dt,vs,fsrf,fdm);
  }
  free(*vp); free(vp);
  free(*vs); free(vs);
  /* sponge abc setup */
  spo = sponge_make(fdm->nb);

  /*------------------------------------------------------------*/
  /* allocate wavefield arrays */
  /* -----------------------------------------------------------*/
  /* particle displacement */
  umz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  uoz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  upz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  uaz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  umx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  uox=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  upx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  uax=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  /* particle velocity */
  vmx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  vox=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  vpx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  vax=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  vmz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  voz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  vpz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  vaz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  /* filtration displacement */
  wmz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  woz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  wpz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  waz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  wmx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  wox=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  wpx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  wax=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  /* filtration velocity */
  qmx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  qox=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  qpx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  qax=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  qmz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  qoz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  qpz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  qaz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  /* stress tensor */
  tzz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  tzx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  txx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  /* pore pressure */
  p  =sf_floatalloc2(fdm->nzpad,fdm->nxpad);

  /* Zero out values */

  zero2d(fdm->nxpad, fdm->nzpad, umz);
  zero2d(fdm->nxpad, fdm->nzpad, uoz);
  zero2d(fdm->nxpad, fdm->nzpad, upz);
  zero2d(fdm->nxpad, fdm->nzpad, uaz);

  zero2d(fdm->nxpad, fdm->nzpad, umx);
  zero2d(fdm->nxpad, fdm->nzpad, uox);
  zero2d(fdm->nxpad, fdm->nzpad, upx);
  zero2d(fdm->nxpad, fdm->nzpad, uaz);

  zero2d(fdm->nxpad, fdm->nzpad, wmz);
  zero2d(fdm->nxpad, fdm->nzpad, woz);
  zero2d(fdm->nxpad, fdm->nzpad, wpz);
  zero2d(fdm->nxpad, fdm->nzpad, waz);

  zero2d(fdm->nxpad, fdm->nzpad, wmx);
  zero2d(fdm->nxpad, fdm->nzpad, wox);
  zero2d(fdm->nxpad, fdm->nzpad, wpx);
  zero2d(fdm->nxpad, fdm->nzpad, waz);

  zero2d(fdm->nxpad, fdm->nzpad, vmz);
  zero2d(fdm->nxpad, fdm->nzpad, voz);
  zero2d(fdm->nxpad, fdm->nzpad, vpz);
  zero2d(fdm->nxpad, fdm->nzpad, vaz);

  zero2d(fdm->nxpad, fdm->nzpad, vmx);
  zero2d(fdm->nxpad, fdm->nzpad, vox);
  zero2d(fdm->nxpad, fdm->nzpad, vpx);
  zero2d(fdm->nxpad, fdm->nzpad, vaz);

  zero2d(fdm->nxpad, fdm->nzpad, qmz);
  zero2d(fdm->nxpad, fdm->nzpad, qoz);
  zero2d(fdm->nxpad, fdm->nzpad, qpz);
  zero2d(fdm->nxpad, fdm->nzpad, qaz);

  zero2d(fdm->nxpad, fdm->nzpad, qmx);
  zero2d(fdm->nxpad, fdm->nzpad, qox);
  zero2d(fdm->nxpad, fdm->nzpad, qpx);
  zero2d(fdm->nxpad, fdm->nzpad, qaz);

  zero2d(fdm->nxpad, fdm->nzpad, tzz);
  zero2d(fdm->nxpad, fdm->nzpad, tzx);
  zero2d(fdm->nxpad, fdm->nzpad, txx);

  zero2d(fdm->nxpad, fdm->nzpad,   p);

  if (opot) {
    qp=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    qs=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
  }

  /* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  /*------------------------------------------------------------*/
  /*
   *  MAIN LOOP - actually run stuff
   */
  /*------------------------------------------------------------*/
  if (verb) fprintf(stderr,"\n");
  for (it=0; it<nt; it++) {
    if (verb) fprintf(stderr,"%d/%d \r",it,nt);

    /*------------------------------------------------------------*/
    /* from displacement to strain and to stress  at t = n        */
    /*------------------------------------------------------------*/
    for  (ix=NOP; ix<fdm->nxpad-NOP-1; ix++) {
      for (iz=NOP; iz<fdm->nzpad-NOP-1; iz++) {
        txx[ix][iz] = (lambdau[ix][iz] + 2.0*shm[ix][iz])*Dx(uox,ix,iz,idx) + lambdau[ix][iz]*Dz(uoz,ix,iz,idz)
                      + alpha[ix][iz]*biotmod[ix][iz]*(Dx(wox,ix,iz,idx) + Dz(woz,ix,iz,idz));
        tzz[ix][iz] = (lambdau[ix][iz] + 2.0*shm[ix][iz])*Dz(uoz,ix,iz,idz) + lambdau[ix][iz]*Dx(uox,ix,iz,idx)
                      + alpha[ix][iz]*biotmod[ix][iz]*(Dx(wox,ix,iz,idx) + Dz(woz,ix,iz,idz));
        p[ix][iz]   = -1.0*alpha[ix][iz]*biotmod[ix][iz]*(Dx(uox,ix,iz,idx) + Dz(uoz,ix,iz,idz)) - biotmod[ix][iz]*(Dx(wox,ix,iz,idx) + Dz(woz,ix,iz,idz));
      }
    }
    /*------------------------------------------------------------*/
    /* free surface */
    /*------------------------------------------------------------*/
    /* Francesco: the z component of the traction must be zero at the free surface */
    if (fsrf) {
      for (ix=0; ix<fdm->nxpad; ix++) {
	       for (iz=0; iz<fdm->nb; iz++) {
	          txx[ix][iz]=0;
	          tzz[ix][iz]=0;
	          tzx[ix][iz]=0;
              p[ix][iz]=0;
          }
      }
      for (ix=0; ix<fdm->nxpad; ix++) {
	       tzz[ix][iz]=0.0;
	      }
    }

	/*------------------------------------------------------------*/
	/* inject stress source                                       */
	/*------------------------------------------------------------*/
	if (srctype == STRESS || srctype == TENSOR) {
	  lint2d_bell(tzz,ww[it][0],cs);
	  lint2d_bell(txx,ww[it][1],cs);
	}
	if (srctype == TENSOR){
	  lint2d_bell(tzx,ww[it][2],cs);
	}
	if (dabc){
	  sponge2d_apply(txx,spo,fdm);
	  sponge2d_apply(tzx,spo,fdm);
	  sponge2d_apply(tzz,spo,fdm);
    sponge2d_apply(p, spo, fdm);
	}

	/*------------------------------------------------------------*/
	/* from stress to acceleration  t = n + 1/2                   */
	/*------------------------------------------------------------*/
	/*
	 * ax = Bx(txx) + Fz(txz)
	 * az = Fx(txz) + Bz(tzz)
	 */
	// #ifdef _OPENMP
	// #pragma omp parallel for			\
	//   schedule(dynamic,fdm->ompchunk)		\
	//   private(iz,ix)				\
	//   shared(fdm,tzz,tzx,txx,uaz,uax,idz,idx)
	// #endif
  for (ix = NOP; ix<fdm->nxpad-NOP; ix++) {
    for (iz = NOP; iz<fdm->nzpad-NOP; iz++) {
      /* filtration velocity, t = n + 1/2 */
      qpz[ix][iz] = coeff_D[ix][iz]*qoz[ix][iz] + coeff_E[ix][iz]*(Dx(tzx,ix,iz,idx) + Dz(tzz,ix,iz,idz)) + coeff_F[ix][iz]*Dz(p,ix,iz,idz);
      qpx[ix][iz] = coeff_D[ix][iz]*qox[ix][iz] + coeff_E[ix][iz]*(Dx(txx,ix,iz,idx) + Dz(tzx,ix,iz,idz)) + coeff_F[ix][iz]*Dx(p,ix,iz,idx);

      /* particle velocity, t = n + 1/2 */
      vpz[ix][iz] = voz[ix][iz] + coeff_A[ix][iz]*(Dx(tzx,ix,iz,idx) + Dz(tzz,ix,iz,idz)) + coeff_B[ix][iz]*Dz(p,ix,iz,idz) + coeff_C[ix][iz]*(qpz[ix][iz] + qoz[ix][iz])*0.5;
      vpx[ix][iz] = vox[ix][iz] + coeff_A[ix][iz]*(Dx(txx,ix,iz,idx) + Dz(tzx,ix,iz,idz)) + coeff_B[ix][iz]*Dx(p,ix,iz,idx) + coeff_C[ix][iz]*(qpx[ix][iz] + qox[ix][iz])*0.5;

      /* particle acceleration, t = n + 1/2 */
      uaz[ix][iz] = (mro[ix][iz]*(Dx(tzx,ix,iz,idx) + Dz(tzz,ix,iz,idz)) + fro[ix][iz]*Dz(p,ix,iz,idz) + fro[ix][iz]*fmo[ix][iz]*qpz[ix][iz])/r_bar[ix][iz];
      uax[ix][iz] = (mro[ix][iz]*(Dx(txx,ix,iz,idx) + Dz(tzx,ix,iz,idz)) + fro[ix][iz]*Dx(p,ix,iz,idx) + fro[ix][iz]*fmo[ix][iz]*qpx[ix][iz])/r_bar[ix][iz];

      /* filtration acceleration, t = n + 1/2 */
      waz[ix][iz] = (-1.0*fro[ix][iz]*(Dx(tzx,ix,iz,idx) + Dz(tzz,ix,iz,idz)) + bro[ix][iz]*Dz(p,ix,iz,idz) + bro[ix][iz]*fmo[ix][iz]*qpz[ix][iz])/r_bar[ix][iz];
      wax[ix][iz] = (-1.0*fro[ix][iz]*(Dx(txx,ix,iz,idx) + Dz(tzx,ix,iz,idz)) + bro[ix][iz]*Dx(p,ix,iz,idx) + bro[ix][iz]*fmo[ix][iz]*qpx[ix][iz])/r_bar[ix][iz];
    }
  }

	/*------------------------------------------------------------*/
	/* inject ACCELERATION source                                 */
	/*------------------------------------------------------------*/
	if (srctype == ACCELERATION) {
	  lint2d_bell(uaz,ww[it][0],cs);
	  lint2d_bell(uax,ww[it][1],cs);
	}

	/*------------------------------------------------------------*/
	/* step forward in time                                       */
	/*------------------------------------------------------------*/
	// #ifdef _OPENMP
	// #pragma omp parallel for				\
	//   schedule(dynamic,fdm->ompchunk)			\
	//   private(iz,ix)					\
	//   shared(fdm,uoz,uox,umz,umx,upz,upx,uaz,uax,ro)
	// #endif
	for   (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	  for  (iz=NOP; iz<fdm->nzpad-NOP; iz++) {
      /* filtration displacement, t = n + 1 */
      wpz[ix][iz] = 2.0*woz[ix][iz] - wmz[ix][iz] + dt*dt*waz[ix][iz];
      wpx[ix][iz] = 2.0*wox[ix][iz] - wmx[ix][iz] + dt*dt*wax[ix][iz];

      /* particle displacement, t = n + 1 */
      upz[ix][iz] = 2.0*uoz[ix][iz] - umz[ix][iz] + dt*dt*uaz[ix][iz];
      upx[ix][iz] = 2.0*uox[ix][iz] - umx[ix][iz] + dt*dt*uax[ix][iz];
	  }
	}

  /*------------------------------------------------------------*/
	/* inject DISPLACEMENT source                                 */
	/*------------------------------------------------------------*/
	if (srctype == DISPLACEMENT) {
	  lint2d_bell(upz,ww[it][0],cs);
	  lint2d_bell(upx,ww[it][1],cs);
	}
	/*------------------------------------------------------------*/
	/* apply the boundary condition                               */
	/* undrained, fixed all around                                */
	/*------------------------------------------------------------*/
	for (ix = 0; ix < fdm->nxpad; ix++){
	  for (iz = 0; iz < NOP; iz++){
	    upz[ix][iz] = 0.0f;
	    upx[ix][iz] = 0.0f;
      wpz[ix][iz] = 0.0f;
      wpx[ix][iz] = 0.0f;
	  }
	}
	for (ix = 0; ix < fdm->nxpad; ix++){
	  for (iz = fdm->nzpad-NOP; iz < fdm->nzpad; iz++){
      upz[ix][iz] = 0.0f;
	    upx[ix][iz] = 0.0f;
      wpz[ix][iz] = 0.0f;
      wpx[ix][iz] = 0.0f;
	  }
	}
	for (ix = 0; ix < NOP; ix++){
	  for (iz = 0; iz < fdm->nzpad; iz++){
      upz[ix][iz] = 0.0f;
	    upx[ix][iz] = 0.0f;
      wpz[ix][iz] = 0.0f;
      wpx[ix][iz] = 0.0f;
	  }
	}
	for (ix = fdm->nxpad-NOP; ix < fdm->nxpad; ix++){
	  for (iz = 0; iz < fdm->nzpad; iz++){
      upz[ix][iz] = 0.0f;
	    upx[ix][iz] = 0.0f;
      wpz[ix][iz] = 0.0f;
      wpx[ix][iz] = 0.0f;
	  }
	}
	/* circulate wavefield arrays */
	/* Change pointers around */

  /* particle displacement */
	utz=umz; utx=umx;
	umz=uoz; umx=uox;
	uoz=upz; uox=upx;
	upz=utz; upx=utx;

  /* relative displacement */
  wtz=wmz; wtx=wmx;
  wmz=woz; wmx=wox;
  woz=wpz; wox=wpx;
  wpz=wtz; wpx=wtx;

  /* partcle velocity */
  vtz=vmz; vtx=vmx;
	vmz=voz; vmx=vox;
	voz=vpz; vox=vpx;
	vpz=vtz; vpx=vtx;

  /* relative velocity */
  qtz=qmz; qtx=qmx;
  qmz=qoz; qmx=qox;
  qoz=qpz; qox=qpx;
  qpz=qtz; qpx=qtx;

	/* Apply the zero-incidence boundary condition */
	if (abcone){
	  abcone2d_apply(uoz,umz,NOP,abcp,fdm);
	  abcone2d_apply(uox,umx,NOP,abcp,fdm);

    abcone2d_apply(woz,wmz,NOP,abcp,fdm);
	  abcone2d_apply(wox,wmx,NOP,abcp,fdm);

	  abcone2d_apply(uoz,umz,NOP,abcs,fdm);
	  abcone2d_apply(uox,umx,NOP,abcs,fdm);

    abcone2d_apply(woz,wmz,NOP,abcs,fdm);
	  abcone2d_apply(wox,wmx,NOP,abcs,fdm);
  }
	/* sponge ABC */
	if (dabc){
	  sponge2d_apply(umz,spo,fdm);
	  sponge2d_apply(uoz,spo,fdm);
	  sponge2d_apply(upz,spo,fdm);

    sponge2d_apply(wmz,spo,fdm);
	  sponge2d_apply(woz,spo,fdm);
	  sponge2d_apply(wpz,spo,fdm);

	  sponge2d_apply(umx,spo,fdm);
	  sponge2d_apply(uox,spo,fdm);
	  sponge2d_apply(upx,spo,fdm);

    sponge2d_apply(wmx,spo,fdm);
	  sponge2d_apply(wox,spo,fdm);
	  sponge2d_apply(wpx,spo,fdm);
	}
	/*------------------------------------------------------------*/
	/* cut wavefield and save */
	/*------------------------------------------------------------*/
	lint2d_extract(uoz,dd[0],cr);
	lint2d_extract(uox,dd[1],cr);

	if (opot){
	  if (snap && it%jsnap==0) {
	    // #ifdef _OPENMP
	    // #pragma omp parallel for			\
	    //   schedule(dynamic,fdm->ompchunk)		\
	    //   private(iz,ix)				\
	    //   shared(fdm,uoz,uox,qp,qs,idz,idx)
	    // #endif

	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	      for (iz=NOP; iz<fdm->nzpad-NOP; iz++) {
          qp[ix][iz] = Dz( uoz,ix,iz,idz ) + Dx( uox,ix,iz,idx );
          qs[ix][iz] = Dz( uox,ix,iz,idz ) - Dx( uoz,ix,iz,idx );
        }
	    }
	    cut2d(qp,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);

	    cut2d(qs,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	  }
	  if (it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);

	}
  else
  {
	  if (snap && it%jsnap==0) {
	    cut2d(uoz,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);

	    cut2d(uox,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	  }
	  if (it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
	}
	/* Save snapshots of wavefield if selected */
	/* Save data snapshots*/
  }

  if (verb) fprintf(stderr,"\n");
  /* ------------------------------------------------------------ */
  /* deallocate arrays */
  /* ------------------------------------------------------------ */

  free(**ww); free(*ww); free(ww);
  free(ss); free(rr);
  free(*dd);  free(dd);
  free(*ro);  free(ro);
  /* Only free double pointers that are not null... */

  switch (type){
  case TRICLINIC:
   free(*c15); free(c15);
   free(*c35); free(c35); //flow through
  case ORTHORHOMBIC:
   free(*c11); free(c11);
   free(*c33); free(c33);
   free(*c55); free(c55);
   free(*c13); free(c13);
   break;
  }
  /* Free all outstanding parameters */
  free(*umz); free(umz);
  free(*uoz); free(uoz);
  free(*upz); free(upz);
  free(*uaz); free(uaz);

  free(*umx); free(umx);
  free(*uox); free(uox);
  free(*upx); free(upx);
  free(*uax); free(uax);

  free(*vmz); free(vmz);
  free(*voz); free(voz);
  free(*vpz); free(vpz);
  free(*vaz); free(vaz);

  free(*vmx); free(vmx);
  free(*vox); free(vox);
  free(*vpx); free(vpx);
  free(*vax); free(vax);

  free(*wmz); free(wmz);
  free(*woz); free(woz);
  free(*wpz); free(wpz);
  free(*waz); free(waz);

  free(*wmx); free(wmx);
  free(*wox); free(wox);
  free(*wpx); free(wpx);
  free(*wax); free(wax);

  free(*vmz); free(qmz);
  free(*voz); free(qoz);
  free(*vpz); free(qpz);
  free(*vaz); free(qaz);

  free(*qmx); free(qmx);
  free(*qox); free(qox);
  free(*qpx); free(qpx);
  free(*qax); free(qax);

  free(*tzz); free(tzz);
  free(*txx); free(txx);
  free(*tzx); free(tzx);

  free(*p);   free(p);

  if (snap){free(*uc);  free(uc);    }
  if (opot) {
    free(*qp);  free(qp);
	  free(*qs);  free(qs);
  }

  exit (0);
/* ------------end 2d ------------------ */
}
