###############################################################
# Created by Robert Walker
# Dec. 2018
# 2D poroelastic finite difference modeling
# Test script
###############################################################
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------
par = {
    'nx':141,  'ox':0, 'dx':0.5,      'lx':'x', 'ux':'m',
    'ny':621,  'oy':0, 'dy':0.5,      'ly':'y', 'uy':'m',
    'nz':621,  'oz':0, 'dz':0.5,      'lz':'z', 'uz':'m',
    'nt':300,   'ot':0, 'dt':1.0e-4,   'lt':'t', 'ut':'s',
    'kt':15,
    'jsnap':100,
    'height':10,
    'nb':20,
    'frq':80.0, # source wavelet frequency
    'nbell':5
    }
#fdmod.param(par)
par['nframe']=5
par['iframe']=4
par['dabc']='y'
#-------------------------------------------------------------
def arr2str(array,sep=' '):
    return sep.join(map(str,array))
#-------------------------------------------------------------
# Create Layers
layers = ((0.01,0.01,0.01,0.01,0.01),
      	  (0.15,0.15,0.15,0.15,0.15),
	      (0.25,0.25,0.25,0.25,0.25))

n1 = len(layers[0])
n2 = len(layers)

layers1 = ' '.join(map(arr2str,layers))

#Flow('layers',None,
#     '''
#     echo %s
#     n1=%g n2=%g o1=%g d1=%g
#     data_format=ascii_float in=$TARGET
#     ''' % (layers1,n1,n2,par['ox'],par['dx'] ))

# ------------------------------------------------------------
# poroelastic parameters

# solid grain density, kg/m**3
sro =(  2250.0,
	      2250.0,
	      2588.0,
	      2588.0)
sro1= arr2str(sro,',')

# fluid density, kg/m**3
fro =(  1040.0,
	      1040.0,
	      952.4,
	      952.4)
fro1= arr2str(fro,',')

# drained bulk modulus, Pa [ kg / (m * s**2) ]
kdr =(  6.21e08,
	      6.21e08,
	      6.21e08,
	      6.21e08)
kdr1= arr2str(kdr,',')

# solid grain bulk modulus, Pa [ kg / (m * s**2) ]
ksg =(  3.6e10,
	      3.6e10,
	      3.6e10,
	      3.6e10)
ksg1= arr2str(ksg,',')

# fluid bulk modulus, Pa [ kg / (m * s**2) ]
kfl =(  2.25e09,
	      2.25e09,
	      2.25e09,
	      2.25e09)
kfl1= arr2str(kfl,',')

# fluid viscosity, Pa*s [ kg / (m * s) ]
fvs =(  0.001,
	      0.001,
	      0.001,
	      0.001)
fvs1= arr2str(fvs,',')

# porosity, %
phi =(  0.1,
	      0.1,
	      0.25,
	      0.25)
phi1= arr2str(phi,',')

# permeability, m**2
prm =(  1.0e-14,
	      1.0e-14,
	      1.0e-14,
	      1.0e-14)
prm1= arr2str(prm,',')

# shear modulus, Pa [ kg / (m * s**2) ]
shm =(  5.25e09,
	      5.25e09,
  	    5.25e09,
	      5.25e09)
shm1= arr2str(shm,',')

# tortuosity, -
tor =(  2.42,
	      2.42,
	      2.49,
	      2.49)
tor1= arr2str(tor,',')

# ==============================================================================
# Testing

dt = par['dt']
dx = par['dx']
dz = par['dz']

nx = par['nx']
nz = par['nz']
nt = par['nt']

fmax = par['frq']*1.5
intervals = 4
frq = par['frq']
F = 0.62  # Formation Factor

# Existing
kdr = np.array(kdr)
kfl = np.array(kfl)
sro = np.array(sro)
fro = np.array(fro)
ksg = np.array(ksg)
shm = np.array(shm)
tor = np.array(tor)
prm = np.array(prm)
phi = np.array(phi)
fvs = np.array(fvs)
safety = 0.5

# Derived
bro = (1.0 - phi)*sro + phi*fro
alpha = 1.0 - kdr/ksg
skem = (1.0/kdr - 1.0/ksg) / (1.0/kdr - 1.0/ksg + phi*(1.0/kfl - 1.0/ksg))
kud = kdr/(1.0 - skem*(1.0 - kdr/ksg))
biotmod = skem * kud/alpha           # biot modulus, Pa
lambdau = kud - 2.0*shm/3.0
mro = tor*fro/phi
r_bar = mro*bro - fro*fro
fmo = fvs/prm
b = fvs/prm
biotfreq = (b*phi)/(tor*fro)
Pwmdr = kdr + 4.0*shm / 3.0          # P Wave modulus, drained, Pa
Pwmud = Pwmdr + alpha**2 * biotmod   # P Wave modulus, undrained, Pa
D = (prm*biotmod*Pwmdr)/(fvs*Pwmud)  # Hydraulic diffusivity, m**2 / s

# Coeffs
coeff_A = (mro*dt)/r_bar
coeff_B = (fro*dt)/r_bar
coeff_C = (fro*bro*b*dt)/r_bar
coeff_D = (1.0 + bro*b*dt) / (2.0*r_bar)
coeff_E = (-1.0*bro*b*dt) / (2.0*r_bar)
coeff_F = (-1.0*fro*dt) / r_bar
coeff_G = bro*dt / r_bar
coeff_H = shm*dt
coeff_I = lambdau * biotmod * dt
coeff_J = alpha * biotmod*dt
coeff_K = -1.0 * biotmod * dt

# Velocities
Vp = np.sqrt((lambdau + 2.0*shm)/bro)
Vs = np.sqrt(shm/bro)

# Stability
tdp = dt * Vp * np.sqrt(2)
tds = dt * Vs * np.sqrt(2)
wplength = safety*Vp / fmax
wslength = safety*Vs / fmax
passcrit = intervals*np.sqrt(dx*dx + dz*dz)
fmaxlim_p = safety*Vp/intervals/np.sqrt(dx*dx+dz*dz)
fmaxlim_s = safety*Vs/intervals/np.sqrt(dx*dx+dz*dz)
lambdad = np.sqrt(D/fmax)
De = fmax / biotfreq
dtmax = 0.606*np.min([dx,dz])/Vp
dtmaxcrewes = np.min([dx,dz]) / np.sqrt(Vp**2 - Vs**2)

pi1 = fro*bro*(F - fro/bro)
pi2 = fro*F*(lambdau + 2*shm) + bro*biotmod - 2*alpha*biotmod*fro
pi3 = biotmod*(lambdau + 2*shm - alpha**2*biotmod)

# ==============================================================================
# Stability Analysis








# ==============================================================================
# Python Test

def Dx(x,dx):
  dx1 = x[1:,:] - x[:-1,:]
  dx2 = x[2:,:] - x[:-2,:]
  x[1:-1,:] = (27.0*dx1[1:,:] - dx2)/(24.0*dx)
  return x

def Dz(z,dz):
  dz1 = z[:,1:] - z[:,:-1]
  dz2 = z[:,2:] - z[:,:-2]
  z[:,1:-1] = (27.0*dz1[:,1:] - dz2)/(24.0*dz)
  return z

def source(t, t0, f, factor):
  a = f*f*np.pi*np.pi
  return factor*np.exp(-a*(t-t0)**2)/(-2.0*a)

Lx = nx*dx  # length x, metres
Lz = nz*dz  # length z, metres
Lt = nt*dt  # length t, sec

# Source Properties
sx = np.int(nx/2)
sz = np.int(nz/2)
f0 = 80.0
t0 = 1.0/f0  # time delay, sec
factor = 1.0e2

# Generate spatial grid
layers = np.array(layers)
grid = np.ones([nx,nz])

qpx = np.zeros([nx,nz])
qmx = np.zeros([nx,nz])
qox = np.zeros([nx,nz])
qpz = np.zeros([nx,nz])
qmz = np.zeros([nx,nz])
qoz = np.zeros([nx,nz])

vpx = np.zeros([nx,nz])
vmx = np.zeros([nx,nz])
vox = np.zeros([nx,nz])
vpz = np.zeros([nx,nz])
vmz = np.zeros([nx,nz])
voz = np.zeros([nx,nz])

txx = np.zeros([nx,nz])
txz = np.zeros([nx,nz])
tzz = np.zeros([nx,nz])
p   = np.zeros([nx,nz])


# Coefficients
A_grid = coeff_A[0]
B_grid = coeff_B[0]
C_grid = coeff_C[0]
D_grid = coeff_D[0]
E_grid = coeff_E[0]
F_grid = coeff_F[0]
G_grid = coeff_G[0]
H_grid = coeff_H[0]
I_grid = coeff_I[0]
J_grid = coeff_J[0]
K_grid = coeff_K[0]


M_grid = biotmod[0]
alpha_grid = alpha[0]
lamu_grid = lambdau[0]
mu_grid = shm[0]
phi_grid = phi[0]
bro_grid = bro[0]
fro_grid = fro[0]
mro_grid = mro[0]
b_grid = b[0]
r_bar_grid = r_bar[0]

# Receivers
nrec = 1
nsrc = 1
rec_vx  = np.zeros([nrec,nt])
rec_vz  = np.zeros([nrec,nt])
rec_qx  = np.zeros([nrec,nt])
rec_qz  = np.zeros([nrec,nt])
rec_p   = np.zeros([nrec,nt])
rec_txx = np.zeros([nrec,nt])
rec_tzz = np.zeros([nrec,nt])
rec_txz = np.zeros([nrec,nt])
rec_src = np.zeros([nsrc,nt])


# Adjust frequency of source
#frq = frq*1000
f0 = 80.

# Adjust oscillating source parameters
osamp = 0.1
osper = 0.0005
osphase = 0.0
osindent = 121

# ==============================================================================
# Time Loop

for i in np.arange(0,nt):

  t = dt*i
  # ----------------------------------------------------------------------------
  # Displacement and strain to stress, t = n

  # Stress
  txx = txx + ( (lamu_grid + 2.0*mu_grid)*Dx(vox,dx) + lamu_grid*Dz(voz,dz) + alpha_grid*M_grid*(Dx(qox,dx) + Dz(qoz,dz)) ) * dt
  tzz = tzz + ( (lamu_grid + 2.0*mu_grid)*Dz(voz,dz) + lamu_grid*Dx(vox,dx) + alpha_grid*M_grid*(Dx(qox,dx) + Dz(qoz,dz)) ) * dt
  txz = txz + ( mu_grid*(Dz(vox,dz) + Dx(voz,dx)) ) * dt

  p   =   p - ( alpha_grid * M_grid * (Dx(vox,dx) + Dz(voz,dz)) + M_grid*(Dx(qox,dx) + Dz(qoz,dz)) ) * dt

  # ----------------------------------------------------------------------------
  # Introduce Stress Source
  # Source
#  txx[sx,sz] = txx[sx,sz] - source(t, t0, frq, factor)*(1.0 - phi_grid)
#  tzz[sx,sz] = tzz[sx,sz] - source(t, t0, frq, factor)*(1.0 - phi_grid)
#  p[sx,sz]   =   p[sx,sz] - source(t, t0, frq, factor)*phi_grid

  # Oscillating Source
  #p[0:,osindent:-osindent] = osamp*np.cos(2.0*np.pi/osper*t + osphase)

  # Gaussian Source
  a = np.pi*np.pi*f0*f0
  source_term = factor * np.exp(-1.0*a*(t-t0)**2)/(-2.0*a)
  
  p[sx,sz] = p[sx,sz] + source_term*M_grid

  # ----------------------------------------------------------------------------
  # Stress to Velocity, t = n + 1/2

  # Relative Velocity
  qpx = qox + ( (-1.0*fro_grid*(Dx(txx,dx) + Dz(txz,dz)) - bro_grid*b_grid*qox - bro_grid*Dx(p,dx) ) / r_bar_grid ) * dt
  qpz = qoz + ( (-1.0*fro_grid*(Dx(txz,dx) + Dz(tzz,dz)) - bro_grid*b_grid*qoz - bro_grid*Dz(p,dz) ) / r_bar_grid ) * dt

  qmx = (qox + qpx) / 2.0
  qmz = (qoz + qpz) / 2.0

  # Matrix Velocity
  vpx = vox + ( ( mro_grid*(Dx(txx,dx) + Dz(txz,dz)) + fro_grid*mro_grid*b_grid*qmx + fro_grid*Dx(p,dx) ) / r_bar_grid ) * dt
  vpz = voz + ( ( mro_grid*(Dx(txz,dx) + Dz(tzz,dz)) + fro_grid*mro_grid*b_grid*qmz + fro_grid*Dz(p,dz) ) / r_bar_grid ) * dt


  # Record Seismograms
#  rec_src[0,i] = source(t, t0, frq, factor)
  rec_src[0,i] = osamp*np.cos(2.0*np.pi/osper*t + osphase)
  rec_vx[0,i]  = vpx[np.int(nx/2), 2]
  rec_vz[0,i]  = vpz[np.int(nx/2), 2]
  rec_qx[0,i]  = qpx[np.int(nx/2), 2]
  rec_p[0,i]   = p[np.int(nx/2), 2]
  rec_txx[0,i] = txx[np.int(nx/2), 2]
  rec_tzz[0,i] = tzz[np.int(nx/2), 2]
  rec_txz[0,i] = txz[np.int(nx/2), 2]


  # Boundary Conditions
  vpx[-1,:] = 0
  vpx[0, :] = 0
  vpz[:,-1] = 0
  vpz[:, 0] = 0
  qpx[-1,:] = 0
  qpx[0, :] = 0
  qpz[:,-1] = 0
  qpz[:, 0] = 0
  
  # Switch pointers
  vox = vpx
  voz = vpz

  qox = qpx
  qoz = qpz
