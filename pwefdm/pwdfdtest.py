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
    'nx':140,  'ox':0, 'dx':3,     'lx':'x', 'ux':'m',
    'ny':750,  'oy':0, 'dy':3,     'ly':'y', 'uy':'m',
    'nz':620,  'oz':0, 'dz':3,     'lz':'z', 'uz':'m',
    'nt':500,  'ot':0, 'dt':2e-4, 'lt':'t', 'ut':'s',
    'kt':100,
    'jsnap':100,
    'height':10,
    'nb':20,
    'frq':50.0, # source wavelet frequency
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
sro =(  2650.0,
	      2650.0,
	      2650.0,
	      2650.0)
sro1= arr2str(sro,',')

# fluid density, kg/m**3
fro =(  1000.0,
	      1000.0,
	      1000.0,
	      1000.0)
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
phi =(  0.3,
	      0.3,
	      0.3,
	      0.3)
phi1= arr2str(phi,',')

# permeability, m**2
prm =(  1.0e-12,
	      1.0e-12,
	      1.0e-12,
	      1.0e-12)
prm1= arr2str(prm,',')

# shear modulus, Pa [ kg / (m * s**2) ]
shm =(  4.55e08,
	      4.55e08,
	      4.55e08,
	      4.55e08)
shm1= arr2str(shm,',')

# tortuosity, -
tor =(  10.0,
	      10.0,
	      10.0,
	      10.0)
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
coeff_C = (fro*fmo*dt)/r_bar
coeff_D = (2.0*r_bar - bro*fmo*dt) / (2.0*r_bar + bro*fmo*dt)
coeff_E = -1.0*(fro*dt) / (r_bar + bro*dt)
coeff_F = -1.0*(bro*dt) / (r_bar + bro*dt)
coeff_G = shm*dt
coeff_H = lambdau*dt
coeff_I = alpha*biotmod*dt
coeff_J = -1.0*biotmod*dt

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
t0 = 0.0  # time delay, sec
factor = 1e9

# Generate spatial grid
layers = np.array(layers)
grid = np.ones([nx+1,nz+1])

qpx = np.zeros([nx+1,nz+1])
qmx = np.zeros([nx+1,nz+1])
qpz = np.zeros([nx+1,nz+1])
qmz = np.zeros([nx+1,nz+1])

vpx = np.zeros([nx+1,nz+1])
vmx = np.zeros([nx+1,nz+1])
vpz = np.zeros([nx+1,nz+1])
vmz = np.zeros([nx+1,nz+1])

txx = np.zeros([nx+1,nz+1])
txz = np.zeros([nx+1,nz+1])
tzz = np.zeros([nx+1,nz+1])
p   = np.zeros([nx+1,nz+1])


# Coefficients
#A_grid = mro[0]/r_bar[0]
#B_grid = (fro[0]*b[0])/(2.0*r_bar[0])
#C_grid = fro[0]/r_bar[0]
#D_grid = (2.0*mro[0]*bro[0] - 2.0*fro[0]**2 - fro[0]*b[0]*dt)/(2.0*mro[0]*bro[0] - 2.0*fro[0]**2 + fro[0]*b[0]*dt)
#E_grid = (2.0*fro[0])/(2.0*mro[0]*bro[0] - 2.0*fro[0]**2 + fro[0]*b[0]*dt)
#F_grid = (2.0*bro[0])/(2.0*mro[0]*bro[0] - 2.0*fro[0]**2 + fro[0]*b[0]*dt)

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

M_grid = biotmod[0]
alpha_grid = alpha[0]
lamu_grid = lambdau[0]
mu_grid = shm[0]
phi_grid = phi[0]


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

for i in np.arange(0,nt):
  
  t = dt*i
  
  
  # Relative Velocity
  qpx = D_grid*qmx  + E_grid*(Dx(txx,dx) + Dz(txz,dz)) + F_grid*Dx(p,dx)
  qpz = D_grid*qmz  + E_grid*(Dx(txz,dx) + Dz(tzz,dz)) + F_grid*Dz(p,dz)
  
  # Matrix Velocity
  vpx = vmx + A_grid*(Dx(txx,dx) + Dz(txz,dz)) + B_grid*Dx(p,dx) + C_grid*(qpx + qmx)/2
  vpz = vmz + A_grid*(Dx(txz,dx) + Dz(tzz,dz)) + B_grid*Dz(p,dz) + C_grid*(qpz + qmz)/2 
  
  # Pressure
  p   =   p + I_grid*(Dx(vpx,dx) + Dz(vpz,dz)) + J_grid*(Dx(qpx,dx) + Dz(qpz,dz))
  
  # Stress
  txx = txx + G_grid*(Dx(vpx,dx) + Dx(vpx,dx)) + ( H_grid*(Dx(vpx,dx) + Dz(vpz,dz)) + I_grid*(Dx(qpx,dx) + Dz(qpz,dz)) )
  txz = txz + G_grid*(Dx(vpz,dx) + Dz(vpx,dz)) 
  tzz = tzz + G_grid*(Dz(vpz,dz) + Dz(vpz,dz)) + ( H_grid*(Dx(vpx,dx) + Dz(vpz,dz)) + I_grid*(Dx(qpx,dx) + Dz(qpz,dz)) )

  # Source
  txx[sx,sz] = txx[sx,sz] + source(t, t0, frq, factor)*(1.0 - phi_grid)
  tzz[sx,sz] = tzz[sx,sz] + source(t, t0, frq, factor)*(1.0 - phi_grid)  
  p[sx,sz]   = p[sx,sz]   + source(t, t0, frq, factor)*phi_grid
  
  # Record Seismograms
  rec_src[0,i] = source(t, t0, frq, factor)
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
  vmx = vpx
  vmz = vpz
  
  qmx = qpx
  qmz = qpz
 
  
  
