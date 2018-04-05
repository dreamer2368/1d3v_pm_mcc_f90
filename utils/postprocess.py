import numpy as np
from wasserstein import *

dir1 = 'Debye'
dir2 = 'Debye_perturbed2'

params = open(dir1+'/record','r')
text = params.read().strip().split()
params.close()
Ns = int(text[0])
Ng = int(text[1])
Nt = int(text[2])
L = float(text[3])
Rmod = int(text[4])
Nr = Nt/Rmod
Np = np.fromfile(dir1+'/Np.bin',dtype='int32')
wp = L/Np[0]

Time = 150.0
dt = Time/Nt
T = np.linspace(0.0,150.0,Nr+1)

dx = L/Ng
xg = np.arange(Ng)*dx + 0.5*dx

E0 = np.fromfile(dir1+'/E.bin',dtype='double')
E0 = np.reshape(E0,(Nr+1,Ng))

phi0 = np.fromfile(dir1+'/phi.bin',dtype='double')
phi0 = np.reshape(phi0,(Nr+1,Ng))

rho0 = np.fromfile(dir1+'/rho.bin',dtype='double')
rho0 = np.reshape(rho0,(Nr+1,Ng))

E1 = np.fromfile(dir2+'/E.bin',dtype='double')
E1 = np.reshape(E1,(Nr+1,Ng))

phi1 = np.fromfile(dir2+'/phi.bin',dtype='double')
phi1 = np.reshape(phi1,(Nr+1,Ng))

rho1 = np.fromfile(dir2+'/rho.bin',dtype='double')
rho1 = np.reshape(rho1,(Nr+1,Ng))

dE = np.sqrt(np.sum( (E1-E0)**2.0, axis=1)*dx)
dphi = np.sqrt(np.sum( (phi1-phi0)**2.0, axis=1)*dx)
drho = np.sqrt(np.sum( (rho1-rho0)**2.0, axis=1)*dx)

i = 144
xp0 = np.fromfile(dir1+'/xp/'+str(i)+'_1.bin',dtype='double')
vp0 = np.fromfile(dir1+'/vp/'+str(i)+'_1.bin',dtype='double')
vp0 = np.reshape(vp0,(3,1e5))

xp1 = np.fromfile(dir2+'/xp/'+str(i)+'_1.bin',dtype='double')
vp1 = np.fromfile(dir2+'/vp/'+str(i)+'_1.bin',dtype='double')
vp1 = np.reshape(vp0,(3,1e5))

Nw = 100
idxT = int(Np[0]/Nw)
Nw = len(xp0[0::idxT])

x0 = np.transpose([xp0[0::idxT],vp0[0,0::idxT]])
x1 = np.transpose([xp1[0::idxT],vp1[0,0::idxT]])
W = wasserstein(Nw,L)
W.measure(x0,x1)
