import numpy as np
import os
from mpi4py import MPI
from wasserstein import *

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

dir1 = '../data/Debye'
dir2 = '../data/Debye_perturbed2'

if( rank == 0 ):    
    params = open(dir1+'/record','r')
    text = params.read().strip().split()
    params.close()
    Np = np.fromfile(dir1+'/Np.bin',dtype='int32')
    Np = Np[0]
else:
    text, Np = None, None

text = comm.bcast(text, root=0)
Np = comm.bcast(Np, root=0)

Ns = int(text[0])
Ng = int(text[1])
Nt = int(text[2])
L = float(text[3])
Rmod = int(text[4])

Nr = Nt/Rmod
wp = L/Np

Time = 150.0
dt = Time/Nt
T = np.linspace(0.0,150.0,Nr+1)

dx = L/Ng
xg = np.arange(Ng)*dx + 0.5*dx

Nw = 100
idxT = int(Np/Nw)
Nw = len(range(Np)[::idxT])

W = wasserstein(Nw,L)

def loadAndMeasure(i,W):
    xp0 = np.fromfile(dir1+'/xp/'+str(i)+'_1.bin',dtype='double')
    vp0 = np.fromfile(dir1+'/vp/'+str(i)+'_1.bin',dtype='double')
    vp0 = np.reshape(vp0,(Np,3))

    xp1 = np.fromfile(dir2+'/xp/'+str(i)+'_1.bin',dtype='double')
    vp1 = np.fromfile(dir2+'/vp/'+str(i)+'_1.bin',dtype='double')
    vp1 = np.reshape(vp0,(Np,3))

    x0 = np.transpose([xp0[0::idxT],vp0[0::idxT,0]])
    x1 = np.transpose([xp1[0::idxT],vp1[0::idxT,0]])
    return W.measure(x0,x1)

amode = MPI.MODE_WRONLY | MPI.MODE_CREATE
filename = 'wasserstein.bin'
fh = MPI.File.Open(comm,filename,amode)
dW = loadAndMeasure(rank,W)
fh.Write_at_all(8*rank,dW)
fh.Close()
