import sys
sys.path.insert(0,'/g/g92/chung28/Programs/python-quartz/lib/python2.7/site-packages')
from scipy.optimize import linprog
from scipy import sparse
import numpy as np

class wasserstein(object):
    dtype = np.double
    
    def __init__(self,N,L):
        self.N = N
        self.L = L
        
        Aeq = np.zeros([2*N,N**2])
        A0 = np.eye(N)
        A1 = np.ones([1,N])
        Aeq[0:N,:] = np.kron(A0,A1)
        Aeq[N:2*N,:] = np.kron(A1,A0)
        Aeq = Aeq[1:2*N,:]

        beq = np.ones([2*N-1,1],dtype='double')/N

        bound = (0.0,None)
        bounds=[]
        for i in range(N**2):
            bounds.append(bound)
            
        self.Aeq = Aeq
        self.beq = beq
        self.bounds = bounds
        
    def measure(self, x0, x1):
        xN, xdim = x0.shape
        if( xN != self.N ):
            print "Error! apply different number of particles."
            print "Setup N: ",self.N
            print "Given N: ",xN
            return
        
        X0 = np.zeros((xN,xN,xdim))
        X1 = np.zeros((xN,xN,xdim))
        for i in range(xdim):
            X0[:,:,i], X1[:,:,i] = np.meshgrid(x0[:,i],x1[:,i])
            
        dX = abs(X0-X1)
        dX[:,:,0] = np.minimum( dX[:,:,0], self.L-dX[:,:,0] )
        dX = np.reshape(np.sqrt(np.sum( dX**2, 2 )),(xN**2,))

#        import timeit
#        start = timeit.default_timer()
        s = linprog(dX,None,None,self.Aeq,self.beq,self.bounds,method='interior-point')
#        stop = timeit.default_timer()
#        print stop-start
        if( s.status != 0 ):
            print s.status, s.message
        return s.fun

class boundedLipschitz(object):
    dtype = np.double
    
    def __init__(self,N,Ng,L):
        self.N = N
        self.Ng = Ng
        self.L = L
        self.dx = L/Ng
        
        #mesh-mesh distance
        Am = sparse.diags([-1.,1.],[-1,0],shape=(Ng,Ng))
        Am[0,Ng-1] = -1.0
        Amm = np.zeros(2*Ng,Ng)
        Amm[0:Ng,:] = Am.todense()
        Amm[Ng:2*Ng,:] = -Am.todense()

        #particle-particle distance
        Ap = sparse.diags([-1.,1.],[-1,0],shape=(2*N,2*N))
        Ap[0,2*N-1] = -1.0
        App = np.zeros(4*N,2*N)
        App[0:2*N,:] = Ap.todense()
        App[2*N:4*N,:] = -Ap.todense()

        #particle-mesh distance
        Apm_ = sparse.diags([-1.,1.],[-1,0],shape=(Ng+N,Ng+N))
        Apm_[0,Ng+N-1] = -1.0
        Apm = np.zeros(2*(Ng+N),Ng+N)
        Apm[0:Ng+N,:] = Apm_.todense()
        Apm[Ng+N:2*(Ng+N),:] = -Apm_.todense()

        bmm = np.ones([2*Ng,1],dtype='double')*self.dx

        bound = (0.0,None)
        bounds=[]
        for i in range(N**2):
            bounds.append(bound)
            
        self.Amm, self.App, self.Apm = Amm, App, Apm
        self.beq = beq
        self.bound = bound
        
    def mm(self, rho0, rho1):
        Ng = rho0.size
        if( (Ng != self.Ng) and (Ng !=rho1.size) ):
            print "Error! apply different number of particles."
            print "Setup Ng: ",self.Ng
            print "rho0: ",Ng
            print "rho1: ",rho1.size
            return
        
        drho = (rho0-rho1)*self.dx

#        import timeit
#        start = timeit.default_timer()
        s = linprog(drho,self.Amm,self.bmm,None,None,None,method='interior-point')
#        stop = timeit.default_timer()
#        print stop-start
        if( s.status != 0 ):
            print s.status, s.message
        return s.fun

