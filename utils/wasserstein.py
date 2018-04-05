import numpy as np
from scipy.optimize import linprog

class wasserstein(object):
    dtype = np.double
    
    def __init__(self,N,L):
        self.N = N
        self.L = L
        
        Aeq = np.zeros([2*N+1,N**2])
        Aeq[0,:] = 1.0
        A0 = np.eye(N)
        A1 = np.ones([1,N])
        Aeq[1:N+1,:] = np.kron(A0,A1)
        Aeq[N+1:2*N+1,:] = np.kron(A1,A0)

        beq = np.ones([2*N+1,1])
        beq[1:2*N+1] *= 1.0/N

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
        
        s = linprog(dX,None,None,self.Aeq,self.beq,self.bounds)
        return s.fun
