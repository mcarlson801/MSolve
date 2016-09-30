import numpy as np
import scipy.integrate as integrate
import scipy.special as special
from scipy.misc import derivative as Dv

# Exact solution for m>=1 and n>=0 over the whole cylinder domain
def ExactSolutionFull(m=1,n=0):
    kap = special.jnp_zeros(m,n+1)[n]
    if m >= 1 and n >= 0:
        return lambda r,theta,z: special.jn(m,kap*r)*np.cos(m*theta)*(np.cosh(kap*(1.0+z))/np.cosh(kap))
    else:
	return lambda r,theta,z: 0

# Exact solution for m>=1 and n>=0 over the surface
def ExactSolutionSurface(m=1,n=0):
    kap = special.jnp_zeros(m,n+1)[n]
    if m >= 1 and n >= 0:
        return lambda r,theta: special.jn(m,kap*r)*np.cos(m*theta)
    else:
	return lambda r,theta: 0

def IntegralTest1(U,p,m,radius=1.0):
    intvals = [0]*p
    for i in range(0,p):
        u = U[i]
        intvals[i] = integrate.dblquad(lambda r,theta: np.cos(m*theta)*u(r,0),0,2*np.pi,
				       lambda r: 0,lambda r: radius)
    return intvals

def ComputeErrors(u,m,n,Nx,Ny):
    uhat = ExactSolutionFull(m,n)
    Uhat = np.zeros((Nx+1)*(Ny+1))
    U = np.zeros((Nx+1)*(Ny+1))
    for i in range(0,Nx+1):
	for j in range(0,Ny+1):
	    Uhat[i+j*(Nx+1)] = uhat((1.0*i)/(1.0*Nx),0,(-1.0*j)/(1.0*Ny))
	    U[i+j*(Nx+1)] = u((1.0*i)/(1.0*Nx),(-1.0*j)/(1.0*Ny))
    scalefactor = np.linalg.norm(Uhat,2)/np.linalg.norm(U,2)
    if (uhat(1.0,0,0) > 0 and u(1.0,0) < 0) or (uhat(1.0,0,0) < 0 and u(1.0,0) > 0):
	scalefactor = -1.0 * scalefactor

    residual = Uhat - scalefactor*U
    Re1 = np.linalg.norm(residual,1)/np.linalg.norm(Uhat,1)    
    Re2 = np.linalg.norm(residual,2)/np.linalg.norm(Uhat,2)
    ReInf = np.linalg.norm(residual,np.inf)/np.linalg.norm(Uhat,np.inf)
    
    return (scalefactor,Re1,Re2,ReInf)
    


