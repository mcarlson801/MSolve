import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from libs import solver
from libs import debug
from libs import plot
from libs import validation

# Knowing the exact solution for the cylindrical tank, the
# following example code runs a converge study for increasing
# granularity of the computational domain

M = 1			# Mode of vibration
P = 1			# Number of eigenvalues to compute
N = 40			# Number of times to refine the mesh
RE1 = np.zeros([P,N])	
RE2 = np.zeros([P,N])
REInf = np.zeros([P,N])


for j in range(N):
    A1 = solver.solver(m=M,p=P,mu=0,order=1,Opt="cyltank "+str(10*(j+1))+" "+str(10*(j+1))+" 1.0 1.0")
    (evals,evecs,U) = A1.solve()
    for i in range(0,P):
        (sf,re1,re2,reinf) = validation.ComputeErrors(U[i],M,0,10*(j+1),10*(j+1))
        RE1[i,j] = re1
        RE2[i,j] = re2
        REInf[i,j] = reinf

plt.figure(1)
plt.loglog(np.linspace(0,N-1,N),RE1[0,:])
plt.figure(2)
plt.loglog(np.linspace(0,N-1,N),RE2[0,:])
plt.figure(3)
plt.loglog(np.linspace(0,N-1,N),REInf[0,:])
plt.show(block=False)
raw_input("Press enter to continue.")
plt.close("all")
