from libs import solver
from libs import debug
from libs import plot
from libs import validation

# Knowing the exact solution for the cylindrical tank, the relative error
# of FEniCS' computation can be computed

M = 1		# Mode of vibration
Nr = 200	# Cylinder: Divisions along the radial direction
Nt = 200	# Cylinder: Divisions along the theta direction

A1 = solver.solver(m=M,p=1,mu=0,order=2,Opt="cyltank "+str(Nr)+" "+str(Nt)+" 1.0 1.0")
(evals,evecs,U) = A1.solve()

(sf,re1,re2,reinf) = validation.ComputeErrors(U[0],M,0,Nx,Ny)
print('Relative Error (1): ' + str(re1))
print('Relative Error (2): ' + str(re2))
print('Relative Error (Infinity): ' + str(reinf))

plot.PlotHeightProfiles(U,M,scalefactor=sf)
plot.PlotHeightProfiles([validation.ExactSolutionSurface()],M)


