from libs import solver
from libs import debug
from libs import plot
from libs import validation

# This example code illustrates what FEniCS computes for our original
# weak formulation when surface tension is introduced

M = 1			# Mode of vibration
P = 5			# Number of eigenvalues to solve for
Mu = 0.01		# Surface tension constant, 1/Bo
Nx = 100		# Cylinder: Divisions along the radial direction
Ny = 100		# Cylinder: Divisions along the theta direction

Order = [1,2,3]		# Run the test for the following list of orders
Options = {"spectral_transform" : "shift-and-invert"
	   "problem_type" : "gen_hermitian"
	   "spectrum" : "smallest magnitude"
	   "spectral_shift" : 1e-6}

for i in Order:
	A1 = solver.solver(m=M,p=P,mu=Mu,order=i,Opt="cyltank "+str(Nx)+" "+str(Ny)+" 1.0 1.0")
	(evals,evecs,U) = A1.solve(Options)
	plot.PlotHeightProfiles(U,M)
