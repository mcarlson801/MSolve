from libs import solver
from libs import debug
from libs import plot
from libs import validation
from matplotlib import pyplot as plt

# The following example code generates the graphics seen in the handout,
# specifically it is a plot of the domain being computed (mirrored around r=0)
# as well as the computed solution

M = 2		# Mode of vibration
Nr = 100	# Cylinder: Number of divisions in the radial direction
Nt = 100	# Cylinder: Number of divisions in the theta direction
H = 500		# Plotting: Number of samples for plotting height profile
P = 6		# Number of eigenvalues to compute and plot

# Dictionary of options to pass to the solver
Options = {"spectral_transform" : "shift-and-invert"
	   "problem_type" : "shift-and-invert"
	   "spectrum" : "smallest real"
	   "spectral_shift" : 1e-6}

A1 = solver.solver(m=M,p=P,mu=0,order=2,Opt="cyltank "+str(Nt)+" "+str(Nt)+" 1.0 1.0")
(evals,evecs,U) = A1.solve(Options)
plot.DualPlots(U,M,H,container='cylinder')

A1 = solver.solver(m=M,p=P,mu=0,order=2,Opt="halfsphere")
(evals,evecs,U) = A1.solve(Options)
plot.DualPlots(U,M,H,container='halfcircle')

A1 = solver.solver(m=M,p=P,mu=0,order=2,Opt="tallglass")
(evals,evecs,U) = A1.solve(Options)
plot.DualPlots(U,M,H,container='tallglass')

A1 = solver.solver(m=M,p=P,mu=0,order=2,Opt="dentedglass")
(evals,evecs,U) = A1.solve(Options)
plot.DualPlots(U,M,H,container='dentedglass')
