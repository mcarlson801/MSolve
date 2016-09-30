from dolfin import *
import numpy as np

class solver:
    # ================== Initialize Solver ======================
    def __init__(self,m=0,p=1,mu=0,order=1,Opt="cyltank 40 40 1.0 1.0"):

	# ----------- Parameters -----------
	self.Opts = Opt.split()
	self.m = m				# Mode of vibration
	self.p = p				# Number of times to plot / Number of eigenvalues to compute
	self.U = [0]*p				# Container of computed eigenfunctions
	self.order = order			# Order of Basis Elements
	self.mu = mu				# Surface Tension Constant
	# ----------------------------------


	# ----------- Check Problem Type, Set Problem-Specific Parameters -----------
	# Parameters specific to Cylinder Tank problem
	if self.Opts[0] == "cyltank":
	    self.hr = int(self.Opts[1])		# Number of divisions along the R axis
	    self.hz = int(self.Opts[2])         # Number of divisions along the Z axis
	    self.s = float(self.Opts[3])	# Radius of tank
	    self.d = float(self.Opts[4])	# Depth of tank
	    self.mesh = RectangleMesh(Point(0.0,-1.0*self.d),Point(self.s,0.0),self.hr,self.hz,"crossed")
	# General mesh
	else:
	    st1 = "meshes/"+Opt+".xml"
	    self.mesh = Mesh(st1)
	# ---------------------------------------------------------------------------
	    
	# ----------- Mark Free Surface Subdomain -----------
	class Steklov(SubDomain):
	    def __init__(self):
		super(Steklov,self).__init__()
            def inside(self,x,on_boundary):
                return on_boundary and near(x[1],0) and x[0] <= 1.0
        self.subdomains = MeshFunction("size_t",self.mesh,1)
        self.subdomains.set_all(0)
    	self.steklov = Steklov()
    	self.steklov.mark(self.subdomains,50)
	# ---------------------------------------------------
    # ===========================================================

    def solve(self,opts):
	# Define Function Space, Trial/Test Functions, boundary region for integration, and expression for r
	self.V = FunctionSpace(self.mesh, "Lagrange", self.order)
        u = TrialFunction(self.V)
        v = TestFunction(self.V)
	ds = Measure("ds")[self.subdomains]
	self.r = Expression('x[0]',cell=triangle,degree=1)
	
	# Weak forms for m = 0 and otherwise
	if self.m == 0:
	    a = (inner(grad(u),grad(v))*self.r**(2*self.m+1) + 2*self.m**2*u*v*self.r**(2*self.m-1) + self.m*(v*Dx(u,0) +
		u*Dx(v,0))*self.r**(2*self.m))*dx
            b = u*v*self.r**(2*self.m+1)*ds(50)
	else:
	    a = (self.r**3*inner(grad(u),grad(v)) + self.r**2*(u*Dx(v,0)+v*Dx(u,0)) + 
		Constant(self.m**2+1)*self.r*u*v)*dx
	    a = a + Constant(self.mu)*( Dx(v,0)*Dx(Dx(u,1),0)*self.r**3 + Constant(self.m**2+1)*self.r*Dx(u,1)*v +
		self.r**2*(Dx(u,1)*Dx(v,0) + Dx(Dx(u,1),0)*v) )*ds(50)
	    b = self.r**3*u*v*ds(50)
	
	# Assemble Matrices
	A = PETScMatrix()
	assemble(a, tensor=A)
	B = PETScMatrix()
	assemble(b, tensor=B)

	eigensolver = SLEPcEigenSolver(A,B)
	for key in opts:
	    eigensolver.parameters[key] = opts[key]
	print("Computing eigenvalues. This can take a minute.")
	eigensolver.solve(self.p)

	# Compute p eigenfunctions and store eigenvalues and eigenvectors
	self.evals = [0]*self.p
	self.evecs = [0]*self.p
	self.C = [0]*self.p
	for i in range(0,self.p):
	    l, c, rx, cx = eigensolver.get_eigenpair(i)
	    self.evals[i] = l
	    self.evecs[i] = rx
	    self.U[i] = Function(self.V)
	    self.U[i].vector()[:] = self.evecs[i]
	    if self.m!=0:
	        self.U[i] = project(self.U[i]*self.r,self.V)
	return [self.evals,self.evecs,self.U]
