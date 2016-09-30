# MSolve
Free Surface Sloshing numerical solution using the FEniCS framework for solving PDEs using the finite element method

Libraries: 
    solver:
        This is the core solver class. It is initialized with a number of problem-specific parameters and then solved using the solve() method. The solve method takes a dictionary of options specific to solving the resultant numerical eigenvalue problem.
    validation:
        Contains the code necessary for validating solutions including the exact solution for the cylindrical tank and any other tests that may be useful.
    debug:
        Contains methods for exposing FEniCS information such as the underlying matrices as well as a map of the degrees of freedom to mesh coordinates.
    plot:
        These are the plotting tools used to generate all the graphics and animations.
    
Examples:
    Ex1-ErrorAnalysis:
        Illustrates how to use the validation class to compute relative error for the cylindrical tank case
    Ex2-DualPlot:
        Example code for generating "dual plots". A dual plot is a plot that shows a computed eigenfunction with its corresponding container
    Ex3-CylinderConvergence:
        Convergence study using the cylindrical tank exact solution. Takes about 3 minutes to run.
    Ex4-SurfaceTensionTests:
        Example of what happens when surface tension is introduced to our original weak formulation
