# MSolve
Free Surface Sloshing numerical solution using the FEniCS framework for solving PDEs using the finite element method

Libraries: <br/>
    solver: <br/>
        This is the core solver class. It is initialized with a number of problem-specific parameters and then solved using the solve() method. The solve method takes a dictionary of options specific to solving the resultant numerical eigenvalue problem. <br/>
    validation: <br/>
        Contains the code necessary for validating solutions including the exact solution for the cylindrical tank and any other tests that may be useful. <br/>
    debug: <br/>
        Contains methods for exposing FEniCS information such as the underlying matrices as well as a map of the degrees of freedom to mesh coordinates. <br/>
    plot: <br/>
        These are the plotting tools used to generate all the graphics and animations. <br/>
    
Examples: <br/>
    Ex1-ErrorAnalysis: <br/>
        Illustrates how to use the validation class to compute relative error for the cylindrical tank case <br/>
    Ex2-DualPlot: <br/>
        Example code for generating "dual plots". A dual plot is a plot that shows a computed eigenfunction with its corresponding container <br/>
    Ex3-CylinderConvergence: <br/>
        Convergence study using the cylindrical tank exact solution. Takes about 3 minutes to run. <br/>
    Ex4-SurfaceTensionTests: <br/>
        Example of what happens when surface tension is introduced to our original weak formulation <br/>
