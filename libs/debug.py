import numpy as np

# Helper method for exporting and viewing matrices
def ExportPETScMatrixFile(A,k):
    retval = np.zeros([A.size(0),A.size(1)])
    M = A.mat()
    for i in range(0,A.size(0)):
        for j in range(0,A.size(1)):
	    retval[i][j] = M[i,j]
    np.savetxt("matrix"+k+".csv", retval, delimiter=',')

def ConvertPETScToNumpy(A):
    M = A.mat()
    retval = np.zeros([M.size(0),M.size(1)])
    for i in range(0,M.size(0)):
        for j in range(0,M.size(1)):
	    try:
	    	retval[i][j] = M[i,j]
	    except:
		pass
    return retval

def PrintDOFs(S):
    dofs_0 = S.dofmap().dofs()
    print(dofs_0)
    print(S.dofmap().tabulate_all_coordinates(mesh))
