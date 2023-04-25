import numpy as np

# We make the estimates of mutation matrices into stochastic matrices.
def normalize_matrix(matrix):
    z = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            z[i,j] = np.real(matrix[i,j])
    matrix = z

    matrix[matrix<0]=0 # Remove negative entries

    # Renormalize (rows sum up to 1)
    matrix[0, :]=matrix[0,:]/np.linalg.norm(matrix[0,:], ord=1)
    matrix[1, :]=matrix[1,:]/np.linalg.norm(matrix[1,:], ord=1)
    matrix[2, :]=matrix[2,:]/np.linalg.norm(matrix[2,:], ord=1)
    matrix[3, :]=matrix[3,:]/np.linalg.norm(matrix[3,:], ord=1)
    return matrix

def normalize_distribution(distr):
    z = np.zeros(4)
    for i in range(4):
        z[i] = np.real(distr[i])
    distr = z

    distr[distr<0]=0 # Remove negative entries

    # Renormalize (rows sum up to 1)
    distr=distr/np.linalg.norm(distr, ord=1)
    return distr
