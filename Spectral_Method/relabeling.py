import itertools
import numpy as np
L=list(itertools.permutations([0, 1, 2, 3]))
#print(L)

def permutation(matrix):

    #print(matrix)
    #print(" ")

    best_score = 0
    best_perm = 0
    current_perm = 0
    for perm in L:
        matrix_perm = np.zeros((4,4), dtype=complex)
        for i in range(4):
            for j in range(4):
                matrix_perm[j, i] = matrix[j, perm[i]]
        Score = sum(np.abs(np.diag(matrix_perm)))
        if (Score > best_score):
            best_score = Score
            best_perm = current_perm
        current_perm += 1

    final_perm = L[best_perm]
    matrix_perm = np.zeros((4,4), dtype=complex)
    for i in range(4):
        for j in range(4):
            matrix_perm[j, i] = matrix[j, final_perm[i]]
    Score = sum(np.abs(np.diag(matrix_perm)))
    #print(matrix_perm)
    #print("Score: ", Score)
    #print("-----------------------------------------------")
    return matrix_perm
