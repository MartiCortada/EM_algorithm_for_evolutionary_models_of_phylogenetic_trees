from FullRecon import FullRecon
from simulate import matrix_generation
import matplotlib.pyplot as plt
import numpy as np
import sys

def init_EM():

    repetitions = 25

    tree = sys.argv[1]
    length = sys.argv[2]

    real_matrices, sequences_file_name = matrix_generation(tree, length)

    for r in range(repetitions):

        ###print(r)
        if r == 0:
            Results, final_nodes_distributions, net = FullRecon(sequences_file_name, tree)

        else:
            Res, distr, _ = FullRecon(sequences_file_name, tree)
            for MutationMatrixA in Res:
                for MutationMatrixB in Results:
                    if (str(MutationMatrixA.source) == str(MutationMatrixB.source) and str(MutationMatrixA.target)==str(MutationMatrixB.target)):
                        MutationMatrixB.matrix = MutationMatrixB.matrix + MutationMatrixA.matrix
            for k in final_nodes_distributions.keys():
                for k_ in distr.keys():
                    if (str(k) == str(k_)):
                        final_nodes_distributions[k] += distr[k_]
    
    length_results = len(Results)
    for i in range(length_results):
        Results[i].matrix = 1/repetitions * Results[i].matrix
        assert(str(real_matrices[i].source) == str(Results[i].source))
        assert(str(real_matrices[i].target) == str(Results[i].target))
        # print(Results[i].source, Results[i].target, "Estimated Matrix")
        # print(" ")
        # print(Results[i].matrix)
        # print(" ")
    
    return Results
