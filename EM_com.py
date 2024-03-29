# Import required libraries
import numpy as np
from Bio import SeqIO, Phylo
from io import StringIO
from collections import Counter, OrderedDict
import itertools
import math
import copy
import sys
import os
import time
# Import required functions from modules we have programmed
from data_simulation import simulate, generate_alignment, generate_sequences
from evaluate import init_EM

newick_input = sys.argv[1]
length = int(sys.argv[2])
repetitions = int(sys.argv[3])

eps = 10**(-3)                  # epsilon for convergence

# Create a directory to store the results
try:
    directory = "./OUTPUTS/"
    os.mkdir(directory)
except FileExistsError:
    pass
try:
    directory += "RANDOM_vs_SPECTRAL_EM_" + str(length) + "_a_b_10/" 
    os.mkdir(directory)
except FileExistsError:
    pass
try:
    results = directory + "RESULTS/" 
    os.mkdir(results)
except FileExistsError:
    pass

######################################## Read newick tree ########################################
newick_tree = open("./" + newick_input, "r")
newick = newick_tree.read()
newick_tree.close()
tree = Phylo.read(StringIO(newick), "newick")
# Change nodes names
for idx, clade in enumerate(tree.get_nonterminals()):
    clade.name = "Int_" + str(idx) if idx > 0 else "Int_0" 
# Change leaves names
for idx, clade in enumerate(tree.get_terminals()):
    clade.name = "Leaf_" + clade.name 
net = Phylo.to_networkx(tree) # to graph

######################################### Simulate data #########################################
edges, node_distr, filename = simulate(net, length, directory) 
# Preprocess leave alignments
leave_alignments = [i for i in SeqIO.parse(directory+filename, "fasta")]
# Codify the alphabet {A,G,C,T} -> {0,1,2,3}
leave_sequences = []
for i in leave_alignments:
    seq = ""
    for j in i.seq:
        if j == "A":
            seq += "0"
        elif j == "G":
            seq += "1"
        elif j == "C":
            seq += "2"
        elif j == "T":
            seq += "3"
        else:
            raise ValueError("Invalid nucleotide")
    leave_sequences.append(seq)
# Store the sequences in a dictionary
real_matrices = dict()
for i in range(len(edges)):
    u = edges[i].edge[0].name.split("_")[1]
    v = edges[i].edge[1].name.split("_")[1]
    matrix_name = "M_" + str(u) + "_to_" + str(v) 
    real_matrices[matrix_name] = np.array(edges[i].transition_matrix)
np.save(directory + "real_matrices", real_matrices) 
np.save(directory + "real_root_distr", node_distr["Int_0"]) 
# We are not interested in this .fasta file anymore
# file_path = directory+filename 
# # Check if the file exists
# try:
#     if os.path.exists(file_path):
#         # Delete the file
#         os.remove(file_path)
# except OSError as error:
#     print("File does not exist.")
#     pass 

###################################### Auxiliary functions ######################################
def init_root_distribution():
    """
    Initialise estimation for the root distribution (all probabilities between 0.2 and 0.3)
    """
    estimated_root_distribution = []
    while True:
        estimated_root_distribution = np.random.dirichlet([1,1,1,1])
        i = 0
        while i<4:
            if estimated_root_distribution[i] > 0.2 and estimated_root_distribution[i] < 0.3:
                i += 1
            else:
                break
        if i == 4:
            assert(np.sum(estimated_root_distribution) < 1.0000001 and np.sum(estimated_root_distribution) > 0.9999999)
            assert(estimated_root_distribution[i] > 0.2 and estimated_root_distribution[i] < 0.3 for i in range(4))
            break
    return estimated_root_distribution

def init_transition_matrix():
    """
    Initialise estimation for the transition matrices (which we estimate as LARGE-DIAGONAL-VALUE matrices)
    """
    M = np.zeros((4,4))
    i=0
    while i<4:
        R = np.random.dirichlet([1,1,1,1])
        if R[i] > 0.6 and R[i] < 0.8:
            M[i,:] = R 
            i = i + 1
    assert (M[i,i] > 0.6 and M[i,i] < 0.8 for i in range(4))
    assert (np.sum(M[i,:]) < 1.0000001 and np.sum(M[i,:]) > 0.9999999 for i in range(4))
    return M

def log_likelihood(states, u_i, params, root_distribution, n_int):
    """
    Compute the log-likehood of the tree
    """
    logL = 0
    for obs in states:
        p_i = 0 
        observed_data = obs[n_int:]
        if observed_data in u_i:
            for state in states:
                if state[n_int:] == observed_data:
                    pi = root_distribution[int(state[0])]
                    for p in params.items():
                        u, v = int(p[1].edge[0].name.split("_")[1]), int(p[1].edge[1].name.split("_")[1])
                        pi *= p[1].transition_matrix[int(state[u]), int(state[v])]
                    p_i += pi
        logL += (u_i[observed_data] * math.log(p_i))
    return logL

def compute_branch_length(params, node_distr, real_params):
    """
    Compute the branch length of the tree
    """
    debg_node_distr = dict()
    branch_length = dict()
    if real_params == True:
        for i in node_distr.items():
            debg_node_distr[i[0].split("_")[1]] = i[1]
        for p in params.items():
            u, v = int(p[0].split("_")[1]), int(p[0].split("_")[3])
            branch_length[f"M_{u}_to_{v}"] = (1/4)*(-math.log((np.sqrt(np.abs(np.linalg.det(np.diag(debg_node_distr[str(u)]))))*np.abs(np.linalg.det(p[1])))/(np.sqrt(np.abs(np.linalg.det(np.diag(debg_node_distr[str(v)])))))))
    else:
        debg_node_distr["0"] = node_distr
        for p in params.items():
            u, v = int(p[0].split("_")[1]), int(p[0].split("_")[3])
            debg_node_distr[str(v)] = np.matmul(debg_node_distr[str(u)], p[1].transition_matrix)
            branch_length[f"M_{u}_to_{v}"] = (1/4)*(-math.log((np.sqrt(np.abs(np.linalg.det(np.diag(debg_node_distr[str(u)]))))*np.abs(np.linalg.det(p[1].transition_matrix)))/(np.sqrt(np.abs(np.linalg.det(np.diag(debg_node_distr[str(v)])))))))
    return branch_length

def get_fasta_alignment(params, node_distr, directory, length): 
    """
    Generate a fasta file with the alignment of the given parameters.
    """          
    node_sequence = dict()
    node_sequence["Int_0"] = generate_alignment(length, node_distr["Int_0"])
    for p in params.items():
        node_sequence[p[1].edge[1].name] = generate_sequences(p[1].transition_matrix, node_sequence[p[1].edge[0].name])
    leaves_seq = {k: v for k, v in node_sequence.items() if k.startswith('L')}
    sequences_in_leaves = list(leaves_seq.values())
    keys_for_sequences = list(leaves_seq.keys())
    iter = 0
    file_name = str(len(sequences_in_leaves))+ "_leaves_" + str(len(sequences_in_leaves[0])) + "length_sequences.fasta"
    file = open(directory+file_name, "w")
    for seq in sequences_in_leaves:
        file.write(">Seq" + str(keys_for_sequences[iter]) + "\n" + seq + "\n")
        iter += 1
    file.close() 
    return file_name 

class Param:
    """
    Class to store the parameters of the tree
    """
    def __init__(self, edge, transition_matrix, alignment=None):
        self.edge = edge
        self.transition_matrix = transition_matrix
        self.alignment = alignment # it will only be placed in the leaves

# Obtain number of leaves and internal nodes
n_leaves, n_int = 0, 0 
for node in net.nodes():
    if node.name.startswith("L"):
        n_leaves += 1
    else:
        n_int += 1

# Define all possible (hidden) states of the internal nodes 
hidden_combinations = itertools.product([0,1,2,3], repeat=n_int)
sorted_hidden_combinations = sorted(hidden_combinations)
for i in range(len(sorted_hidden_combinations)):
    sorted_hidden_combinations[i] = ''.join([str(s) for s in sorted_hidden_combinations[i]])

# Branch length of the real tree
bl = compute_branch_length(real_matrices, node_distr, True)
np.save(directory + "real_branch_lengths", bl) 

PARAMS_init = dict() # Dictionary to store the estimation matrices for each edge
for edge in net.edges():
    # If we are in a leaf, we need to specify the alignment
    iter_leaves = 0
    u = str(edge[0].name.split("_")[1])
    v = str(edge[1].name.split("_")[1])
    matrix_name = "M_" + str(u) + "_to_" + str(v)
    if edge[1].name.startswith("L"):
        new_edge = Param(edge, real_matrices[matrix_name], leave_alignments[iter_leaves].seq)
        iter_leaves += 1
    # Otherwise, are in internal nodes
    else:
        new_edge = Param(edge, real_matrices[matrix_name])
    PARAMS_init[matrix_name] = new_edge

################################# EXPECTATION-MAXIMIZATION ALGORITHM #################################
print("running EM...")
print("---"*10)
for r in range(repetitions):
    try:
        rep_directory = results + "repetition_" + str(r+1) + "/"
        os.mkdir(rep_directory)
    except FileExistsError:
        pass
    start_time = time.time()
    ################################################################################################
    ################################ STEP 0: Initialise parameters #################################
    ################################################################################################
    estimated_root_distribution = init_root_distribution() # First estimation for the root distribution

    PARAMS = dict() # Dictionary to store the estimation matrices for each edge
    for edge in net.edges():
        # If we are in a leaf, we need to specify the alignment
        iter_leaves = 0
        if edge[1].name.startswith("L"):
            new_edge = Param(edge, init_transition_matrix(), leave_alignments[iter_leaves].seq)
            iter_leaves += 1
        # Otherwise, are in internal nodes
        else:
            new_edge = Param(edge, init_transition_matrix())
        u = str(edge[0].name.split("_")[1])
        v = str(edge[1].name.split("_")[1])
        name = "M_" + u + "_to_" + v
        PARAMS[name] = new_edge

    params = copy.copy(PARAMS)      # copy of the parameters (will be updated at each iteration)
    
    # Count all the leaves ocurrences 
    sequence_length = len(leave_sequences[0])
    number_of_sequences = len(leave_sequences)
    occurrences = []
    for i in range(sequence_length):
        _ = ''
        for j in range(number_of_sequences):
            _ += leave_sequences[j][i]
        occurrences.append(_)
    c = Counter(occurrences)
    u_i = OrderedDict(sorted(c.items()))

    file_name = get_fasta_alignment(PARAMS_init, node_distr, rep_directory, length)
    # Preprocess leave alignments
    leave_alignments = [i for i in SeqIO.parse(rep_directory+filename, "fasta")]
    # Codify the alphabet {A,G,C,T} -> {0,1,2,3}
    leave_sequences = []
    for i in leave_alignments:
        seq = ""
        for j in i.seq:
            if j == "A":
                seq += "0"
            elif j == "G":
                seq += "1"
            elif j == "C":
                seq += "2"
            elif j == "T":
                seq += "3"
            else:
                raise ValueError("Invalid nucleotide")
        leave_sequences.append(seq)

    # Define all possible fully observed states (i.e. combining u_i + internal nodes)
    # '000000' would mean to have an 'A' in all nodes: Int_0, Int_1, Leaf_2, Leaf_3, Leaf_4, Leaf_5 
    states = list(itertools.product(list(sorted_hidden_combinations), list(u_i.keys())))
    states = [i[0]+i[1] for i in states]
    
    iter = 0                        # iteration counter 
    logL = log_likelihood(states, u_i, params, estimated_root_distribution, n_int)
    logL_ = 0
    loglikeligood = []
    loglikeligood.append(logL)

    # Log-likelihood of real parameters
    logL_real_parameters = log_likelihood(states, u_i, PARAMS_init, estimated_root_distribution, n_int)
    if r == 0:
        out = open(directory+"Loglikelihood_real_param.txt", "w", encoding='latin-1')
        out.write(f"{logL_real_parameters}")
        out.close() 

    while np.abs(logL_ - logL) > eps and iter < 100:
        if iter > 0:
            logL = logL_

        ##########################################################################################
        ######################################### E-step #########################################
        ##########################################################################################
        # Define the expected hidden data matrix
        U = dict()
        root_distr = estimated_root_distribution 

        for obs in states:
            U[obs] = u_i[obs[n_int:]] # initialise by computing u_i of the expected hidden data matrix formula
            f_ij = estimated_root_distribution[int(obs[0])] 
            for p in params.items():
                u, v = int(p[1].edge[0].name.split("_")[1]), int(p[1].edge[1].name.split("_")[1])
                f_ij *= p[1].transition_matrix[int(obs[u]), int(obs[v])]
            # Now f_ij := numerator of the expected hidden data matrix formula 

            # We need to marginalise to obtain f_i := denominator of the expected hidden data matrix formula
            f_i = 0 # we initialise the partial data model
            observed_data = obs[n_int:]
            for state in states:
                if state[n_int:] == observed_data:
                    pi = estimated_root_distribution[int(state[0])]
                    for p in params.items():
                        u, v = int(p[1].edge[0].name.split("_")[1]), int(p[1].edge[1].name.split("_")[1])
                        pi *= p[1].transition_matrix[int(state[u]), int(state[v])]
                    f_i += pi

            U[obs] *= f_ij/f_i

        assert(len(U) <= 4**(n_int+n_leaves))
        sum = 0
        sum2 = 0
        for i in U.items():
            sum += i[1]
        for i in u_i.items():
            sum2 += i[1]
        assert(sum < sequence_length+0.001 and sum > sequence_length-0.001)
        assert(sum2 < sequence_length+0.001 and sum2 > sequence_length-0.001)

        ########################################################################################
        ######################################## M-step ########################################
        ########################################################################################
        estimated_parameters = dict()
        root_distr = np.zeros(4)
        iter_root = 0 # counter for the root distribution, we only want to count once

        for p in params.items():
            u, v = int(p[1].edge[0].name.split("_")[1]), int(p[1].edge[1].name.split("_")[1])
            name = "M_" + str(u) + "_to_" + str(v)
            estimated_parameters[name] = np.zeros((4,4))
            for i in range(4):
                for j in range(4):
                    ui = 0
                    for obs in U.items():
                        # For the root distribution, we only need to consider the first internal node
                        if iter_root == 0 and obs[0][0] == str(i):
                            root_distr[i] += obs[1]
                        # For the transition matrices:
                        if obs[0][u] == str(i) and obs[0][v] == str(j):
                            estimated_parameters[name][i][j] += obs[1]
                        if obs[0][u] == str(i):
                            ui += obs[1]
                    # Update root distribution
                    if iter_root == 0:
                        root_distr[i] /= sequence_length
                    iter_root = 1
                    estimated_parameters[name][i][j] /= ui
                iter_root = 0
                
            # Update parameters
            params[name].transition_matrix = estimated_parameters[name]

        ########################################################################################
        ################## Compute log-likelihood of the estimated parameters ##################
        ########################################################################################
        logL_ = log_likelihood(states, u_i, params, root_distr, n_int)
        loglikeligood.append(logL_)

        iter += 1
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Repetition {r+1} done!")
    print("Number of iterations: ", iter)
    print("Execution time: ", elapsed_time)
    print("---"*10)

    # Save the results
    out = open(rep_directory+"Niter.txt", "w", encoding='latin-1')
    out.write(f"{iter}")
    out.close() 
    out = open(rep_directory+"TExec.txt", "w", encoding='latin-1')
    out.write(f"{elapsed_time}")
    out.close() 
    M_estimation = dict()
    for p in params.items():
        M_estimation[p[0]] = p[1].transition_matrix
    np.save(rep_directory+"M_estimation.npy", M_estimation)
    np.save(rep_directory+"root_estimation.npy", root_distr)
    np.save(rep_directory+"loglikelihood.npy", loglikeligood)
    bl = compute_branch_length(params, root_distr, False)
    np.save(rep_directory + "estimated_branch_lengths", bl) 

print("Done!")

################################# EXPECTATION-MAXIMIZATION ALGORITHM #################################
# Create a directory to store the results
try:
    directory = "./OUTPUTS/"
    os.mkdir(directory)
except FileExistsError:
    pass
try:
    directory += "SPECTRAL_vs_RANDOM_EM_" + str(length) + "_a_b_10/" 
    os.mkdir(directory)
except FileExistsError:
    pass
try:
    results = directory + "RESULTS/" 
    os.mkdir(results)
except FileExistsError:
    pass
print("running EM...")
print("---"*10)
for r in range(repetitions):
    try:
        rep_directory = results + "repetition_" + str(r+1) + "/"
        os.mkdir(rep_directory)
    except FileExistsError:
        pass
    start_time = time.time()
    ################################################################################################
    ################################ STEP 0: Initialise parameters #################################
    ################################################################################################
    estimated_root_distribution = init_root_distribution() # First estimation for the root distribution

    Results = init_EM(newick, length, save=True, directory=directory)
    i=0
    PARAMS = dict() # Dictionary to store the estimation matrices for each edge
    for edge in net.edges():
        # If we are in a leaf, we need to specify the alignment
        iter_leaves = 0
        if edge[1].name.startswith("L"):
            new_edge = Param(edge, Results[i].matrix, leave_alignments[iter_leaves].seq)
            iter_leaves += 1
        # Otherwise, are in internal nodes
        else:
            new_edge = Param(edge, Results[i].matrix)
        u = str(edge[0].name.split("_")[1])
        v = str(edge[1].name.split("_")[1])
        name = "M_" + u + "_to_" + v
        PARAMS[name] = new_edge
        i += 1
    # # Check if the file exists
    # try:
    #     file_name = f"{n_leaves}_leaves_{length}length_sequences.fasta"
    #     if os.path.exists(directory+file_name):
    #         # Delete the file
    #         os.remove(file_path)
    # except OSError as error:
    #     print("File does not exist.")
    #     pass 
    params = copy.copy(PARAMS)      # copy of the parameters (will be updated at each iteration)
    
    # Count all the leaves ocurrences 
    sequence_length = len(leave_sequences[0])
    number_of_sequences = len(leave_sequences)
    occurrences = []
    for i in range(sequence_length):
        _ = ''
        for j in range(number_of_sequences):
            _ += leave_sequences[j][i]
        occurrences.append(_)
    c = Counter(occurrences)
    u_i = OrderedDict(sorted(c.items()))

    file_name = get_fasta_alignment(PARAMS_init, node_distr, rep_directory, length)
    # Preprocess leave alignments
    leave_alignments = [i for i in SeqIO.parse(rep_directory+filename, "fasta")]
    # Codify the alphabet {A,G,C,T} -> {0,1,2,3}
    leave_sequences = []
    for i in leave_alignments:
        seq = ""
        for j in i.seq:
            if j == "A":
                seq += "0"
            elif j == "G":
                seq += "1"
            elif j == "C":
                seq += "2"
            elif j == "T":
                seq += "3"
            else:
                raise ValueError("Invalid nucleotide")
        leave_sequences.append(seq)

    # Define all possible fully observed states (i.e. combining u_i + internal nodes)
    # '000000' would mean to have an 'A' in all nodes: Int_0, Int_1, Leaf_2, Leaf_3, Leaf_4, Leaf_5 
    states = list(itertools.product(list(sorted_hidden_combinations), list(u_i.keys())))
    states = [i[0]+i[1] for i in states]
    
    iter = 0                        # iteration counter 
    logL = log_likelihood(states, u_i, params, estimated_root_distribution, n_int)
    logL_ = 0
    loglikeligood = []
    loglikeligood.append(logL)

    # Log-likelihood of real parameters
    logL_real_parameters = log_likelihood(states, u_i, PARAMS_init, estimated_root_distribution, n_int)
    if r == 0:
        out = open(directory+"Loglikelihood_real_param.txt", "w", encoding='latin-1')
        out.write(f"{logL_real_parameters}")
        out.close() 

    while np.abs(logL_ - logL) > eps and iter < 100:
        if iter > 0:
            logL = logL_

        ##########################################################################################
        ######################################### E-step #########################################
        ##########################################################################################
        # Define the expected hidden data matrix
        U = dict()
        root_distr = estimated_root_distribution 

        for obs in states:
            U[obs] = u_i[obs[n_int:]] # initialise by computing u_i of the expected hidden data matrix formula
            f_ij = estimated_root_distribution[int(obs[0])] 
            for p in params.items():
                u, v = int(p[1].edge[0].name.split("_")[1]), int(p[1].edge[1].name.split("_")[1])
                f_ij *= p[1].transition_matrix[int(obs[u]), int(obs[v])]
            # Now f_ij := numerator of the expected hidden data matrix formula 

            # We need to marginalise to obtain f_i := denominator of the expected hidden data matrix formula
            f_i = 0 # we initialise the partial data model
            observed_data = obs[n_int:]
            for state in states:
                if state[n_int:] == observed_data:
                    pi = estimated_root_distribution[int(state[0])]
                    for p in params.items():
                        u, v = int(p[1].edge[0].name.split("_")[1]), int(p[1].edge[1].name.split("_")[1])
                        pi *= p[1].transition_matrix[int(state[u]), int(state[v])]
                    f_i += pi

            U[obs] *= f_ij/f_i

        assert(len(U) <= 4**(n_int+n_leaves))
        sum = 0
        sum2 = 0
        for i in U.items():
            sum += i[1]
        for i in u_i.items():
            sum2 += i[1]
        assert(sum < sequence_length+0.001 and sum > sequence_length-0.001)
        assert(sum2 < sequence_length+0.001 and sum2 > sequence_length-0.001)

        ########################################################################################
        ######################################## M-step ########################################
        ########################################################################################
        estimated_parameters = dict()
        root_distr = np.zeros(4)
        iter_root = 0 # counter for the root distribution, we only want to count once

        for p in params.items():
            u, v = int(p[1].edge[0].name.split("_")[1]), int(p[1].edge[1].name.split("_")[1])
            name = "M_" + str(u) + "_to_" + str(v)
            estimated_parameters[name] = np.zeros((4,4))
            for i in range(4):
                for j in range(4):
                    ui = 0
                    for obs in U.items():
                        # For the root distribution, we only need to consider the first internal node
                        if iter_root == 0 and obs[0][0] == str(i):
                            root_distr[i] += obs[1]
                        # For the transition matrices:
                        if obs[0][u] == str(i) and obs[0][v] == str(j):
                            estimated_parameters[name][i][j] += obs[1]
                        if obs[0][u] == str(i):
                            ui += obs[1]
                    # Update root distribution
                    if iter_root == 0:
                        root_distr[i] /= sequence_length
                    iter_root = 1
                    estimated_parameters[name][i][j] /= ui
                iter_root = 0
                
            # Update parameters
            params[name].transition_matrix = estimated_parameters[name]

        ########################################################################################
        ################## Compute log-likelihood of the estimated parameters ##################
        ########################################################################################
        logL_ = log_likelihood(states, u_i, params, root_distr, n_int)
        loglikeligood.append(logL_)

        iter += 1
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Repetition {r+1} done!")
    print("Number of iterations: ", iter)
    print("Execution time: ", elapsed_time)
    print("---"*10)

    # Save the results
    out = open(rep_directory+"Niter.txt", "w", encoding='latin-1')
    out.write(f"{iter}")
    out.close() 
    out = open(rep_directory+"TExec.txt", "w", encoding='latin-1')
    out.write(f"{elapsed_time}")
    out.close() 
    M_estimation = dict()
    for p in params.items():
        M_estimation[p[0]] = p[1].transition_matrix
    np.save(rep_directory+"M_estimation.npy", M_estimation)
    np.save(rep_directory+"root_estimation.npy", root_distr)
    np.save(rep_directory+"loglikelihood.npy", loglikeligood)
    bl = compute_branch_length(params, root_distr, False)
    np.save(rep_directory + "estimated_branch_lengths", bl) 

print("Done!")