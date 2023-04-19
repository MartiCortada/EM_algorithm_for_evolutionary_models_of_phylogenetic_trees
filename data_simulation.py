# Required packages
from Bio import Phylo
from io import StringIO
from Bio import SeqIO
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve
from scipy.linalg import expm

class Edge:
    def __init__(self, edge, transition_matrix=None):
        self.edge = edge
        self.transition_matrix = transition_matrix

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

def generate_alignment(length, distribution):
    """
    Generates an alignment of length `length` using a given `distribution`.
    """
    seq = ''
    for i in range(length):
        # generate a random sample from the multinomial distribution
        nucleotide = np.random.choice(['A', 'G', 'C', 'T'], p=distribution)
        seq += nucleotide
    return seq

# Return a matrix from a dicrionary
def get_matrix_from_dict(d):
    """
    Return a matrix from a dictionary
    """
    Q2 = np.zeros((4,4))
    coefficients = list(d.values())
    for i in range(4):
        for j in range(4):
            Q2[i,j] = coefficients[i*4+j]
    return Q2

def alpha(new_distribution, Q, i, k):
    return min(1, (new_distribution[k]*Q[k,i])/(new_distribution[i]*Q[i,k]))

def get_M2(new_distribution,d2,l):
    P = np.zeros((4,4))
    iter = True
    while iter:
        Q = np.zeros((4,4))
        i=0
        while i<4:
            # posar el paràmere de la diagoal a 10 i els altres 1
            # 10 * e**(-l)
            dir = np.ones(4)
            # dir[i] = 50*np.exp(-l)
            dir[i] = (50/np.sqrt(l))*np.exp(-l)
            R = np.random.dirichlet(dir)
            if R[i] > 0.3:
                Q[i,:] = R 
                i = i + 1   
        for i in range(4):
            for j in range(4):
                if i == j:
                    sum = 0
                    for k in range(4):
                        if k != i:
                            sum += (Q[i,k] * (1 - alpha(new_distribution,Q,i,k)))
                    P[i,j] = Q[i,i] + sum
                else:
                    P[i,j] = Q[i,j]*alpha(new_distribution,Q,i,j)
        assert (np.abs(np.sum(new_distribution - np.matmul(new_distribution,P)))) < 10**-6
        # Adjust the matrix diagonalising
        vaps, _ = np.linalg.eig(P)
        vaps = sorted(vaps, reverse=True)
        A = symbols('A')
        eq = Eq(-d2+(((1-A)*vaps[1]+A)*((1-A)*vaps[2]+A)*((1-A)*vaps[3]+A)),0)
        sol = solve(eq, A)
        # We only want the real solution between 0 and 1
        # sol = [s for s in sol if s.is_real and s > 0 and s < 1]  ## MARTA, POSEM TOTES LES ENTRE 0 I 1???
        # 10*-22 I com ho considerem? Fem try except o ens el prenem real
        res = 0
        for s in sol:
            if s.is_real and s > 0 and s < 1:
                res = s
                res = np.float64(res)
                P = (1-res)*P + res*np.identity(4)
                iter = False
                break
            elif s.is_complex:
                b = np.imag(s)
                a = sp.re(s)
                if np.abs(b) < 10**-20 and a > 0 and a < 1:
                    res = sp.re(s)
                    res = np.float64(res)
                    P = (1-res)*P + res*np.identity(4)
                    iter = False
                    break
    return P

def generate_random_matrix(distribution, l):
    res = 1
    while res >= 1:
        M1 = np.zeros((4,4))
        i=0
        while i<4:
            # posar el paràmere de la diagoal a 10 i els altres 1
            # 10 * e**(-l)
            dir = np.ones(4)
            # dir[i] = 50*np.exp(-l)
            dir[i] = (50/np.sqrt(l))*np.exp(-l)
            R = np.random.dirichlet(dir)
            if R[i] > 0.3:
                M1[i,:] = R 
                i = i + 1
        #print("M1", M1)
            
        new_distribution = np.matmul(distribution,M1)
        D = np.diag(distribution)
        D_ = np.diag(new_distribution)
        res = np.exp(-l)*np.sqrt(np.linalg.det(D_))/np.sqrt(np.linalg.det(D))
        detM1 = np.linalg.det(M1)
        if detM1 > np.exp(-l)*np.sqrt(np.linalg.det(D_))/np.sqrt(np.linalg.det(D)):
            pass
        else:
            res = 1
    d2 = np.exp(-l)*np.sqrt(np.linalg.det(D_))/(detM1*np.sqrt(np.linalg.det(D)))
    M2 = get_M2(new_distribution,d2,l) ##
    detM2 = np.linalg.det(M2) 
    assert(np.abs(detM2 - d2) < 10**-6)
    M = np.matmul(M1,M2)
    return M

def generate_sequences(M, seq):
    new_seq = ""
    for s in seq:
        if s == "A":
            new_seq += np.random.choice(['A', 'G', 'C', 'T'], p=M[0,:])
        elif s == "G":
            new_seq += np.random.choice(['A', 'G', 'C', 'T'], p=M[1,:])
        elif s == "C":
            new_seq += np.random.choice(['A', 'G', 'C', 'T'], p=M[2,:])
        else:
            new_seq += np.random.choice(['A', 'G', 'C', 'T'], p=M[3,:])
    return new_seq

def simulate(net):
    node_distribution = dict()
    node_distribution["Int_0"] = init_root_distribution()
    node_sequence = dict()
    node_sequence["Int_0"] = generate_alignment(1000, node_distribution["Int_0"])
    iter = 0
    edges = []
    for edge in net.edges():
        # Extract branch length
        l = edge[1].branch_length
        new_edge = Edge(edge, generate_random_matrix(node_distribution[edge[0].name], l))
        edges.append(new_edge)
        node_distribution[edge[1].name] = np.matmul(node_distribution[edge[0].name],new_edge.transition_matrix)
        for i in range(4):
            assert(np.sum(new_edge.transition_matrix[i,:])<1.000000001 and np.sum(new_edge.transition_matrix[i,:])>0.999999999)
        # create alignment for the node
        node_sequence[edge[1].name] = generate_sequences(new_edge.transition_matrix, node_sequence[edge[0].name])
        iter += 1
    assert(iter == len(net.edges()))
    leaves_seq = {k: v for k, v in node_sequence.items() if k.startswith('L')}
    sequences_in_leaves = list(leaves_seq.values())
    keys_for_sequences = list(leaves_seq.keys())
    iter = 0
    file = open("test_4_leaves.fasta", "w")
    for seq in sequences_in_leaves:
        file.write(">Seq" + str(keys_for_sequences[iter]) + "\n" + seq + "\n")
        iter += 1
    file.close()
    return edges