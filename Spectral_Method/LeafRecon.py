import copy
import numpy as np
import networkx as nx
import operator
import random
from eigenvalue_decomposition import eigenvalue_decomposition
from leaves_data import leaf_distribution
from leaves_data import matrix_from_leaves
from my_classes import MM
from my_classes import SepReconInput
from typing import NamedTuple

# Recover a mutation matrix previously computed
def obtain_matrix(Results, source, target):
    for matrix in Results:
        if matrix.source == source and matrix.target==target:
            return matrix.matrix

# Algorithm: LeafRecon
def LeafRecon(leaves_sequences, tree, a, U, S, net, Results, SR_input, Internal_Nodes_Distr):
    #print("a: ", a)
    W = {a} #Current leaf
    #print("W: ", W)
    U = U.union({a}) # Set of uncovered leaves

    #Tree's depth
    Delta = 0
    for clade in tree.get_terminals():
        p=tree.get_path(clade)
        Delta=max(Delta, len(p))
    #print("Delta: ", Delta)

    # Subset of nodes at distance at most Delta from leaf a
    B = set()
    for node in net.nodes():
        sp=nx.shortest_path_length(net, source=a, target=node)
        if (sp > 0 and sp <= Delta):
            B.add(node)
    #print("B :", B)

    # Outside neighbors of W without using edges in S
    N = set()
    aux_net = nx.Graph(net)
    for edge in S:
        aux_net.remove_edge(edge[0], edge[1])
    for node in nx.node_boundary(aux_net, W):
        N.add(node)
    #print("N :", N)

    while B.intersection(N) != set():
        intersect = B.intersection(N)
        #print("Intersection: ", intersect)

        # Pick any vertex in the intersection computed
        r = random.choice(tuple(intersect))
        #print("r: ", r)

        # Let (ro, r) be the edge with endpoint r in the path from a to r
        path = nx.shortest_path(net,a,r)
        r0 = path[-2]
        #print("r0: ", r0)

        # If r0 = a, we set P_a_ro to be the identity matrix
        if (r0 == a):
            Results.append(MM(a, r0, np.identity(4)))

        # If r is not a leaf
        if (net.degree(r) > 1):
            #print("r no es una fulla")

            # let (r,r1) and (r, r2) be the other two edges out of r
            aux_net = nx.Graph(net)
            aux_net.remove_edge(r,r0)
            r1 = list(aux_net.edges(r))[0][1]
            r2 = list(aux_net.edges(r))[1][1]
            #print ("r1: ", r1)
            #print ("r2: ", r2)

            # Find the closest leaf b connected to r1 by a path not using (r0, r)
            aux_net = nx.Graph(net)
            aux_net.remove_edge(r0,r)
            sp = nx.shortest_path_length(aux_net, r1)
            sp = {key:val for key, val in sp.items() if aux_net.degree(key) == 1}
            sp = sorted(sp.items(), key=operator.itemgetter(1))
            b = sp[0][0]
            #print("b: ", b)

            # Find the closest leaf c connected to r2 by a path not using (r0, r)
            aux_net = nx.Graph(net)
            aux_net.remove_edge(r0,r)
            sp = nx.shortest_path_length(aux_net, r2)
            sp = {key:val for key, val in sp.items() if aux_net.degree(key) == 1}
            sp = {key:val for key, val in sp.items() if key != b}
            sp = sorted(sp.items(), key=operator.itemgetter(1))
            c = sp[0][0]
            #print("c: ", c)

            # Eigenvalue decomposition to obtain to P_a_r
            P_a_r = eigenvalue_decomposition(leaves_sequences, a, b, c)
            Results.append(MM(a, r, P_a_r))
            P_a_r = obtain_matrix(Results, a, r)

        # Otherwise r is a leaf
        else:
            #print("r es una fulla")
            P_a_r=matrix_from_leaves(leaves_sequences, a, r)
            #Results.append(MM(a, r, P_a_r))
            U = U.union({r})
            #print("U: ", U)

        # Compute P_ro_r
        P_a_r0 = obtain_matrix(Results, a, r0)
        P_r0_r = np.matmul(np.linalg.inv(P_a_r0), P_a_r)
        Results.append(MM(r0, r, P_r0_r))

        # Compute  π_r
        π_a = leaf_distribution(leaves_sequences, a)[0]
        π_r = np.matmul(π_a, P_a_r)
        if (net.degree(r) > 1):
            Internal_Nodes_Distr[r] = π_r
        #print("π_a", π_a)
        #print("π_r", π_r)

        # Obtain P_r_ro
        if (net.degree(r0) == 1):
            π_r0 = leaf_distribution(leaves_sequences, r0)[0]
        else:
            π_r0 = Internal_Nodes_Distr[r0]
        P_r_r0 = np.matmul(np.matmul(np.linalg.inv(np.diag(π_r)),np.matrix.transpose(P_r0_r)), np.diag(π_r0))
        Results.append(MM(r, r0, P_r_r0))

        # Obtain P_r_a (if r0 != a)
        if (r0 != a):
            P_r_a = np.matmul(np.matmul(np.linalg.inv(np.diag(π_r)),np.matrix.transpose(P_a_r)), np.diag(π_a))
            Results.append(MM(r, a, P_r_a))

        # Update W
        W = W.union({r})
        #print(W)

        # If r is not a leaf and the distance between a and r is Delta
        if (net.degree(r) > 1 and nx.shortest_path_length(net, source=a, target=r) == Delta):
            if ((r, r1) not in S):
                t_prime = 1 + len(S)
                SR_input.append(SepReconInput(t_prime, b, a, r1, r))
                S.add((r,r1))
            if ((r, r2) not in S):
                t_prime = 1 + len(S)
                SR_input.append(SepReconInput(t_prime, c, a, r2, r))
                S.add((r,r2))

        #print("S: ", S)

        # Recompute B (subset of nodes at distance at most Delta from leaf a)
        B = set()
        for node in net.nodes():
            sp=nx.shortest_path_length(net, source=a, target=node)
            if (sp > 0 and sp <= Delta):
                B.add(node)
        #print("B :", B)

        # Recompute N (Outside neighbors of W without using edges in S)
        N = set()
        aux_net = nx.Graph(net)
        for edge in S:
            aux_net.remove_edge(edge[0], edge[1])
        for node in nx.node_boundary(aux_net, W):
            N.add(node)
        #print("N :", N)

        for MutationMatrix in Results:
            pass
            #print(MutationMatrix.source, MutationMatrix.target)
            #print(MutationMatrix.matrix)
            #print("--------------------------------------------------")

        #print("LeafRecon loop finished")

    return U, Results, SR_input, Internal_Nodes_Distr
