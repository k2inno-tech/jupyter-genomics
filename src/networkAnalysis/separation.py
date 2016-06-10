#! /usr/bin/env python

"""
 -----------------------------------------------------------------------

Based on code by: Joerg Menche

Substantial modifications by: Brin Rosenthal (sbrosenthal@ucsd.edu)

 This code is based on the separation.py created by Joerg Menche in the paper:
    Uncovering Disease-Disease Relationships Through The Human Interactome

    by Joerg Menche, Amitabh Sharma, Maksim Kitsak, Susan Dina
    Ghiassian, Marc Vidal, Joseph Loscalzo & Albert-Laszlo Barabasi
 
 
 

 -----------------------------------------------------------------------
"""


import networkx as nx
import numpy as np
import optparse
import sys


def read_network(network_file):
    """
    Reads a network from an external file.

    * The edgelist must be provided as a tab-separated table. The
    first two columns of the table will be interpreted as an
    interaction gene1 <==> gene2

    * Lines that start with '#' will be ignored
    """

    G = nx.Graph()
    for line in open(network_file,'r'):
        # lines starting with '#' will be ignored
        if line[0]=='#':
            continue
        # The first two columns in the line will be interpreted as an
        # interaction gene1 <=> gene2
        line_data   = line.strip().split('\t')
        node1 = line_data[0]
        node2 = line_data[1]
        G.add_edge(node1,node2)

    print "\n> done loading network:"
    print "> network contains %s nodes and %s links" %(G.number_of_nodes(),
                                                       G.number_of_edges())
    
    return G


# =============================================================================
def read_gene_list(gene_file):
    """
    Reads a list genes from an external file.

    * The genes must be provided as a table. If the table has more
    than one column, they must be tab-separated. The first column will
    be used only.

    * Lines that start with '#' will be ignored
    """

    genes_set = set()
    for line in open(gene_file,'r'):
        # lines starting with '#' will be ignored
        if line[0]=='#':
            continue
        # the first column in the line will be interpreted as a seed
        # gene:
        line_data = line.strip().split('\t')
        gene      = line_data[0]
        genes_set.add(gene)

    print "\n> done reading genes:"
    print "> %s genes found in %s" %(len(genes_set),gene_file)

    return genes_set


# =============================================================================
def remove_self_links(G):

    sl = G.selfloop_edges()
    G.remove_edges_from(sl)





# =============================================================================
def calc_single_set_distance(G,given_gene_set):

    """
    Calculates the mean shortest distance for a set of genes on a
    given network    
    

    PARAMETERS:
    -----------
        - G: network
        - gene_set: gene set for which distance will be computed 

    RETURNS:
    --------
         - mean shortest distance 

    """


    # remove all nodes that are not in the network, just to be safe
    all_genes_in_network = set(G.nodes())
    gene_set = given_gene_set & all_genes_in_network

    # get the network distances for all gene pairs:
    all_path_lenghts = get_pathlengths_for_single_set(G,gene_set)

    all_distances = []

    # going through all gene pairs
    for geneA in gene_set:

        all_distances_A = []
        for geneB in gene_set:

            # I have to check which gene is 'smaller' in order to know
            # where to look up the distance of that pair in the
            # all_path_lengths dict
            if geneA < geneB:
                if all_path_lenghts[geneA].has_key(geneB):
                    all_distances_A.append(all_path_lenghts[geneA][geneB])
            else:
                if all_path_lenghts[geneB].has_key(geneA):
                    all_distances_A.append(all_path_lenghts[geneB][geneA])

        if len(all_distances_A) > 0:
            l_min = min(all_distances_A)
            all_distances.append(l_min)

    # calculate mean shortest distance
    mean_shortest_distance = np.mean(all_distances)

    return all_distances


import time


def calc_pwd_seed_random(G,gene_set_1,gene_set_2,num_sims):
    # SBR added: compare the overlap in gene_set_1, gene_set_2 to random overlap
    
    all_genes = G.nodes()
    
    MSD_seed = calc_set_pair_distances(G,gene_set_1,gene_set_2)
    print('seed mean shortest distance= ' + str(np.mean(MSD_seed)))
    
    num_genes_1 = len(np.intersect1d(list(gene_set_1),all_genes))
    num_genes_2 = len(np.intersect1d(list(gene_set_2),all_genes))
    
    MSD_rand = []
    for s in range(num_sims):

        print('calculating random set ' + str(s) + ' out of ' + str(num_sims))
        G_rand = nx.configuration_model(G.degree().values())
        G_rand = nx.relabel_nodes(G_rand,dict(zip(range(len(G_rand.nodes())),G.degree().keys())))
        
        MSD_rand_temp = calc_set_pair_distances(G_rand,gene_set_1,gene_set_2)
        print('rand mean shortest distance= ' + str(np.mean(MSD_rand_temp)))
        MSD_rand.extend(MSD_rand_temp)
        
    
    return MSD_seed, MSD_rand

def calc_sAB(G,gene_set_1,gene_set_2):
    # SBR: calculate the set distance for two gene sets
    dAA = np.mean(calc_single_set_distance(G,gene_set_1))
    dBB = np.mean(calc_single_set_distance(G,gene_set_2))
    
    genes_1_2 = gene_set_1.union(gene_set_2)
    dAB = np.mean(calc_set_pair_distances(G,gene_set_1,gene_set_2))
    
    sAB = dAB-(dAA+dBB)/2
    return sAB

def calc_sAB_seed_random(G,gene_set_1,gene_set_2,num_sims,printflag=True):
    # SBR added: compare the set distance in gene_set_1, gene_set_2 to random set distances
    
    all_genes = G.nodes()
    
    sAB_seed = calc_sAB(G,gene_set_1,gene_set_2)
    if printflag:
        print('seed gene sets separation = ' + str(sAB_seed))
    
    num_genes_1 = len(np.intersect1d(list(gene_set_1),all_genes))
    num_genes_2 = len(np.intersect1d(list(gene_set_2),all_genes))
    
    sAB_rand = []
    
    for s in range(num_sims):
        t0 = time.time()
        if printflag:
            print(s)
        # get two random gene sets of the same size as num_genes_1 and num_genes_2
        rand_seeds_1 = set(random.sample(all_genes,num_genes_1))
        rand_seeds_2 = set(random.sample(all_genes,num_genes_2))
        
        sAB_rand_temp = calc_sAB(G,rand_seeds_1,rand_seeds_2)
        if printflag:
            print('rand gene sets separation = ' + str(sAB_rand_temp))
            print(time.time()-t0)
        sAB_rand.append(sAB_rand_temp)
        
        
        
    
    return sAB_seed, sAB_rand


