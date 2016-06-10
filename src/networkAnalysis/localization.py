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
import random 
import numpy as np
import optparse
import sys

import separation as tools




# =================================================================================
def get_lcc_size(G,seed_nodes):
    """
    return the lcc size
    """

    # getting subgraph that only consists of the black_nodes
    g = nx.subgraph(G,list(seed_nodes))

    if g.number_of_nodes() != 0:
        # get all components 
        max_CC = max(nx.connected_component_subgraphs(g), key=len)
        return len(max_CC.nodes())  # size of largest connected component"

    else:
        return 0
    
# test connectivity of focal node set compared to random node set

def measure_connectivity(G,focal_nodes,method='jaccard',num_reps=10):
    avg_focal_degree = np.mean(G.degree(focal_nodes).values())
    var_focal_degree = np.std(G.degree(focal_nodes).values())
    
    
    
    if method=='jaccard':
        
        #focal_sim = jaccard_sim_array(G,focal_nodes)
        
        # select random node set
        rand_sim=[]
        focal_sim=[]
        for r in range(num_reps):
            
            # bootstrap a sample from focal_nodes
            focal_subset = np.random.choice(list(focal_nodes),size=len(focal_nodes),replace=True)
            
            focal_sim.append(jaccard_sim_array(G,focal_subset))
            
            print('calculating random set ' + str(r) + ' out of ' + str(num_reps))
            G_temp = nx.configuration_model(G.degree().values())
            G_rand = nx.Graph()  # switch from multigraph to digraph
            G_rand.add_edges_from(G_temp.edges())
            # remove self-loops
            #G_rand.remove_edges_from(G_rand.selfloop_edges())
            G_rand = nx.relabel_nodes(G_rand,dict(zip(range(len(G_rand.nodes())),G.degree().keys())))
            
            rand_sim.append(jaccard_sim_array(G_rand,focal_subset)) # measure the jaccard similarity of degree-preserving edge shuffled network
            
    elif method=='edge_overlap':
        
        focal_sim = num_shared_neighbors(G,focal_nodes)
        # select random node set
        rand_sim=[]
        for r in range(num_reps):
            print('calculating random set ' + str(r) + ' out of ' + str(num_reps))
            G_rand = nx.configuration_model(G.degree().values())
            G_rand = nx.relabel_nodes(G_rand,dict(zip(range(len(G_rand.nodes())),G.degree().keys())))
            
            rand_sim.append(num_shared_neighbors(G_rand,focal_nodes)) # rand_sim is array of length num_reps
        
        
    
    return focal_sim,rand_sim
    
def rand_node_set(G,num_nodes,avg_degree,var_degree,tolerance=.1):
    node_list = G.nodes()
    node_shuf = list(node_list)  # make a copy
    
    frac_diff_avg,frac_diff_var = 10,10
    rand_nodes = []
    
    # loop over random nodes- look for ones with similar avg and var degree to focal node set
    while (frac_diff_avg>tolerance) or (frac_diff_var>tolerance):

    
        np.random.shuffle(node_shuf)
        rand_nodes = node_shuf[0:num_nodes]
        avg_degree_temp = np.mean(G.degree(rand_nodes).values())
        var_degree_temp = np.std(G.degree(rand_nodes).values())

        frac_diff_avg = np.abs(1-np.abs(avg_degree_temp-avg_degree)/avg_degree)
        frac_diff_var = np.abs(1-np.abs(var_degree_temp-var_degree)/var_degree)

    return rand_nodes
    

def jaccard_sim_array(G,node_list):
    jsim_list = []
    num_nodes = len(node_list)
    for i in range(num_nodes-1):
        for j in range(i+1,num_nodes):
            nA = node_list[i]
            nB = node_list[j]
            
            jsim = jaccard_sim_AB(G,nA,nB)
            if jsim>-1:
                # only append if node has nonzero num neighbors
                jsim_list.append(jsim)
    
    jsim_avg = np.mean(jsim_list)
    return jsim_avg
        
    
def jaccard_sim_AB(G,nodeA,nodeB):
    neighA = G.neighbors(nodeA)
    neighB = G.neighbors(nodeB)
    if (len(neighA)+len(neighB))>0:
        jsim = len(np.intersect1d(neighA,neighB))/float(len(np.union1d(neighA,neighB)))
    else:
        # set to -1 if node has no neighbors
        jsim=-1
    
    return jsim

def num_shared_neighbors(G,node_list):
    
    edge_count = 0
    for node in node_list:
        edge_count += len(np.intersect1d(G.neighbors(node),node_list))
        
    return edge_count
    
        
