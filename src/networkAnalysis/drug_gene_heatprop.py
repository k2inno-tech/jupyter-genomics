# import some useful packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn
import networkx as nx
import pandas as pd
import random
import json


import sys
code_path = 'source'
sys.path.append(code_path)
import network_prop

def load_DB_data(fname):
    '''
    Function to load drug bank data (in format of this file 'drugbank.0.json.new')
    '''

    with open(fname, 'r') as f:
        read_data = f.read()
    f.closed

    si = read_data.find('\'\n{\n\t"source":')
    sf = read_data.find('\ncurl')

    DBdict = dict()

    # fill in DBdict
    while si > 0:
        
        db_temp = json.loads(read_data[si+2:sf-2])
        DBdict[db_temp['drugbank_id']]=db_temp

        # update read_data
        read_data = read_data[sf+10:]
        
        si = read_data.find('\'\n{\n\t"source":')
        sf = read_data.find('\ncurl')
        
    return DBdict

def load_cluster_data(fname):
    # NOTE: THIS ONLY WORKS FOR SYMMETRIC CLUSTERS FOR NOW
    
    # load a cluster in which to run heat propagation
    sample_mat = pd.read_csv(fname,index_col=0)

    idx_to_node = dict(zip(range(len(sample_mat)),list(sample_mat.index)))

    sample_mat = np.array(sample_mat)
    sample_mat = sample_mat[::-1,0:-1] # reverse the indices for use in graph creation
    
    G_cluster = nx.Graph()
    G_cluster = nx.from_numpy_matrix(np.abs(sample_mat))
    G_cluster = nx.relabel_nodes(G_cluster,idx_to_node)
    
    return G_cluster


def find_drugs_from_hot_genes(Fnew,G_DB,seed_genes,keep_seed_genes =True):
    
    '''
    Function to find drugs associated with hot genes, return a dictionary containing drug info and heat rank
    inputs:
        - Fnew: heat vector from network_prop.network_propagation function
        - G_DB: Bipartite graph with drugs and drug targets (genes) as nodes
        - seed_genes:  seed genes network propagation was started from
        - keep_seed_genes:  decide whether to include seed genes in output drug list (default True)
    
    '''
    
    # make sure Fnew is sorted
    Fnew.sort(ascending=False)
    
    ranked_genes = list(Fnew.index)
    
    if not keep_seed_genes:
        # (Should we only keep non-seed genes??)
        ranked_genes = list(np.setdiff1d(ranked_genes,seed_genes)) 
        ranked_genes = Fnew[ranked_genes]
        ranked_genes.sort(ascending=False)
        ranked_genes = list(ranked_genes.index)

    gene_drug_dict = dict() # build up a list of genes and drugs that may be related to input list
    for g in ranked_genes:
        if g in G_DB.nodes():  # check if g is in drugbank graph
            drug_neighs_temp = list(nx.neighbors(G_DB,g))
            
            # add drug neighbors and ranked score to gene_drug_dict
            gene_drug_dict[g] = {'drugs':drug_neighs_temp,'heat_rank':ranked_genes.index(g)}
            
        else:
            # fill in dictionary when there are no drugs related to focal gene
            gene_drug_dict[g] = {'drugs':[],'heat_rank':ranked_genes.index(g)}
            
    # return a sorted dataframe- more useful
    gene_drug_df = pd.DataFrame(gene_drug_dict).transpose()  # note we have to transpose so index is genes
    gene_drug_df['heat_value'] = Fnew[ranked_genes]
    gene_drug_df = gene_drug_df.sort(columns='heat_rank')
            
    return gene_drug_df
            
            

def drug_gene_heatprop(seed_genes,path_to_DB_file,path_to_cluster_file,plot_flag=False):
    
    '''
    Function to establish drugs potentially related to an input gene list, using network propagation methods
    
    inputs:
        - seed_genes:  genes from which to initiate heat propagation simulation
        - path_to_DB_file:  path to drug bank file, including filename
        - path_to_cluster_file: path to cluster file, including filename
        - plot_flag: should we plot the subnetwork with heat overlaid? Default False
        
    '''
    
    
    # load and parse the drug-bank file into a dict ()
    DBdict = load_DB_data(path_to_DB_file)
    
    # make a network out of drug-gene interactions
    DB_el = []
    for d in DBdict.keys():
        node_list = DBdict[d]['node_list']
        for n in node_list:
            DB_el.append((DBdict[d]['drugbank_id'],n['name']))
            
            
    G_DB = nx.Graph()
    G_DB.add_edges_from(DB_el)
    
    G_cluster = load_cluster_data(path_to_cluster_file)
    
    # calculate the degree-normalized adjacency matrix
    Wprime = network_prop.normalized_adj_matrix(G_cluster,weighted=True)
    
    # run the network_propagation simulation starting from the seed genes
    Fnew = network_prop.network_propagation(G_cluster,Wprime,seed_genes)
    
    # sort heat vector Fnew
    Fnew.sort(ascending=False)
    
    # if plot_flag is on plot the cluster genes with heat overlaid
    if plot_flag:
        pos = nx.spring_layout(G_cluster)

        plt.figure(figsize=(10,10))
        nx.draw_networkx_edges(G_cluster,pos=pos,alpha=.03)
        nx.draw_networkx_nodes(G_cluster,pos=pos,node_size=20,alpha=.8,node_color=Fnew[G_cluster.nodes()],cmap='jet',
                               vmin=0,vmax=np.max(Fnew)/10)
        nx.draw_networkx_nodes(G_cluster,pos=pos,nodelist=seed_genes,node_size=50,alpha=.7,node_color='red',linewidths=2)

        plt.grid('off')
        plt.title('Sample subnetwork: post-heat propagation',fontsize=16)
    
    # find the drugs related to hot genes
    gene_drug_df = find_drugs_from_hot_genes(Fnew,G_DB,seed_genes,keep_seed_genes =True)
    
    return gene_drug_df
    
    