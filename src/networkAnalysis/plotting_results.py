"""
 -----------------------------------------------------------------------

Author: Brin Rosenthal (sbrosenthal@ucsd.edu)

 -----------------------------------------------------------------------
"""


import scipy
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import community
import seaborn
import separation
import localization
import random
import mygene

# latex rendering of text in graphs
import matplotlib as mpl
mpl.rc('text', usetex = False)
mpl.rc('font', family = 'serif')


# compare the largest connected components of individual diseases to num_reps randomly selected gene sets
def compare_LCC(disease_nodes, Gint, num_reps=10):
    G_disease = nx.subgraph(Gint,list(disease_nodes))
    LCC_disease = max(nx.connected_component_subgraphs(G_disease), key=len)
    LCC_size_disease = len(LCC_disease.nodes())
    
    # get random distribution
    LCC_size_rand = []
    for i in range(num_reps):
        print('calculating random set ' + str(i) + ' out of ' + str(num_reps))
        G_temp = nx.configuration_model(Gint.degree().values())
        G_rand = nx.Graph()  # switch from multigraph to digraph
        G_rand.add_edges_from(G_temp.edges())
        # remove self-loops
        #G_rand.remove_edges_from(G_rand.selfloop_edges())
        G_rand = nx.relabel_nodes(G_rand,dict(zip(range(len(G_rand.nodes())),Gint.degree().keys())))
        
        rand_seeds = disease_nodes  #set(random.sample(Gint.nodes(),len(G_disease.nodes())))
        G_rand_small = nx.subgraph(G_rand,list(rand_seeds))

        # get rand lcc
        LCC_rand = max(nx.connected_component_subgraphs(G_rand_small), key=len)
        LCC_size_rand.append(len(LCC_rand.nodes()))
        
    return LCC_size_disease, LCC_size_rand


# compare the average shortest distances between gene pairs
def compare_SD(disease_nodes, Gint, num_reps=10):

    SD_disease = separation.calc_single_set_distance(Gint,disease_nodes)

    
    # get random distribution
    SD_rand = []
    for i in range(num_reps):
        print('calculating random set ' + str(i) + ' out of ' + str(num_reps))
        
        G_temp = nx.configuration_model(Gint.degree().values())
        G_rand = nx.Graph()  # switch from multigraph to digraph
        G_rand.add_edges_from(G_temp.edges())
        # remove self-loops
        #G_rand.remove_edges_from(G_rand.selfloop_edges())
        G_rand = nx.relabel_nodes(G_rand,dict(zip(range(len(G_rand.nodes())),Gint.degree().keys())))
        
        rand_seeds = disease_nodes #set(random.sample(Gint.nodes(),len(disease_nodes)))
        
        # get random shortest distances
        SD_rand.extend(separation.calc_single_set_distance(G_rand,rand_seeds))

        
    return SD_disease, SD_rand

# plot the localization of a disease
def plot_disease_localization(SD_disease,SD_rand,disease_name='temp'):
    #n_rand,x_rand = np.histogram(SD_rand,bins=(np.unique(SD_rand)))
    n_rand = np.bincount(SD_rand)
    x_rand = range(np.amax(SD_rand)+1)
    x_rand = list(x_rand)
    n_rand_norm = [float(n)/np.sum(n_rand) for n in n_rand]

    #n_disease,x_disease = np.histogram(SD_disease,bins=(np.unique(SD_disease)))
    n_disease = np.bincount(SD_disease)
    x_disease = range(np.amax(SD_disease)+1) #np.unique(SD_disease)
    x_disease = list(x_disease)
    n_disease_norm = [float(n)/np.sum(n_disease) for n in n_disease]

    # calculate p-value that two diseases are significantly different
    U,p = scipy.stats.mannwhitneyu(SD_disease,SD_rand)
    psig = nsf(p,2)
    

    plt.plot(x_rand,n_rand_norm,'o-',label='random')
    plt.plot(x_disease,n_disease_norm,'ro-',label=disease_name + ': \n p = ' + str(psig))
    # also plot vertical lines for means
    plt.axvline(x=np.mean(SD_disease),ymin=0,ymax=1,color='r',linestyle='--')
    plt.axvline(x=np.mean(SD_rand),ymin=0,ymax=1,color='b',linestyle='--')
    plt.xlabel('Shortest distance ($d_s$)',fontsize=14)
    plt.ylabel('$P(d_s)$',fontsize=14)
    plt.legend(fontsize=12,loc='best')
    
# plot comparison of rand LCC to disease LCC
def plot_disease_LCC(LCC_size_disease, LCC_size_rand,disease_name='temp'):
    n_rand,x_rand = np.histogram(LCC_size_rand,bins=15)
    n_rand = list(n_rand)
    x_rand = list(x_rand[0:-1])
    n_rand_norm = [float(n)/np.sum(n_rand) for n in n_rand]
    
    # get the z-score
    z_LCC = (LCC_size_disease-np.mean(LCC_size_rand))/np.std(LCC_size_rand)
    z_2digits = nsf(z_LCC,n=2)
    
    x_total = LCC_size_rand
    x_total.append(LCC_size_disease)
    width = (max(x_total)-min(x_total))/30.0
    
    plt.bar(x_rand,n_rand_norm,width=width,label='random')
    plt.bar(LCC_size_disease+width,1,width=width/1.5,color='r',label=disease_name + ': z-score = ' + str(z_2digits))
    plt.xlabel('largest connected component',fontsize=16)
    plt.ylabel('frequency',fontsize=16)
    
    plt.ylim([0,max(n_rand_norm)+.1])
    
    plt.legend(fontsize=12,loc='best')

# function to return significant digits in exp form
def nsf(num, n=1):
    """n-Significant Figures"""
    numstr = ("{0:.%ie}" % (n-1)).format(num)
    return float(numstr)


# plot the distributions of within disease distances and between disease distances
def plot_d1_d2_distributions(MSD_both,MSD_1,MSD_2,d1name='d1',d2name='d2',saveflag=False,savefilename='temp.png'):
    plt.figure(figsize=(7,5))
    n_2 = plt.hist(MSD_2,bins = [1,2,3,4,5])
    n_1 = plt.hist(MSD_1,bins=[1,2,3,4,5])
    n_both = plt.hist(MSD_both,bins=[0,1,2,3,4,5])
    x_2 = n_2[1][0:-1]
    x_1 = n_1[1][0:-1]
    x_both = n_both[1][0:-1]
    n_2 = list(n_2[0])
    n_1 = list(n_1[0])
    n_both = list(n_both[0])

    n_both_norm = [n/sum(n_both) for n in n_both]
    n_2_norm = [n/sum(n_2) for n in n_2]
    n_1_norm = [n/sum(n_1) for n in n_1]

    # clear current axes
    plt.cla()

    plt.plot(x_2,n_2_norm,'.-',color='firebrick',markersize=16,label=d2name)
    plt.plot(x_1,n_1_norm,'.-',color='darkorange',markersize=16,label=d1name)
    plt.plot(x_both,n_both_norm,'.-',color='black',markersize=16,label='Pairwise')
    plt.legend()
    plt.xlim(-.1,5)
    plt.ylim(-.1,.75)
    plt.xlabel('Shortest distance $d$ \n ($d_{AA}, d_{BB}$ or $d_{AB}$)',fontsize=16)
    plt.ylabel('Distribution $P(d)$',fontsize=16)

    if saveflag:
        plt.savefig(savefilename,dpi=300,bbox_inches='tight')

        
# this function plots the separation of seed node set pairs compared to random distribution
def plot_sAB_seed_rand(sAB_seed,sAB_rand,dname='temp',loc='best'):
    
    # calculate the z-score following supplemental methods:
    z_sAB = (sAB_seed - np.mean(sAB_rand))/np.std(sAB_rand)
    zsig = nsf(z_sAB,2)
    
    seaborn.boxplot(x=None,y=sAB_rand,)
    #plt.boxplot(MSD_rand)
    plt.plot(0,sAB_seed,'r*',markersize=16,label=dname + ':\n z-score = ' + str(zsig))
    plt.ylim(-.5,.5)
    plt.xlim(-1,1)
    plt.ylabel('separation $s_{AB}$',fontsize=16)
    plt.legend(fontsize=12,loc=loc)

    

# write some useful functions for clustering 
def create_G_neigh(Gint,genes_D1,genes_D2):
    # create the subgraph composed of the genes in disease 1, disease 2, and neighbors of at least one gene in both groups
    
    G_D1 = nx.subgraph(Gint,list(genes_D1))
    G_D2 = nx.subgraph(Gint,list(genes_D2))
    genes_D1_D2 = []
    genes_D1_D2.extend(list(np.intersect1d(list(genes_D1),Gint.nodes())))
    genes_D1_D2.extend(list(np.intersect1d(list(genes_D2),Gint.nodes())))
    G_D1_D2 = nx.subgraph(Gint,list(genes_D1_D2))

    num_D1 = len(np.intersect1d(list(genes_D1),Gint.nodes()))
    num_D2 = len(np.intersect1d(list(genes_D2),Gint.nodes()))

    # find genes which are neighbors of seed genes
    D1_in_graph = list(np.intersect1d(Gint.nodes(),list(genes_D1)))
    D2_in_graph = list(np.intersect1d(Gint.nodes(),list(genes_D2)))

    neigh_list = []
    for gS in list(D1_in_graph):
        neigh_temp_S = Gint.neighbors(gS)
        for gB in list(D2_in_graph):
            neigh_temp_B = Gint.neighbors(gB)
            neigh_overlap = list(np.intersect1d(neigh_temp_S,neigh_temp_B))  # common neighbors to both
            neigh_list.extend(neigh_overlap)

    neigh_list.extend(D1_in_graph)
    neigh_list.extend(D2_in_graph)
    neigh_list = list(np.unique(neigh_list))
    neigh_list_non_SB = np.setdiff1d(neigh_list,D2_in_graph)
    neigh_list_non_SB = np.setdiff1d(neigh_list_non_SB,D1_in_graph)
    G_neigh = nx.subgraph(Gint,neigh_list)

    node_colors = np.ones(len(G_neigh.nodes()))
    node_colors = pd.Series(node_colors,index = G_neigh.nodes())
    node_colors[D1_in_graph] = .5
    node_colors[D2_in_graph] = 0
    
    return G_neigh,D1_in_graph,D2_in_graph

def get_partition_no_small(G_neigh):
     # find communities of neighbor subgraph


    partition = community.best_partition(G_neigh)
    print('modularity of partition = ' + str(community.modularity(dict(partition),G_neigh)))
    
    # get rid of small groups
    partition = pd.Series(partition)
    part_VC = partition.value_counts()

    small_groups = list(part_VC[part_VC<5].index)
    for i in list(partition.index):
        if partition[i] in small_groups:
            partition[i] = -1
            
    return partition


def find_paths_to_LCC(G_disease,G_neigh):
    # this function finds the largest connected component, and returns an edge_list containing the shortest paths from 
    # unconnected nodes to LCC
    
    # find largest connected components for both diseases
    largest_CC = max(nx.connected_component_subgraphs(G_disease),key=len)
    LCC = largest_CC.nodes()
    
    nodes_in_disease = G_disease.nodes()

    not_in_LCC = np.setdiff1d(nodes_in_disease,LCC)

    # loop over nodes not in LCC
    Plen_all = []
    edge_list = []
    for nsource in not_in_LCC:
        # initialize path length at a high value
        Plen = 1000
        shortest_path = []
        for ntarget in LCC:
            # check if path exists
            if nx.has_path(G_neigh,nsource,ntarget):
                shortest_path_temp = nx.shortest_path(G_neigh,nsource,ntarget)
                PL_temp = len(shortest_path_temp)
                if PL_temp<Plen:
                    Plen = PL_temp
                    shortest_path = shortest_path_temp
            Plen_all.append(Plen)
            if len(shortest_path)>1:
                edge_list_temp = [(shortest_path[i],shortest_path[i+1]) for i in range(len(shortest_path)-1)]
            else:
                edge_list_temp = shortest_path
            edge_list.extend(edge_list_temp)

    # get rid of duplicate edges
    el_series = pd.Series(edge_list)
    edge_list = list(pd.unique(el_series))
    return edge_list




def plot_network_2_diseases(G_neigh,pos,G_1,G_2=None,d1name='d1',d2name='d2',saveflag=False,
                            savefilename='temp.png'):
    # this function plots the network, with orange up triangles for disease 1, blue down triangles for disease 2,
    # and purple diamonds for genes common to both diseases.
    # background neighborhood genes are also plotted in transparent small gray circles.  
    # edges are drawn between nodes in G_1 and/or G_2 (but no background edges are drawn for visual purposes)

    fig = plt.figure(figsize=(20,20),facecolor='#333236')
    axes = fig.add_subplot(1, 1, 1, axisbg='#333236')
    
    # add functionality for only 1 disease
    if G_2==None:
        d_label=d1name
        G_2=G_1
    else:
        d_label=d1name + ' and ' +  d2name

    # create lists of nodes for plotting
    nodes_both = list(np.intersect1d(G_1.nodes(),G_2.nodes()))
    nodes_only_1 = list(np.setdiff1d(G_1.nodes(),nodes_both))
    nodes_only_2 = list(np.setdiff1d(G_2.nodes(),nodes_both))
    
    


    # uncomment the following three lines to plot the shortest path lengths from nodes which are disconnected from LCC
    #nx.draw_networkx_edges(G_neigh,edgelist=G_SCZ_BIP.edges(),alpha=.3,pos=pos)
    #nx.draw_networkx_edges(G_neigh,edgelist = edge_list_SCZ,alpha=.2,pos=pos,edge_color='firebrick',style='--')
    #nx.draw_networkx_edges(G_neigh,edgelist = edge_list_BIP,alpha=.2,pos=pos,edge_color='royalblue',style='--')
    nx.draw_networkx_edges(G_neigh,edgelist = G_1.edges(),alpha=.4,pos=pos,edge_color='salmon',width=2.0)
    nx.draw_networkx_edges(G_neigh,edgelist = G_2.edges(),alpha=.4,pos=pos,edge_color='lightskyblue',width=2.0)
    nx.draw_networkx_nodes(G_neigh, alpha=.4,node_size=30,pos=pos,node_color='whitesmoke',cmap='Set1',with_labels=False,
                          label='hot neighbors of '+d1name)
    nx.draw_networkx_nodes(G_neigh, nodelist=nodes_only_1,alpha=.7,
                           linewidths=1.5,node_shape='^',node_size=160,pos=pos,
                           node_color='salmon',with_labels=False,label=d1name)
    nx.draw_networkx_nodes(G_neigh, nodelist=nodes_only_2,alpha=.7,
                           linewidths=1.5,node_shape='v',node_size=160,pos=pos,
                           node_color='lightskyblue',with_labels=False,label=d2name)

    nx.draw_networkx_nodes(G_neigh, nodelist=nodes_both,alpha=.7,
                           linewidths=1.5,node_shape='d',node_size=160,pos=pos,
                           node_color='violet',with_labels=False,label=d_label)
    leg = plt.legend(fontsize=16)
    for text in leg.get_texts():
        plt.setp(text, color = 'w')
    plt.grid('off')

    if saveflag:
        plt.savefig(savefilename,dpi=300)
        
# function to plot nodes color-coded by community
def plot_nodes_by_community(G_neigh,G_1,G_2,G_1_2,partition,pos,d1name='d1',d2name='d2',saveflag=False,
                           savefilename='temp_community.png'):
    fig = plt.figure(figsize=(22,22),facecolor='#E2E0E5')
    axes = fig.add_subplot(1, 1, 1, axisbg='#E2E0E5')
    partition = pd.Series(partition)

    # create lists of nodes for plotting
    nodes_both = list(np.intersect1d(G_1.nodes(),G_2.nodes()))
    nodes_only_1 = list(np.setdiff1d(G_1.nodes(),nodes_both))
    nodes_only_2 = list(np.setdiff1d(G_2.nodes(),nodes_both))

    nx.draw_networkx_edges(G_neigh,edgelist=G_1_2.edges(),alpha=.3,pos=pos)


    nx.draw_networkx_nodes(G_neigh,nodelist=G_neigh.nodes(),alpha=.6,node_shape='o',node_size=30,pos=pos,node_color=partition[G_neigh.nodes()],
                           vmin = min(partition),vmax=max(partition),cmap='Set1',with_labels=False,
                           label='neighbors of ' + d1name + ' and ' + d2name)

    nx.draw_networkx_nodes(G_neigh, nodelist=list(nodes_only_1),alpha=.9,node_shape='v',node_size=120,pos=pos,node_color=partition[nodes_only_1],
                           vmin = min(partition),vmax=max(partition),cmap='Set1',with_labels=False,
                           linewidths=.6,label=d1name)
    nx.draw_networkx_nodes(G_neigh, nodelist=list(nodes_only_2),alpha=.9,node_shape='^',node_size=120,pos=pos,node_color=partition[nodes_only_2],
                           vmin = min(partition),vmax=max(partition),cmap='Set1',with_labels=False,
                           linewidths=.6,label=d2name)
    nx.draw_networkx_nodes(G_neigh, nodelist=list(nodes_both),alpha=.9,node_shape='d',node_size=120,pos=pos,node_color=partition[nodes_both],
                           vmin = min(partition),vmax=max(partition),cmap='Set1',with_labels=False,
                           linewidths=.6,label=d1name + ' and ' + d2name)

    plt.grid('off')
    if saveflag:
        plt.savefig(savefilename,dpi=300)
        

def plot_focal_community(G_neigh,G_focal,G_1,G_2,G_1_2,partition,pos,d1name='d1',d2name='d2',
                         saveflag=False,savefilename='temp_comm.png'):
    # plot the focal group, with all edges drawn between nodes in this group
    fig = plt.figure(figsize=(20,20),facecolor='#E2E0E5')
    axes = fig.add_subplot(1, 1, 1, axisbg='#E2E0E5')
    partition = pd.Series(partition)

    focal_group = G_focal.nodes()
    
    nodes_not_focal = np.setdiff1d(G_neigh.nodes(),G_1_2.nodes())
    g1_focal = list(np.intersect1d(G_1.nodes(),focal_group))
    g2_focal = list(np.intersect1d(G_2.nodes(),focal_group))
    nodes_both_focal = list(np.intersect1d(g1_focal,g2_focal))
    nodes_only_1_focal = list(np.setdiff1d(g1_focal,nodes_both_focal))
    nodes_only_2_focal = list(np.setdiff1d(g2_focal,nodes_both_focal))

    nx.draw_networkx_edges(G_neigh,edgelist=G_focal.edges(),alpha=.1,pos=pos)
    nx.draw_networkx_nodes(G_neigh, node_list = nodes_not_focal,alpha=.3,node_size=30,pos=pos,
                           node_color='gray',cmap='Set1',with_labels=False)

    nx.draw_networkx_nodes(G_neigh, nodelist=list(focal_group),alpha=.6,node_shape='o',node_size=30,pos=pos,node_color=partition[focal_group],
                           vmin = min(partition),vmax=max(partition),cmap='Set1',with_labels=False,
                           linewidths=1.1,label='neighbors of ' + d1name + ' and ' + d2name)
    nx.draw_networkx_nodes(G_neigh, nodelist=list(nodes_only_1_focal),alpha=.9,node_shape='^',node_size=70,pos=pos,node_color=partition[nodes_only_1_focal],
                           vmin = min(partition),vmax=max(partition),cmap='Set1',with_labels=False,
                           linewidths=1.3,label=d1name)
    nx.draw_networkx_nodes(G_neigh, nodelist=list(nodes_only_2_focal),alpha=.9,node_shape='v',node_size=70,pos=pos,node_color=partition[nodes_only_2_focal],
                           vmin = min(partition),vmax=max(partition),cmap='Set1',with_labels=False,
                           linewidths=1.3,label=d2name)
    nx.draw_networkx_nodes(G_neigh, nodelist=list(nodes_both_focal),alpha=.9,node_shape='d',node_size=70,pos=pos,node_color=partition[nodes_both_focal],
                           vmin = min(partition),vmax=max(partition),cmap='Set1',with_labels=False,
                           linewidths=1.3,label=d1name + ' and ' + d2name)

    leg = plt.legend(fontsize=16)
    plt.grid('off')
    if saveflag:
        plt.savefig(savefilename,dpi=300)
        
        






    
    