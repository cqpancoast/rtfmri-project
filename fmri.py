# fMRI data processing file
# Takes in all of the rtfMRI data, processes it all into a dictionary.
# Before you run this, put all the actual data into a constituent directory
# - (I called mine "data").
# Also contains functions for doing things we might want to do with whatever.

#TODO (THINGS TO POSSIBLY DO) ((mostly for binary networks))
# - delete zero nodes for every participant (possibly — figure this shit out)
# - run minimum spanning tree, build this up like leaves
#   - this is because of the possiblity of having disconnected nodes

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio
import networkx as nx
import os
import math


## INITIALIZATION

''' CASEY'S FXN INPUTS'''
rtfmridirpath = '/Users/caseypancoast/Documents/_School/Complex Networks/rtfmri/'
datadir = 'data/'

# Import the necessary files, put them all into a massive dictionary,
# Concurrently create a dictionary of graphs with the same structure.
# Also makes a keyed atlas: that is, {0:this brain region, 1: that br...,}

def initialize_from_folder(rtfmridirpath, datadir):
    datadirpath = rtfmridirpath + datadir
    triallist = os.listdir(datadirpath) #assumes that data are in constituent dir
    
    all_data = {}
    G_all = {}
    for i in triallist:
        filepath = datadirpath + i
        all_data[i] = {}
        G_all[i] = {}
        for j in range(spio.loadmat(filepath)['Z'].shape[2]): # num of subjects
            #Do we need to turn all nan's into zero?
            all_data[i][j] = (spio.loadmat(filepath)['Z'])[:,:,j]
            # [i] is the trial, [j] is the subject number.
            G_all[i][j] = nx.convert_matrix.from_numpy_array(all_data[i][j])
            # CREATES WEIGHTED GRAPH, INCLUDING NEGATIVE WEIGHTS AND NAN SELF-LOOPS
    
    # Convert atlas into dictionary keyed by node #
    atlasfile = 'atlas.txt'
    atlaspath = rtfmridirpath + atlasfile
    file = open(atlaspath, mode='r')
    atlas = np.asarray(file.read().split(sep='\n'))
    keyed_atlas = {}
    for i in range(atlas.shape[0]):
        keyed_atlas[i] = atlas[i]
        
    return G_all, keyed_atlas 
    
# Splits G_all into two dicts of the same structure holding sz and ct subjects.
def separate_diagnoses(rtfmridirpath, G_all):
    
    sz_patients = [0, 1, 2, 3, 4, 5, 8, 9, 12, 13, 14]
    nt_patients = [6, 7, 10, 11]
    # Didn't feel like actually doing data processing on the excel file.
    # Excel is one indexed, these are zero indexed.
    
    G_nt = {}
    for i in G_all:
        G_nt[i] = {}
        for j in nt_patients:
            G_nt[i][j] = G_all[i][j]
    G_sz = {}
    for i in G_all:
        G_sz[i] = {}
        for j in sz_patients:
            G_sz[i][j] = G_all[i][j]
    G_split = {'NT Graphs': G_nt, 'SZ Graphs': G_sz}
    return G_split

def load():
    G_all, keyed_atlas = initialize_from_folder(rtfmridirpath, datadir)
    G_split = separate_diagnoses(rtfmridirpath, G_all)
    return G_all, G_split, keyed_atlas

#G_all, G_split, keyed_atlas = load()
    
# DMN and CEN nodes
dmn_nodes = [26, 27, 28, 29, 30, 31, 48, 55, 56, 79, 80, 91, 92, 114, 115, 116, 117, 118, 119, 120, 121]
cen_nodes = [0, 1, 4, 5, 34, 35]


## PROCESSING

# process_diagnoses : G_split, ProcessingFxn(Graph -> X) -> G_split_processed(X)
# Process NT and SZ connectivity graphs using some metric that is passed in.
def process_diagnoses(G_split, process_graph):
    
    NT_processed = {}
    SZ_processed = {}

    # For each subject in a trial, process their graph in a certain way.
    for trial in G_split['NT Graphs']:
        NT_processed[trial] = {k: process_graph(v) for k, v in G_split['NT Graphs'][trial].items()}
    for trial in G_split['SZ Graphs']:
        SZ_processed[trial] = {k: process_graph(v) for k, v in G_split['SZ Graphs'][trial].items()}

    G_split_processed = {'NT_processed': NT_processed, 'SZ_processed': SZ_processed}
    return G_split_processed

### FUNCTIONS TO USE WITH THE ABOVE FUNCTION
#   filter_edges_by_weight   
    # this allows you to work with functions on non-weighted graphs
#   strength
#   attack_graph
#   performance/coverage (lambda'd with DMN/CEN nodes)
        
# average_diagnoses : G_split_processed(X), AveragingFxn(List-of X -> X) 
#                         -> G_diagnosis_comparison(X)
# Takes in a processed G_split, averages the NT diagnoses and the SZ diagnoses
# according to a provided function.
# Used to compare NT/SZ graphs on network measures.
def average_diagnoses(G_split_processed):
    
    # list_average : [List-of Float] -> Float
    def list_average(l):
        if len(l) == 0:
            return 0
        else:
            return sum(l)/len(l)
    
    # dict_average : [list-of Dict(Int, Float)] -> Dict(Int, Float)
    # Constructs a single dict whose keys are all keys from the original dicts
    # and whose values are the averages of all values at those keys.
    def dict_average(dict_list):
        
        keys = range(max(map(lambda d: len(d), dict_list)))
        
        averaged_dict = {}
        for i in keys:
            values_at_index = []
            for dict_ in dict_list:
                if i in dict_.keys():
                    values_at_index.append(dict_[i])
                else: continue
            averaged_dict[i] = list_average(values_at_index)
        
        return averaged_dict
    
    # Invoke the correct averaging function depending on the processing type
    processed_item = list(list(list(G_split_processed.values())[0].values())[0].values())[0]
    if isinstance(processed_item, dict):
        average_processing = dict_average
    elif isinstance(processed_item, float or int):
        average_processing = lambda l: sum(l)/len(l)
    
    NT_averages = []
    SZ_averages = []
    for trial in G_split_processed['NT_processed']:
        NT_averages.append(average_processing(list(G_split_processed['NT_processed'][trial].values())))
    for trial in G_split_processed['SZ_processed']:
        SZ_averages.append(average_processing(list(G_split_processed['SZ_processed'][trial].values())))
    
    NT_average = average_processing(NT_averages)
    SZ_average = average_processing(SZ_averages)
    
    G_diagnosis_comparison = {'NT': NT_average, 'SZ': SZ_average}
    return G_diagnosis_comparison

# filter_edges_by_weight : Graph(Weighted), Float or Int -> Graph
# Removes all edges from a graph below a given weight then removes all weights.
# kwargs: weighted = Boolean (default is False)
#TODO filter by density instead of weight.
def filter_edges_by_weight(G_in, weight_limit, **kwargs):
    
    if kwargs is None:
        weighted = False
    else:
        weighted = kwargs['weighted']
    
    G = G_in.copy() #don't alter the original nodes!
    edge_list = list(G.edges()) #to prevent weird stuff
    
    for u, v in edge_list:
        if math.isnan(G[u][v]['weight']) or G[u][v]['weight'] <= weight_limit:
            G.remove_edge(u, v)    
        elif not weighted: 
            del G[u][v]['weight']
        
    return G

# strength : Graph -> Dict(Node, Strength)
# Plots the strength dict of a given graph. Strength values can be negative.
# *args:
# - default: no condition on which edges are considered
# - pos_only: only positive edges are considered
# - neg_only: only negative edges are considered
# - range: values between the next two inputted numbers will be considered
def strength(G, *args):
    
    if len(args) == 0:
        cond = lambda x: True
    elif args[0] == 'pos_only':
        cond = lambda x: x > 0
    elif args[0] == 'neg_only':
        cond = lambda x: x < 0
    elif args[0] == 'range':
        cond = lambda x: x > args[1] and x < args[2]   
    
    G_strength = {}
    for node in G.nodes():
        node_strength = 0
        for neighbor in G.neighbors(node):
            neighbor_strength = G[node][neighbor]['weight']
            if math.isnan(neighbor_strength):
                continue
            elif cond(neighbor_strength):
                node_strength += neighbor_strength
                
        G_strength[node] = node_strength
        
    return G_strength
            
# regional_connectivity : Graph [List-of Node] -> Float
# Takes in a list of nodes and uses a coverage function to measure the 
# effectiveness of that partition.
# (NOTE: Treating the rest of the network as a single community is okay in this
# instance, as we are only comparing the relative connectivities of the regions
# across trials, not making some sort of non-relative claim about them.)
def regional_connectivity(G, region):
    
    # partition G
    partition = [region, list(set(G.nodes())^set(region))]
    
    # All edges in the given region
    region_edges = []
    for u in region:
        for v in region:
            region_edges.append((u, v))
            
    print(region_edges)
    
    if nx.is_weighted(G):
        in_edge_weight = 0
        all_edge_weight = 0
        for u, v in G.edges():
            if u == v:
                continue
            elif (u, v) in region_edges:
                in_edge_weight += G[u][v]['weight']
                all_edge_weight += G[u][v]['weight']
            else:
                all_edge_weight += G[u][v]['weight']
        connectivity = in_edge_weight/all_edge_weight
    elif not nx.is_weighted(G):
        connectivity = nx.coverage(G, partition)
        
    if connectivity < 0:
        print("There's a problem.")
    
    return connectivity

# attack_graph : Graph NodePropertyFxn(Graph -> Dict) -> Dict
# Attack the graph by removing the nodes with the highest cof a given property
# Returns a dictionary where the keys are the fraction of nodes that have been
# removed and the values are largest connected component sizes.
def attack_graph(G, node_property):
    
    def find_max_value_index(dict_):
        
        max_value = 0
        n_to_rm = 0
        for index in dict_:
            if dict_[index] > max_value:
                max_value = dict_[index]
                n_to_rm = index
                
        # If the nodes are all the same size, just choose the first one.
        if n_to_rm == 0:
            n_to_rm = list(dict_.keys())[0]
    
        return n_to_rm
    
    recalc_t = 3 #How many nodes do you go through before recalculating?
    
    G_init_size = G.number_of_nodes()
    largest_CC = {}
    
    for i in range(math.floor(G_init_size/recalc_t) - 1):
            
        for j in range(recalc_t):
            G.remove_node(find_max_value_index(node_property(G)))
        
        # Fraction of nodes removed
        frac_rm = i*recalc_t/G_init_size
        # Largest connected component
        largest_CC[frac_rm] = max(list(map(lambda G: G.number_of_nodes(), nx.connected_component_subgraphs(G))))

    return largest_CC


## VISUALIZATION

# plot_dicts : Dict(Number, Number) or [List/Dict of Dict(Number, Number)] -> Image
# Takes in a Dictionary and plots the values against the sorted keys.
# Cited: https://stackoverflow.com/questions/37266341/plotting-a-python-dict-in-order-of-key-values
def plot_dicts(dicts_):
    #Just a single dict?
    if isinstance(dicts_, dict) and isinstance(list(dicts_.values())[0], float or int):
        plt.plot(*zip(*sorted(dicts_.items())))
        plt.show()
    #List of dicts?
    elif isinstance(dicts_, list) and isinstance(dicts_[0], dict):
        for dict_ in dicts_:
            plt.plot(*zip(*sorted(dict_.items())))
        plt.show()
    #Dict of dicts?
    elif isinstance(dicts_, dict) and isinstance(list(dicts_.values())[0], dict):
        for dict_ in list(dicts_.values()):
            plt.plot(*zip(*sorted(dict_.items())))
        plt.show()
    else:
        print("plot_dicts: given unrecognized type")

# I'm not sure how, but this should create a network in the shape of a brain.
# - SHOULD IT BE 3D
# - Not totally sure whether this function is necessary, but it probably is.
# - Really, I'm not even sure where to begin with this one.
def visualize_graph(G):
    nx.draw(G)