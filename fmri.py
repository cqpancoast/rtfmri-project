# fMRI data processing file
# Takes in all of the rtfMRI data, processes it all into a dictionary.
# Before you run this, put all the actual data into a constituent directory
# - (I called mine "data").
# Also contains functions for doing things we might want to do with whatever.

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

def load():
    G_all, keyed_atlas = initialize_from_folder(rtfmridirpath, datadir)
    G_split = separate_diagnoses(rtfmridirpath, G_all)
    return G_all, G_split, keyed_atlas

#G_all, G_split, keyed_atlas = load()

# Import the necessary files, put them all into a massive dictionary,
# Concurrently create a dictionary of graphs with the same structure.
# Also makes a keyed atlas: that is, {0:this brain region, 1: that br...,}

# IMPORTANT NOTE: nx.convert_matrix.from_numpy_array is the thing that we are
# considering changing if we want to make a different kind of graph for each 
# trial.

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
            # CREATES WEIGHTED GRAPH, INCLUDING NEGATIVE WEIGHTS
    
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

# process_diagnoses : G_split, ProcessingFxn(Graph -> ???) -> G_split_processed
# Process NT and SZ connectivity graphs using some metric that is passed in.
def process_diagnoses(G_split, process_graph):

    NT_processed = {}
    SZ_processed = {}

    # For each subject in a trial, process their graph in a certain way.
    for trial in G_split['NT Graphs']:
        NT_processed[trial] = {k: process_graph(v) for k, v in G_split['NT Graphs'][trial].items()}
    for trial in G_split['SZ Graphs']:
        SZ_processed[trial] = {k: process_graph(v) for k, v in G_split['SZ Graphs'][trial].items()}

    G_split_processed = {'NT Processed': NT_processed, 'SZ Processed': SZ_processed}
    return G_split_processed

### FUNCTIONS TO USE WITH THE ABOVE FUNCTION
#   nx.average_clustering
#   strength
#   attack_graph
#TODO function that computes DMN/CEN connectivity for a graph.
    # probably uses community detection
    # manual vs. finding something online?
    # manual:
        #
    

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

# I'm not sure how, but this should create a network in the shape of a brain.
# - IT SHOULD BE 3D
# - Not totally sure whether this function is necessary, but it probably is.
# - Really, I'm not even sure where to begin with this one.
def visualize_graph(G):
    nx.draw(G)