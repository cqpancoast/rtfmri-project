# fMRI data processing file
# Takes in all of the rtfMRI data, processes it all into a dictionary.
# Before you run this, put all the actual data into a constituent directory
# - (I called mine "data").
# Also contains functions for doing things we might want to do with whatever.

import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd
# import scipy as sp
import collections
import scipy.io as spio
import networkx as nx
import os


## INITIALIZATION

''' CASEY'S FXN INPUTS'''
rtfmridirpath = '/Users/caseypancoast/Documents/_School/Complex Networks/rtfmri/'
datadir = 'data/'

# Import the necessary files, put them all into a massive dictionary,
# Concurrently create a dictionary of graphs with the same structure.
# Also makes a keyed atlas: that is, {0:[this brain region], 1: [that br...],]}

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
 
    
## FUNCTIONS    
    
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

# I'm not sure how, but this should create a network in the shape of a brain.
# - IT SHOULD BE 3D
# - Not totally sure whether this function is necessary, but it probably is.
def visualize_graph(G):
    ...
    
# Plots the strength distribution of a given graph. Strength values can be negative.        
# THIS DOES NOT WORK YEET
# https://networkx.github.io/documentation/stable/auto_examples/drawing/plot_degree_histogram.html
def strength_distribution(G):
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    # print "Degree sequence", degree_sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())
    
    fig, ax = plt.subplots()
    # I am trying to make this a scatter plot instead of a bar chart.
    plt.bar(deg, cnt, width=0.80, color='b')
    
    plt.title("Degree Distribution")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    ax.set_xticklabels(deg)
    
    '''
    # draw graph in inset
    plt.axes([0.4, 0.4, 0.5, 0.5])
    Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[0]
    pos = nx.spring_layout(G)
    plt.axis('off')
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos, alpha=0.4)
    '''
    
    # Show the plot. Disable this if necessary.
    plt.show()