# fMRI data processing file
# Takes in all of the rtfMRI data, processes it all into a dictionary.
# Before you run this, put all the actual data into a constituent directory
# - (I called mine "data").
# Also contains functions for doing things we might want to do with whatever.

import numpy as np
# import matplotlib as mpl
# import pandas as pd
# import scipy as sp
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
            # THIS GRAPH IS A WEIGHTED, FULLY CONNECTED GRAPH OF ALL 128 REGIONS
    
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
    
# Splits G_all into a dict of the same structure holding sz and ct subjects.
def separate_diagnoses(rtfmridirpath, G_all):
    
    sz_patients = [0, 1, 2, 3, 4, 5, 8, 9, 12, 13, 14]
    ct_patients = [6, 7, 10, 11]
    # Didn't feel like actually doing data processing on the excel file.
    # Excel is one indexed, these are zero indexd.
    
    G_ct = {}
    for i in G_all:
        G_ct[i] = {}
        for j in ct_patients:
            G_ct[i][j] = G_all[i][j]
    G_sz = {}
    for i in G_all:
        G_sz[i] = {}
        for j in sz_patients:
            G_sz[i][j] = G_all[i][j]
    G_split = {'Neurotypical Graphs': G_ct, 'Schizophrenia Graphs': G_sz}
    return G_split

