# fMRI data processing file
# Takes in all of the rtfMRI data, processes it all into a dictionary.
# Before you run this, put all the actual data into a constituent directory
# - (I called mine "data").
# Also contains functions for doing things we might want to do with whatever.

import numpy as np
# import matplotlib as mpl
import pandas as pd
import scipy as sp
import networkx as nx
import os


## INITIALIZATION

''' CASEY'S FXN INPUTS
rtfmridirpath = '/Users/caseypancoast/Documents/MATLAB/rtfmri/'
datadir = 'data/'
'''

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
        for j in range(sp.io.loadmat(filepath)['Z'].shape[2]): # num of subjects
            #Do we need to turn all nan's into zero?
            all_data[i][j] = (sp.io.loadmat(filepath)['Z'])[:,:,j]
            # [i] is the trial, [j] is the subject number.
            G_all[i][j] = nx.convert_matrix.from_numpy_array(all_data[i][j])
    
    # Convert atlas into dictionary keyed by node #
    atlasfile = 'atlas.txt'
    atlaspath = rtfmridirpath + atlasfile
    file = open(atlaspath, mode='r')
    atlas = np.asarray(file.read().split(sep='\n'))
    keyed_atlas = {}
    for i in range(atlas.shape[0]):
        keyed_atlas[i] = atlas[i]
    
    
## FUNCTIONS    
    
# Produces two lists of subject numbers from excel diagnosis data.
# Splits the alldata graph G_all up into dictionaries of schiz/non-schiz subjects.
        
### AS OF RIGHT NOW ThIS FUNCTION DOES NOT WORK ###         
        
def separate_diagnoses(rtfmridirpath, G_all):
    
    excelfile = 'rtfmri_info.xlsx'
    excelpath = rtfmridirpath + excelfile
    diagdata = pd.read_excel(excelpath)['sz_patient']
    # The above is a dictionary wth 0 - 14 mapping to 0 for control and 1 for schz.
    
    G_ctr = {}
    for i in G_all:
        for j in diagdata:
            G_ctr[i][j] = G_all[i][j]
    G_sch = {}
    for i in G_all:
        for j in diagdata['Schizophrenia']:
            G_sch[i][j] = G_all[i][j]
    G_split = {'Neurotypical Graphs': G_ctr, 'Schizophrenia Graphs': G_sch}
    return G_split

