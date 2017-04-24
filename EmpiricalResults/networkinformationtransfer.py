# Demo code for network-to-network information transfer mapping (empirical data),  Ito et al., 2017
# Author: Takuya Ito (takuya.ito@rutgers.edu)
# Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017). Cognitive task information is transferred between brain regions via resting-state netw    ork topology. bioRxiv. https://doi.org/10.1101/101782

## Import modules
import sys
# append path to utils directory
sys.path.append('../utils/')
# append path to information transfer mapping module
sys.path.append('../')
import numpy as np
import scipy.stats as stats
import multregressionconnectivity as mreg
import informationtransfermapping as itm
import os

# Import global variable network partitions
networkdef = np.loadtxt('./data/network_array.csv',delimiter=',')
networkmappings = {'fpn':7, 'vis':1, 'smn':2, 'con':3, 'dmn':6, 'aud1':8, 'aud2':9, 'dan':11}
# Force aud2 key to be the same as aud1 (merging two auditory networks)
aud2_ind = np.where(networkdef==networkmappings['aud2'])[0]
networkdef[aud2_ind] = networkmappings['aud1']
# Redefine new network mappings with no aud1/aud2 distinction
networkmappings = {'fpn':7, 'vis':1, 'smn':2, 'con':3, 'dmn':6, 'aud':8, 'dan':11}

def computeRestFC(subj):
    """
    Computes the Glasser parcellated restFC of a given subject using multiple linear regression
    Takes in a subject number as a string, and outputs a FC matrix
    """
    filename = './data/FC_Estimates/' + subj + '_rest_nuisanceResids_Glasser.csv'
    timeseries = np.loadtxt(filename, delimiter=',')
    fcmat = mreg.multregressionconnectivity(timeseries)
    return fcmat

def loadActivations(subj, net='all'):
    """
    Loads in task betas
    
    Parameters: 
        subj = subject ID as a string
        net = network affiliation as a 3 letter string (see above in key values for networkmappings dictionary)
    Returns: 
        miniblock activations (128) X regions matrix
    """

    datafile = './data/MiniblockActivations/' + subj + '_miniblock_taskBetas_Glasser.csv'
    betas = np.loadtxt(datafile, delimiter=',')
    # The first 17 parameters nuisance regressors
    betas = betas[:,17:]
    # Isolate parcels belonging to specific networks
    if net == 'all':
        return betas.T
    else:
        net_ind = np.where(networkdef==net)[0]
        return betas[net_ind,:].T
    
def networkToNetworkInformationTransferMapping(subj,fcmat,null=False):
    """
    Given a subject number and a resting-state FC Matrix, compute the information transfer mapping estimates for all pairs of networks

    Parameters:
        subj = subject ID as a string
        fcmat = full 360x360 region-by-region matrix (using Glasser atlas)
        null = if null is true, run the null model (i.e., permuting FC topology)

    Returns:
        ite_matrix = a dictionary (for all 3 rule dimensions); each value contains a 7x7 matrix for all network-to-network information transfer estimates
    """
    # Use this network key ordering for computing matrix
    netkeys = {0:'vis',1:'smn',2:'con',3:'dmn',4:'fpn', 5:'aud', 6:'dan'}

    # Load in activation pattern
    betas  = loadActivations(subj, net='all')

    # Load in rule dimension identifiers (i.e., rule labels for each miniblock
    # Also instantiate empty information transfer estimate matrix
    ruledims = ['logic','sensory','motor']
    rule_labels = {}
    ite_matrix = {}
    for ruledim in ruledims:
        # Subtract miniblock indices by 1, since they were originally created for Matlab array indices
        rule_labels[ruledim] = np.loadtxt('./data/CPROTaskIdentifiers/' + subj + '_' + ruledim + '_miniblocksPerRule_MatlabIndices.csv',delimiter=',') - 1 
        # Instantiate empty network-to-network information transfer mapping matrix
        ite_matrix[ruledim] = np.zeros((len(networkmappings),len(networkmappings)))

    for net1 in netkeys.keys():
        for net2 in netkeys.keys():
            if net1==net2: continue
            # Find parcels for net1 and reshape
            net1_ind = np.where(networkdef==networkmappings[netkeys[net1]])[0]
            net1_ind.shape = (len(net1_ind),1)
            # Find parcels for net2 and reshape
            net2_ind = np.where(networkdef==networkmappings[netkeys[net2]])[0]
            net2_ind.shape = (len(net2_ind),1)
            
            # Compute betas for net1 and net2
            net1actual = np.squeeze(betas[:,net1_ind.T]) # Reformat to miniblocks x parcels
            net2actual = np.squeeze(betas[:,net2_ind.T]) # Reformat to miniblocks x parcels
            
            # Compute FC for net1 to net2
            w12 = fcmat[net1_ind,net2_ind.T].copy()
            w21 = fcmat[net2_ind,net1_ind.T].copy()
            
            # Randomize topology if null keyword is True
            if null:
                 np.random.shuffle(w12)
                 np.random.shuffle(w21)
            
            # Now perform network-to-network actflow
            net1tonet2predictions = itm.activityFlowXtoY(net1actual,w12)
            net2tonet1predictions = itm.activityFlowXtoY(net2actual,w21)

            # Do predicted-to-actual similarity analysis for all rule dimensions
            for ruledim in ruledims:
                ite_matrix[ruledim][net1, net2] = itm.predictedToActualSimilarity(net1tonet2predictions,net2actual,rule_labels[ruledim])        
                ite_matrix[ruledim][net2, net1] = itm.predictedToActualSimilarity(net2tonet1predictions,net1actual,rule_labels[ruledim])        

    return ite_matrix




            

