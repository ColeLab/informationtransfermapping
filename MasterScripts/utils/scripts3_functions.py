# Taku Ito
# 2/16/16

# Set of Functions required for MC 18 analyses
# Some functions were taken from previous pipelines i.e., mc16, etc.

import numpy as np
import time
import multiprocessing as mp
import nibabel as nib
import scipy.stats as stats
import sys
import matplotlib.pyplot as plt
import statsmodels.sandbox.stats.multicomp as mc
import random
import os
import statsmodels.tools.tools as tools
import scipy.stats.mstats as mstats
import sklearn.svm as svm
import sklearn.feature_selection as feature_selection
import sklearn.preprocessing as preprocessing
from sklearn import cross_validation

# Load Gordon Parcels on Surface
GordonParcels_filename= '/projects/IndivRITL/data/Parcels/Parcels_LR.dlabel.nii'
GordonParcels = nib.load(GordonParcels_filename)
GordonParcels = GordonParcels.get_data()
gordon_array = np.squeeze(GordonParcels)

def importRuleTimingsV3(subj, ruledim, hrfconv=True):
    """ 
    Import the rule timings for entire miniblocks for IC connectivity 
    Input: 
        subject number (as a string)
        rule dimension: either 'logic', 'sensory', 'motor', as a string
    Output:
        rules, rules_miniblock
        rules: dict ranging from 1-4 (for each rule), and the associated binary rule array for each rule
        rules_miniblock: dict ranging from 1-4 (for each rule), and then a sub-dict for each miniblock of that rule containing a 1xtime array of the whole timeseries indicating when that miniblock is occurring
    """
    
    rules = {}
    if hrfconv:
        stimfiledir = '/projects2/ModalityControl2/data/stimfiles_cprorules_v3/'
        rulemat = np.loadtxt(stimfiledir + subj + '_cprorules_hrflagged_stims.csv', delimiter=',')
        if ruledim == 'logic':
            rules[1] = rulemat[:,4]
            rules[2] = rulemat[:,5]
            rules[3] = rulemat[:,6]
            rules[4] = rulemat[:,7]
        elif ruledim == 'sensory':
            rules[1] = rulemat[:,0]
            rules[2] = rulemat[:,1]
            rules[3] = rulemat[:,2]
            rules[4] = rulemat[:,3]
        elif ruledim == 'motor':
            rules[1] = rulemat[:,8]
            rules[2] = rulemat[:,9]
            rules[3] = rulemat[:,10]
            rules[4] = rulemat[:,11]
        
    else:
        stimfiledir = '/projects/ModalityControl/data/stimfiles_flexhubpattern/'
        if ruledim == 'logic': 
            rules[1] = subj + '_stimfile_CPRORules_ModalityControlv3_EV5_BOTH.1D'
            rules[2] = subj + '_stimfile_CPRORules_ModalityControlv3_EV6_NOTBOTH.1D'
            rules[3] = subj + '_stimfile_CPRORules_ModalityControlv3_EV7_EITHER.1D'
            rules[4] = subj + '_stimfile_CPRORules_ModalityControlv3_EV8_NEITHER.1D'
            
        elif ruledim == 'sensory':
            rules[1] = subj + '_stimfile_CPRORules_ModalityControlv3_EV1_CONSTANT.1D'
            rules[2] = subj + '_stimfile_CPRORules_ModalityControlv3_EV2_HIPITCH.1D'
            rules[3] = subj + '_stimfile_CPRORules_ModalityControlv3_EV3_VERTICAL.1D'
            rules[4] = subj + '_stimfile_CPRORules_ModalityControlv3_EV4_RED.1D'
        elif ruledim == 'motor':
            rules[1] = subj + '_stimfile_CPRORules_ModalityControlv3_EV9_LMID.1D'
            rules[2] = subj + '_stimfile_CPRORules_ModalityControlv3_EV10_LINDEX.1D'
            rules[3] = subj + '_stimfile_CPRORules_ModalityControlv3_EV11_RMID.1D'
            rules[4] = subj + '_stimfile_CPRORules_ModalityControlv3_EV12_RINDEX.1D'
            
        rules[1] = np.loadtxt(stimfiledir + rules[1])
        rules[2] = np.loadtxt(stimfiledir + rules[2])
        rules[3] = np.loadtxt(stimfiledir + rules[3])
        rules[4] = np.loadtxt(stimfiledir + rules[4])
    
    # Now separate out miniblocks per for each rule
    # iterate through each TR
    rules_miniblock = {}
    for rule in rules:
        rules_miniblock[rule] = {}
        
    tr = 0
    miniblockcount = 0 # count miniblocks
    while tr <= len(rules[rule]):
        # if next 5 TRs all are stimulus, then we know it's the start of a new miniblock
        # encodings are a block of 5 TRs
        for rule in rules:
            if np.sum(rules[rule][tr:(tr+5)])==5:
                # if wholemb option selected, include whole task stims for entire miniblock
    #                if wholemb:
    #                    # Make all trs in that miniblock 1
    #                    rules[rule][tr:(tr+36)] = 1 # currently this will not be in use since only measuring miniblock connections
    #                
    #                    # create a separate binary timeseries for this miniblock
    #                    rules_miniblock[rule][miniblockcount] = np.zeros(shape=(4648,))
    #                    # Return TR indices for this miniblock
    #                    if tr + 36 > 4648: # Set this case for the last miniblock of the run
    #                        rules_miniblock[rule][miniblockcount] = range(tr,4648)
    #                    else:
    #                        rules_miniblock[rule][miniblockcount] = range(tr,(tr+36))
    #                    miniblockcount += 1
    #                else:
                #isolate the block to only include task times within that window
                # 0 out any other TRs
                tmp = np.zeros(shape=(rules[rule].shape))
                tmp[tr:tr+36] = rules[rule][tr:(tr+36)]
                block_ind = np.where(tmp==1)[0]
                rules_miniblock[rule][miniblockcount] = list(block_ind)
                # Go to next miniblock
                miniblockcount += 1
                # If true, we want to just skip ahead to the next miniblock
                tr += 36
          #  else:

        # Check if current TR is the onset of a new block... otherwise keep checking next TR
        newblock = False
        for rule in rules:
            if np.sum(rules[rule][tr:(tr+5)])==5:
                newblock=True
        if newblock==False:
            tr += 1
    
    return rules, rules_miniblock

#####################################################################################
# ModalityControl18a

def getBetaMatrixForNetwork(subj, network, ruledim):
    """
    Given a set of indices and a subject, find the task regressors of those set of ROI indices
    """
    resultdir = '/projects/ModalityControl/data/results/miniblock_beta_series_glm_results/'

    if ruledim =='vis' or ruledim=='aud':
        taskbetas = np.loadtxt(resultdir + subj + '_miniblock_taskBetas_sensory_GordonSurface.csv', delimiter=',')
    else:
        taskbetas = np.loadtxt(resultdir + subj + '_miniblock_taskBetas_' + ruledim + '_GordonSurface.csv', delimiter=',')
    taskbetas = taskbetas[:,17:]
    # Remove the first 17 betas since those are nuisance regressors
    if network=='all':
        output = taskbetas
    else:
        net_ind = np.where(order['oldlabels']==network)[0]
        output = np.zeros(shape=(1,len(taskbetas[0,:])))
        for roi in net_ind:
            roi_indices = np.where(gordon_array==roi)[0]
            output = np.vstack((output,taskbetas[roi_indices,:]))
        output = np.delete(output,0,axis=0) # Delete the first row
        
    rules, rulesmb = importRuleTimings(subj, 'sensory', hrfconv=True)
    if ruledim == 'vis':
        visrule_start = len(rulesmb[1].keys()) + len(rulesmb[2].keys())
        output = output[:,visrule_start:]
    elif ruledim == 'aud':
        audrule_end = len(rulesmb[1].keys()) + len(rulesmb[2].keys())
        output = output[:,:audrule_end]    
    
    return output

def setUpSVMMatrix(subj, ruledim, network):
    """
    Given a subject, rule dimension, and a set of Gordon indices (presumably of a CCN) obtain an SVM matrix and
    associated labels
    """
    labels = []
    if ruledim == 'vis':
        tmp, tmpmb = importRuleTimings(subj, 'sensory', hrfconv=True)
        rules = {}
        rulesmb = {}
        rules[1] = tmp[3]
        rules[2] = tmp[4]
        rulesmb[1] = tmpmb[3]
        rulesmb[2] = tmpmb[4]
    elif ruledim == 'aud':
        rules, rulesmb = importRuleTimings(subj, 'sensory', hrfconv=True)
        del rules[3], rules[4], rulesmb[3], rulesmb[4] # Delete visual rules if Auditory
    else:
        rules, rulesmb = importRuleTimings(subj, ruledim, hrfconv=True)

    # Import task betas
    taskbetas = getBetaMatrixForNetwork(subj, network, ruledim)
    
    # Count number of mini blocks per rule type 
    num_mbs = {}
    mbcounter = 0
    total_mbs = 0
    for rule in rulesmb:
        # Compute number of good miniblocks for this particular rule
        num_mbs[rule] = len(rulesmb[rule].keys())
        total_mbs += len(rulesmb[rule].keys())
    svm_mat = np.zeros(shape=(total_mbs,taskbetas.shape[0]))
    
    for rule in rulesmb:
        # Create labels for each rule
        for i in range(num_mbs[rule]): labels.append(rule)
        # Extract the betas from beta regressor matrix
        if rule == 1:
            svm_mat[0:num_mbs[rule],:] = taskbetas[:,0:num_mbs[rule]].T # Transpose so it's mbs x ROIs
            mbcounter += num_mbs[rule]
        else:
            mbstart = mbcounter
            mbend = mbcounter + num_mbs[rule]
            svm_mat[mbstart:mbend, :] = taskbetas[:,mbstart:mbend].T
            mbcounter += num_mbs[rule]

    return svm_mat, labels

def bootstrapTrainset(svm_mat, labels, n=100):
    """ 
    Given an unblanced SVM matrix and associated set of labels, will balance a data set
    Take the minimum number of samples for a given rule, and each rule will have that number of samples  
    Balances the trainset by sampling WITH REPLACEMENT
    """

    # Find number of features
    num_features = svm_mat.shape[1]

#     # First find the minimum number of miniblocks for each rule type
#     # Will determine the balancing number based on the fewest miniblock of a particular rule
#     rules_ind = {}
#     for rule in np.unique(labels):
#         # Find the indices of a particular rule
#         rules_ind[rule] = np.where(np.asarray(labels)==rule)[0]
#         # If it's the first rule, set the number to the number of miniblocks for the first rule
#         if rule==1:
#             num = len(np.where(np.asarray(labels)==rule)[0])
#         else:
#             # Otherwise, if the new number of miniblocks is fewer than the fewest number of miniblocks, reset variable
#             newnum = len(np.where(np.asarray(labels)==rule)[0])
#             if newnum < num:
#                 num = newnum

    # For each rule, randomly sample the minimum (i.e., 'num') number of samples of each rule type for a balanced matrix
    svm_mat2 = np.zeros(shape=(n*len(np.unique(labels)), num_features)) 
    labels2 = []
    samplecount = 0
    for rule in np.unique(labels):
        rule_ind = np.where(labels==rule)[0]
        ind = np.random.choice(rule_ind, size=n, replace=True)
        # Create balanced svm_mat with `n` samples per rule
        svm_mat2[samplecount:(samplecount+n),:] = svm_mat[ind,:]
        # Create balanced label array
        for i in range(n): labels2.append(rule)
        samplecount += n
    
    return svm_mat2, labels2

def balanceSamples(svm_mat, labels):
    min_samples = len(labels)
    for cond in np.unique(labels):
        # If the current condition has the fewest number of samples per condition, set this to the min
        num_samples = len(np.where(np.asarray(labels)==cond)[0])
        if num_samples < min_samples:
            min_samples = num_samples
    # Instantiate empty new SVM matrices and labels
    svm_mat2 = np.zeros(shape=(min_samples*len(np.unique(labels)), svm_mat.shape[1]))
    labels2 = np.zeros(shape=(min_samples*len(np.unique(labels)),))
    
    # Actual balancing loop
    start = 0
    for cond in np.unique(labels):
        rule_ind = np.where(labels==cond)[0]
        ind = np.random.choice(rule_ind, size=min_samples, replace=False)
        svm_mat2[start:(start+min_samples),:] = svm_mat[ind,:]
        labels2[start:(start+min_samples)] = np.asarray(labels)[ind]
        start += min_samples
    
    return svm_mat2, labels2

def runCV((svm_mat, labels, clf, permutation)):
    """ 
    Leave 4 out cross validation, leaving one of the rules out per CV
    """
    # Normalize data
    svm_mat = preprocessing.scale(svm_mat,axis=0)
    
    # For each rule, obtain an array of its indices
    rule_ind = {}
    for rule in np.unique(labels):
        rule_ind[rule] = np.where(np.asarray(labels)==rule)[0]
        if rule == 1:
            num_cvs = len(rule_ind[1])    # Compute number of CVs
        else:
            if len(rule_ind[rule]) < num_cvs: num_cvs = len(rule_ind[rule])

    accuracy = []
    for cv in range(num_cvs):
        # Specify trainset indices
        # Indexing array for separating train and test sets
        notcv = np.ones(shape=(num_cvs,),dtype=bool)
        notcv[cv] = 0
        # Create lists of trainset and testsets
        trainset_ind = []
        testset_ind = []
        for rule in rule_ind:
            testset_ind.append(rule_ind[rule][cv])
            trainset_ind.extend(list(rule_ind[rule][notcv]))
        trainset = svm_mat[trainset_ind,:]
        trainset_labels = np.asarray(labels)[trainset_ind]
        
#         ### ADDED FEATURE SELECTION
#         f, p = feature_selection.f_classif(trainset, trainset_labels)
#         p_ind = np.where(p < 0.1)[0]
#         trainset = trainset[:,p_ind]
#         ### END
        
        # balance trainset or bootstrap trainset
#         trainset, trainset_labels = bootstrapTrainset(trainset, trainset_labels)
        trainset, trainset_labels = balanceSamples(trainset, trainset_labels)
        # Permute training labels if 'permutation' is True
        if permutation:
            trainset_labels = np.random.permutation(trainset_labels)

        testset = svm_mat[testset_ind,:]
        
#         ### FEATURE SELECTION add
#         testset = testset[:,p_ind] 
#         ### END
        
        testset_labels = np.asarray(labels)[testset_ind]

        clf.fit(trainset, trainset_labels)
        accuracy.append(clf.score(testset, testset_labels))
    
    return np.mean(accuracy)

#####################################################################################
# ModalityControl18b

def computeActFlow(v,timeseries, beta_mat):
    # Compute fc from vertex to vertices
    num_vertices = timeseries.shape[0]
    fcmat = np.zeros(shape=(num_vertices,1))
    for i in range(num_vertices):
        fcmat[i] = stats.pearsonr(v, timeseries[i,:])[0]
    
    actflow = beta_mat * fcmat
    actflow = np.mean(actflow)
    return actflow

def computeVerticesInNetworks(network):
    """
    For a give network, outputs all the vertices within that network as an array
    """
    # identify the ROIs in this network
    net_ind = np.where(order['oldlabels']==network)[0]
    vertex_ind = np.zeros(shape=(1,)) # number of rest TRs
    for roi in net_ind:
        roi_indices = np.where(gordon_array==roi)[0]
        vertex_ind = np.hstack((vertex_ind,roi_indices))
    vertex_ind = np.delete(vertex_ind,0,axis=0) # Delete the first row
    vertex_ind = vertex_ind.astype(int)
    return vertex_ind

def correlateBetweenNetworks(timeseries, net1str, net2str):
    """
    Parameters:
    Cross-correlate vertex-wise resting state FC from net1str to net2str
    timeseries : timeseries for a particular subject's resting state
    net1str : vertices from first network to cross-correlate with
    net2str : vertices from second network to cross-corrleate with
    
    Output:
    output : fcmatrix with m x n dimensions, where m is the number of vertices from net1str
             and n is the number of vertices from net2str
    """
    # Get vetices from each of the networks
    net1 = computeVerticesInNetworks(net1str)
    net2 = computeVerticesInNetworks(net2str)
    # Create empty fc cross correlation mat
    fcmat = np.zeros(shape=(len(net1), len(net2)))
    timestart = time.time()
    icount = 0
    for i in net1:
        if np.mod(icount,500)==0:
#             print 'Completed vertex', icount, 'out of', len(net1), 'of', net1str
            timeend = time.time()
#             print 'Elapsed time (s) after 500 vertices', (timeend - timestart)
            timestart = timeend
        jcount = 0
        for j in net2: 
            fcmat[icount,jcount] = stats.pearsonr(timeseries[i,:], timeseries[j,:])[0]
            jcount += 1
        icount += 1
        
    # If net1 and net2 are the same, 0 the diagonal
    if net1str==net2str: np.fill_diagonal(fcmat, 0)

    return fcmat

def actFlowBetweenNetworks(fcmat, net1str, net2str, ruledim, subj, posOnly=True, outdir=None):
    """
    Run ActFlow between two specified networks.
    fcmat : vertex-wise crosscorrelation between net1str and net2str, where net1 is rows and net2 is columns
    net1str : must correspond to the network vertices on the rows axis of fcmat
    net2str : must correspond to the network vertices on the cols axis of cmat
    ruledim : rule dimension to simulate actflow
    subj : subject number
    """
    # If Pos only, use positive correlations only
    if posOnly:
        tmp = fcmat > 0
        fcmat = np.multiply(fcmat, tmp)
    # Get network vertices
    net1 = computeVerticesInNetworks(net1str)
    net2 = computeVerticesInNetworks(net2str)
    # Get betas
    betas = getBetaMatrixForNetwork(subj, 'all', ruledim)
    num_mbs = betas.shape[1]
    # Instantiate actflow of net1
    actflownet1 = np.zeros(shape=(len(net1),num_mbs))
    icount = 0
    for i in net1:
        for mb in range(num_mbs): 
            actflownet1[icount,mb] = np.mean(betas[net2,mb] * fcmat[icount,:]) # activity flow from net2 to net1
        icount += 1
        
    actflownet2 = np.zeros(shape=(len(net2),num_mbs))
    jcount = 0
    for j in net2:
        for mb in range(num_mbs):
            actflownet2[jcount,mb] = np.mean(betas[net1,mb] * fcmat[:,jcount])
        jcount += 1
    
    if outdir==None:
        outdir = '/projects/ModalityControl/data/results/actflow_network2network/'
    
    np.savetxt(outdir + subj + '_' + net1str + '_actflow_from_' + net2str + '_' + ruledim + '.txt', actflownet1, delimiter=',')
    np.savetxt(outdir + subj + '_' + net2str + '_actflow_from_' + net1str + '_' + ruledim + '.txt', actflownet2, delimiter=',')
    return actflownet1, actflownet2

def getActFlowActivity(subj, net, fromnet, ruledim, basedir=None):
    if basedir==None:
        basedir = '/projects/ModalityControl/data/results/actflow_network2network/'
    
    if ruledim =='vis' or ruledim=='aud':
        actflow = np.loadtxt(basedir + subj + '_' + net + '_actflow_from_' + fromnet + '_sensory.txt', delimiter=',')
    else:
        actflow = np.loadtxt(basedir + subj + '_' + net + '_actflow_from_' + fromnet + '_' + ruledim + '.txt', delimiter=',')
        
    rules, rulesmb = importRuleTimings(subj, 'sensory', hrfconv=True)
    if ruledim == 'vis':
        visrule_start = len(rulesmb[1].keys()) + len(rulesmb[2].keys())
        actflow = actflow[:,visrule_start:]
    elif ruledim == 'aud':
        audrule_end = len(rulesmb[1].keys()) + len(rulesmb[2].keys())
        actflow = actflow[:,:audrule_end]    
   
    return actflow

# SVM matrix set up for all subjects
def setUpSVMMatrixActFlow(subj, ruledim, net, fromnet, basedir=None):
    """
    Given a subject, rule dimension, and a set of Gordon indices (presumably of a CCN) obtain an SVM matrix and
    associated labels
    """

    taskbetas = getActFlowActivity(subj, net, fromnet, ruledim, basedir=basedir)

    svm_mat = np.zeros(shape=taskbetas.T.shape)
    labels = []
    if ruledim == 'vis':
        tmp, tmpmb = importRuleTimings(subj, 'sensory', hrfconv=True)
        rules = {}
        rulesmb = {}
        rules[1] = tmp[3]
        rules[2] = tmp[4]
        rulesmb[1] = tmpmb[3]
        rulesmb[2] = tmpmb[4]
    elif ruledim == 'aud':
        rules, rulesmb = importRuleTimings(subj, 'sensory', hrfconv=True)
        del rules[3], rules[4], rulesmb[3], rulesmb[4] # Delete visual rules if Auditory
    else:
        rules, rulesmb = importRuleTimings(subj, ruledim, hrfconv=True)


    # Count number of mini blocks per rule type 
    num_mbs = {}
    mbcounter = 0
    total_mbs = 0
    for rule in rulesmb:
        # Compute number of good miniblocks for this particular rule
        num_mbs[rule] = len(rulesmb[rule].keys())
        total_mbs += len(rulesmb[rule].keys())
    svm_mat = np.zeros(shape=(total_mbs,taskbetas.shape[0]))
    
    for rule in rulesmb:
        # Compute number of good miniblocks for this particular rule
        num_mbs[rule] = len(rulesmb[rule].keys())
        # Create labels for each rule
        for i in range(num_mbs[rule]): labels.append(rule)
        # Extract the betas from beta regressor matrix
        if rule == 1:
            svm_mat[0:num_mbs[rule],:] = taskbetas[:,0:num_mbs[rule]].T # Transpose so it's mbs x ROIs
            mbcounter += num_mbs[rule]
        else:
            mbstart = mbcounter
            mbend = mbcounter + num_mbs[rule]
            svm_mat[mbstart:mbend, :] = taskbetas[:,mbstart:mbend].T
            mbcounter += num_mbs[rule]
    return svm_mat, labels


def actflowSVM((svm_train, svm_test, trainset_labels, testset_labels, clf, permutation)):
    """ 
    Leave 4 out cross validation, leaving one of the rules out per CV
    """
    # Normalize data together
    num_train = svm_train.shape[0]
#     scaled = np.vstack((svm_train, svm_test))
#     scaled = preprocessing.scale(scaled, axis=0)
#     svm_train = scaled[0:num_train,:]
#     svm_test = scaled[num_train:,:]

    svm_train = preprocessing.scale(svm_train, axis=0)
    svm_test = preprocessing.scale(svm_test, axis=0)
    
#         ### ADDED FEATURE SELECTION
#         f, p = feature_selection.f_classif(trainset, trainset_labels)
#         p_ind = np.where(p < 0.1)[0]
#         trainset = trainset[:,p_ind]
#         ### END
        
    
    trainset = svm_train
    testset = svm_test
    
    trainset, trainset_labels = balanceSamples(trainset, trainset_labels)
    testset, testset_labels = balanceSamples(testset, testset_labels)
    
    # Permute training labels if 'permutation' is True
    if permutation:
        trainset_labels = np.random.permutation(trainset_labels)

#         ### FEATURE SELECTION add
#         testset = testset[:,p_ind] 
#         ### END

    clf.fit(trainset, trainset_labels)
    accuracy = clf.score(testset, testset_labels)

    return accuracy


def actflowSVMleaveBlockOut((svm_train, svm_test, trainset_labels, testset_labels, clf, permutation)):
    """ 
    Leave 4 out cross validation, leaving one of the rules out per CV
    """
    # Normalize data together
    num_train = svm_train.shape[0]
#     scaled = np.vstack((svm_train, svm_test))
#     scaled = preprocessing.scale(scaled, axis=0)
#     svm_train = scaled[0:num_train,:]
#     svm_test = scaled[num_train:,:]

    svm_train = preprocessing.scale(svm_train, axis=0)
    svm_test = preprocessing.scale(svm_test, axis=0)
    
    # trainset labels and testset labels should be the same for these purposes
    labels = trainset_labels

    ## TEST CODE 
    # For each rule, obtain an array of its indices
    rule_ind = {}
    for rule in np.unique(labels):
        rule_ind[rule] = np.where(np.asarray(labels)==rule)[0]
        if rule == 1:
            num_cvs = len(rule_ind[1])    # Compute number of CVs
        else:
            if len(rule_ind[rule]) < num_cvs: num_cvs = len(rule_ind[rule])

    accuracy = []
    for cv in range(num_cvs):
        # Specify trainset indices
        # Indexing array for separating train and test sets
        notcv = np.ones(shape=(num_cvs,),dtype=bool)
        notcv[cv] = 0
        # Create lists of trainset and testsets
        trainset_ind = []
        testset_ind = []
        for rule in rule_ind:
            testset_ind.append(rule_ind[rule][cv])
            trainset_ind.extend(list(rule_ind[rule][notcv]))
        trainset = svm_train[trainset_ind,:]
        trainset_labels = np.asarray(labels)[trainset_ind]
        
#         ### ADDED FEATURE SELECTION
#         f, p = feature_selection.f_classif(trainset, trainset_labels)
#         p_ind = np.where(p < 0.1)[0]
#         trainset = trainset[:,p_ind]
#         ### END
        
        # balance trainset or bootstrap trainset
#         trainset, trainset_labels = bootstrapTrainset(trainset, trainset_labels)
        trainset, trainset_labels = balanceSamples(trainset, trainset_labels)
        # Permute training labels if 'permutation' is True
        if permutation:
            trainset_labels = np.random.permutation(trainset_labels)

        testset = svm_test[testset_ind,:]
        
#         ### FEATURE SELECTION add
#         testset = testset[:,p_ind] 
#         ### END
        
        testset_labels = np.asarray(labels)[testset_ind]

        clf.fit(trainset, trainset_labels)
        accuracy.append(clf.score(testset, testset_labels))
#         ### ADDED FEATURE SELECTION
#         f, p = feature_selection.f_classif(trainset, trainset_labels)
#         p_ind = np.where(p < 0.1)[0]
#         trainset = trainset[:,p_ind]
#         ### END
        
    ## TEST CODE END

    return np.mean(accuracy)


def actflowSVMleaveBlockOut_FeatSel((svm_train, svm_test, trainset_labels, testset_labels, clf, permutation, featsel)):
    """ 
    Leave 4 out cross validation, leaving one of the rules out per CV
    With univariate feature selection
    """
    # Normalize data together
    num_train = svm_train.shape[0]
#     scaled = np.vstack((svm_train, svm_test))
#     scaled = preprocessing.scale(scaled, axis=0)
#     svm_train = scaled[0:num_train,:]
#     svm_test = scaled[num_train:,:]

    svm_train = preprocessing.scale(svm_train, axis=0)
    svm_test = preprocessing.scale(svm_test, axis=0)
    
    # trainset labels and testset labels should be the same for these purposes
    labels = trainset_labels

    ## TEST CODE 
    # For each rule, obtain an array of its indices
    rule_ind = {}
    for rule in np.unique(labels):
        rule_ind[rule] = np.where(np.asarray(labels)==rule)[0]
        if rule == 1:
            num_cvs = len(rule_ind[1])    # Compute number of CVs
        else:
            if len(rule_ind[rule]) < num_cvs: num_cvs = len(rule_ind[rule])

    accuracy = []
    for cv in range(num_cvs):
        # Specify trainset indices
        # Indexing array for separating train and test sets
        notcv = np.ones(shape=(num_cvs,),dtype=bool)
        notcv[cv] = 0
        # Create lists of trainset and testsets
        trainset_ind = []
        testset_ind = []
        for rule in rule_ind:
            testset_ind.append(rule_ind[rule][cv])
            trainset_ind.extend(list(rule_ind[rule][notcv]))
        trainset = svm_train[trainset_ind,:]
        trainset_labels = np.asarray(labels)[trainset_ind]
        
        ### ADDED FEATURE SELECTION
        f, p = feature_selection.f_classif(trainset, trainset_labels)
        p_ind = np.where(p < featsel)[0]
        trainset = trainset[:,p_ind]
        ### END
        
        # balance trainset or bootstrap trainset
#         trainset, trainset_labels = bootstrapTrainset(trainset, trainset_labels)
        trainset, trainset_labels = balanceSamples(trainset, trainset_labels)
        # Permute training labels if 'permutation' is True
        if permutation:
            trainset_labels = np.random.permutation(trainset_labels)

        testset = svm_test[testset_ind,:]
        
        ### FEATURE SELECTION add
        testset = testset[:,p_ind] 
        ### END
        
        testset_labels = np.asarray(labels)[testset_ind]

        clf.fit(trainset, trainset_labels)
        accuracy.append(clf.score(testset, testset_labels))
        

    return np.mean(accuracy)

