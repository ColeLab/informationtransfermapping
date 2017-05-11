# Demo code for statistical analysis of computational model (simulated data) Ito et al., 2017
# Takuya Ito (takuya.ito@rutgers.edu)
# Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017). Cognitive task information is transferred between brain regions via resting-state netw    ork topology. bioRxiv. https://doi.org/10.1101/101782

## Import modules
import numpy as np
import os
import scipy.stats as stats
import statsmodels.api as sm
from scipy.stats import gamma
import sys 
import multregressionconnectivity as mreg
import sklearn.svm as svm

def setUpActFlowMat(subj,net,fromnet,tasks,nodespercommunity,nblocks=20,datadir='./'):
    """
    Retrieves a predicted target network's activation pattern for a given source network.
    Organizes prediction data as sample (block) x feature (regions) matrix

    Parameters:
        subj = subjID (numeric)
        net = the target network (0 through 5; numeric)
        fromnet = the source network (0 through 5; numeric)
        tasks = taskIDs (array). tasks 1-4 are for the top down only task. tasks 5-9 are for the top down AND bottom up task stimulation.
        nodespercommunity = number of nodes per community; default 50
        nblocks = number of blocks per task condition. Default = 20 per condition, 80 in total.
        datadir = directory that holds the simulation data

    Returns
        svmmat = a sample (block) X feature (regional activations) matrix
        labels = task condition labels associated for each block
    """
    nsamples = len(tasks)*nblocks
    nfeatures = nodespercommunity # regions per network
    svm_mat = np.zeros((nsamples,nfeatures))
    labels = np.zeros((nsamples,))

    indir = datadir + '/actflow_predictions/'
    indcount = 0
    for task in tasks:
        filename = 'subj'+str(subj)+'_task'+str(task)+'_net'+str(fromnet)+'tonet'+str(net)+'_multregFC.txt'
            
        actflowdat = np.loadtxt(indir+filename,delimiter=',')
        svm_mat[indcount:(indcount+nblocks),:] = actflowdat.T
        labels[indcount:(indcount+nblocks)] = task
        
        indcount += nblocks
        
    return svm_mat, labels

def setUpBetasMat(subj,net,tasks,Ci,nodespercommunity,nblocks=20,datadir='./'):

    """
    Retrieves a network's *actual* activation pattern.
    Organizes data as sample (block) x feature (regions) matrix

    Parameters:
        subj = subjID (numeric)
        net = the network's activation pattern to retreive (0 through 5; numeric)
        tasks = taskIDs (array). tasks 1-4 are for the top down only task. tasks 5-9 are for the top down AND bottom up task stimulation.
        nodespercommunity = number of nodes per community; default 50
        nblocks = number of blocks per task condition. Default = 20 per condition, 80 in total.
        datadir = directory that holds the simulation data

    Returns
        svmmat = a sample (block) X feature (regional activations) matrix
        labels = task condition labels associated for each block
    """
    nfeatures = nodespercommunity # Number of regions for each network
    nsamples = len(tasks)*nblocks
    svm_mat = np.zeros((nsamples,nfeatures))
    labels =np.zeros((nsamples,))
    net_ind = np.where(Ci==net)[0]
    #net_ind = np.arange(net*nodespercommunity,net*nodespercommunity+nodespercommunity)

    indir = datadir + '/task_betas/'
    indcount = 0
    for task in tasks:
        filename = 'subj'+str(subj)+'_task'+str(task)+'_allblocks.txt'
        betas = np.loadtxt(indir + filename, delimiter=',')
        # Get relevant network data
        svm_mat[indcount:(indcount+nblocks),:] = betas[net_ind,:].T # get all trials
        labels[indcount:(indcount+nblocks)] = task
        
        indcount += nblocks
        
    # Demean svm_mat across features for each sample
    mean_tmp =  np.mean(svm_mat,axis=1) 
    mean_tmp.shape = (len(mean_tmp),1)
    svm_mat = svm_mat - mean_tmp
       
    return svm_mat, labels

def predictedToActualRSA((subj,net,fromnet,tasks,nblocks,Ci,nodespercommunity,datadir)):
    """
    Runs a leave-block-out CV style RSA analysis (leaving 4 blocks out per CV)
    This code performs the 2nd and 3rd steps in the 'information transfer mapping' procedure
    Trains on predicted ActFlow data
    Tests on real data (betas)

    Parameters:
        subj = subjID (numeric)
        net = the target network (0 through 5; numeric)
        fromnet = the source network (0 through 5; numeric)
        tasks = taskIDs (array). tasks 1-4 are for the top down only task. tasks 5-9 are for the top down AND bottom up task stimulation.
        nblocks = number of blocks per task condition. Default = 20 per condition, 80 in total.
        Ci = community affiliation vector
        nodespercommunity = number of nodes per community; default 50
        datadir = directory that holds the simulation data

    Returns
        info_transfer_est = the average information transfer estimate 
    """

    actflow_mat, labels = setUpActFlowMat(subj,net,fromnet,tasks,nodespercommunity,nblocks=nblocks,datadir=datadir)
    real_mat, labels = setUpBetasMat(subj,net,tasks,Ci,nodespercommunity,nblocks=nblocks,datadir=datadir)
    
    ncvs = nblocks
    indices = np.arange(actflow_mat.shape[0])
    matched_rhos = []
    mismatch_rhos = []
    for cv in range(ncvs):
        task_ind = {}
        prototype = {}
        # Construct prototypes of each task
        for task in tasks: 
            # Get indices for this particular task
            task_ind[task] = np.where(labels==task)[0]
            # Decide which one is your 'comparison test trial' will be
            test_ind = task_ind[task][cv]
            # Find the indices for the prototypes
            train_ind = np.setxor1d(test_ind,task_ind[task])
            prototype[task] = np.mean(real_mat[train_ind,:],axis=0)
            
        # Now compare each pair of tasks with the prototype
        for task_a in tasks:
            for task_b in tasks:
                test_ind = task_ind[task_a][cv] # Compare task a
                rho_tmp = stats.spearmanr(prototype[task_b].T,actflow_mat[test_ind,:].T)[0] # With task b
                if task_a==task_b:
                    # Match!
                    matched_rhos.append(rho_tmp)
                else:
                    mismatch_rhos.append(rho_tmp)

    
    # Get averages 
    matched_rhos_avg = np.arctanh(np.mean(matched_rhos))
    mismatch_rhos_avg = np.arctanh(np.mean(mismatch_rhos))
    info_transfer_est = matched_rhos_avg - mismatch_rhos_avg
    return info_transfer_est
        

def predictedToActualSVM((subj,net,fromnet,tasks,nblocks,Ci,nodespercommunity,datadir)):
    """
    Runs a leave-block-out CV style SVM analysis (leaving 4 blocks out per CV)
    This analysis corresponds to the SVM analysis in Supplementary Fig. 3
    Trains on predicted ActFlow data
    Tests on real data (betas)

    Parameters:
        subj = subjID (numeric)
        net = the target network (0 through 5; numeric)
        fromnet = the source network (0 through 5; numeric)
        tasks = taskIDs (array). tasks 1-4 are for the top down only task. tasks 5-9 are for the top down AND bottom up task stimulation.
        nblocks = number of blocks per task condition. Default = 20 per condition, 80 in total.
        Ci = community affiliation vector
        nodespercommunity = number of nodes per community; default 50
        datadir = directory that holds the simulation data

    Returns
        The mean decoding accuracy on the target network's actual activation pattern
    """
    
    actflow_mat, labels = setUpActFlowMat(subj,net,fromnet,tasks,nodespercommunity,nblocks=nblocks,datadir=datadir)
    real_mat, labels = setUpBetasMat(subj,net,tasks,Ci,nodespercommunity,nblocks=nblocks,datadir=datadir)
    
    # Normalize
    #actflow_mat = preprocessing.scale(actflow_mat, axis=0)
    #real_mat = preprocessing.scale(real_mat, axis=0)
    
    ncvs = nblocks
    indices = np.arange(actflow_mat.shape[0])
    clf = svm.SVC(C=1.0,kernel='linear')
    accuracy = []
    # Starting test indices
    test_ind = np.arange(0,80,20) 
    for cv in range(ncvs):
        train_ind = np.setxor1d(test_ind,indices)
        
        trainset = actflow_mat[train_ind,:]
        testset = real_mat[test_ind,:]
        trainset_labels = labels[train_ind]
        testset_labels = labels[test_ind]
        
        # Normalize the train set across samples (for each feature)
        trainset_mean = np.mean(trainset,axis=0)
        trainset_mean.shape = (1,len(trainset_mean))
        trainset_std = np.std(trainset, axis=0)
        trainset_std.shape = (1,len(trainset_std))
        trainset = np.divide((trainset - trainset_mean),trainset_std)
       # Normalize test set using the mean and std from the trainset (to avoid non-circularity)
        testset = np.divide((testset - trainset_mean),trainset_std)
        
        clf.fit(trainset, trainset_labels)
        accuracy.append(clf.score(testset, testset_labels))
        
        # Change new testing group
        test_ind += 1
        
    return np.mean(accuracy)
        
