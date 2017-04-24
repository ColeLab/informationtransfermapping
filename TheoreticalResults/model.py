# Demo code for Ito et al., 2017
# Author: Takuya Ito (takuya.ito@rutgers.edu)
# Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017). Cognitive task information is transferred between brain regions via resting-state network topology. bioRxiv. https://doi.org/10.1101/101782

## Import modules
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.api as sm
from scipy.stats import gamma
import sys
sys.path.append('../utils/')
import multregressionconnectivity as mreg


# Define transfer function (hyperbolic tangent)
phi = lambda x: np.tanh(x)

def generateStructuralNetwork(ncommunities=5, innetwork_dsity=.35, outnetwork_dsity=.05, 
                              hubnetwork_dsity=.2, nodespercommunity=50, showplot=False):
    """
    Randomly generates a structural network with a single hub network

    Parameters:
        ncommunities = number of communities within the network (one will automatically be a hub-network
        innetwork_dsity = connectivity density of within-network connections
        outnetwork_dsity = connectivity density of out-of-network connections
        hubnetwork_dsity = out-of-network connectivity density for the hub-network
        showplot = if set to True, will automatically display the structural matrix using matplotlib.pyplot

    Returns: 
        Unweighted structural connectivity matrix (with 1s indicating edges and 0s otherwise)
    """

    totalnodes = nodespercommunity * ncommunities

    W = np.zeros((totalnodes,totalnodes))
    # Construct structural matrix
    nodecount = 0
    for i in range(ncommunities):
        for j in range(ncommunities):
            for node in range(nodespercommunity):
                # Set within network community connections
                if i==j:
                    tmp_a = np.random.rand(nodespercommunity,nodespercommunity)<innetwork_dsity
                    indstart = i*nodespercommunity
                    indend = i*nodespercommunity+nodespercommunity
                    W[indstart:indend,indstart:indend] = tmp_a
                else:
                    tmp_b = np.random.rand(nodespercommunity,nodespercommunity)<outnetwork_dsity
                    indstart_i = i*nodespercommunity
                    indend_i = i*nodespercommunity + nodespercommunity
                    indstart_j = j*nodespercommunity
                    indend_j = j*nodespercommunity + nodespercommunity
                    W[indstart_i:indend_i, indstart_j:indend_j] = tmp_b

    # Redo a community as a hub-network
    hubnetwork = 0
    if hubnetwork_dsity>0: # Only do if we want a hub network
        for i in range(ncommunities):
            for j in range(ncommunities):
                if (i==hubnetwork or j==hubnetwork) and i!=j:
                    tmp_b = np.random.rand(nodespercommunity,nodespercommunity)<hubnetwork_dsity
                    indstart_i = i*nodespercommunity
                    indend_i = i*nodespercommunity + nodespercommunity
                    indstart_j = j*nodespercommunity
                    indend_j = j*nodespercommunity + nodespercommunity
                    W[indstart_i:indend_i, indstart_j:indend_j] = tmp_b

    # Make sure self-connections exist
    np.fill_diagonal(W, 1)

    if showplot:
	plt.figure()
        plt.imshow(W, origin='lower',cmap='bwr')
        plt.title('Structural Matrix with 5 communities\nand 1 hub-network', y=1.08)
        plt.xlabel('Regions')
        plt.ylabel('Regions')
        plt.colorbar()
        plt.tight_layout()
    
    return W

def generateSynapticNetwork(W, showplot=False):
    """
    Generate synaptic matrix over structural matrix with randomized gaussian weighs with
    mean = 1.0 and standard deviation of 0.2 (so all weights are positive)
    
    Parameters:
        W = structural connectivity matrix
        showplot = if set to True, will automatically display the structural matrix using matplotlib.pyplot

    Returns:
        Synaptic matrix with Gaussian weights on top of structural matrix
    """
    # Find non-zero connections
    G = np.zeros((W.shape))
    totalnodes = G.shape[0]
    connect_ind = np.where(W!=0)
    nconnects = len(connect_ind[0])
    weights = np.random.normal(loc=1.0,scale=0.2, size=(nconnects,))
    G[connect_ind] = weights
    
    # Find num connections per node
    nodeDeg = np.sum(W,axis=1)

    # Synaptic scaling according to number of incoming connections
    np.fill_diagonal(G,0)
    for col in range(G.shape[1]):
        G[:,col] = np.divide(G[:,col],np.sqrt(nodeDeg))
    #G = G/np.sqrt(totalnodes)

    if showplot:
	plt.figure()
        plt.imshow(G, origin='lower')#, vmin=0, vmax=20)
        plt.colorbar()
        plt.title('Synaptic Weight Matrix -- Coupling Matrix', y=1.08)
        plt.xlabel('Regions')
        plt.ylabel('Regions')
        plt.tight_layout()
        
    return G

def networkModel(G, Tmax=100,dt=.1,g=1.0,s=1.0,tau=1,I=None, noise=None):
    """
    G = Synaptic Weight Matrix
    Tmax = 100      (1sec / 1000ms)
    dt = .1         (1ms)
    g = 1.0         Coupling 
    s = 1.0         Self connection
    tau = 1.0       Time constant 
    I = 0.0         Stimulation/Task
    
    
    """
    T = np.arange(0, Tmax, dt)
    totalnodes = G.shape[0]

    # External input (or task-evoked input) && noise input
    if I==None: I = np.zeros((totalnodes,len(T)))
    # Noise parameter
    if noise == None: noise = np.zeros((totalnodes,len(T)))
    if noise == 1: noise = np.random.normal(size=(totalnodes,len(T)))

    # Initial conditions and empty arrays
    Enodes = np.zeros((totalnodes,len(T)))
    # Initial conditions
    Einit = np.random.rand(totalnodes,)
    Enodes[:,0] = Einit

    spont_act = np.zeros((totalnodes,))
    for t in range(len(T)-1):

        ## Solve using Runge-Kutta Order 2 Method
        # With auto-correlation
        spont_act = (noise[:,t] + I[:,t])
        k1e = -Enodes[:,t] + g*np.dot(G,phi(spont_act)) # Coupling
        k1e += s*phi(Enodes[:,t]) + spont_act# Local processing
        k1e = k1e/tau
        # 
        ave = Enodes[:,t] + k1e*dt
        #
        # With auto-correlation
        spont_act = (noise[:,t+1] + I[:,t+1])
        k2e = -ave + g*np.dot(G,phi(spont_act)) # Coupling
        k2e += s*phi(ave) + spont_act # Local processing
        k2e = k2e/tau

        Enodes[:,t+1] = Enodes[:,t] + (.5*(k1e+k2e))*dt

    return Enodes

def hrf(times):
    """ 
    Return values for HRF at given times 
    """
    # Gamma pdf for the peak
    peak_values = gamma.pdf(times, 6.0)
    # Gamma pdf for the undershoot
    undershoot_values = gamma.pdf(times, 12.0)
    # Combine them
    values = peak_values - 0.35 * undershoot_values
    # Scale max to 0.6
    return values / np.max(values) * 0.6

def convolveTimeseries(timeseries, samplingrate=1.0, TRLength=100):
    """
    Convolve the timeseries with the canonical  hemodynamic response function
    Also downsample the timeseries to an appropriate fMRI sampling rate
    Default sampling rate = 1.0 (10ms)
    Default TRLength = 100 (1 second)

    Parameter:
        timeseries = region x time activation time series
        samplingrate = The sampling rate of the model
        TRLength = The TR length to downsample the HRF-convolved time series too (since fMRI has low-sampling rate relative to neural processes)
    """

    simsample_rate=samplingrate/TRLength
    simsample_times = np.arange(0, 30, simsample_rate, dtype=float)
    hrf_at_simsample = hrf(simsample_times)

    # Convolve simulated time series with HRF
    nregions = timeseries.shape[0]
    timeseries_convolved=np.ones(np.shape(timeseries),dtype=float)
    for region in range(nregions):
        convolved = np.convolve(timeseries[region,:], hrf_at_simsample)
        n_to_remove = len(hrf_at_simsample) - 1
        convolved = convolved[:-n_to_remove]
        timeseries_convolved[region,:]=convolved

    # Downsample fMRI time series
    TR=TRLength
    dt_rec=samplingrate
    n_skip_BOLD = int(TR/dt_rec)
    BOLD_rec = timeseries_convolved[:,::n_skip_BOLD]

    return BOLD_rec

def runTaskGLM(taskdata, timing):
    """
    Runs a Task GLM to extract task regressors from simulated BOLD data
    
    Parameters:
        task_data = region x time matrix for simulated fMRI task data
        stimtimes = stimulus times (convolved with HRF and downsampled to fMRI TRs)

    Returns:
        task_betas = betas extracted from the single task regressor
        task_resids = residuals from the GLM
    """
    X = timing.T
    # Our model needs an intercept so we add a column of 1s:
    X = sm.add_constant(X)
    # Define parameters
    ntimepoints = taskdata.shape[1]
    nregions = taskdata.shape[0]
    # Define empty variables
    task_resids = np.zeros(taskdata.shape)
    task_betas = np.zeros((nregions,))
    # Performing a GLM for each region
    for region in range(0,nregions):
        y=taskdata[region,:]
        model = sm.OLS(y, X)
        results = model.fit()
        # Save betas
        task_betas[region]=results.params[1]
        # Save residuals
        task_resids[region,:]=results.resid
    
    return task_betas, task_resids

def computeActFlowBetweenNetworks(fc_mat, task_betas, ncommunities, nodespercommunity, posOnly=False):
    """
    Computes predicted task activations between all networks for a single task block
    
    Parameters:
        fc_mat = functional connectivity matrix 
        task_betas = GLM coefficients (betas) for task regressor for a single block
        ncommunities = number of network communities (default = 5)
        nodespercommunity = number of nodes per community (default = 50)
    """
    if posOnly: 
        tmp = fc_mat > 0
        fc_mat = np.multiply(tmp,fc_mat)
    actflowpredict = {}
    for i in range(ncommunities):
        actflowpredict[i] = {}
        for j in range(ncommunities):
            comm_ind1 = np.arange(i*nodespercommunity,(i*nodespercommunity)+nodespercommunity)
            comm_ind1.shape = (nodespercommunity,1)
            comm_ind2 = np.arange(j*nodespercommunity,(j*nodespercommunity)+nodespercommunity)
            comm_ind2.shape = (1,nodespercommunity)
            activity = task_betas[comm_ind1].T
            # Predicting act flow from network i to network j
            actflowpredict[i][j] = np.dot(activity,fc_mat[comm_ind1,comm_ind2])/float(len(comm_ind2))
            
    return actflowpredict

def ActFlowPipeline(fcmat_mreg, task_ts, tasktiming, ncommunities, nodespercommunity, samplingrate=1.0, TRLength=100):
    """
    Function that runs network-to-network activity flow mapping with simulated fMRI data (from simulated raw neural signals)
    Default sampling rate = 1.0 (10ms)
    Default TRLength = 100 (1 second)

    Parameters: 
        fcmat_mreg = resting-state functional connectivity matrix estimated with multiple linear regression
        task_ts = simulated task time series (unconvolved); convolution performed within this function
        tasktiming = task timing sampled at original model sampling rate (i.e., 10ms)
        ncommunities = number of network communities (default = 5)
        nodespercommunity = number of regions/nodes per community (default = 50)
        samplingrate = sampling rate (default = 1.0 = 10ms)
        TRLength = length of TR to downsample simulated fMRI data (default = 100 = 1sec = 10ms * 100)

    Returns:
       actflow_mreg = activity flow estimate predictions using multiple regression FC
       task_betas_postfmri = task regressor coefficients (betas) from simulated fMRI data

    """
    
    ####
    ## Task Analysis
    # Downsample and then convolve stimulus timing
    n_skip_BOLD = int(TRLength)
    # Temporally downsample task timing to seconds (before convolving)
    timingconv = tasktiming[:,::n_skip_BOLD]
    timingconv.shape = (timingconv.shape[1],)
    hrfsample_rate=1.0 # HRF is sampled at seconds
    hrfsample_times = np.arange(0, 30, hrfsample_rate, dtype=float)
    hrf_at_simsample = hrf(hrfsample_times)
    hrfconvtime = np.convolve(timingconv, hrf_at_simsample)
    n_to_remove = len(hrf_at_simsample) - 1
    convolved = hrfconvtime[:-n_to_remove]
    # Output
    timing_convolved_downsampled = convolved
    timing_convolved = timingconv

    # Convolve simulated task data into fMRI signal
    task_bold = convolveTimeseries(task_ts, 
                                   samplingrate=samplingrate, 
                                   TRLength=TRLength)

    ## Run GLM on simulated task fMRI data
    task_betas_postfmri, task_resids_postfmri = runTaskGLM(task_bold, timing_convolved_downsampled)
        
    ## Run ActFlow between every pair of networks do only postfmri for now using MultReg FC
    actflow_mreg = computeActFlowBetweenNetworks(fcmat_mreg, 
                                                 task_betas_postfmri, 
                                                 ncommunities, 
                                                 nodespercommunity, posOnly=False)
        
    return actflow_mreg, task_betas_postfmri


def runSubjectRuns(subj, s=1.0, g=1.0):
    """
    MAIN COMPUTATION METHOD

    Procedures performed in this function:
    1. Construct structural and synaptic connectivity matrices for single subject
    2. Run resting-state simulation for single subject; transform to fMRI data via convolution
    3. Estimate simulated resting-state FC
    4. Run task simulations for single subject; transform to fMRI data via convolution
    5. Perform network-to-network activity flow mapping
    """
    ######################
    #### 0. SET UP PARAMETERS
    ## Set up simulation parameters
    nblocks = 20
    ## First four tasks are top-down stimulation only (mimicking rule-encoding processes)
    topdown_only = range(1,5)
    ## Second four tasks are a combination of top-down and bottom-up stimulation (mimicking rule-encoding processes with subsequent stimuli presentatiosn)
    topdown_and_bottomup = range(5,9)

    ntasks = len(topdown_and_bottomup) + len(topdown_only)

    # Set number of time points for each task
    Tmax = 10000
    # Set number of time points for each rest
    Tmaxrest = 60000
  
    samplingrate = 1.0 # 10ms
    TRLength=100 # 1 second, i.e., 10ms * 100
    ######################


    ######################
    #### 1. CONSTRUCT STRUCTURAL AND SYNAPTIC CONNECTIVITY MATRICES
    # Parameters for subject's networks
    ncommunities = 5 # Number of communities
    innetwork_dsity = .35 # within network connectivity density
    outnetwork_dsity = .05 # out of network connectivity density
    hubnetwork_dsity = .20 # out of network connectivity density for hub network only

    nodespercommunity = 50 # number of nodes per community
    Ci = np.repeat(np.arange(ncommunities),nodespercommunity) # Construct a community affiliation vector
    hub_ind = np.where(Ci==0)[0] # Identify the regions associated with the hub network (hub network is by default the 0th network)
    totalnodes = nodespercommunity*ncommunities
   
    ##########
    # Construct structural matrix
    W = generateStructuralNetwork(ncommunities=ncommunities, innetwork_dsity=innetwork_dsity,
                                  outnetwork_dsity=outnetwork_dsity, hubnetwork_dsity=hubnetwork_dsity, 
                                  nodespercommunity=nodespercommunity, showplot=False)
    # Construct synaptic matrix
    G = generateSynapticNetwork(W, showplot=False)
    ######################

    
    ######################
    #### 2. RUN RESTING-STATE SIMULATION AND CONVOLVE TO FMRI DATA
    # Basic static parameters for simulation
    dt = 1.0 # Sampled @ 10ms
    T = np.arange(0,Tmax,dt)
    tau = 1.0
    
    # Instantiate task data dictionary
    taskdata = {}
    
    # Run task for 4 seconds every 10 seconds (aside from first 10 seconds)
    stimtimes = {}
    
    # Construct timing array for convolution -- this timing is irrespective of the task being performed
    # Tasks are only determined by which nodes are stimulated
    tasktiming = np.zeros((1,len(T)))
    for t in range(len(T)):
        if t%2000>500 and t%2000<1000:
            tasktiming[0,t] = 1.0
    

    #### Compute rest portion of analysis
    ## Run rest simulation prior to task simulations
    restdata = networkModel(G, Tmax=Tmaxrest,
                            dt=dt,g=g,s=s,tau=tau,I=None,noise=1)

    #print 'Convolving resting state data with HRF and generating intrinsic FC Matrix'
    # Convolve and downsample simulated timeseries
    rest_bold = convolveTimeseries(restdata,
                                   samplingrate=samplingrate,
                                   TRLength=TRLength)
    ######################


    ######################
    #### 3. ESTIMATE RESTING-STATE FC 
    # Pre-fMRI
    print 'Computing Resting-state FC Estimates'
    fcmat_pearson = np.corrcoef(rest_bold)
    fcmat_mreg = mreg.multregressionconnectivity(rest_bold)
    ######################

    
    ######################
    #### 4. RUN TASK SIMULATION AND CONVOLVE TO FMRI DATA
    # Generate output variables
    actflow_mreg2 = {}
    task_betas = {}
    for task in range(1,ntasks+1):
        print 'Running subject', subj, 'task', task
        # Generate empty output variables
        actflow_mreg2[task] = {}
        for i in range(ncommunities):
            actflow_mreg2[task][i] = {}
            for j in range(ncommunities):
                if i==j: continue
                actflow_mreg2[task][i][j] = np.zeros((nodespercommunity,nblocks))
        task_betas[task] = np.zeros((totalnodes,nblocks))
        
        #### Task Identification:  Selection stimulation nodes -- Stimulate 1/4 nodes out of nodespercommunity
        if task in topdown_only:
            taskcount = task-np.min(topdown_only)
            stimsize = np.floor(nodespercommunity/4.0)
            stim_nodes = np.arange((taskcount)*stimsize,(taskcount)*stimsize+stimsize,dtype=int)
            stimtimes[task] = np.zeros((totalnodes,len(T)))
        
        ## Runs simulations for stimulation to the hub network and ALL local network
        if task in topdown_and_bottomup:
            # Identify task number for this type of stimulation
            taskcount = task-np.min(topdown_and_bottomup) # Starts @ 0; iterates for len(flexandlocalnets)
            # Find stimulation nodes for the hub network
            networkstim = taskcount # since number of stimuli depend on number of tasks
            # Split half the stim size for each network for an equal number of stimulated nodes across networks
            stimsize = np.floor(nodespercommunity/4.0)
            stim_nodes = np.arange((networkstim)*stimsize,(networkstim)*stimsize+stimsize,dtype=int)
            
            # Find stimulation nodes of local communities
            # This starts @ 1 for task 1. Then iterates to community 2-4 for each subsequent task
            community = taskcount + 1  
            local_ind = np.where(Ci==community)[0] # Identify indices for the local community associated for this task (bottom-up stimulation regions)
            # Identify efferent connections from local network to hub network
            W_mask = np.zeros((W.shape))
            W_mask[local_ind,hub_ind] = 1.0
            local2hub_connects = np.multiply(W,W_mask)
            local_regions = np.where(local2hub_connects==1)[0]  # only look for local regions that innervate hub regions

            # Identify regions in the hub network and this local network that have connections
            stimsize = np.floor(nodespercommunity/4.0)
            stim_nodes2 = local_regions[0:stimsize] 

            # Stack hub nodes with local network nodes
            stim_nodes = np.hstack((stim_nodes,stim_nodes2))
            stimtimes[task] = np.zeros((totalnodes,len(T)))

        # Generate stim times
        for t in range(len(T)):
            if tasktiming[0,t] == 1:
                # Task stimulation every 10 seconds for 4 seconds, excluding the first 10 seconds
                stimtimes[task][stim_nodes,t] = .5
               # # Now account for tasks in which double as many nodes are stimulated
               # if task in topdown_and_bottomup:
               #     stimtimes[task][stim_nodes,t] = .25

        
        # Run every block for each task
        for block in range(nblocks):
            taskdata[task] = networkModel(G,Tmax=Tmax,dt=dt,g=g,s=s,tau=tau,
                                          I=stimtimes[task], noise=1)
            ######################

            
            ######################
            #### 5. RUN NETWORK-TO-NETWORK ACTIVITY FLOW MAPPING
            #### Run ActFlow analysis ####
            out = ActFlowPipeline(fcmat_mreg, taskdata[task], tasktiming, ncommunities, 
                                  nodespercommunity, samplingrate=samplingrate, TRLength=TRLength)
            # Collect outputs
            actflow_mreg, betas = out
            # Store actflow variables
            for i in range(ncommunities):
                for j in range(ncommunities):
                    if i==j: continue
                    actflow_mreg2[task][i][j][:,block] = actflow_mreg[i][j]
        
            task_betas[task][:,block] = betas

    return actflow_mreg2, task_betas, fcmat_pearson, fcmat_mreg


def subjectSimulationAndSaveToFile(subj,s=1.0,g=1.0,outdir='.'):
    """
    Run entire ActFlow procedure for single subject and saves to file specified as the 'outdir' parameter.
    Saves all files to an output directory named ItoEtAl2017_Simulations/

    Parameters:
        subj = subject number (as an int)
        s = self coupling parameter (default = 1)
        g = global coupling parameter (default = 1)
        outdir = Output directory
    """

    ## First write out to regular directory for regular subjects
    actflow_mreg, task_betas, fcmat_pearson, fcmat_mreg = runSubjectRuns(subj,s,g)
    
    ## Save files
    # On Linux server
    outdir = outdir + '/ItoEtAl2017_Simulations/'
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    ncommunities = 5

    # Make directories
    if not os.path.exists(outdir + 'actflow_predictions'): os.makedirs(outdir + 'actflow_predictions')
    if not os.path.exists(outdir + 'task_betas'): os.makedirs(outdir + 'task_betas')
    if not os.path.exists(outdir + 'restfc'): os.makedirs(outdir + 'restfc')

    # Save actflow predictions
    for i in range(ncommunities):
        for j in range(ncommunities):
            if i==j: continue
            for task in actflow_mreg:
                actflowconfig = '_net'+str(i)+'tonet'+str(j)
                file1b = 'actflow_predictions/subj'+str(subj)+'_task'+str(task)+actflowconfig+'_multregFC.txt'
                np.savetxt(outdir + file1b, actflow_mreg[task][i][j], delimiter=',')
    
    # Save task betas
    for task in task_betas:
        file4 = 'task_betas/subj'+str(subj)+'_task'+str(task)+'_allblocks.txt'
        np.savetxt(outdir + file4, task_betas[task], delimiter=',')

    # Save FC Matrices
    file6a = 'restfc/subj'+str(subj)+'_restfc_pearson.txt'
    file6b = 'restfc/subj'+str(subj)+'_restfc_multreg.txt'
    np.savetxt(outdir + file6a, fcmat_pearson,delimiter=',')
    np.savetxt(outdir + file6b, fcmat_mreg, delimiter=',')
    
