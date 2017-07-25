# Taku Ito
# 07/14/2017
##Modified by MWC June 20, 2017
##Modified by TI June 20, 2017

# Code to perform permutation testing to control for family-wise error (FWE)
# Using max-T approach as described in Nichols & Holmes (2002)
# Nichols TE, Holmes AP. (2002). Nonparametric permutation tests for functional neuroimaging: A primer with Examples. Hum. Brain Mapp., 15: 1-25. doi:10.1002/hbm.1058

import numpy as np
import scipy.stats as stats
import multiprocessing as mp
from statsmodels.distributions.empirical_distribution import ECDF


def maxT(diff_arr, nullmean=0, alpha=.05, tail=1, permutations=1000, nproc=1, pvals=False):
    """
    Performs family-wise error correction using permutation testing (Nichols & Holmes 2002)
    Note! Assumes a one-sided t-test (specify tail of test by tail parameter) 

    Citation: 
        Nichols TE, Holmes AP. (2002). Nonparametric permutation tests for functional neuroimaging: A primer with Examples. Hum. Brain Mapp., 15: 1-25. doi:10.1002/hbm.1058
    Required Parameters:
        diff_arr    =   MxN matrix of set of M independent tests for condition 1 minus condition 2 across N subjects
                        diff_arr can also be an array of multiple values (or tests) compared against the nullmean (or null mean)
    Optional Parameters:
        nullmean    =   Expected value of the null hypothesis {default = 0, for a t-test against 0}
        alpha       =   alpha value to return the maxT threshold {default = .05}
        tail        =   [0, 1, or -1] 
                        If tail = 1, reject the null hypothesis if the statistic is greater than the null dist (upper tailed test).  
                        If tail = -1, reject the null hypothesis if the statistic is less than the null dist (lower tailed test). 
                        If tail = 0, reject the null hypothesis for a two-tailed test
                        {default : 0} 
        permutations =  Number of permutations to perform {default = 1000}
        nproc       =   number of processes to run in parallel {default = 1}
        pvals       =   if True, returns equivalent p-value distribution for all t-values {default = True}

    Returns:
        t: Array of T-values of correct contrast map (Mx1 vector, for M tests)
        maxTThreshold   : The t-value threshold corresponding to the corrected alpha value. If a two-tailed test is specified, the maxR is provided as an absolute value
        p (optional)    : Array of FWE-corrected p-values (Mx1 vector, for M tests); 

    """
    # Focus on difference matrix -- more computationally feasible (and less data to feed into parallel processing)

    # Prepare inputs for multiprocessing
    inputs = []
    for i in range(permutations):
        seed = np.random.randint(0,100000,1)[0]
        inputs.append((diff_arr,nullmean,tail,seed))

    pool = mp.Pool(processes=nproc)
    result = pool.map_async(_maxTpermutation,inputs).get()
    pool.close()
    pool.join()

    # Returns an array of T-values distributions (our null distribution of "max-T" values)
    maxT_dist = np.asarray(result)

    #Find threshold for alpha
    maxT_dist_sorted = np.sort(maxT_dist)
    # Specify which tail we want
    if tail == 1:
        topPercVal_maxT_inx = int(len(maxT_dist_sorted)*(1-alpha))
        maxT_thresh = maxT_dist_sorted[topPercVal_maxT_inx]
    elif tail == -1:
        topPercVal_maxT_inx = int(len(maxT_dist_sorted)*(alpha))
        maxT_thresh = maxT_dist_sorted[topPercVal_maxT_inx]
    elif tail == 0:
        topPercVal_maxT_inx = int(len(maxT_dist_sorted)*(1-alpha))
        maxT_thresh = maxT_dist_sorted[topPercVal_maxT_inx]
#    elif tail == 0: 
#        topPercVal_maxT_inx = int(len(maxT_dist_sorted)*(alpha/2.0))
#        botPercVal_maxT_inx = int(len(maxT_dist_sorted)*(1-alpha/2.0))
#        # Provide two r thresholds for two-tailed test 
#        topT_thresh = maxT_dist_sorted[topPercVal_maxT_inx]
#        botT_thresh = maxT_dist_sorted[botPercVal_maxT_inx]


    # Obtain real t-values 
    t = stats.ttest_1samp(diff_arr, nullmean, axis=1)[0]

    if pvals:
        # Construct ECDF from maxT_dist
        ecdf = ECDF(maxT_dist)

        # Return p-values from maxT_dist using our empirical CDF (FWE-corrected p-values)
        p_fwe = ecdf(t)

        if tail == 1:
            p_fwe = 1.0 - p_fwe
        
        return t, maxT_thresh, p_fwe

    else:
        return t, maxT_thresh


def _maxTpermutation((diff_arr,nullmean,tail,seed)):
    """
    Helper function to perform a single permutation
    """

    np.random.seed(seed)

    # Create a random matrix to shuffle conditions (randomly multiply contrasts by 1 or -1)
    shufflemat = np.random.normal(0,1,diff_arr.shape)
    pos = shufflemat > 0
    neg = shufflemat < 0
    # matrix of 1 and -1
    shufflemat = pos + neg*(-1)

    # Shuffle raw values
    diff_arr = np.multiply(diff_arr, shufflemat)

    # Take t-test against 0 for each independent test 
    t_matrix = stats.ttest_1samp(diff_arr,nullmean,axis=1)[0] 

    if tail==1:
        maxT = np.max(t_matrix)
    elif tail==-1:
        maxT = np.min(t_matrix)
    elif tail==0:
        maxT = np.max(np.abs(t_matrix))
    
    return maxT


### 
# maxR approach (for individual difference correlations)
def maxR(data_arr, behav_arr, alpha=.05, tail=0, permutations=1000, nproc=1, pvals=False):
    """
    Performs family-wise error correction using permutation testing (Nichols & Holmes 2002)
    Using correlations as test-statistic (as opposed to T-Stat. Can be used for RSA type analysis or individual difference correlations.
    Note! Assumes a one-tailed test (specify tail of test by tail parameter) 

    Citation: 
        Nichols TE, Holmes AP. (2002). Nonparametric permutation tests for functional neuroimaging: A primer with Examples. Hum. Brain Mapp., 15: 1-25. doi:10.1002/hbm.1058
    Required Parameters:
        data_arr    =   MxN matrix of set of M independent measurements (e.g., FC-values) across N subjects
        behav_arr   =   Nx1 array of behavioral measures for N subjects
    Optional Parameters:
        alpha       =   alpha value to return the maxT threshold {default = .05}
        tail        =   [0,1, or -1] 
                        If tail = 1, reject the null hypothesis if the correlation is greater than the null dist (upper tailed test).  
                        If tail = -1, reject the null hypothesis if the correlation is less than the null dist (lower tailed test). 
                        If tail = 0, reject the null hypothesis for a two-tailed test
                        {default : 0} 
        permutations =  Number of permutations to perform {default = 1000}
        nproc       =   number of processes to run in parallel {default = 1}
        pvals       =   if True, returns equivalent p-value distribution for all t-values {default = True}

    Returns:
        r               : Array of Pearson-r values of the true correlations map (Mx1 vector, for M tests)
        maxRThreshold   : The Pearson-r value corresponding to the corrected alpha value. If a two-tailed test is specified, the absolute value of the maxR threshold is provided.
        p (optional)    : Array of FWE-corrected p-values (Mx1 vector, for M tests); 

    N.B.: Only works for paired one-sample t-tests
    """
    # Calculating the TRUE Pearson correlations in a vectorized format (increasing speed)
    data_normed = stats.zscore(data_arr,axis=1)
    behav_normed = stats.zscore(behav_arr)
    trueR = np.mean(np.multiply(behav_normed,data_normed),axis=1)

    # Prepare inputs for multiprocessing
    inputs = []
    for i in range(permutations):
        seed = np.random.randint(0,100000,1)[0]
        inputs.append((data_normed,behav_normed,tail,seed))

    pool = mp.Pool(processes=nproc)
    result = pool.map_async(_maxRpermutation,inputs).get()
    pool.close()
    pool.join()

    # Returns an array of T-values distributions (our null distribution of "max-T" values)
    maxR_dist = np.asarray(result)

    #Find threshold for alpha
    maxR_dist_sorted = np.sort(maxR_dist)
    # Specify which tail we want
    if tail == 1:
        topPercVal_maxR_inx = int(len(maxR_dist_sorted)*(1-alpha))
        maxR_thresh = maxR_dist_sorted[topPercVal_maxR_inx]
    elif tail == -1:
        topPercVal_maxR_inx = int(len(maxR_dist_sorted)*(alpha))
        maxR_thresh = maxR_dist_sorted[topPercVal_maxR_inx]
    elif tail == 0:
        topPercVal_maxR_inx = int(len(maxR_dist_sorted)*(1-alpha))
        maxR_thresh = maxR_dist_sorted[topPercVal_maxR_inx]
#    elif tail == 0: 
#        topPercVal_maxR_inx = int(len(maxR_dist_sorted)*(alpha/2.0))
#        botPercVal_maxR_inx = int(len(maxR_dist_sorted)*(1-alpha/2.0))
#        # Provide two r thresholds for two-tailed test 
#        topR_thresh = maxR_dist_sorted[topPercVal_maxR_inx]
#        botR_thresh = maxR_dist_sorted[botPercVal_maxR_inx]
    

    if pvals:
        # Construct ECDF from maxT_dist
        ecdf = ECDF(maxR_dist)

        # Return p-values from maxT_dist using our empirical CDF (FWE-corrected p-values)
        p_fwe = ecdf(trueR)

        if tail == 1:
            p_fwe = 1.0 - p_fwe
        
        #if tail!=0:
        return trueR, maxR_thresh, p_fwe

    else:
        #if tail!=0:
        return trueR, maxR_thresh


def _maxRpermutation((data_normed,behav_normed,tail,seed)):
    """
    Helper function to perform a single permutation
    Assumes the first row are the labels (or behavioral measures)
    """

    np.random.seed(seed)

    # Randomly permute behavioral data along 2nd dimension (subjects). Note: np.random.shuffle() requires transposes
    np.take(behav_normed,np.random.rand(len(behav_normed)).argsort(),out=behav_normed)
    # Randomly permute measurement data along 2nd dimension (subjects). Note: np.random.shuffle() requires transposes
    #np.take(data_normed,np.random.rand(data_normed.shape[1]).argsort(),axis=1,out=data_normed)
    # Calculating Pearson correlations in a vectorized format (increasing speed)
    r_values = np.mean(np.multiply(behav_normed,data_normed),axis=1)

    if tail==1:
        maxR = np.max(r_values)
    elif tail==-1:
        maxR = np.min(r_values)
    elif tail==0:
        maxR = np.max(np.abs(r_values))

    return maxR


