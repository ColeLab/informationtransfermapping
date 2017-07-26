# Taku Ito
# 07/14/2017

# Code to perform permutation testing to control for family-wise error (FWE)
# Using max-T approach as described in Nichols & Holmes (2002)


import numpy as np
import scipy.stats as stats
import multiprocessing as mp
from statsmodels.distributions.empirical_distribution import ECDF


def permutationFWE(diff_arr, nullmean=0, permutations=1000, nproc=1):
    """
    Performs family-wise error correction using permutation testing (Nichols & Holmes 2002)
    
    Parameters:
        diff_arr = MxN matrix of set of M independent tests for condition 1 minus condition 2 across N subjects
        permutations = Number of permutations to perform
        nproc = number of processes to run in parallel

    Returns:
        t: Array of T-values of correct contrast map (Mx1 vector, for M tests)
        p: Array of FWE-corrected p-values (Mx1 vector, for M tests); 
           Note, p-values correspond to values on the CDF. One-sided or or two-sided p-values can be computed accordingly.

    N.B.: Only works for paired one-sample t-tests
    """
    # Focus on difference matrix -- more computationally feasible (and less data to feed into parallel processing)

    # Prepare inputs for multiprocessing
    inputs = []
    for i in range(permutations):
        seed = np.random.randint(0,100000,1)[0]
        inputs.append((diff_arr,nullmean,seed))

    pool = mp.Pool(processes=nproc)
    result = pool.map_async(_permutation,inputs).get()
    pool.close()
    pool.join()

    # Returns an array of T-values distributions (our null distribution of "max-T" values)
    maxT_dist = np.asarray(result)


    # Obtain real t-values 
    t = stats.ttest_1samp(diff_arr, nullmean, axis=1)[0]
    #t = np.mean(diff_arr,axis=1)

    # Construct ECDF from maxT_dist
    ecdf = ECDF(maxT_dist)

    # Return p-values from maxT_dist using our empirical CDF (FWE-corrected p-values)
    p_fwe = ecdf(t)
    
    return t, p_fwe




def _permutation((diff_arr,nullmean,seed)):
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

    # Shuffle information estimates
    diff_arr = np.multiply(diff_arr, shufflemat)

    # Take t-test against 0 for each independent test 
    t_matrix = stats.ttest_1samp(diff_arr,nullmean,axis=1)[0] 
    #t_matrix = np.mean(diff_arr,axis=1)

    maxT = np.max(t_matrix)
    
    return maxT



