# Taku Ito
# 05/29/16

#Multiple linear regression for FC approximation in python (templated off Mike's MATLAB script)

import numpy as np
import statsmodels.api as sm

def multregressionconnectivity(activityMatrix):
    """
    Activity matrix should be region/voxel X time
    Assumes all time series are de-meaned
    """

    nregions = activityMatrix.shape[0]
    timepoints = activityMatrix.shape[1]
    if nregions > timepoints:
         raise Exception('More regions (regressors) than timepoints! Use Principle component regression')

    interaction_mat = np.zeros((nregions,nregions))

    for targetregion in range(nregions):
        otherregions = range(nregions)
        otherregions = np.delete(otherregions, targetregion) # Delete target region from 'other regiosn'
        X = activityMatrix[otherregions,:].T
        # Add 'constant' regressor
        X = sm.add_constant(X)
        y = activityMatrix[targetregion,:]
        model = sm.OLS(y, X)
        results = model.fit()
        interaction_mat[otherregions, targetregion] = results.params[1:] # all betas except for constant betas
        
    return interaction_mat 


