# Taku Ito
# 05/29/16

#Multiple linear regression for FC approximation in python (templated off Mike's MATLAB script)

import numpy as np
import statsmodels.api as sm
from sklearn.decomposition import PCA

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


def pcaRegressionConnectivity(activityMatrix, ncomponents=500, outdir='/projects/ModalityControl/data/results/multregconn_restfc_vertices/',filename='default.csv'):
    """
    Inputs:
    activityMatrix - regions X time matrix
    ncomponents - number components to run regression on

    Returns:
    region X region multiple regression connectivity matrix
    """

    nregions = activityMatrix.shape[0]
    interaction_mat = np.zeros((ncomponents,ncomponents));

    pca = PCA(n_components=ncomponents)
    reduced_mat = pca.fit_transform(activityMatrix.T)
    reduced_mat = reduced_mat.T

    for targetregion in range(ncomponents):
        otherregions = range(ncomponents)
        otherregions = np.delete(otherregions, targetregion) # Delete target region from 'other regiosn'
        X = reduced_mat[otherregions,:].T
        # Add 'constant' regressor
        X = sm.add_constant(X)
        y = reduced_mat[targetregion,:]
         ##*** STARTHERE
        model = sm.OLS(y, X)
        results = model.fit()
        interaction_mat[otherregions, targetregion] = results.params[1:] # all betas except for constant betas

    # Now tranpose matrix back to original space, i.e., ncomponents -> nvertices/nregions
    # Since it's not in samples x component space and just run each 'sample' or component separately
    multreg_mat = np.zeros((nregions,nregions))
    for sample in range(ncomponents):
        feat_array = pca.inverse_transform(interaction_mat[sample,:])
        multreg_mat[sample,:] = feat_array

    outfile = outdir + filename
    # Output a csv file
    np.savetxt(outfile, multreg_mat, delimiter=',')
        
    return multreg_mat

