# Information transfer mapping -- network-to-network information transfer results (based on empirical fMRI data)

#### Author: Takuya Ito (takuya.ito@rutgers.edu)
#### Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017). Cognitive task information is transferred between brain regions via resting-state network topology. bioRxiv. https://doi.org/10.1101/101782
#### Last update: 04/24/2017
#### Refer to methods and supplementary materials for full description of method

## Directory organization:
**File:** data.zip
* A zip file containing all empirical preprocessed data necessary to generate panels in Supplementary Fig. 1
* Contains MiniblockActivations: Whole brain, miniblock beta activation coefficients (for all 128 miniblocks). Obtained using a beta series regression on each miniblock in the experimental paradigm, for every brain region in the Glasser et al. (2016) atlas. One file for every subject.
* FC_Estimate: Whole brain functional connectivity matrix obtained during resting-state fMRI. Estimated using multiple linear regression (see Methods).
* CPROTaskIdentifiers: Task-rule condition file. Every subject has three files, one for each rule domain (logic, sensory, and motor). Each csv file is organized as a sample (i.e., miniblock) X condition (task-rule) matrix. For each condition (or task-rule), there are 32 samples. Since there are four task-rules for each rule domain, we have 128 elements in total (corresponding to the 128 miniblocks).
* network_array.csv: Network assignments for the Glasser et al. (2016) atlas. Organized as a 1D array.

**File:** Demo_Network2NetworkInformationTransfer.ipynb
* Jupyter notebook that performs network-to-network information transfer mapping on the provided empirical data. Replicates panels shown in Supplementary Fig. 1

**File:** networkInformationTransferMapping.py
* Main python module containing functions to perform empirical network-to-network information transfer mapping.

