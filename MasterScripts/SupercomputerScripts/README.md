# Supercomputer scripts - Scripts that required use of a compute cluster 
## Raw code used to generate the computations and results/figures for Fig. 6, Fig. 7, and Supplementary Fig. 4-7 required use of a compute cluster.

#### Takuya Ito (takuya.ito@rutgers.edu)
#### Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017). Cognitive task information is transferred between brain regions via resting-state network topology. bioRxiv. https://doi.org/10.1101/101782
#### Last update: 06/26/2017

## Directory organization:
**Directory:** Fig6_RegionToRegionInformationTransferMapping
* **File:** Fig6_RegionToRegionInformationTransferMapping/ActFlow_ITE_DecodePerformance_LogRegression_v2.m
    * The main script that contains MATLAB functions (and helper functions) used to perform region-to-region information transfer mapping, as well as decode task performance using the block-wise region-to-region information transfer estimates.
* **File:** Fig6_RegionToRegionInformationTransferMapping/ActFlow_ITE_DecodePerformance_LogRegression_v2.sh
    * The batch script to submit jobs to the cluster, calling the main function of the MATLAB script.


**Directory:** Fig4_FCPermutationTest
* **File:** Fig4_FCPermutationTest/runPermutations_withinNetPerm.py
    * The main script that contains Python code functions (and helper functions) used to perform network-to-network information transfer mapping permutation tests. Here, each permutation permuted the FC patterns within a pair of networks, ensuring that the successful transfer of information depended on the FC pattern/topology between networks.
* **File:** Fig4_FCPermutationTest/runPermutations_withinNetPerm.sh
    * The batch script to submit jobs to the cluster, calling the main function of the Python script.

**N.B.: All other files are 'helper functions' referenced in the main scripts**
