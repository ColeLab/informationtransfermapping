# Information transfer mapping
#### Public code release for Ito et al. (2017)
#### Takuya Ito (takuya.ito@rutgers.edu)
#### Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017) “Cognitive task information is transferred between brain regions via resting-state network topology.”. Nature Communications. 8:1027. PMID: 29044112 http://doi.org/10.1038/s41467-017-01000-w 
#### Last update: 07/26/2017
#### Refer to methods and supplementary materials for full descriptions

## Directory organization:
**Module:** informationtransfermapping.py
* Contains demo code and demo functions for performing information transfer mapping. See EmpiricalResults demo for example implementation. 

**Directory:** TheoreticalResults/
* Contains demos for our computational (simulation) model results. Specifically, replicates Fig. 4 and Supplementary Fig. 3. See subdirectory for other details. 

**Directory:** EmpiricalResults/
* Contains demos for our empirical network-to-network information transfer mapping results. Specifically, replicates Supplementary Fig. 1. See subdirectory for other details.

**Directory:** MasterScripts/
* Contains all raw code used to run analyses and generate figures from HCP minimally preprocessed material. Code is not organized as "tutorial" code. See subdirectory for more details.

**Directory:** utils/
* **File: multregressionconnectivity.py** Contains a function to estimate functional connectivity using multiple linear regression.
* **File: permutationTesting.py** Contains a function to run FWE-correction using permutation testing.
