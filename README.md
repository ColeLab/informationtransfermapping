# Information transfer mapping
#### Public code release for Ito et al. (2017)
#### Author: Takuya Ito (takuya.ito@rutgers.edu)
#### Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017). Cognitive task information is transferred between brain regions via resting-state network topology. bioRxiv. https://doi.org/10.1101/101782
#### Last update: 04/24/2017
#### Refer to methods and supplemental materials for full descriptions

## Directory organization:
**File:** ItoEtAl2017_Simulations.zip
* A zip file containing simulated data generated from computational model
* Unzip prior to running ItoEtAl2017_ComputationalModelGroupAnalysis.ipynb

**File:** ItoEtAl2017_ComputationalModelTutorial.ipynb
* Contains a introductory tutorial to the formulation of the computational model used in the study
* Generates simple and illustrative figure of the simulation data

**Module:** informationtransfermapping.py
* Contains demo code and demo functions for performing information transfer mapping. See EmpiricalResults demo for example implementation. 

**Directory:** TheoreticalResults/
* Contains demos for our computational (simulation) model results. Specifically, replicates Fig. 4 and Supplementary Fig. 3. See subdirectory for other details. 

**Directory:** EmpiricalResults/
* Contains demos for our empirical network-to-network information transfer mapping results. Specifically, replicates Supplementary Fig. 1. See subdirectory for other details. 

**Directory:** utils/
* **File: multregressionconnectivity.py** Contains a function to estimate functional connectivity using multiple linear regression
