# Information transfer mapping -- Master scripts
## Scripts included in this directory were used to generate all reported results (including statistics and figures) in the main manuscript

#### Takuya Ito (taku.ito1@gmail.com)
#### Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017). Cognitive task information is transferred between brain regions via resting-state network topology. bioRxiv. https://doi.org/10.1101/101782
#### Last update: 04/24/2017

## Directory organization:
**Directory:** glmPostProcessingScripts/
* Directory containing all scripts (mostly in MATLAB) that takes HCP-minimally preprocessed data and performs-post processing using regression 
* Includes scripts for resting-state nuisance regression
* Includes scripts for running standard fMRI GLM analysis (for task-evoked activity estimation) using MATLAB
* Performs both types of regression on both parcellated time series (using the Glasser parcels) and 64k-vertex-level time series
* See sub-directory for more details

**Directory:** SupercomputerScripts/
* Directory containing scripts (shell-scripts, Python, MATLAB) that performs computations on a supercomputer
* Any usage of this code is mentioned/reference in the corresponding jupyter notebooks
* Includes:
  * Code for region-to-region information transfer mapping (as well as principal components regression for vertex-wise FC estimation and decoding task performance with information transfer estimates))
  * Code for network-to-network FC permutation testing (SFig. 4)

**Directory:** utils/
* Directory that contains miscellaneous code necessary to run analyses in corresponding jupyter notebooks (Manuscript\*ipynb)
 

**Files:** Manuscript\*.ipynb
* Contains the primary jupyter notebooks used to report statistical results and create figures for the manuscript
* File names are organized by figures. For example, Manuscript3\*.ipynb corresponds to the notebook that generates results related to Figure 3. ManuscriptS1\*.ipynb corresponds to the notebook that generates results/figures related to Supplementary Figure 1.
* Note: These notebooks contain raw code used to generate results/figures, and will fail/break if run on any other system. They were written to perform on a specific lab server.
