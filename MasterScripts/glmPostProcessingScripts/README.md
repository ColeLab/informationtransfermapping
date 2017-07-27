# Post-processing scripts - nuisance regression for resting-state fMRI, and standard fMRI GLM activation scripts for obtaining task-evoked activation estimates
## Scripts included in this directory were used process minimally preprocessed data, and to prepare the data to be analyzed. Two task-evoked activation GLMs are run: (1) at the 64k vertex-level (on the CIFTI surface) and (2) at the parcellated region-level (using the Glasser et al. (2016) atlas). Two nuisance regressions are run for resting-state fMRI: (1) at the 64k vertex-level (on the CIFTI surface) and (2) at the parcellated region-level using the Glasser et al. (2016) atlas).

#### Takuya Ito (takuya.ito@rutgers.edu)
#### Citation: Ito T, Kulkarni KR, Schultz DH, Mill RD, Chen RH, Solomyak LI, Cole MW (2017). Cognitive task information is transferred between brain regions via resting-state network topology. bioRxiv. https://doi.org/10.1101/101782
#### Last update: 06/26/2017

## Directory organization:
**File:** MainScript.m
* The main script that contains calls to all the major/main functions, running either a standard fMRI task GLM post-processing, or functions that perform nuisance regression on resting-state fMRI data.

**Directory:** regionGlasserLevel/
* Directory containing scripts for running regressions on the Glasser et al. (2016) parcellated time series. Data is read in as a 360Region X Time matrix, and regressions are performed on each of the 360 time series individually.
* **File:** GlasserGLM_miniblock_betaseries.m - function that performs a task-evoked GLM analysis, regressing out the activity of each miniblock separately.
* **File:** GlasserGLM_restdata.m - function that performs nuisance regression on parcellated resting-state fMRI data.
* **N.B.** - all other files are helper functions for the main functions (e.g., loading in regressors/task timing files).

**Directory:** vertex64kLevel
* Directory containing scripts for running regressions on the dense surface with 64k vertices. Data is read in from CIFTI format, organized as a 64k Vertex X Time matrix, and regressions are performed on each of the 64k vertices individually.
* **File:** SurfaceGLM_miniblock_betaseries.m - function that performs a task-evoked GLM analysis, regressing out the activity of each miniblock separately.
* **File:** SurfaceGLM_rest.m - function that performs nuisance regression on each of the 64k vertices during resting-state fMRI data.
* **N.B.** - all other files are helper functions for the main functions (e.g., loading in regressors/task timing files).

**File:** parcellate_Glasser2016Parcels.sh
* A bash script that uses the command line workbench command to generate parcellated time series using the Glasser et al. (2016) atlas. Reads in minimally preprocessed data, preprocessed by the HCP Pipelines.
* This function averages across all vertices within a region, and provides a downsampled region X time parcellated time series.

