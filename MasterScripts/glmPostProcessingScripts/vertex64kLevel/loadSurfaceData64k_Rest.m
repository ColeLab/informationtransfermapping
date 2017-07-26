function output = loadSurfaceData64k_Rest(subj)
% Taku Ito
% 08/19/2016
%
% Parameters: subj ( must be in put with single quotations, i.e., as a string)

    addpath('/projects/AnalysisTools/gifti-1.6/')

    % Get data directory based on the assumption that script is in projects/IndivRITL/docs/scripts/modalControl/
    datadir = ['/projects2/ModalityControl2/data/'];
    numTasks = 8;

    % Load in rest surface data first
    inFile = [datadir subj '/analysis/Rest1_Atlas_64k.dtseries.nii'];

    % Now import the data to MATLAB using ft_read_cifti in AnalysisTools
    disp(['Importing surface (64k) cifti file for rest data'])
    tmp = ciftiopen(inFile, 'wb_command');
    rest = tmp.cdata;

    output.rest = rest;

end
