function output = loadGlasserData(subj)
% Taku Ito
% 3/12/15
%
% This script loads in the surface data using the Glasser et al. 2016 380 ROI surface parcellation scheme
% 
% Parameters: subj ( must be in put with single quotations, i.e., as a string)

    % Get data directory based on the assumption that script is in projects/IndivRITL/docs/scripts/modalControl/
    datadir = ['/projects2/ModalityControl2/data/'];
    analysisdir = [datadir '/' subj '/analysis/'];
    numTasks = 8;

    % Each dtseries is 380x581, so we will want to create an empty matrix of 380x4648 (8*581 = 4648)
    % We will create a 3d Matrix first, and reshape it into a 2d roi x dtseries after
    nrois = 360;
    output_tmp = zeros(nrois,581,8);
    
    for task=1:numTasks

        % First need to parcellate each surface's Task${num}.dtseries.nii using workbench command
        inFile = [datadir subj '/analysis/Task' num2str(task) '_Atlas.LR.Glasser2016Parcels.32k_fs_LR.ptseries.nii'];

        % Now, import the data to MATLAB (for a particular task) using ft_read_cifti in AnalysisTools
        disp(['Importing parcellated cifti file for task ' num2str(task)])
        tmp = ciftiopen(inFile,'wb_command');
        data{task} = tmp.cdata;

        % Now concatenate the time series across all 8 tasks. 
        % The variable data is a cell with 8 structs. To get the dtseries of each struct, we acceess them via data{task}.dtseries
        output_tmp(:,:,task) = data{task};
        % Demean each run
        disp(['Demeaning each run separately on task prior to concatenation...'])
        for roi=1:nrois
            output_tmp(roi,:,task) = output_tmp(roi,:,task) - mean(output_tmp(roi,:,task));

        end
    end

    % Parcellate and import rest data too
    inFile = [datadir subj '/analysis/Rest1_Atlas.LR.Glasser2016Parcels.32k_fs_LR.ptseries.nii'];

    % Now import the data to MATLAB using ft_read_cifti in AnalysisTools
    disp(['Importing parcellated cifti file for rest data'])
    tmp = ciftiopen(inFile, 'wb_command');
    rest = tmp.cdata;

    % Reshape into two dimensions, nrois x 4648
    output.task = reshape(output_tmp, [nrois,4648]);
    output.rest = rest;

end
