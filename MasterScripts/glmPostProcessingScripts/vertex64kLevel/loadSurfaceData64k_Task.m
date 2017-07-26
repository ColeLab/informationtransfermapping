function output = loadSurfaceData64k_Task(subj)
% Taku Ito
% 08/19/2016
%
% Parameters: subj ( must be in put with single quotations, i.e., as a string)

    addpath('/projects/AnalysisTools/gifti-1.6/')

    % Get data directory based on the assumption that script is in projects/IndivRITL/docs/scripts/modalControl/
    datadir = ['/projects2/ModalityControl2/data/'];
    numTasks = 8;

    % Load in task data
    for task=1:numTasks

        inFile = [datadir subj '/analysis/Task' num2str(task) '_Atlas_downsampled64k.dtseries.nii'];
        % Now, import the data to MATLAB (for a particular task) using ft_read_cifti in AnalysisTools
        disp(['Importing 64k Surface cifti file for task ' num2str(task)])
        tmp = ciftiopen(inFile,'wb_command');
        data{task} = tmp.cdata;

    end

    % Now concatenate as a matrix
    nrois = size(data{1},1);
    nTRs = size(data{1},2);
    output_tmp = zeros(nrois,nTRs,numTasks);

    for task=1:numTasks
        % Now demean each task
        output_tmp(:,:,task) = data{task};
        % Demean each run
        disp(['Demeaning each run separately on task prior to concatenation...'])
        mean_ts = mean(output_tmp(:,:,task),2);
        output_tmp(:,:,task) = output_tmp(:,:,task) - mean_ts(:,ones(1,nTRs));
    end


    % Reshape into two dimensions, nrois x 4648
    output.task = reshape(output_tmp, [nrois,4648]);
end
