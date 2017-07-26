function output = GlasserGLM_restdata(subj, gsr, nproc)
% This function runs a GLM on the Glasser Parcels (380) on a single subject
% This will only regress out noise parameters and HRF convolved miniblock onsets 
%

% Input parameter:
%   subj = subject number as a string
    
    numTasks = 8;
    numTRs = numTasks*581;
    basedir = ['/projects/ModalityControl/data/' subj '/analysis/'];

    data = loadGlasserData(subj); 
    data = data.rest;
    numROIs = size(data,1);
    

    % Load only noise regressors and task regressors for data
    X = loadStimFiles_byMiniblockV3(subj,gsr);
    taskRegs = X.restRegressors; % These include the two binary regressors

    parfor (regionNum=1:numROIs,nproc)
        ROITimeseries = data(regionNum,:);
        %disp(['Regressing out region number ' num2str(regionNum) ' out of ' num2str(numROIs)])
        stats = regstats(ROITimeseries', taskRegs, 'linear', {'r', 'beta', 'rsquare'});

        % Collect regression results
        residual_dtseries(regionNum, :) = stats.r';
        betas_dtseries(regionNum, :) = stats.beta;

    end

    % Write out 2D Residual dtseries to CSV
    if gsr==0
        outname1 = ['/projects2/ModalityControl2/data/resultsGlasser/glmRest_GlasserParcels/' subj '_rest_nuisanceResids_Glasser.csv'];
        
    elseif gsr==1
        outname1 = ['/projects2/ModalityControl2/data/resultsGlasser/glmRest_GlasserParcels/' subj '_rest_nuisanceResids_Glasser_noGSR.csv'];
    end
    outname2 = ['/projects2/ModalityControl2/data/resultsGlasser/glmRest_GlasserParcels/' subj '_rest_taskBetas_Glasser.csv'];

    csvwrite(outname1, residual_dtseries)
    csvwrite(outname2, betas_dtseries)
    output.residual_dtseries = residual_dtseries;
    output.betas_dtseries = betas_dtseries;
end



