% Taku Ito
% 08/19/2016 

% Group analysis (sanity check)


addpath('/projects/AnalysisTools/BCT/2014_04_05 BCT/')
addpath('/projects/AnalysisTools/gifti-1.6/')
subjNums = '032 033 037 038 039 045 013 014 016 017 018 021 023 024 025 026 027 031 035 046 042 028 048 053 040 049 057 062 050 030 047 034';
subjNumStr = strread(subjNums, '%s', 'delimiter', ' ');

addpath('regionLevel/')
addpath('vertexLevel/')


%%
% Runs subject-wise GLM on Glasser 360 parcels for noise regression (rest) and task-evoked regression (task) in parallel
execute = 0;
if execute == 1
    nproc = 3; 
    gsr = 0;
    for (i=1:32)
        disp(['Running GLMs on Glasser Parcels for each subject ' num2str(i)])
        GlasserGLM_restdata(subjNumStr{i},gsr,nproc);
        GlasserGLM_miniblock_betaseries(subjNumStr{i},gsr,nproc);
    end
end

%% 
% Runs subject-wise 64k surface GLM for task and rest
execute = 0;
if execute == 1
    nproc = 12; 
    gsr = 0;
    for (i=1:32)
        disp(['Running 64k surface glm (miniblock) on subject ' num2str(i)])
        SurfaceGLM_miniblock_betaseries(subjNumStr{i},gsr,nproc);
        surfaceGLM_rest.m(subjNumStr{i},gsr,nproc);
    end
end
%%
% Save binary task timing files (with HRF convolution) (for task FC windows) to CSVs
execute = 0;
if execute == 1
    for (i=1:32)
        gsr = 0;
        
        % First run for 64 unique tasks
        X = loadStimFiles_by64TaskSet_mc16(subjNumStr{i},gsr);
        taskStims_HRFBinary = X.taskStims_HRF > 0.5;
        outputfile = ['/projects2/ModalityControl2/data/stimfiles_unique64taskset/' subjNumStr{i} '_64taskset_hrflagged_stims.csv'];
        csvwrite(outputfile, taskStims_HRFBinary);

        % Next run for 128 unique miniblocks (beta series)
        X = loadStimFiles_byMiniblockV3(subjNumStr{i},gsr);
        taskStims_HRFBinary = X.taskStims_HRF > 0.5;
        outputfile = ['/projects2/ModalityControl2/data/stimfiles_miniblockAll_v3/' subjNumStr{i} '_Miniblock_hrflagged_stims.csv'];
        csvwrite(outputfile, taskStims_HRFBinary);

    end
end




