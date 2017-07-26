function output = loadStimFiles_cprorules(subj, gsr)
% Taku Ito
% 3/12/15
%
% This script imports the stim files for IndivRITL modality control into a regressor matrix including wm/ventricle/wholebrain timeseries, 12 motion parameters, and 8 stimulus time series
% Regressing out WM, Ventricles (and corresponding derivatives), motion parameters, and auditory versus visual rule contrasts (i.e., constant and high pitch versus vertical and red)
% 
% Parameters: 
%   subj (must be input with single quotations, i.e., as a string!)
%   gsr - 1 if you would like to include a GSR regressor in the matrix, 0 if not   

    % Set up basic parameters
    basedir = ['/projects/ModalityControl/data/indivritl/' subj];
    datadir = [basedir '/MNINonLinear/Results'];
    analysisdir = [basedir '/analysis'];
    totalRestTRs = 1070; 
    trsPerRun = 581;
    numTaskStims = 12; % Sensory task rules (4)
    numMotionParams = 12; % HCP Pipe outputs 12
    trLength = .785; 
   
    %%
    % Need to create the derivative time series for ventricle, white matter, and whole brain signal
    disp(['Creating derivative time series for ventricle, white matter, and whole brain signal for subject ' subj])
    eval(['!1d_tool.py -overwrite -infile ' analysisdir '/' subj '_WM_timeseries_task.1D -derivative -write ' analysisdir '/' subj '_WM_timeseries_deriv_task.1D'])
    eval(['!1d_tool.py -overwrite -infile ' analysisdir '/' subj '_ventricles_timeseries_task.1D -derivative -write ' analysisdir '/' subj '_ventricles_timeseries_deriv_task.1D'])
    eval(['!1d_tool.py -overwrite -infile ' analysisdir '/' subj '_wholebrainsignal_timeseries_task.1D -derivative -write ' analysisdir '/' subj '_wholebrainsignal_timeseries_deriv_task.1D'])
    eval(['!1d_tool.py -overwrite -infile ' analysisdir '/' subj '_WM_timeseries_rest.1D -derivative -write ' analysisdir '/' subj '_WM_timeseries_deriv_rest.1D'])
    eval(['!1d_tool.py -overwrite -infile ' analysisdir '/' subj '_ventricles_timeseries_rest.1D -derivative -write ' analysisdir '/' subj '_ventricles_timeseries_deriv_rest.1D'])
    eval(['!1d_tool.py -overwrite -infile ' analysisdir '/' subj '_wholebrainsignal_timeseries_rest.1D -derivative -write ' analysisdir '/' subj '_wholebrainsignal_timeseries_deriv_rest.1D'])


    %% 
    % Import Rest noise parameters
    disp(['Importing wm, ventricle and global brain time series into MATLAB for subj ' subj])
    % First 2 columns will be white matter
    timeseriesRegressors_rest(:,1) = importdata([analysisdir '/' subj '_WM_timeseries_rest.1D']);
    timeseriesRegressors_rest(:,2) = importdata([analysisdir '/' subj '_WM_timeseries_deriv_rest.1D']);
    % Columns 3 and 4 will be ventricles
    timeseriesRegressors_rest(:,3) = importdata([analysisdir '/' subj '_ventricles_timeseries_rest.1D']);
    timeseriesRegressors_rest(:,4) = importdata([analysisdir '/' subj '_ventricles_timeseries_deriv_rest.1D']);
    % Columns 5 and 6 will be whole brain signal, though we will opt to remove these out later if we do not want to perform GSR
    timeseriesRegressors_rest(:,5) = importdata([analysisdir '/' subj '_wholebrainsignal_timeseries_rest.1D']);
    timeseriesRegressors_rest(:,6) = importdata([analysisdir '/' subj '_wholebrainsignal_timeseries_deriv_rest.1D']);
    % Import rest movement regressors
    movementReg = [datadir '/Rest1/Movement_Regressors.txt'];
    movementRegressors_rest = importdata(movementReg); 
    
    %% 
    % if GSR is not selected (i.e., 0), then remove it from the timeseriesRegressor
    if gsr == 0 
        timeseriesRegressors_rest = timeseriesRegressors_rest(:,1:4);
    end


    %%
    % If we want to regress out 
    noiseRegressors = [movementRegressors_rest timeseriesRegressors_rest];

    output.noiseRegressors = noiseRegressors;
end
    

    
