function output = ActFlow_ITE_DecodePerformance_LogRegression(subj)
% Taku Ito
% 08/31/2016
% Runs region-to-region ActFlow using the Glasser Parcels
    
    % Load in basic parameters and atlases
    datadir = ['/work/ti61/projects/ModalityControl/data/'];
    pathtoparcel = [datadir 'GlasserKKPartition/'];
    % Define parcel csv array file (in vertex space)
    dlabels = csvread([pathtoparcel 'Q1-Q6_RelatedParcellation210.LR.CorticalAreas_dil_Colors.32k_fs_LR.csv']);
    dilationMatrix = csvread([pathtoparcel 'ParcelLabels/GlasserParcelsAll_Dilated.csv']);

    % Define paths to data
    pathtorest = [datadir 'results/rest_glm_64k/'];
    pathtobetas = [datadir 'results/glmMiniblockBetaSeries/'];
    % Define paths to output directory

    % Load resting state data
    disp('Loading in rest data')
    restdata = csvread([pathtorest subj '_rest_nuisanceResids_64kSurface.csv']);
    betas = csvread([pathtobetas subj '_miniblock_taskBetas_Surface64k.csv']);
    % First 17 betas are nuisance regressors, so remove those
    betas = betas(:,18:end);
    ind_nan = isnan(betas);
    betas(ind_nan) = 0;


    % Compute rest FC via principal component regression
    disp('Starting computation of FC rest matrix')
    fcmat = VertexToVertexPCMultRegFC(restdata, dlabels, dilationMatrix);

    % Force NaNs to be 0, since they shouldn't influence activity flow anyway
    ind_nan = isnan(fcmat);
    fcmat(ind_nan) = 0;

    acc_array = csvread(['/work/ti61/projects/ModalityControl/data/behavresults/' subj '_acc_by_mb.txt']);


    % Compute Region to Region ActFlow
    [ITEMatrix, BehavioralAccMatrix] = RegionToRegionActFlow(betas, fcmat, dlabels, subj, acc_array);

    % Write out ITE matrix to output directory
    ruledims = {'logic','sensory','motor'};
    for i=1:length(ruledims)
        outputdir = [datadir 'results/RegionToRegionActFlowGlasserParcelsNM3_Spearman/'];
        filename = [subj '_' ruledims{i} '_RegionToRegionActFlowGlasserParcels.csv'];
        csvwrite([outputdir filename], ITEMatrix(:,:,i));
    end

    filename2 = [subj '_RegionToRegionActFlowGlasserParcels_BehavioralAcc.csv'];
    csvwrite([outputdir filename2], BehavioralAccMatrix);

    output = 0;

end



function [RSADiffMatrix, BehavioralAccMatrix]  = RegionToRegionActFlow(betas, fcmat, dlabels, subj, acc_array)
% Runs region-to-region ActFlow at the vertex level
% Inputs
%   betas = vertex X miniblock matrix for whole brain
%   fcmat = fcmat -- should be 64k x 64k
%   dlabels = 64kx1 array indicating a parcel per vertex

    datadir = ['/work/ti61/projects/ModalityControl/data/'];
    nmbs = size(betas,2) % Get number of miniblocks (or number of samples)
    % Define number of parcels    
    nParcels = 360;
    nVertices = size(betas,1);
    % Instantiate empty region-to-region RSA representational quality matrix
    RSADiffMatrix = zeros(nParcels,nParcels,3); % region X region X rule dimension
    BehavioralAccMatrix = zeros(nParcels,nParcels);
    
    % Start ActFlow for loop
    for roi1=1:nParcels
        disp(['Running All ActFlow from region ' num2str(roi1)])
        % Find indices for roi1
        roi1_ind = dlabels == roi1;
        % get betas for roi1 and transpose to make miniblocks x vertices
        betas1 = betas(roi1_ind,:)';

        for roi2=1:nParcels
            % Skip if roi1 and 2 are the same
            if roi1==roi2 
                continue
            end

            % Find indices for roi2
            roi2_ind = dlabels == roi2;
            % get betas for roi2 and transpose to make miniblocks x vertices
            betas2 = betas(roi2_ind,:)';
            
            w12 = fcmat(roi1_ind, roi2_ind);
            w21 = fcmat(roi2_ind, roi1_ind);
    
            % Compute region to region actflow using dot product
            actflow12 = betas1 * w12; % output is in miniblock x vertex space
            actflow21 = betas2 * w21; % output is in miniblock x vertex space

            %% Run predicted-to-actual similarity analysis for all 3 rule types
            ruledims = {'logic','sensory','motor'};
            ITE12_allruledims = zeros(nmbs, length(ruledims));
            ITE21_allruledims = zeros(nmbs, length(ruledims));
            rulecount = 1;
            for (i=1:length(ruledims))
                disp(['Running Region-To-Region ActFlow for ' ruledims{i}])
                % Load in miniblock to rule assignments
                rule_ind = csvread([datadir 'CPROTaskIdentifiers/' subj '_' ruledims{i} '_miniblocksPerRule_MatlabIndices.csv']);

                % Predicted-to-actual similarity analysis
                % Run RSA CV of these two regions' actflow representations
                [ITE12] = RSACV(actflow12, betas2, rule_ind); 
                [ITE21] = RSACV(actflow21, betas1, rule_ind);

                ITE12_allruledims(:,rulecount) = ITE12;
                ITE21_allruledims(:,rulecount) = ITE21;
                rulecount = rulecount + 1;
            end

            % Classify performance with all three information transfer estimates
            % Use logistic regression model to decode (information transfer estimate (independent var) and performance (dependent var))
            % Instantiate logistic func
            logistic = @(x) 1./(1+exp(-(x)));
            %
            accuracy12 = zeros(nmbs,1);
            accuracy21 = zeros(nmbs,1);
            for mb=1:nmbs
                indices = 1:nmbs;
                indices = indices == indices; % make binary
                test_ind = mb;
                train_ind = indices;
                train_ind(test_ind) = 0; % exclude test index from train indices

                % Train model 12
                trainset12 = ITE12_allruledims(train_ind,:);
                testset12 = ITE12_allruledims(test_ind,:);
                [coefs] = glmfit(trainset12, acc_array(train_ind),'binomial','link','logit');
                % don't include the first beta, i.e., B0, since that biases according to the trainset
                % beta2: coef for logic; beta3: coef for sensory information; beta4: coef for motor information
                % testset12(1:3): information transfer estimates for each of the rule domains for this mb
                y_predict = logistic(coefs(2)*testset12(1) + coefs(3)*testset12(2) + coefs(4)*testset12(3));
                accuracy12(mb) = y_predict;

                %Train model 21
                trainset21 = ITE21_allruledims(train_ind,:);
                testset21 = ITE21_allruledims(test_ind,:);
                [coefs] = glmfit(trainset21, acc_array(train_ind),'binomial','link','logit');
                % predict (with logistic)
                logistic = @(x) 1./(1+exp(-(x)));
                % don't include the first beta, i.e., B0, since that biases according to the trainset
                % beta2: coef for logic; beta3: coef for sensory information; beta4: coef for motor information
                % testset12(1:3): information transfer estimates for each of the rule domains for this mb
                y_predict = logistic(coefs(2)*testset21(1) + coefs(3)*testset21(2) + coefs(4)*testset21(3));
                accuracy21(mb) = y_predict;
            end

            % Binarize prediction
            predictions12 = accuracy12 > .5;
            predictions21 = accuracy21 > .5;
            % Binarize accuracies per miniblock (if < 50% make 0; if > 50% make 1);
            actual_accuracy = acc_array > .5;

            % Compute decoding accuracy
            decodingaccuracy12 = mean(predictions12==actual_accuracy);
            decodingaccuracy21 = mean(predictions21==actual_accuracy);

            % Store ITEs into matrix
            for i=1:length(ruledims)
                RSADiffMatrix(roi1,roi2,i) = mean(ITE12_allruledims(:,i));
                RSADiffMatrix(roi2,roi1,i) = mean(ITE21_allruledims(:,i));
            end

            % Store behavioral performance decodings to a regionXregion matrix
            BehavioralAccMatrix(roi1,roi2) = decodingaccuracy12;
            BehavioralAccMatrix(roi2,roi1) = decodingaccuracy21;
        end
        
   end

end



function [ITE] = RSACV(actflowmat, realmat, rule_ind)
% Compute RSA similarity of predicted activity patterns to matched real activity patterns
% Also compute RSA similarity of predictec activity patterns to mismatched real activity patterns

    ncvs = size(rule_ind,1);
    nrules = size(rule_ind,2);
    % Make sure each CV fold is randomized
    for i=1:nrules
        idx=randperm(ncvs);
        rule_ind(:,i) = rule_ind(idx,i);
    end
    % Perform feature-wise normalization for each column, across samples
    %actflowmat = zscore(actflowmat);
    %realmat = zscore(realmat);

    % Spatially demean
    actflowmat = bsxfun(@minus, actflowmat, mean(actflowmat, 2));
    realmat = bsxfun(@minus, realmat, mean(realmat,2));

    % Number of vertices for this particular ROI
    nfeatures = size(actflowmat,2);
    

    % Instantiate empty vector array for matches and mismatches per cv
    matches = zeros(ncvs,1);
    mismatches = zeros(ncvs,1);
    % Trial variable for keeping track of accuracies
    nmbs = size(actflowmat,1); % number of miniblocks (i.e., samples)
    ITE = zeros(nmbs,1);
    for cv=1:ncvs
        % Find the indices for the prototyped blocks
        ind_prototype = ((1:ncvs) ~= cv);
        
        % Instantiate empty matrix for all training blocks
        prototypes = zeros(nrules, nfeatures); % number of rules x number of vertices in this ROI 
        % Construct prototypes for each rule 
        for rule=1:nrules
            % Find the miniblocks for this particular rule
            rulemb = rule_ind(ind_prototype,rule);
            prototypes(rule,:) = mean(realmat(rulemb,:), 1); % Average across miniblocks for this rule type
        end


        % Correlate each rule type of this current CV to each of the prototypes
        rsamat = zeros(nrules,nrules);
        for rule1=1:nrules
            block_ind = rule_ind(cv, rule1);
            count = 1;
            for rule2=1:nrules
                % Get similarity measure
                r = corr(actflowmat(block_ind,:)',prototypes(rule2,:)');
                r = atanh(r);
                %r = r(1,2); % Get off diagonal for real pearson r
                
                if rule1==rule2
                    correct = r;
                else
                    errortrials(count) = r;
                    count = count + 1;
                end

                rsamat(rule1,rule2) = r;
            end
            % assign this information transfer estimate back to its original miniblock number
            ITE(rule_ind(cv,rule1)) = correct - mean(errortrials);
            clear errortrials;
        end


    end


end



function fcmat = VertexToVertexPCMultRegFC(activityMatrix, dlabels, dilationMatrix)
% Runs principal componenet regression on vertex-wise data
% Runs according to the Glasser parcel set
% Runs a PCA for every vertex outside a dilated parcel and computes multiple regression from all vertices outside that dilation to a single vertex in that parcel
% Inputs
%   activityMatrix - activity matrix to derive FC estimates from
%   pathtoparcel - path to the Glasser parcellation (vertex-wise)
%   dilationMatrix - a vertex X parcel matrix, where for each column (parcel), every row indicates the dilated parcel from which to not compute FC from

    % Define number of parcels
    nParcels = 360;

    % Define number of principal components to downsample to
    numComponents = 500;
    % Define number of vertices in original surface space
    nVertices = 64984;
    % Define empty vertex-to-vetex FC matrix
    interactionMatrix = zeros(nVertices,nVertices);

    for parcel=1:nParcels

        disp(['Running multiple regression on all vertices in parcel ' num2str(parcel)])
        tic
        % Select vertices to be held out and to be predicted (within a parcel)
        targetVertices = dlabels==parcel;
        % Select vertices to be the source -- all vertices outside the dilated parcel
        sourceVertices = dilationMatrix(:,parcel)~=parcel;

        disp('Running PCA')
        % Run PCA
        [PCALoadings, PCAScores] = pca(activityMatrix(sourceVertices,:)', 'NumComponents', numComponents);

        disp('Beginning Principal Component Regression')
        for vertex=find(targetVertices)'
            % Pad matrix with 1s for B0s
            regmat = [ones(length(activityMatrix(vertex,:)),1), PCAScores(:,1:numComponents)]; 
            pcabetas = regress(activityMatrix(vertex,:)', regmat);
            betaPCR = PCALoadings*pcabetas(2:end);

            interactionMatrix(sourceVertices, vertex) = betaPCR;
        end
        toc

    end

    fcmat = interactionMatrix;

end

