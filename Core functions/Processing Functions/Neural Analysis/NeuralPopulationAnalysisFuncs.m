classdef NeuralPopulationAnalysisFuncs
    %NEURALPOPULATIONANALYSISFUNCS
    % class of functions to do neural population analysis this will include
    % functions that handle PCA, Classfiers, Neural representational
    % geometry
    %   Detailed explanation goes here

    properties
        Property1
    end
    properties (Access=private)
        ManData=ManipulateData;
        TrialFunc=TrialFuncs;
        FigParams=fig_params;
        ArtSimFunc=ArtificialSimFuncs;
    end

    methods
        function  obj=NeuralPopulationAnalysisFuncs(varargin)
            global AnalysisOpts
            if nargin~=0 % initialize vars
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
        end
        function  obj=ParseParams(obj,InputArgs)
            %Process optional inputs
            if mod(length(InputArgs), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(InputArgs)
                try
                    obj.(InputArgs{i}) = InputArgs{i+1};
                catch
                    error('Couldn''t set option ''%s''.', InputArgs{2*i-1});
                end
            end
        end

        function [features, Observed, Shuffled] = SVMClassifier_Binary(Data,varargin)
            % Camden MacDowell 2019
            %
            % See MATLAB documentation for fitcsvm for more details.
            % Citation for decoding neural data using classifiers:
            % Pereira, Mitchell, Botvinick, 2008 NeuroImage. Machine Learning
            % Classifiers and fMRI: a tutorial overview
            %
            % @SYNPOSIS
            % Performs non-linear binary classification on matrix Data according to the
            % class identities in the last column of Data (see @Data). Automatically
            % tunes classifier hyperparameters kernalscale and boxcontrain using pre-built
            % baysian optimization (baysopt fnc). Optionally performs feature selection
            % (currently only supports ANOVA feature selection).
            %
            % @INPUTS
            % Data: observations x features matrix. (i.e. trials x pixels for imaging data).
            % Final column contains the response variable with the class identify of
            % each row.
            %
            % @varargin: see opts structure below
            %
            % @Optimized Hyperparameters:
            % @Box contraint: is the regularization parameter, where
            % increase = hard margin and thus high cost on missclassified points.
            %
            % @Kernal sclae is 1/gamma. This is the radius of influence of the sample
            % selected by the model as support vectors. If small then you end up with
            % lots of support vectors that really only apply to those specific points
            % (low radius of influence). If large than the opposite and the model is to
            % 'rigid'/contrained to capture the shape fo the data.
            %
            % @Kernel
            % RBF: sn ~Gaussian kernel. Good for representing highly varying terrains.

            %% Set optional paramters
            opts.holdout = 0.25; %fraction hold out data for validation
            opts.nshuf = 1000; %Number of shuffles for shuffle test
            opts.featureselect  = 'anova'; %Options: 'none', 'anova':
            opts.numberfeatures = 50; %Number of features to select
            opts.pca = 1;  %First perform dimensionality reduction. Uses #PCs to explain 99% of variance.
            opts.solver = 1; %SVM solver to use: 1=SMO,2=ISDA,3=L1QP

            %Process optional inputs
            if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(varargin)
                try
                    opts.(varargin{i}) = varargin{i+1};
                catch
                    error('Couldn''t set option ''%s''.', varargin{2*i-1});
                end
            end

            %%
            %set appropriate solver;
            solver = {'SMO','ISDA','L1QP'};

            Observed = struct(); %Observed results structure
            Shuffled = struct(); %Shuffled results structure
            predictors = Data(:, 1:end-1);
            response = Data(:,end);

            %Hold out data for validation. This is not used in tuning of classifier
            cvp = cvpartition(response, 'Holdout', opts.holdout);
            trainingPredictors = predictors(cvp.training, :);
            trainingResponse = response(cvp.training, :);
            validationPredictors = predictors(cvp.test, :);
            validationResponse = response(cvp.test, :);

            %preform pca on the testing data, project validation data into pca space
            if opts.pca
                [coef, score, ~, ~, explain, mu] = pca(trainingPredictors);
                trainingPredictors = score(:,1:find(cumsum(explain)>99,1));
                validationPredictors = (coef(:,1:find(cumsum(explain)>99,1))'*(validationPredictors-repmat(mu,size(validationPredictors,1),1))')';
            end

            %Optional feature selection on training data;
            [features] = FeatureSelection(trainingPredictors,trainingResponse,opts.featureselect,opts.numberfeatures);
            trainingPredictors(:,~features)=[];
            trainingResponse(~features) = [];
            validationPredictors(:,~features)=[];
            validationResponse(~features) = [];

            % Tune classifier hyperparameters using cross validation: Only use the training data for this.
            c = cvpartition(trainingResponse, 'kfold', 10);
            optstune = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
                'AcquisitionFunctionName','expected-improvement-plus','MaxObjectiveEvaluations',100);
            svmmod = fitcsvm(trainingPredictors,trainingResponse,'KernelFunction','rbf',...
                'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',optstune,'Solver',solver{opts.solver});

            % Build a classifier using tuned parameters and all training data
            classificationSVM = fitcsvm(trainingPredictors,trainingResponse,'KernelFunction','rbf','BoxConstraint',...
                svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
                'KernelScale',svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelScale,'Solver',solver{opts.solver});

            % Create the result struct with predict function
            svmPredictFcn = @(x) predict(classificationSVM, x);
            validationPredictFcn = @(x) svmPredictFcn(x);

            % Compute validation predictions
            [validationPredictions, validationScores] = validationPredictFcn(validationPredictors);
            Outstats = [];
            for i = 1:numel(unique(response))
                [Outstats(i).X,Outstats(i).Y,Outstats(i).T,Outstats(i).AUC,Outstats(i).OPTROCPT,Outstats(i).SUBY] = perfcurve(validationResponse,validationScores(:,i),i);
            end

            % Compute validation accuracy
            correctPredictions = (validationPredictions == validationResponse);
            isMissing = isnan(validationResponse);
            correctPredictions = correctPredictions(~isMissing);

            %Save off all desired information
            Observed.AUC = mean([Outstats(:).AUC]);
            Observed.X = [cat(1,Outstats(:).X),cat(1,Outstats(:).Y)];
            Observed.Y = [cat(1,Outstats(:).Y),cat(1,Outstats(:).Y)];
            Observed.T = [cat(1,Outstats(:).T),cat(1,Outstats(:).T)];
            Observed.SubY = [cat(1,Outstats(:).X),cat(1,Outstats(:).Y)];
            Observed.Optrocpt = mean([Outstats(:).AUC]);
            Observed.Accurary = sum(correctPredictions)/length(correctPredictions);
            Observed.Predictions = validationPredictions;
            Observed.CorrectResponse = validationResponse;
            Observed.Scores = validationScores;
            Observed.Classifier = classificationSVM;
            Observed.trainingResponse = trainingResponse;
            Observed.trainingPredictors = trainingPredictors;
            Observed.validationPredictors = validationPredictors;
            Observed.validationResponse = validationResponse;


            % Now randomly shuffle the held out labels and use the same classifier
            % Compute validation predictions
            for shuf = 1:opts.nshuf
                validationResponse_shuf = validationResponse(randperm(numel(validationResponse)));
                [validationPredictions, validationScores] = validationPredictFcn(validationPredictors);
                Outstats = [];

                for i = 1:numel(unique(response))
                    [Outstats(i).X,Outstats(i).Y,Outstats(i).T,Outstats(i).AUC,Outstats(i).OPTROCPT,Outstats(i).SUBY] =...
                        perfcurve( validationResponse_shuf,validationScores(:,i),i);
                end

                % Compute validation accuracy
                correctPredictions = (validationPredictions ==  validationResponse_shuf);
                isMissing = isnan( validationResponse_shuf);
                correctPredictions = correctPredictions(~isMissing);
                validationAccuracy = sum(correctPredictions)/length(correctPredictions);

                Shuffled(shuf).AUC = mean([Outstats(:).AUC]);
                Shuffled(shuf).X = [cat(1,Outstats(:).X),cat(1,Outstats(:).Y)];
                Shuffled(shuf).Y = [cat(1,Outstats(:).Y),cat(1,Outstats(:).Y)];
                Shuffled(shuf).T = [cat(1,Outstats(:).T),cat(1,Outstats(:).T)];
                Shuffled(shuf).SubY = [cat(1,Outstats(:).X),cat(1,Outstats(:).Y)];
                Shuffled(shuf).Optrocpt = mean([Outstats(:).AUC]);
                Shuffled(shuf).Accurary = validationAccuracy;
                Shuffled(shuf).Predictions = validationPredictions;
                Shuffled(shuf).CorrectResponse =  validationResponse_shuf;
                Shuffled(shuf).Scores = validationScores;
            end
        end
    end
end

