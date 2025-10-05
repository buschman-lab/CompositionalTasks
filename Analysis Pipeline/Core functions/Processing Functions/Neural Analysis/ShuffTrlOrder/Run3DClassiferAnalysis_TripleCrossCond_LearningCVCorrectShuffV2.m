 function [ClassifierResults,ClassifierOpts]=Run3DClassiferAnalysis_TripleCrossCond_LearningCVCorrectShuffV2(obj,FactorData,FactorLevelComb,ClassifierOpts,varargin) % With correct shuffling procidure, runs cross temporal classifier anlaysis for 3 simultanous features in three different cross conditions during learning with cross validation
%             corrected so that the shuffling is applied before any trial range 
%             shuffles the trial order
%             running a classifier analysis for three simultaneous features in three different cross conditions during learning with cross-validation.
%             Here are the main steps and functionalities of this code:
% 
%            * Setting up variables and options based on input parameters and global variables.
%            * Determining the target field (e.g., 'StimCongruency') and cross-validation specifications.
%            * Looping through different conditions and trial ranges.
%            * Performing data preprocessing, including applying a neuron dropping algorithm.
%            * Running a classifier analysis using cross-validation and collecting the results.
%            * Calculating a shuffle distribution if specified.
%            * Applying some conditions and filtering on the data based on various criteria.
%            * Continuing with learning analysis and cross-validation for the first condition.
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            ClassifierOpts.ConcatinateMethod='SepTimePoint'; % for PCA analysis
            ClassifierOpts.MainTargetFactor={'Quadrants'};% can be Quadrants ColorMLComb
            ClassifierOpts.StimCongruencyFactorName=obj.StimCongruencyFactorName;
            TrialNumInd=strcmp(AnalysisOpts.factornames,'TrialNum');

            TargetFactorsInd=obj.GetFactornameInd(ClassifierOpts.TargetFactors{1});
            TargetFactors_2ndDInd=obj.GetFactornameInd(ClassifierOpts.TargetFactors_2ndD{1});
            TargetFactors_3ndDInd=obj.GetFactornameInd(ClassifierOpts.TargetFactors_3ndD{1});
            
            if strcmp(ClassifierOpts.TargetFactors_2ndD{1},ClassifierOpts.TargetFactors{1})
                DimMethod_2ndD=1;else;DimMethod_2ndD=2;
            end
            if strcmp(ClassifierOpts.TargetFactors_3ndD{1},ClassifierOpts.TargetFactors{1})
                DimMethod_3ndD=1;else;DimMethod_3ndD=2;
            end
            % if our factors are based on the rule then adjust the shuffleing procidure so the end factor data is based on the rule as well
            if strcmp(ClassifierOpts.TargetFactors_2ndD{1},'Rule')
                DimMethod_2ndD=3;
            end
            if strcmp(ClassifierOpts.TargetFactors_3ndD{1},'Rule')
                DimMethod_3ndD=3;
            end
            DimTxt={'','2','3'};     
            AnalysisOpts.DimTxt={'','2','3'}; 
            TargetFields=obj.DetermineStimCongruency(ClassifierOpts,1,FactorData,FactorData);
            ClassifierOpts=obj.DetermineDimCVspecs(ClassifierOpts); % determine cross validation specs 
          
            [ShuffTxt,CondSet,TrlRngSet,RepSet]=obj.ReadClassifierSpockRunConds(ClassifierOpts);
            RepSetBuff=RepSet;

            ExtraStrSaveShuffLabel=[ClassifierOpts.Name '_' AnalysisOpts.CurrentCh_AreaName '_' num2str(AnalysisOpts.DividSpockClassifier_Cond) '_ShuffleLabels'];

            % now create a mesh of time points and run classifier train and test repectively
            [ClassifierOpts,TrainTimInd,TestTimInd,nXtimePnt,TimeMatrixSize]=obj.GetTimeRangeforThisCond(ClassifierOpts);

            % apply neuron dropping algorithm
            [FactorDataBuff,IncludedNeurons,ClassifierOpts]=obj.ApplyNeuronDroppingAlgorithm4LearningAnlaysis(FactorData,ClassifierOpts,FactorLevelComb);
            
            FactorDataBuff=rmfield(FactorDataBuff,'dataMean');
            ClassifierOpts.typebuff=ClassifierOpts.type; % take a copy of type
            % get the ind of factors of test we are keeping
            ClassifierOpts.FactorInds2Keep=AnalysisOpts.FactorInds2Keep;
            FactorInds2Keep=cellfun(@(x) find(ClassifierOpts.(x).Ind),ClassifierOpts.FactorInds2Keep);
            AnalysisOpts.FactorInds2KeepInd=FactorInds2Keep; % save the indices
            ClassifierResults=[];ClassifierResults_Shuffled=[];

            % now get the classifer performance per crosstime point
            for nCond=1:length(CondSet) % loop on conditions
                Cond=CondSet(nCond);
                ClassifierOpts.CurrCond=Cond;
                % if there is SeqHist you should consider them first
                if sum(ClassifierOpts.StimCongruency==9) || sum(ClassifierOpts.StimCongruency==10) || sum(ClassifierOpts.StimCongruency==13) || sum(ClassifierOpts.StimCongruency==14) || sum(ClassifierOpts.StimCongruency==15) || sum(ClassifierOpts.StimCongruency==17)
                    obj.SeqHistCond=ClassifierOpts.SeqHistCond{Cond};end
                if iscell(ClassifierOpts.type);type=ClassifierOpts.type{Cond};else;type=ClassifierOpts.type;end
                
                % create a shuffle distribution for the subsequenct shufflings for each dimension
               if obj.CalShuff==1
                    ClassifierOpts=CreateShuffleDist4ClassifierCorrectShuffV5(obj,FactorDataBuff,ClassifierOpts,FactorLevelComb,1,Cond,1);
                    ClassifierOpts=CreateShuffleDist4ClassifierCorrectShuffV5(obj,FactorDataBuff,ClassifierOpts,FactorLevelComb,2,Cond,1);
                    ClassifierOpts=CreateShuffleDist4ClassifierCorrectShuffV5(obj,FactorDataBuff,ClassifierOpts,FactorLevelComb,3,Cond,1);
               end

                if ~AnalysisOpts.GetOnlyShuffLabelsClassifier % we are only getting classifier label then skip the rest
                    %% run the code
                    % create learning trials for this condition
                    [TrainTrlRange,TestTrlRange,TrainTrlInd,TestTrlInd]=obj.GetTrialRangeforThisCond(ClassifierOpts,Cond);

                    % if we have the whole range of TrlRngSet then choose the one for this condition
                    if iscell(TrlRngSet);ThisTrlRngSet=TrlRngSet{Cond};else;ThisTrlRngSet=TrlRngSet;end
                   if AnalysisOpts.UseRep4Cluster==1;ThisTrlRngSet=unique(ThisTrlRngSet);end

                    % determine what condition we are computing the data for
                    if obj.CalShuff==1 % if we are running shuffle
                        if RepSetBuff==ClassifierOpts.NrepShufperFold+1 % this is observed
                            ShuffObservedRound=0;obj.ShuffObservedRoundObj=0;
                            repshuffSet=1;
                            fprintf(2,'\n*****Running Observed reps')
                        else   % this is shuffle
                            ShuffObservedRound=1;obj.ShuffObservedRoundObj=1;
                            if AnalysisOpts.UseRep4Cluster==1
                                repshuffSet=RepSetBuff;
                            else
                                repshuffSet=RepSet;
                            end
                            fprintf(2,'\n*****Running Shuffle reps')
                        end
                        RepSet=1:ClassifierOpts.NrepShuf;
                    else % this is for mean and STD
                        fprintf(2,'\n*****Running Mean STD reps')
                        ShuffObservedRound=0;obj.ShuffObservedRoundObj=0;
                        repshuffSet=1;
                    end

                    for repshuff=1:length(repshuffSet) % shuffle classifier repetition
                        fprintf(2,'\nRunning classifer, Rep %i',repshuff)

                        if  obj.CalShuff & ShuffObservedRound
                            if length(repshuffSet)>1;error('This code does not work here');end
                            NewFactorDataDim1Buff=obj.ManData.SwapShuffleTrialsFactorData(ClassifierOpts,FactorDataBuff,1,1,FactorLevelComb,TargetFactorsInd,2,TargetFactorsInd,Cond,1);
                            NewFactorDataDim2Buff=obj.ManData.SwapShuffleTrialsFactorData(ClassifierOpts,FactorDataBuff,1,2,FactorLevelComb,TargetFactors_2ndDInd,2,TargetFactorsInd,Cond,1);
                            NewFactorDataDim3Buff=obj.ManData.SwapShuffleTrialsFactorData(ClassifierOpts,FactorDataBuff,1,3,FactorLevelComb,TargetFactors_3ndDInd,2,TargetFactorsInd,Cond,1);
                        else
                            NewFactorDataDim1Buff=FactorDataBuff;
                            NewFactorDataDim2Buff=FactorDataBuff;
                            NewFactorDataDim3Buff=FactorDataBuff;
                        end

                        % if we are skipping applying SeqHist to traning data
                        NewFactorDataDim1Buff=obj.SkipSeqHistTrain(NewFactorDataDim1Buff,ClassifierOpts,Cond);
                        NewFactorDataDim2Buff=obj.SkipSeqHistTrain(NewFactorDataDim2Buff,ClassifierOpts,Cond);
                        NewFactorDataDim3Buff=obj.SkipSeqHistTrain(NewFactorDataDim3Buff,ClassifierOpts,Cond);

                        for nTrlRng=1:length(ThisTrlRngSet)
                            TrlRng=ThisTrlRngSet(nTrlRng);

                            ThisTrainTrialRange=TrainTrlRange{TrainTrlInd(TrlRng)};
                            ThisTestTrialRange=TestTrlRange{TestTrlInd(TrlRng)};
                           
                            % now limit factor data trials to these trial range we want generate factor data for training
                            NewFactorDataDim1=obj.LimitTrialsBasedonFactor(NewFactorDataDim1Buff,'TrialNum',ThisTrainTrialRange);
                            NewFactorDataDim2=obj.LimitTrialsBasedonFactor(NewFactorDataDim2Buff,'TrialNum',ThisTrainTrialRange);
                            NewFactorDataDim3=obj.LimitTrialsBasedonFactor(NewFactorDataDim3Buff,'TrialNum',ThisTrainTrialRange);

                            % generate factordata for test
                            FactorDataTest=obj.LimitTrialsBasedonFactor(FactorDataBuff,'TrialNum',ThisTestTrialRange);

                            % if we are limitign based on block performance then
                            if isfield(ClassifierOpts,'LimitFromSwitchPerf')
                                % generate factor data for test
                                NewFactorDataDim1=obj.LimitTrialsBasedonFactor(NewFactorDataDim1,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(1),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{1});
                                NewFactorDataDim2=obj.LimitTrialsBasedonFactor(NewFactorDataDim2,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(1),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{1});
                                NewFactorDataDim3=obj.LimitTrialsBasedonFactor(NewFactorDataDim3,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(1),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{1});

                                % generate factordata for training
                                FactorDataTest=obj.LimitTrialsBasedonFactor(FactorDataTest,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(2),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{2});
                            end
                           
                            % if we are limiting based on an arbitrary factor
                            if isfield(ClassifierOpts,'LimitTrainTrialsBasedonFactor') | isfield(ClassifierOpts,'LimitTestTrialsBasedonFactor')
                                 % generate factor data for test
                                [NewFactorDataDim1,NewFactorDataDim2,NewFactorDataDim3]=obj.LimitAllDimsTrialsBasedonFactor(ClassifierOpts,NewFactorDataDim1,NewFactorDataDim2,NewFactorDataDim3,Cond,1);
                                % generate factordata for training
                                 FactorDataTest=obj.LimitAllDimsTrialsBasedonFactor(ClassifierOpts,FactorDataTest,FactorDataTest,FactorDataTest,Cond,2);
                            end
                            
                            for rep=1:length(RepSet)       % loop on repetition per condition
                                %  CVfoldResults=[]; % initialize this rep
                                FoldTimeTic=tic;
                                for cvf=1:ClassifierOpts.nCVfold(Cond) % loop on cross validation fold

                                    fprintf('\nTrain/Testing Train1:%s Test1:%s | Train2:%s Test2:%s | Train3:%s Test3:%s Cond:%i RepSuff:%i TrlRng:%i Rep:%i CV:%i',obj.ManData.ConvMat2Char(ClassifierOpts.TrainCond{Cond}),...
                                        obj.ManData.ConvMat2Char(ClassifierOpts.TestCond{Cond}),obj.ManData.ConvMat2Char(ClassifierOpts.TrainCond2{Cond}),obj.ManData.ConvMat2Char(ClassifierOpts.TestCond2{Cond}),...
                                        obj.ManData.ConvMat2Char(ClassifierOpts.TrainCond3{Cond}),obj.ManData.ConvMat2Char(ClassifierOpts.TestCond3{Cond}),Cond,repshuff,TrlRng,rep,cvf);

                                    ClassifierOpts.CurrCVfold=cvf;

                                    PrepTimeTic=tic;
                                    isok=0;
                                    %% sample train and test trials
                                    while isok==0 % make sure we have enough trials for both categories in both dimensions
                                        % grab a second set of random data we will use the test data of this set of the second dimension but keep the training data from first set
                                        if strcmp(type,'TripleCrossCond')
                                            try
                                                % grab random set of data for Dim 1
                                                [predictors ,response , TrainDataAllFactors,TestDataAllFactors,TrainStimInds,~,predictorsSpkCnt,TrainDataAllFactors3D,TestDataAllFactors3D]=...
                                                    obj.GrabThisRepClassifierData(NewFactorDataDim1,ClassifierOpts,ClassifierOpts.TrainCond{Cond},...
                                                    ClassifierOpts.TestCond{Cond},FactorLevelComb,Cond,FactorDataTest,ThisTrainTrialRange,ThisTestTrialRange,'LookatDim2',0);                                             
                                               
                                                % grab random set of data for Dim 2 
                                                [predictors2,response2,Train2DataAllFactors, Test2DataAllFactors,Train2StimInds,~,predictors2SpkCnt,Train2DataAllFactors3D,Test2DataAllFactors3D]=...
                                                    obj.GrabThisRepClassifierData(NewFactorDataDim2,ClassifierOpts,ClassifierOpts.TrainCond2{Cond},...
                                                    ClassifierOpts.TestCond2{Cond},FactorLevelComb,Cond,FactorDataTest,ThisTrainTrialRange,ThisTestTrialRange,'LookatDim2',2);
                                               
                                                % grab random set of data for Dim 3
                                                [predictors3,response3,Train3DataAllFactors, Test3DataAllFactors,Train3StimInds,~,predictors3SpkCnt,Train3DataAllFactors3D,Test3DataAllFactors3D]=...
                                                    obj.GrabThisRepClassifierData(NewFactorDataDim3,ClassifierOpts,ClassifierOpts.TrainCond3{Cond},...
                                                    ClassifierOpts.TestCond3{Cond},FactorLevelComb,Cond,FactorDataTest,ThisTrainTrialRange,ThisTestTrialRange,'LookatDim2',3);
                                              
                                                % euqalize trials across conditions
                                                [predictors,predictors2,predictors3,response,response2,response3,TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,...
                                                    TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors,...
                                                    predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt,ClassifierOpts]=obj.EqualizeTrialsAcrossConds(ClassifierOpts,Cond,predictors,predictors2,predictors3,...
                                                    response,response2,response3,TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,...
                                                    TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors,...
                                                    predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt);
                                            catch ME
                                                error(ME.message);
                                            end
                                            if cvf==1
                                                % copy some of the data into ClassifierOpts for future use
                                                ClassifierOpts=obj.ManData.CopyVars2Struct(ClassifierOpts,'TrainStimInds',TrainStimInds,'Train2StimInds',Train2StimInds,'Train3StimInds',Train3StimInds,...
                                                    'TrainDataAllFactors',TrainDataAllFactors,'TestDataAllFactors',TestDataAllFactors,...
                                                    'Train2DataAllFactors',Train2DataAllFactors,'Test2DataAllFactors',Test2DataAllFactors,...
                                                    'Train3DataAllFactors',Train3DataAllFactors,'Test3DataAllFactors',Test3DataAllFactors);
                                            end
                                        end
                                        % if we are changing the distribution of train
                                        % and test trials for example for  congruency test train on congruent  test on incongruent or vice versa)
                                        [TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,TestDataAllFactors,Test2DataAllFactors,...
                                            Test3DataAllFactors,predictors,predictors2,predictors3,response,response2,response3,predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt]=...
                                            obj.LimitLearningTrainTestTrialsWithFactor(ClassifierOpts,Cond,...
                                            TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors,...
                                            predictors,predictors2,predictors3,response,response2,response3,1,predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt);

                                        [TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,TestDataAllFactors,Test2DataAllFactors,...
                                            Test3DataAllFactors,predictors,predictors2,predictors3,response,response2,response3,predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt]=...
                                            obj.LimitLearningTrainTestTrialsWithFactor(ClassifierOpts,Cond,...
                                            TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors,...
                                            predictors,predictors2,predictors3,response,response2,response3,2,predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt);

                                        [TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,TestDataAllFactors,Test2DataAllFactors,...
                                            Test3DataAllFactors,predictors,predictors2,predictors3,response,response2,response3,predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt]=...
                                            obj.LimitLearningTrainTestTrialsWithFactor(ClassifierOpts,Cond,...
                                            TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors,...
                                            predictors,predictors2,predictors3,response,response2,response3,3,predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt);
                                    
                                        % adjust the labels
                                        [DataLabels_Train,DataLabels_Train_2ndD,DataLabels_Train_3ndD, DataLabels_Test,DataLabels_Test_2ndD,...
                                            DataLabels_Test_3ndD, Levels_Train,Levels_2ndD_Train,Levels_3ndD_Train, ...
                                            Levels_Test,Levels_2ndD_Test,Levels_3ndD_Test,response,response_2ndD,response_3ndD]=obj.AdjustClassifierLabels(ClassifierOpts,Cond,TrainDataAllFactors,Train2DataAllFactors,...
                                            Train3DataAllFactors,TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors,response,response2,response3,TargetFactorsInd,TargetFactors_2ndDInd,TargetFactors_3ndDInd);

                                        if ((length(unique(DataLabels_Train_2ndD))==1 | length(unique(DataLabels_Test_2ndD))==1 |... % test dim 2
                                                ~isempty(setdiff(unique(DataLabels_Train_2ndD)',Levels_2ndD_Train)) | ...
                                                ~isempty(setdiff(unique(DataLabels_Test_2ndD)',Levels_2ndD_Test))) & ...
                                                isnan(ClassifierOpts.One_Class_ResponseLbl2{Cond})) | ...
                                                ((length(unique(DataLabels_Train_3ndD))==1 | length(unique(DataLabels_Test_3ndD))==1 |... % test dim 3
                                                ~isempty(setdiff(unique(DataLabels_Train_3ndD)',Levels_3ndD_Train)) | ...
                                                ~isempty(setdiff(unique(DataLabels_Test_3ndD)',Levels_3ndD_Test))) & ...
                                                isnan(ClassifierOpts.One_Class_ResponseLbl3{Cond}))
                                            isok=0;
                                        else
                                            isok=1;
                                        end
                                    end

                                    %% if one of the dimensions has rule information then add mean activity to that
                                    [predictors,predictors2,predictors3,predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt]=obj.AdjustMeanActivity4Rule(ClassifierOpts,Cond,IncludedNeurons,predictors,predictors2,predictors3,...
                                        predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt,...
                                        TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,...
                                        TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors);

                                     %% if we are adjusting the mean activity of each iteration for train and test seperately
                                    [predictors,predictors2,predictors3,predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt]=obj.AdjustMeanActivity4Iteration(ClassifierOpts,Cond,IncludedNeurons,predictors,predictors2,predictors3,...
                                        predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt,...
                                        TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,...
                                        TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors);

                                    PrepTimeToc=toc(PrepTimeTic);

                                    %% now train and test classifier(no shuffle in this stage at all)
                                    if ShuffObservedRound==0
                                        [ClassifierResults,ClassifierOpts]=obj.RunClassifier3D_Optimized(ClassifierResults,Cond,nCond,nTrlRng,rep,ClassifierOpts,predictors,response,predictors2,response_2ndD,predictors3,response_3ndD,...
                                            predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt,TrainTimInd,TestTimInd,nXtimePnt,TimeMatrixSize,0,NaN);

                                        %% run some preliminary subspace analysis as well in case we want to correlate with results
                                        SubspaceTimeTic=tic;
                                        if sum(strcmp(TargetFields,'Congruency')) & ~sum(ClassifierOpts.StimCongruency==12) & ~sum(ClassifierOpts.StimCongruency==1) & ClassifierOpts.RunPrelimSubspace==1
                                            for Dim=1:3
                                                % Run subspace analysis on different time points on the all dimensions
                                                eval(sprintf('[Subspace%s,~,~,~,FactorDataTiled] =obj.DiscoverSubspaces(cat(1,RunData.predictors%s{1},RunData.predictors%s{2}),ClassifierOpts,[Train%sDataAllFactors;Test%sDataAllFactors],[ones(1,size(RunData.predictors%s{1},1)) 2*ones(1,size(RunData.predictors%s{2},1))]);',DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim}));
                                                eval(sprintf('[ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).PairAngled%s,~,ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).CosTheta%s] =arrayfun(@(x) obj.CalculateAnglesBetSubspaces({Subspace%s.CondScore{1}{x} Subspace%s.CondScore{2}{x}},[]),1:AnalysisOpts.NTim);',DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim})); % calculates angles between the planes that are fit to each condition
                                                eval(sprintf('[ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).Subspace%sCompression1]=arrayfun(@(x) obj.CalculateCompressionBetSubspaces(FactorDataTiled{1}(x)),1:AnalysisOpts.NTim);',DimTxt{Dim}));  % Calculates compression for each subspace
                                                eval(sprintf('[ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).Subspace%sCompression2]=arrayfun(@(x) obj.CalculateCompressionBetSubspaces(FactorDataTiled{2}(x)),1:AnalysisOpts.NTim);',DimTxt{Dim}));  % Calculates compression for each subspace

                                                % if we have any spike count classification that we need to do then do that
                                                for SpkCnt=1:size(ClassifierOpts.SpkCountPeriod,1)
                                                    % run geometry analysis for the third dimension
                                                    eval(sprintf('[Subspace%sSpkCnt{SpkCnt},~,~,~,FactorDataTiledSpkCnt{SpkCnt}] =obj.DiscoverSubspaces(cat(1,RunData.predictors%sSpkCnt{1}{SpkCnt},RunData.predictors%sSpkCnt{2}{SpkCnt}),ClassifierOpts,[Train%sDataAllFactors;Test%sDataAllFactors],[ones(1,size(RunData.predictors%s{1},1)) 2*ones(1,size(RunData.predictors%s{2},1))]);',DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim},DimTxt{Dim}));
                                                    eval(sprintf('[ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).PairAngled%s_SpkCntPrd%i,~,ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).CosTheta%s_SpkCntPrd%i] =obj.CalculateAnglesBetSubspaces({Subspace%sSpkCnt{SpkCnt}.CondScore{1}{1} Subspace%sSpkCnt{SpkCnt}.CondScore{2}{1}},[]);',DimTxt{Dim},SpkCnt,DimTxt{Dim},SpkCnt,DimTxt{Dim},DimTxt{Dim})); % calculates angles between the planes that are fit to each condition
                                                    eval(sprintf('[ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).Subspace%sCompression1_SpkCntPrd%i]= obj.CalculateCompressionBetSubspaces(FactorDataTiledSpkCnt{SpkCnt}{1}(1));',DimTxt{Dim},SpkCnt));  % Calculates compression for each subspace
                                                    eval(sprintf('[ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).Subspace%sCompression2_SpkCntPrd%i]= obj.CalculateCompressionBetSubspaces(FactorDataTiledSpkCnt{SpkCnt}{2}(1));',DimTxt{Dim},SpkCnt));  % Calculates compression for each subspace
                                                end
                                            end
                                        end
                                        SubspaceTimeToc=toc(SubspaceTimeTic);

                                        %% save test data
                                        if ~obj.CalShuff
                                            % save the factor data just keep ColorML and ShapeML reward Rule and response location
                                            ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).TestDataAllFactors  =TestDataAllFactors;%(:,FactorInds2Keep);%[CVTestDataAllFactors(:,1:5) mean(CVTestDataAllFactors(:,6:9,:),3)];
                                            ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).Test2DataAllFactors =Test2DataAllFactors;%(:,FactorInds2Keep);%[CVTest2DataAllFactors(:,1:5) mean(CVTest2DataAllFactors(:,6:9,:),3)];
                                            ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).Test3DataAllFactors =Test3DataAllFactors;%(:,FactorInds2Keep);%[CVTest3DataAllFactors(:,1:5) mean(CVTest3DataAllFactors(:,6:9,:),3)];
                                        end
                                        %% save this classifier result folds
                                        CVfoldResults(cvf)= ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep);
                                    end
                                    %% train and test of shuffled data
                                    if ShuffObservedRound==1 & obj.CalShuff % then shuffle training labels
                                        ShuffleRun=1; % are we shuffleing the trining data
                                        % RunData uses the same data for this run for the shuffle as well
                                        % [ClassifierResults_Shuffled]=obj.RunClassifier3D(ClassifierResults_Shuffled,Cond,nCond,nTrlRng,repshuff,ClassifierOpts,RunData.predictors,RunData.response,RunData.predictors2,RunData.response_2ndD,RunData.predictors3,RunData.response_3ndD,...
                                        % RunData.predictorsSpkCnt,RunData.predictors2SpkCnt,RunData.predictors3SpkCnt,TrainTimInd,TestTimInd,nXtimePnt,TimeMatrixSize,ShuffleRun,repshuff);
                                       
                                        [ClassifierResults_Shuffled,ClassifierOpts]=obj.RunClassifier3D_Optimized(ClassifierResults_Shuffled,Cond,nCond,nTrlRng,rep,ClassifierOpts,predictors,response,predictors2,response_2ndD,predictors3,response_3ndD,...
                                            predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt,TrainTimInd,TestTimInd,nXtimePnt,TimeMatrixSize,ShuffleRun,repshuff);

                                        CVfoldResults_Shuffled(repshuff,cvf)= ClassifierResults_Shuffled(nCond).TrialRange(nTrlRng).Rep(rep);
                                        SubspaceTimeToc=0;
                                    end
                                    TocRepInnerTime=toc(PrepTimeTic);
                                    fprintf('** %s RepTime:%0.2f PrepTime:%0.2f SubSpcTime:%0.2f ',ShuffTxt,TocRepInnerTime,PrepTimeToc,SubspaceTimeToc);
                                end % loop on fold
                                FoldTimeToc=toc(FoldTimeTic);

                                %% take average of the values across CV folds
                                if ShuffObservedRound==0
                                    fprintf(2,'\n*********Total Fold Time:%0.2f Saving Observed*******',FoldTimeToc);
                                    ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep)=obj.AverageCVfoldData(CVfoldResults(1:ClassifierOpts.nCVfold(Cond)),ClassifierOpts);
                                end

                                if ShuffObservedRound==1 & obj.CalShuff
                                    fprintf(2,'\n*********Total Fold Time:%0.2f Saving Shuffle*******',FoldTimeToc);
                                    %  for repshuff=1:ClassifierOpts.NrepShufperFold
                                    ClassifierResults_Shuffled(nCond).TrialRange(nTrlRng).Rep(repshuff,rep)=obj.AverageCVfoldData(CVfoldResults_Shuffled(repshuff,1:ClassifierOpts.nCVfold(Cond)),ClassifierOpts);
                                    %  end
                                end

                            end  % loop on repetition
                        end % loop on shuffle repetition
                    end % loop on ntrl range
                end
            end % loop on Cond

            if AnalysisOpts.GetOnlyShuffLabelsClassifier % if we are only saving shuffle label classifier
                ClassifierResults=[];
                ClassifierOpts=obj.ManData.rmfieldExept(ClassifierOpts,{'ClassifierShuffleLabel','ClassifierShuffleTrialIndex','ClassifierShuffleTrainCondIndex'});
                obj.ManData.SaveVar('Classifier',ClassifierOpts,'ClassifierOpts',ExtraStrSaveShuffLabel,'WantedDate','ALL');
                return
            end

            ClassifierOpts.AnalysisOpts=AnalysisOpts;
            if ~obj.CalShuff
                % save off the variables( this is
                ExtraStrSave= ['_' ShuffTxt ClassifierOpts.Name '_' AnalysisOpts.CurrentCh_AreaName '_' AnalysisOpts.SpkCntStartFieldName '_' AnalysisOpts.TrlSpkTimeFieldName '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
                obj.ManData.DeleteFile('Classifier',ExtraStrSave,1,'WantedDate','ALL');   % first delete existing file
                obj.ManData.SaveVar('Classifier',ClassifierResults,'ClassifierResults',ExtraStrSave,'WantedDate','ALL','SaveAnalysisOpts',0);
                if strcmp(ShuffTxt,'C1_') || isempty(ShuffTxt) % if this is the first condition then save
                    obj.ManData.SaveVar('Classifier',ClassifierOpts,'ClassifierOpts',ExtraStrSave,'WantedDate','ALL','SaveAnalysisOpts',0);
                end
            elseif obj.CalShuff
                ExtraStrSaveShuff=['_' ShuffTxt ClassifierOpts.Name '_' AnalysisOpts.CurrentCh_AreaName '_' AnalysisOpts.SpkCntStartFieldName '_' AnalysisOpts.TrlSpkTimeFieldName '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
                % if we have also calculated the shuffle save the shuffle data too
                % save both observed and shuffle in one fule
                obj.ManData.DeleteFile('Classifier',ExtraStrSaveShuff,1,'WantedDate','ALL');
                if ShuffObservedRound==1
                    obj.ManData.SaveVar('Classifier',ClassifierResults_Shuffled,'ClassifierResults_Shuffled',ExtraStrSaveShuff,'WantedDate','ALL','SaveAnalysisOpts',0); % save shuffled results
                else % save observed
                    obj.ManData.SaveVar('Classifier',ClassifierResults,'ClassifierResults',ExtraStrSaveShuff,'WantedDate','ALL','SaveAnalysisOpts',0); % save observed results
                end
                if strcmp(ShuffTxt,'ShLb_C1_') || isempty(ShuffTxt)    % if this is the first condition then save
                    obj.ManData.SaveVar('Classifier',ClassifierOpts,'ClassifierOpts',ExtraStrSaveShuff,'WantedDate','ALL','SaveAnalysisOpts',0);
                end
            end
        end
