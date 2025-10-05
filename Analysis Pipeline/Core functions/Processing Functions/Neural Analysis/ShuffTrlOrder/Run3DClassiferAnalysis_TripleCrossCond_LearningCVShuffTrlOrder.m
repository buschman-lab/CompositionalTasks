function [ClassifierResults,ClassifierOpts]=Run3DClassiferAnalysis_TripleCrossCond_LearningCVShuffTrlOrder(obj,FactorData,FactorLevelComb,ClassifierOpts,varargin) % With correct shuffling procidure, runs cross temporal classifier anlaysis for 3 simultanous features in three different cross conditions during learning with cross validation
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
            if obj.CalShuffTrlOrder==1;obj.CalShuff=1;end % if we are shuffling the trial order then we are  shuffleing 

            DimTxt={'','2','3'};                      
            TargetFields=obj.DetermineStimCongruency(ClassifierOpts,1,FactorData,FactorData);
            ClassifierOpts=obj.DetermineDimCVspecs(ClassifierOpts); % determine cross validation specs 
          
            [ShuffTxt,CondSet,TrlRngSet,RepSet]=obj.ReadClassifierSpockRunConds(ClassifierOpts);
            ExtraStrSaveShuffLabel=[ClassifierOpts.Name '_' AnalysisOpts.CurrentCh_AreaName '_ShuffleLabelsTrialOrder'];
            RepSetBuff=RepSet;
            % now create a mesh of time points and run classifier train and test repectively
            [ClassifierOpts,TrainTimInd,TestTimInd,nXtimePnt,TimeMatrixSize]=obj.GetTimeRangeforThisCond(ClassifierOpts);

              % if we are cutting the shuffling only to two episodes
            if contains(ClassifierOpts.Name,'CutShuf')
                ClassifierOpts.TestTrlNumRange ={[1 125 50 5]}; 
            end
            % apply neuron dropping algorithm
            [FactorDataBuff,IncludedNeurons,ClassifierOpts]=obj.ApplyNeuronDroppingAlgorithm4LearningAnlaysis(FactorData,ClassifierOpts,FactorLevelComb);
            
            % if we are cutting the shuffling only to two episodes
            if contains(ClassifierOpts.Name,'CutShuf')
                ClassifierOpts.TestTrlNumRange ={[1 125 50 75]}; 
            end
          

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

                % if we are shuffling the data at the level of trial number order then shuffle it here

                %  if AnalysisOpts.GetOnlyShuffLabelsClassifier % generate shuffle labels once for each repetition
                ClassifierOpts=CreateShuffleDist4ClassifierCorrectShuffV5(obj,FactorDataBuff,ClassifierOpts,FactorLevelComb,1,Cond,1);
                ClassifierOpts=CreateShuffleDist4ClassifierCorrectShuffV5(obj,FactorDataBuff,ClassifierOpts,FactorLevelComb,2,Cond,1);
                ClassifierOpts=CreateShuffleDist4ClassifierCorrectShuffV5(obj,FactorDataBuff,ClassifierOpts,FactorLevelComb,3,Cond,1);
                %                elseif obj.CalShuff==1 % generate shuffle labels once for each repetition
                %                     ClassifierOpts4Shuff=obj.ManData.LoadVar('Classifier','ClassifierOpts',ExtraStrSaveShuffLabel,0,'WantedDate','ALL');
                %                     ClassifierOpts=obj.ManData.CopyVars2Struct(ClassifierOpts,'ClassifierShuffleLabel',ClassifierOpts4Shuff.ClassifierShuffleLabel,...
                %                         'ClassifierShuffleTrialIndex',ClassifierOpts4Shuff.ClassifierShuffleTrialIndex,'ClassifierShuffleTrainCondIndex',ClassifierOpts4Shuff.ClassifierShuffleTrainCondIndex);
                %                 end

                if ~AnalysisOpts.GetOnlyShuffLabelsClassifier % we are only getting classifier label then skip the rest
                    %% run the code
                    % create learning trials for this condition
                    [TrainTrlRange,TestTrlRange,TrainTrlInd,TestTrlInd]=obj.GetTrialRangeforThisCond(ClassifierOpts,Cond);

                    % if we have the whole range of TrlRngSet then choose the one for this condition
                    if iscell(TrlRngSet);ThisTrlRngSet=TrlRngSet{Cond};else;ThisTrlRngSet=TrlRngSet;end

                    % determine what condition we are computing the data for
                    if obj.CalShuffTrlOrder==1 % if we are running shuffle
                        if RepSetBuff==ClassifierOpts.NrepShufperFold+1 % this is observed
                            ShuffObservedRound=0;
                            repshuffSet=1;
                        else   % this is shuffle
                            ShuffObservedRound=1;
                            repshuffSet=RepSet;
                        end
                        RepSet=1:ClassifierOpts.NrepShuf;
                    else % this is for mean and STD
                        ShuffObservedRound=0;
                        repshuffSet=1;
                    end

                    for repshuff=1:length(repshuffSet) % shuffle classifier repetition
                        fprintf(2,'\nRunning permutation test on the data, Rep %i/%i',repshuff,length(repshuffSet))

                        if obj.CalShuff & ShuffObservedRound % switch the trialnums and use only the Method 2
                            if length(repshuffSet)>1; error('This Code would not work here');end
                            NewFactorDataTestDim1Buff=obj.ManData.SwapShuffleTrialsFactorData(ClassifierOpts,FactorDataBuff,1,1,FactorLevelComb,TrialNumInd,2,TargetFactorsInd,Cond,1);
                            NewFactorDataTestDim2Buff=obj.ManData.SwapShuffleTrialsFactorData(ClassifierOpts,FactorDataBuff,1,2,FactorLevelComb,TrialNumInd,2,TargetFactorsInd,Cond,1);
                            NewFactorDataTestDim3Buff=obj.ManData.SwapShuffleTrialsFactorData(ClassifierOpts,FactorDataBuff,1,3,FactorLevelComb,TrialNumInd,2,TargetFactorsInd,Cond,1);
                        else
                            NewFactorDataTestDim1Buff=FactorDataBuff;
                            NewFactorDataTestDim2Buff=FactorDataBuff;
                            NewFactorDataTestDim3Buff=FactorDataBuff;
                        end

                        for nTrlRng=1:length(ThisTrlRngSet)
                            TrlRng=ThisTrlRngSet(nTrlRng);

                            ThisTrainTrialRange=TrainTrlRange{TrainTrlInd(TrlRng)};
                            ThisTestTrialRange=TestTrlRange{TestTrlInd(TrlRng)};

                            % now limit factor data trials to these trial range we want generate factor data for training
                            NewFactorDataTestDim1=obj.LimitTrialsBasedonFactor(NewFactorDataTestDim1Buff,'TrialNum',ThisTestTrialRange);
                            NewFactorDataTestDim2=obj.LimitTrialsBasedonFactor(NewFactorDataTestDim2Buff,'TrialNum',ThisTestTrialRange);
                            NewFactorDataTestDim3=obj.LimitTrialsBasedonFactor(NewFactorDataTestDim3Buff,'TrialNum',ThisTestTrialRange);

                            % generate factordata for training
                            FactorData=obj.LimitTrialsBasedonFactor(FactorDataBuff,'TrialNum',ThisTrainTrialRange);

                            % if we are limitign based on block performance then
                            if isfield(ClassifierOpts,'LimitFromSwitchPerf')
                                % generate factor data for test
                                NewFactorDataTestDim1=obj.LimitTrialsBasedonFactor(NewFactorDataTestDim1,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(2),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{2});
                                NewFactorDataTestDim2=obj.LimitTrialsBasedonFactor(NewFactorDataTestDim2,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(2),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{2});
                                NewFactorDataTestDim3=obj.LimitTrialsBasedonFactor(NewFactorDataTestDim3,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(2),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{2});

                                % generate factordata for training
                                FactorData=obj.LimitTrialsBasedonFactor(FactorData,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(1),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{1});
                            end

                            for rep=1:length(RepSet)       % loop on repetition per condition
                                %  CVfoldResults=[]; % initialize this rep
                                FoldTimeTic=tic;
                                for cvf=1:ClassifierOpts.nCVfold(Cond) % loop on cross validation fold

                                    fprintf('\nTrain/Testing TrlShuff Train1:%s Test1:%s | Train2:%s Test2:%s | Train3:%s Test3:%s TrlRng:%i Rep:%i CV:%i',obj.ManData.ConvMat2Char(ClassifierOpts.TrainCond{Cond}),...
                                        obj.ManData.ConvMat2Char(ClassifierOpts.TestCond{Cond}),obj.ManData.ConvMat2Char(ClassifierOpts.TrainCond2{Cond}),obj.ManData.ConvMat2Char(ClassifierOpts.TestCond2{Cond}),...
                                        obj.ManData.ConvMat2Char(ClassifierOpts.TrainCond3{Cond}),obj.ManData.ConvMat2Char(ClassifierOpts.TestCond3{Cond}),TrlRng,rep,cvf);

                                    ClassifierOpts.CurrCVfold=cvf;

                                    PrepTimeTic=tic;
                                    isok=0;
                                    %% sample train and test trials
                                    while isok==0 % make sure we have enough trials for both categories in both dimensions
                                        % grab a second set of random data we will use the test data of this set of the second dimension but keep the training data from first set
                                        if strcmp(type,'TripleCrossCond')
                                            % grab random set of data for Dim 1
                                            [predictors ,response , TrainDataAllFactors,TestDataAllFactors,TrainStimInds,~,predictorsSpkCnt,TrainDataAllFactors3D,TestDataAllFactors3D]=...
                                                obj.GrabThisRepClassifierData(FactorData,ClassifierOpts,ClassifierOpts.TrainCond{Cond},...
                                                ClassifierOpts.TestCond{Cond},FactorLevelComb,Cond,NewFactorDataTestDim1,ThisTrainTrialRange,ThisTestTrialRange,'LookatDim2',0);
                                            % grab random set of data for Dim 2
                                            [predictors2,response2,Train2DataAllFactors, Test2DataAllFactors,Train2StimInds,~,predictors2SpkCnt,Train2DataAllFactors3D,Test2DataAllFactors3D]=...
                                                obj.GrabThisRepClassifierData(FactorData,ClassifierOpts,ClassifierOpts.TrainCond2{Cond},...
                                                ClassifierOpts.TestCond2{Cond},FactorLevelComb,Cond,NewFactorDataTestDim2,ThisTrainTrialRange,ThisTestTrialRange,'LookatDim2',2);
                                            % grab random set of data for Dim 3
                                            [predictors3,response3,Train3DataAllFactors, Test3DataAllFactors,Train3StimInds,~,predictors3SpkCnt,Train3DataAllFactors3D,Test3DataAllFactors3D]=...
                                                obj.GrabThisRepClassifierData(FactorData,ClassifierOpts,ClassifierOpts.TrainCond3{Cond},...
                                                ClassifierOpts.TestCond3{Cond},FactorLevelComb,Cond,NewFactorDataTestDim3,ThisTrainTrialRange,ThisTestTrialRange,'LookatDim2',3);
                                            % euqalize trials across conditions
                                            [predictors,predictors2,predictors3,response,response2,response3,TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,...
                                                TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors,...
                                                predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt,ClassifierOpts]=obj.EqualizeTrialsAcrossConds(ClassifierOpts,Cond,predictors,predictors2,predictors3,...
                                                response,response2,response3,TrainDataAllFactors,Train2DataAllFactors,Train3DataAllFactors,...
                                                TestDataAllFactors,Test2DataAllFactors,Test3DataAllFactors,...
                                                predictorsSpkCnt,predictors2SpkCnt,predictors3SpkCnt);
                                           
                                            if cvf==1 % if this is the first CVf then save the data since this is the data that has been copied
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
                                     %   if ~obj.CalShuff
                                            % save the factor data just keep ColorML and ShapeML reward Rule and response location
                                            ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).TestDataAllFactors  =TestDataAllFactors;%(:,FactorInds2Keep);%[CVTestDataAllFactors(:,1:5) mean(CVTestDataAllFactors(:,6:9,:),3)];
                                     %       ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).Test2DataAllFactors =Test2DataAllFactors;%(:,FactorInds2Keep);%[CVTest2DataAllFactors(:,1:5) mean(CVTest2DataAllFactors(:,6:9,:),3)];
                                     %       ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep).Test3DataAllFactors =Test3DataAllFactors;%(:,FactorInds2Keep);%[CVTest3DataAllFactors(:,1:5) mean(CVTest3DataAllFactors(:,6:9,:),3)];
                                     %   end
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
                                        
                                        % because in learning conditions are usually equalized we don't save all of the dimensions for shuffle
                                        ClassifierResults_Shuffled(nCond).TrialRange(nTrlRng).Rep(rep).TestDataAllFactors  =TestDataAllFactors;%(:,FactorInds2Keep);%[CVTestDataAllFactors(:,1:5) mean(CVTestDataAllFactors(:,6:9,:),3)];
                                    %    ClassifierResults_Shuffled(nCond).TrialRange(nTrlRng).Rep(rep).Test2DataAllFactors =Test2DataAllFactors;%(:,FactorInds2Keep);%[CVTest2DataAllFactors(:,1:5) mean(CVTest2DataAllFactors(:,6:9,:),3)];
                                    %    ClassifierResults_Shuffled(nCond).TrialRange(nTrlRng).Rep(rep).Test3DataAllFactors =Test3DataAllFactors;%(:,FactorInds2Keep);%[CVTest3DataAllFactors(:,1:5) mean(CVTest3DataAllFactors(:,6:9,:),3)];
                                  
                                        CVfoldResults_Shuffled(repshuff,cvf)= ClassifierResults_Shuffled(nCond).TrialRange(nTrlRng).Rep(rep);
                                        SubspaceTimeToc=0;
                                    end
                                    TocRepInnerTime=toc(PrepTimeTic);
                                    fprintf('** %s RepTime:%0.2f PrepTime:%0.2f SubSpcTime:%0.2f ',ShuffTxt,TocRepInnerTime,PrepTimeToc,SubspaceTimeToc);
                                end % loop on fold
                                FoldTimeToc=toc(FoldTimeTic);
                                fprintf(2,'\n*********Total Fold Time:%0.2f *******',FoldTimeToc);

                                %% take average of the values across CV folds
                                if ShuffObservedRound==0
                                    ClassifierResults(nCond).TrialRange(nTrlRng).Rep(rep)=obj.AverageCVfoldData(CVfoldResults(1:ClassifierOpts.nCVfold(Cond)),ClassifierOpts);
                                end

                                if ShuffObservedRound==1 & obj.CalShuff
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
                if strcmp(ShuffTxt,'Shuf_C1_') || isempty(ShuffTxt)    % if this is the first condition then save
                    obj.ManData.SaveVar('Classifier',ClassifierOpts,'ClassifierOpts',ExtraStrSaveShuff,'WantedDate','ALL','SaveAnalysisOpts',0);
                end
            end
        end
 