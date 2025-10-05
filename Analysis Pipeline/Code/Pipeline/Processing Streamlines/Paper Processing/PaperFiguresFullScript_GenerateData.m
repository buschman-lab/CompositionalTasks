%% this script is to generate the main data for the paper
%% for classifier 
%% for classifier 
% @RunMode:0 dont run shuffle, run main and plot; 1= run missing files main and shuffle and then plot
% 2=plot only 3=Run main, no shuffle no plot 4=run missing main only no shuffle or plot (in subspace it is run shuffle only)
% 5: plot comparision only 6: collect comparision results
% 7: GLM: plot summery for each area 8: concatinate remaining files 9: run shuffle only  
% all of these figure configurations go together with

%% for GLM 
% Runmode 0: fit each neuron without shuffle (PS11g)
% Runmode 1: fit each neuron with it's shuffle (PS11), concatinate (PS15) it and create a summery(PS14) for each neuron(No Plots) 
% Runmode 2: Plot compact singlecell cahractristics with GLM(PS10)
% Runmode 5: run model comparision(PS12)
% Runmode 6 plot model comparision results(PS13)
% Runmode 7 plot all GLM result summery for each area(PS16)
% Runmode 8  concatinate GLM result summery for each area(PS14)
% Runmode 9  concatinate GLM result summery for each neuron(PS15)

global AnalysisOpts

% FigureSummery072523.ppt saved in OneDrive folder of the project
CL=ClusterFuncs;
CL.GenClusterStatusFile;
NeuAna=NeuralAnalysisFuncs;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%<<<<<<<<<<<Section Double check>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double check if we have all of the files and generate what we don't have
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001};%,'Enforce_nCondsRun_Shuff',{[1001]}};

ColorShapeEncoding.Conditions={'3D_Color_Color_BalRespDir'};%  '3D_Color_Response_XgenBalInCongV6',
ColorShapeEncoding.AreaXTemp=1;
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=[1:5]; 

% if we need to find the shuffle classifiers that are failed
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,2,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    'RunMode2',1,ColorShapeEncoding.ExtraConds{:},'ReWriteClassifierData',1);

% run PFC stuff
ColorShapeEncoding.Conditions={'3D_Response_Xgen_BalInCong'};%'3D_Response_Xgen_BalInCongRBC' '3D_Color_Response_XgenBalInCongV5',
ColorShapeEncoding.ConditionsXTemp={'3D_Color_Response_XgenBalInCongV2'}; 
ColorShapeEncoding.AreaXTemp=1;
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1; 

% if we need to find the shuffle classifiers that are failed
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,2,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    'RunMode2',2,ColorShapeEncoding.ExtraConds{:},'ReWriteClassifierData',0);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Learning conditionsLearning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,'ReWriteClassifierData',0,'ThisIsSinaPC',1};
Learning.ExtraConds={'Classifier_Nrep',100,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysisTemp'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','Learning3D_Color_Color_Color_Xgen_AltRule_RB'};%,'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB',...
%    'Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Color_Response_Rule_Xgen_AltRule'}; 

%   { 'Learning3D_Shape_Color_Rule_Xgen_AltRuleSameRule_RB'};'Learning3D_ShapeRB_ColorRB_Rule_Xgen_AltRule_RB'     'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB_EQ' 
%'Learning3D_ShapeRB_ColorRB_Rule_Xgen_AltRule_RB'};%'Learning3D_Shape_Color_Rule_Xgen_AltRuleSameRule_RB'};
%{'Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};

% {'Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Color_Response_Rule_Xgen_AltRule'};
% 'Learning3D_Shape_Color_Rule_Xgen_AltRuleR1_RB','Learning3D_Shape_Color_Rule_Xgen_SameRuleR1_RB' 
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;

% geneCrate and plot data
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',0,'RunMode2',1,...
    'ntrlPerCondTestArea',nan,'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',nan,'NTrlStpLearningArea',5);

% Parameter Sweep
GenerateAllFiguresPipeline('Classifier','ALL', Learning.Area,0,100,[{'Decoding'}],Learning.Conditions,Learning.StrTime,'SweepClassifierConds',1,Learning.ExtraConds{:});
NeuAna.CollectTableClassifierSweepParam(Learning.Conditions);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Learning analysis with Shuffle trial order %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysisTemp'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB'};%,'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};%,...
  %  {'Learning3D_Color_Response_Rule_Xgen_AltRule','Learning3D_Color_Response_Rule_Xgen_SameRule',...
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;

% geneCrate and plot data
% generate shuffles for First Last based on Feature
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',1);

Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysis'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

% do full shuffling with balanced trials 
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB'};

% generate shuffles for First Last based on shuffling all data
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',0);

Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RBN'};
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',1);

% generate shuffle with balancing but don't euqalize across dimensions
% uses the previous classifier specs generated
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RBN'};
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'ClassifierNameExtra','NoEq');

GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'ClassifierNameExtra','NoEq');

% generate shuffle with balancing and random subtract the mean but don't euqalize across dimensions
% uses the previous classifier specs generated
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RBN'};

GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'ClassifierNameExtra','NoEq');

GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'ClassifierNameExtra','NoEq');

Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB'};

% run shuffle across all trials random without equalizing 
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,'Subtract4TrialRange',0);

% run shuffle across all trials balanced without equalizing 
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,'Subtract4TrialRange',0);

% run shuffle across all trials random without equalizing ; only for the first block 1-50
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',2,...
    'Subtract4TrialRange',0);

% run shuffle across all trials random without equalizing ; Do for trials 1-63 and 63-126
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysis'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

Learning.Conditions={'Learning3D_Color_Color_Color_Xgen_AltRule_RB'};
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',1,'EqualizeTrialsXConds',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',3,...
    'Subtract4TrialRange',0);

Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysis'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

% generate classifier for each step 
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,3,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',0,'RunMode2',2,...
    'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,'NTrlRngTestLearningArea',nan,...
    'NTrlStpLearningArea',1,'ClassifierNameExtra','ES');

% run shuffle across all trials random without equalizing ; Do for trials  and 63-126
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysis'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB'};
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',4,...
    'Subtract4TrialRange',0);

GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',4,...
    'Subtract4TrialRange',0);

%% run with all trials and equalize across conditions for feature 
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysisTemp'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;

Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB'};

% our normal but with RND
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',0,'Classifier_NrepShufperFold',100,'Classifier_NrepShuf',25,...
    'ntrlPerCondTestArea',nan,'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',nan,'NTrlStpLearningArea',nan);

% 60 trial of test and 2 trls 
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',1,'EqualizeTrialsXConds',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',0,'Classifier_NrepShufperFold',100,'Classifier_NrepShuf',50,...
    'ntrlPerCondTestArea',2,'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',50,'NTrlStpLearningArea',nan);

% 60 trial of test and 3 trls 
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',0,'Classifier_NrepShufperFold',100,'Classifier_NrepShuf',25,...
    'ntrlPerCondTestArea',3,'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',60,'NTrlStpLearningArea',nan,'ClassifierNameExtra','TS3');

% 70 trial of test and 2 trls Noshuffle
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',0,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',0,'Classifier_NrepShufperFold',100,'Classifier_NrepShuf',25,...
    'ntrlPerCondTestArea',2,'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',70,'NTrlStpLearningArea',nan);

% other conditions
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};%'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};%,,'Learning3D_Color_Response_Rule_Xgen_AltRule'};
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,8,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0);%;,'Enforce_nCondsRun_Shuff',{[481,482]});

%% generate shuffle for learning Axis decoding 
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysis'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

Learning.StrTime={'SAMPLE_ON'};
 Learning.Area=1:5;
Learning.Conditions={'Learning3D_AxisDecoding'};
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',0,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',0);


for Area=1:5
    Learning.Area=Area;
    PopulationAnalysisTemp(0, 0, [], Learning.Area, Learning.Conditions{1}, 'ProcessingStep',4, ...
        'IncludeAllClassifierInfo',0, ...
        'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',0,'RunMode2',2,...
        'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
        'Classifier_TrlShuff_GenClassifierSpecs',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',5,...
        'Subtract4TrialRange',0,'IncludeAllClassifierInfo',0);
end

 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%<<<<<<<<<<<Section 0:Behavior>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generates figure for all animals 
% the results will be saved in BasePath\Results\Behavior
GenerateBhvPaperFigs([], 'ALL')
GenerateBhvPaperFigs([], 'Silas')
GenerateBhvPaperFigs([], 'Chico')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%<<<<<<<<<<<Section 1:GLM>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run GLM model
GenerateAllFiguresPipeline('GLM','ALL', [1:5],0,100,[],[],{'SAMPLE_ON'},'ThisisSinaPC',1)
% fit GLM for each neuron with 1000 shuffle 
GenerateAllFiguresPipeline('GLM','ALL', [1:5],1,100,[],[],{'SAMPLE_ON'})
% Plot compact singlecell charactristics with GLM
GenerateAllFiguresPipeline('GLM','ALL', [1:5],2,100,[],[],{'SAMPLE_ON'})
% generate data for model comparision
GenerateAllFiguresPipeline('GLM','ALL', [1:5],5,100,[],[],{'SAMPLE_ON'})
% plot model comparision results
GenerateAllFiguresPipeline('GLM','ALL', [1:5],6,100,[],[],{'SAMPLE_ON'})
% concatinate GLM result summery for each neuron and then for area(PS15)
GenerateAllFiguresPipeline('GLM','ALL', [1],9,100,[],[],{'SAMPLE_ON'})
% concatinate GLM result summery for each area(PS14)
GenerateAllFiguresPipeline('GLM','ALL', [1:5],8,100,[],[],{'SAMPLE_ON'},'ClassifierFunctiononClust','PopulationAnalysis')

%% plotting single cell results
% plot all GLM result summery for each area(PS16)
GenerateAllFiguresPipeline('GLM','ALL', [1:5],7,100,[],[],{'SAMPLE_ON'})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%<<<<<<<<<<<Section 2:Compositionality of Reperesentations>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%<<<<<<<<<<<Section 2-1:Color & Shape Information>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
% runs color and shape information for balanced and unbalanced conditions
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1};%,'Enforce_nCondsRun_Shuff',{1001}};
ColorShapeEncoding.Conditions={'3D_Color_Response_XgenBalInCongV5'};%,'3D_Shape_Cat_Xgen'};%,'3D_Response_Xgen_BalInCong'};%3D_Color_Response_XgenBalInCongV2 '3D_Color_Response_XgenBalInCongV2'};%,'3D_Color_Response_XgenBalInCong','3D_Color_Cat_Xgen','3D_Shape_Cat_Xgen'};%'3D_Shape_Response_Xgen',,'3D_Color_Response_Xgen'}; 
ColorShapeEncoding.ConditionsXTemp={'3D_Color_Response_XgenBalInCongV2'};%,'3D_Color_Cat_Xgen'};
ColorShapeEncoding.AreaXTemp=1;
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1:5; 

% step 1 sweep parameters
GenerateAllFiguresPipeline('Classifier','ALL', ColorShapeEncoding.Area,0,100,[{'Decoding'}],ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,'SweepClassifierConds',1);

% find the best number of trials for each area
NeuAna.CollectTableClassifierSweepParam(ColorShapeEncoding.Conditions);

% step 2 run main conditions
GenerateAllFiguresPipeline('Classifier','ALL', ColorShapeEncoding.Area,3,100,[{'Decoding'}],ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:});

% step 3 run cross temporal color classifier for these same conditions for now run this without shuffle since it takes a long time
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.AreaXTemp,0,100,{'Decoding'},ColorShapeEncoding.ConditionsXTemp,ColorShapeEncoding.StrTime,...
    'RunCrossTemporalClassifer',1,'DividSpockClassifier',2,'RunMode2',3,'RunDummyFile',0);AnalysisOpts.RunCrossTemporalClassifer=0;
AnalysisOpts.RunCrossTemporalClassifer=0;

% Run shuffle classifiers  
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,9,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
   ColorShapeEncoding.ExtraConds{:},'ReWriteClassifierData',1);

% if we need to find the shuffle classifiers that are failed
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,9,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    'RunMode2',2,ColorShapeEncoding.ExtraConds{:},'ReWriteClassifierData',0);

% if we need to concatinate shuffle classifier results
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,8,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    'RunMode2',2,ColorShapeEncoding.ExtraConds{:});

% Collect the figures generated by this analysis 
GenerateAllFiguresPipeline('Classifier','ALL', ColorShapeEncoding.Area,7,100,[{'Decoding'}],ColorShapeEncoding.Conditions,...
    {'SAMPLE_ON'},'Page2SaveNum',{[1 4],[1 4],[1 4]});

%% Comparision of encoding of color and response across areas
% runs color and shape information for balanced and unbalanced conditions
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'RunMode2',1,'ThisIsSinaPC',1};
ColorShapeEncoding.Conditions={'3D_Color_Response_XgenBalInCongV5_AreaComp'};%,'3D_Shape_Cat_Xgen'};%,'3D_Response_Xgen_BalInCong'};%3D_Color_Response_XgenBalInCongV2 '3D_Color_Response_XgenBalInCongV2'};%,'3D_Color_Response_XgenBalInCong','3D_Color_Cat_Xgen','3D_Shape_Cat_Xgen'};%'3D_Shape_Response_Xgen',,'3D_Color_Response_Xgen'}; 
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1; 
global AnalysisOpts
AnalysisOpts.PopulationAna.Classifier_TaskConcateSpockTime=1500*ones(1,100);

% Run resampling of classifier conditions classifiers 
% Compare PFC with STR
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,8,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:},'NeuResamplingSubPopSize',110);
% Compare PFC with IT
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,8,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:},'NeuResamplingSubPopSize',195);
% Compare PFC with FEF
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,8,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:},'NeuResamplingSubPopSize',116);
% Compare PFC with LIP
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,8,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:},'NeuResamplingSubPopSize',54);

AnalysisOpts.PopulationAna.Classifier_TaskConcateSpockTime=800*ones(1,100);

% to plot the area comparision
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'RunMode2',1,'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysisTemp'};
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1; 
ColorShapeEncoding.Conditions={'3D_Xgen_Compare'};
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,5,100,{'Comparision'},...
    ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:},'NeuResamplingSubPopSize',x),[110 195 116 54]);% 

%% Run Beta values for area comparision
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'RunMode2',1,'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysisTemp'};

ColorShapeEncoding.Conditions={'3D_Color_Shape_AxisV2'};%,'3D_Color_Shape_Axis'};
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1; 
 
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,1,100,{'Decoding'}, ...
    ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:});

%% Plotting for paper
%plot Main figure (figure 2)
ColorShapeEncoding.DataGenDate=datetime(2023,12,21);
ColorShapeEncoding.Conditions={'3D_Color_Response_XgenBalInCongV5'};%'3D_Response_Xgen_BalInCong',,'3D_Color_Cat_Xgen','3D_Shape_Cat_Xgen','3D_Response_Xgen_BalInCong'};%,'3D_Color_Response_XgenBalInCongV2'};
ColorShapeEncoding.Area=[1:5]; % only plot PFC for this
ColorShapeEncoding.StrTime={'SAMPLE_ON'}; % for now it is only sample on
ColorShapeEncoding.ExtraConds={'DividSpockClassifier',3,'RunMode2',3,'Classifier_FileDateTimeTh',ColorShapeEncoding.DataGenDate,'ReWriteClassifierData',0};

GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,2,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:});

% run on PC
cellfun(@(x) PopulationAnalysisTemp(0, 0, [], [ColorShapeEncoding.Area], x, 'ProcessingStep',4,'PopulationAna.PSTHbin',100,'DividSpockClassifier_TrlRng',[1],'DividSpockClassifier_Cond',2,'CalShuffleClassifier',1,'CalShuffTrlOrderClassifier',0,...
    'SpkCntStartFieldName',ColorShapeEncoding.StrTime{1},'NTrlRngTrainLearningArea',nan,'NTrlRngTestLearningArea',nan,'ntrlPerCondArea',nan),ColorShapeEncoding.Conditions,'UniformOutput',0)

% plot comparision of timing 
ColorShapeEncoding.Conditions={'3D_Cat_Color_Response_Compare'};
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,5,100,{'Comparision'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:});

% comparision of timing between generalization of color for other areas
ColorShapeEncoding.Conditions={'3D_Cat_Color_XgenCol_Compare'};
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,5,100,{'Comparision'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:});

% comparision of timing between generalization of response for other areas
ColorShapeEncoding.Conditions={'3D_Response_XgenResponse_Compare'};
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,5,100,{'Comparision'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:});

% run comparision on PC
ColorShapeEncoding.ComparisionConditions={'3D_Cat_Color_Response_Compare','3D_Cat_Color_XgenCol_Compare','3D_Response_XgenResponse_Compare'};
cellfun(@(x) PopulationAnalysisTemp(0, 0, [], [ColorShapeEncoding.Area], x, 'ProcessingStep',5,'PopulationAna.PSTHbin',100,'DividSpockClassifier_TrlRng',[1],'DividSpockClassifier_Cond',2,'CalShuffleClassifier',1,'CalShuffTrlOrderClassifier',0,...
    'SpkCntStartFieldName',ColorShapeEncoding.StrTime{1},'NTrlRngTrainLearningArea',nan,'NTrlRngTestLearningArea',nan,'ntrlPerCondArea',nan),ColorShapeEncoding.ComparisionConditions,'UniformOutput',0)

% plot cross temporal analaysis
ColorShapeEncoding.ConditionsXTemp={'3D_Color_Response_XgenBalInCongV2'};%'3D_Color_Cat_Xgen'};%
ColorShapeEncoding.AreaXTemp=1;
GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.AreaXTemp,2,100,{'Decoding'},ColorShapeEncoding.ConditionsXTemp,{'SAMPLE_ON'},...
    'RunCrossTemporalClassifer',1,'DividSpockClassifier',2,'RunMode2',3,'RunDummyFile',0,'Classifier_FileDateTimeTh',datetime(2023,8,4));AnalysisOpts.RunCrossTemporalClassifer=0;
AnalysisOpts.RunCrossTemporalClassifer=0;

% SUPP figures
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'pvalClassifierAnalysis_ClusterCorrect',0.001};
ColorShapeEncoding.Conditions={'3D_Shape_Cat_Xgen'};%,'3D_Color_Response_XgenBalInCongV5'};%,,'3D_Color_Response_XgenBalInCong'};%'3D_Color_Cat_Xgen','3D_Shape_Cat_Xgen'};%'3D_Shape_Response_Xgen',,'3D_Color_Response_Xgen'}; 
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=[1:5];
GenerateAllFiguresPipeline('Classifier','ALL', ColorShapeEncoding.Area,2,100,[{'Decoding'}],ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    ColorShapeEncoding.ExtraConds{:});

%% %%%%%%<<<<<<<<<<<Section 2-2:Response Information>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
% runs response information
ResponseEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,'ReWriteClassifierData',1};
ResponseEncoding.Conditions={'3D_Response_Xgen_BalInCong'};%'3D_Response_Xgen_BalCorr','3D_Response_Xgen_BalCong',
ResponseEncoding.Area=1;
% step 1 sweep parameters
GenerateAllFiguresPipeline('Classifier','ALL', ResponseEncoding.Area,0,100,[{'Decoding'}],ResponseEncoding.Conditions,{'SAMPLE_ON'},'SweepClassifierConds',1);
% find the best number of trials for each area
NeuAna.CollectTableClassifierSweepParam(ResponseEncoding.Conditions);
% step 2 run conditions
GenerateAllFiguresPipeline('Classifier','ALL', ResponseEncoding.Area,1,100,[{'Decoding'}],ResponseEncoding.Conditions,{'SAMPLE_ON'},ResponseEncoding.ExtraConds{:});
% step 3 run cross temporal color classifier for these same conditions
% for now run this without shuffle since it takes a long time
GenerateAllFiguresPipeline('Classifier','ALL', ResponseEncoding.Area,1,100,[{'Decoding'}],ResponseEncoding.Conditions,{'SAMPLE_ON'},'RunCrossTemporalClassifer',1);

%% %%%%%%<<<<<<<<<<<Section 2-3:Entropy Analysis>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
% runs response information
Entropy.ExtraConds={'Classifier_Nrep',500,'DividSpockClassifier',2,'RunDummyFile',0,'SweepClassifierConds',0};
Entropy.Conditions={'3D_Shared_Color_Response_EntropyR3','3D_Shared_Color_Response_EntropyInCongR3','3D_Shared_Color_Response_EntropyR2','3D_Shared_Color_Response_EntropyInCongR2'};
Entropy.Area=1;Entropy.StrTime={'SACCADE_START','SAMPLE_ON'};
% step 1 sweep parameters
GenerateAllFiguresPipeline('Classifier','ALL',Entropy.Area,0,100,{'Decoding'},Entropy.Conditions,{'SAMPLE_ON'},'SweepClassifierConds',1);
% find the best number of trials for each area
NeuAna.CollectTableClassifierSweepParam(Entropy.Conditions);
% shared response and shared color category (color 2/3 -> resp (2/3) and resp (1/3))
% shared color transfering to same axis or the other axis (color 3/2 -> resp (2/2) and rep (1/2)
GenerateAllFiguresPipeline('Classifier','ALL',Entropy.Area,0,100,{'Decoding'},Entropy.Conditions(1:2),Entropy.StrTime,Entropy.ExtraConds{:});
% same analysis but take only incongruent trials 
GenerateAllFiguresPipeline('Classifier','ALL',Entropy.Area,0,100,{'Decoding'},Entropy.Conditions(3:4),Entropy.StrTime,Entropy.ExtraConds{:});

% run old entropy analysis
Entropy.ExtraConds={'Classifier_Nrep',500,'DividSpockClassifier',2,'RunDummyFile',0,'SweepClassifierConds',0};
Entropy.Conditions={'3D_Shared_Color_Response_EntropyOld','3D_Color_Response_XgenR2R3Entropy'};
Entropy.Area=1;Entropy.StrTime={'SACCADE_START'};
GenerateAllFiguresPipeline('Classifier','ALL',Entropy.Area,2,100,{'Decoding'},Entropy.Conditions(1:2),Entropy.StrTime,Entropy.ExtraConds{:});

%% exploratory entropy analysis
Entropy.ExtraConds={'Classifier_Nrep',1000,'DividSpockClassifier',3,'RunDummyFile',0,...
    'SweepClassifierConds',0,'ClassifierFunctiononClust','PopulationAnalysisTempCongEntropy','ThisIsSinaPC',1};
Entropy.Conditions={'3D_Shared_Color_Response_EntropyR3'};%'3D_Shared_Color_Response_EntropyR3Bal',,'3D_Shared_Color_Response_EntropyR2Bal','3D_Shared_Color_Response_EntropyR2'}; %3D_Shared_Color_Response_EntropyR3Bal_BalCongTest
Entropy.Area=1;Entropy.StrTime={'SACCADE_START'};
NTrlRngTrainLearningAreaSet=[-50]; % -75 800 number of trials to take from end of the block for training per area
NTrlRngTestLearningAreaSet=[-50]; % -75 800 number of trials shift for Test
ntrlPerCondAreaSet=[4]; % number of training trails per area

arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Entropy.Area,2,100,{'Decoding'},Entropy.Conditions,Entropy.StrTime,...
    Entropy.ExtraConds{:},'NTrlRngTrainLearningArea',NTrlRngTrainLearningAreaSet(x),...
    'NTrlRngTestLearningArea',NTrlRngTestLearningAreaSet(x),'ntrlPerCondArea',ntrlPerCondAreaSet(x)),1:length(ntrlPerCondAreaSet));

%% plot entropy results 
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Entropy.Area,2,100,{'Decoding'},Entropy.Conditions,Entropy.StrTime,...
    Entropy.ExtraConds{:},'NTrlRngTrainLearningArea',NTrlRngTrainLearningAreaSet(x),...
    'NTrlRngTestLearningArea',NTrlRngTestLearningAreaSet(x),'ntrlPerCondArea',ntrlPerCondAreaSet(x)),1:length(ntrlPerCondAreaSet));

PopulationAnalysisTemp(0, 0, [], [1], [11],'ProcessingStep',4,'SpkCntStartFieldName','SACCADE_START','CalShuffTrlOrderClassifier',0,'RunCrossTemporalClassifer',0,'CalShuffleClassifier',0,...
'PopulationAna.PSTHbin',100,'SpkCntStartFieldName','SACCADE_START','DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',16*ones(1,7),'DividSpockClassifier',3,'NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',4)

% plot entropy for congruency results '3D_Shared_Color_Response_EntropyR3Bal'
PopulationAnalysisTempCongEntropy(0, 0, [], [1], [11],'ProcessingStep',4,'SpkCntStartFieldName','SACCADE_START','CalShuffTrlOrderClassifier',0,'RunCrossTemporalClassifer',0,'CalShuffleClassifier',0,...
'PopulationAna.PSTHbin',100,'SpkCntStartFieldName','SACCADE_START','DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',16*ones(1,7),'DividSpockClassifier',3,'NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',4)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%<<<<<Section 3:Learning analysis>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,'ReWriteClassifierData',1};
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};
% 'Learning3D_Color_Response_Rule_Xgen_AltRule',...
%  'Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};
%'Learning3D_Shape_Color_Rule_Xgen_AltRuleR1_RB','Learning3D_Shape_Color_Rule_Xgen_SameRuleR1_RB'}, Learning3D_Shape_Color_Rule_Xgen_AltRule_RB_CutShuff;

Learning.ConditionsShuffleTrlOrder=Learning.Conditions;
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;

% step 1 sweep parameters
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',1);
NeuAna.CollectTableClassifierSweepParam(Learning.Conditions); % collect the results 

% generate and plot data
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,2,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',0);

% generate shuffle trial order data
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.ConditionsShuffleTrlOrder,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1);

% generate remaining shuffle trial order data
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.ConditionsShuffleTrlOrder,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'ReWriteClassifierData',0);

% concatinate shuffle files
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,8,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',0,'RunMode2',2);

% generate extra conditions with limited threshold for perfromance
Learning.Conditions={'Learning3D_Color_Shape_Rule_Xgen_AltRule'};
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'LimitFromSwitchPerfArea',x),[0.55 0.6]);

% generate extra conditions with ntrlPerCondArea
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'ntrlPerCondArea',x),[4]);

% generate extra conditions with steps
Learning.Conditions={'Learning3D_Color_Shape_Rule_Xgen_AltRule'};
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'NTrlStpLearningArea',x),[1]);

%% plot learning analysis
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB', 'Learning3D_Color_Response_Rule_Xgen_AltRule',...
    'Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};
%'Learning3D_Shape_Color_Rule_Xgen_AltRuleR1_RB','Learning3D_Shape_Color_Rule_Xgen_SameRuleR1_RB'}, Learning3D_Shape_Color_Rule_Xgen_AltRule_RB_CutShuff;

GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,2,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0);

% plot learning extra conditions
Learning.Conditions={'Learning3D_Color_Shape_Rule_Xgen_AltRule'};
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,2,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'LimitFromSwitchPerfArea',x),[0.45 40]);

%% SUPP learning analysis 
% learning for other areas
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0};
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB'};
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=2:5;
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0);

% look at learning for R1 same and Altrule for PFC
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0};
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRuleR1_RB','Learning3D_Shape_Color_Rule_Xgen_SameRuleR1_RB','Learning3D_Shape_Color_Rule_Xgen_AltRule_HiTh_RB'};
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,1,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0);

%look at learning from SACCADE START for PFC
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0};
Learning.Conditions={'Learning3D_Color_Shape_Rule_Xgen_AltRule'};
Learning.StrTime={'SACCADE_START'};
Learning.Area=1;
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%<<<<<Section 4:Compression analysis>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,'ThisIsSinaPC',1};
Learning.Conditions={'Learning3D_Shape_Color_Color_Compression_RB'};
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'ntrlPerCondArea',9);

% look at other areas
Learning.Area=2:5;
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,2,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'ntrlPerCondArea',x),[12 9 5],'UniformOutput',0);

%% %%%%%%%%%%%%%%%%%%%%%%%%<<<<<Section 5:Orthogonality analysis>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1};
ColorShapeEncoding.Conditions={'Learning3D_Shape_Color_Response_Xgen_Orhto'};
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=[1]; 

 GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,0,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    ColorShapeEncoding.ExtraConds{:},'ReWriteClassifierData',1);

%% Angle analysis 
% run Learning3D_Shape_Color_Color_Xgen_Angle for angle analysis during learning
% run Learning3D_Color_Color_Shape_Xgen_Angle_BlkEnd for angle analysis end of the block 
% run Learning3D_Color_ColorRB_Response_Xgen_Angle_BlkEnd for angle analysis end of the block and using balanced trials

ColorShapeEncoding.ExtraConds={'Classifier_Nrep',100,'DividSpockClassifier',3,'RunDummyFile',0,'ReWriteClassifierData',1,...
    'ClassifierFunctiononClust','PopulationAnalysisTemp','ThisIsSinaPC',1};
ColorShapeEncoding.Conditions={'Learning3D_Color_Color_Shape_Xgen_Angle_BlkEnd_Final'}; %'Learning3D_Color_ColorRB_Response_Xgen_Angle_BlkEnd_Final',
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1; 

NTrlRngTrainLearningAreaSet=[-50];%,-75,-50,-75];  
ntrlPerCondAreaSet         =[5];%,5,4,4]; 
NTrlRngTestLearningAreaSet =[50];

arrayfun (@(x) GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,2,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    ColorShapeEncoding.ExtraConds{:},'NTrlRngTrainLearningArea',NTrlRngTrainLearningAreaSet(x),...
    'ntrlPerCondArea',ntrlPerCondAreaSet(x),'NTrlRngTestLearningArea',NTrlRngTestLearningAreaSet(x)),1:length(ntrlPerCondAreaSet));

%% BlkSeq Analysis for Rule 2
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysis',...
    'Classifier_NrepShuf',100,'Classifier_NrepShufperFold',250,'CalShuffTrlOrderClassifier',0,...
    'EnforceSpockTime',200};%,'Enforce_nCondsRun',{2}};

ColorShapeEncoding.Conditions={'Learning3D_RuleR2SeqAnalysisLearningCtrlR2',...
    'Learning3D_RuleR2SeqAnalysisLearningCtrl',...
    'Learning3D_RuleR2SeqAnalysisLearningCtrlR2onR2'};

%'Learning3D_RuleR2SeqAnalysisLearningCtrlMS'
%'Learning3D_RuleR2SeqAnalysisLearningV2','Learning3D_RuleR2SeqAnalysisLearning','Learning3D_RuleR2SeqAnalysisLearningCtrl'};
%'Learning3D_RuleR2SeqAnalysisXgen','Learning3D_RuleR2SeqAnalysisLearningSameRule','Learning3D_RuleR2SeqAnalysisLearningAltRule'}; 

ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1; 

GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,3,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    ColorShapeEncoding.ExtraConds{:});%,'ntrlPerCondArea',5);

%% Axis suppression analysis
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysis'};%,'Enforce_nCondsRun',{2}};

ColorShapeEncoding.Conditions={'Learning3D_AxisDecoding'};%'3D_AxisDecodingSupp'};%''};%'3D_AxisDecodingSupp'};'3D_Response_Xgen_BalInCongAxisSupp' 
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1:4; 

GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,9,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
    ColorShapeEncoding.ExtraConds{:}); 

%% Run Targeted Dim reduction analysis 
TDRana.ExtraConds={'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysisTemp','DividSpockClassifier',3,'DividSpockGLM_Cond',1};
TDRana.Conditions={''};
TDRana.StrTime={'SACCADE_START','SAMPLE_ON'};%
TDRana.Area=[1]; 

GenerateAllFiguresPipeline('Classifier','ALL',TDRana.Area,10,100,{'Decoding'},TDRana.Conditions,TDRana.StrTime,...
    TDRana.ExtraConds{:});

%% Run Angle Analysis Targeted Dim reduction analysis 
% generate GLM models 

% fit data
GenerateAllFiguresPipeline('GLM','ALL', [1],0,100,[],[],{'SAMPLE_ON'},...
    'ClassifierFunctiononClust','PopulationAnalysisTemp','ThisisSinaPC',1,'RunThisGLMmodel','SensoryMotorIntegerInteractMdlToSwitch');

GenerateAllFiguresPipeline('GLM','ALL', [1],0,100,[],[],{'SAMPLE_ON'},...
    'ClassifierFunctiononClust','PopulationAnalysisTemp','ThisisSinaPC',1,'RunThisGLMmodel','SensoryMotorIntegerInteractMdlFromSwitch');

% plot results
TDRana.ExtraConds={'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysisTemp','DividSpockClassifier',3,'DividSpockGLM_Cond',1};
TDRana.Conditions={''};
TDRana.StrTime={'SAMPLE_ON'};%
TDRana.Area=[1]; 
TDRana.GLMmlds={'SensoryMotorIntegerInteractMdl','SensoryMotorIntegerInteractMdlToSwitch','SensoryMotorIntegerInteractMdlFromSwitch'};

cellfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',TDRana.Area,10,100,{'Decoding'},TDRana.Conditions,TDRana.StrTime,...
    TDRana.ExtraConds{:},'RunThisGLMmodel',x),TDRana.GLMmlds,'UniformOutput',0);

%% Run fitting ARIMA model on data 
% ARIMA.ExtraConds={'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysisTemp','DividSpockClassifier',3,'DividSpockGLM_Cond',1};
% ARIMA.Conditions={''};
% ARIMA.StrTime={'SAMPLE_ON'};%
% ARIMA.Area=[1]; 
%  
% GenerateAllFiguresPipeline('Classifier','ALL',ARIMA.Area,10,100,{'Decoding'},ARIMA.Conditions,ARIMA.StrTime,...
%     ARIMA.ExtraConds{:},'DividSpockClassifier_Cond',1:12);

%% run statistical tests with ARIMA 
ARIMA.Conditions={'Learning3D_Color_Response_Rule_Xgen_AltRule'};%'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB',...
    'Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Color_Response_Rule_Xgen_AltRule'}; 

GenerateAllFiguresPipeline('Classifier','ALL',ARIMA.Area,10,100,{'Learning'},ARIMA.Conditions,Learning.StrTime,...
    ARIMA.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_TrlShuff_ShuffMethod','FEA','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',0,'Classifier_NrepShufperFold',150,'Classifier_NrepShuf',50,...
    'ntrlPerCondTestArea',nan,'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',nan,'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',1);

%% Mixed Selectivity Analysis
GenerateAllFiguresPipeline('GLM','ALL', [1],0,100,[],[],{'SAMPLE_ON'},...
    'ClassifierFunctiononClust','PopulationAnalysisTemp','ThisisSinaPC',1,'RunThisGLMmodel','SensoryCatBalancedReward');

% generate area summery 
GenerateAllFiguresPipeline('GLM','ALL', [1],8,100,[],[],{'SAMPLE_ON'},...
    'ClassifierFunctiononClust','PopulationAnalysisTemp','ThisisSinaPC',1,'RunThisGLMmodel','SensoryCatBalancedReward')

%% Control for response  
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,'ReWriteClassifierData',1,...
    'ClassifierFunctiononClust','PopulationAnalysis'};
ColorShapeEncoding.Conditions={'3D_Color_Shape_Response_BalRespDir'};
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=2:5; 

% step 2 run main conditions
GenerateAllFiguresPipeline('Classifier','ALL', ColorShapeEncoding.Area,0,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,ColorShapeEncoding.ExtraConds{:});

%% Run Linear Non linear analysis 
ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ThisisSinaPC',1,'SweepClassifierConds',0,'ClassifierFunctiononClust','PopulationAnalysisTemp'};

GenerateAllFiguresPipeline('Classifier','ALL', [1],0,100,{'FineTune'},{'3D_FineTuneColorCategoryClassifier'},...
    {'SAMPLE_ON'},ExtraConds{:});


%% %%%%%%%%%%%%%%%%%%%%%%%%%<<<<<<<<<<<SUPP secion:Classifier Specific analysis>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%

%% run Response Control analysis(checks for timing of response across two axis of response
ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,'SweepClassifierConds',0,'RunMode2',3};
ResponseControl.Conditions={'3D_Color_Response_XgenR2R3','3D_Color_Response_XgenR1R2R3BalInCong','3D_Color_Response_XgenR1R2R3BalInCongV2'};
ResponseControl.Area=1:5;
ResponseControl.StrTime={'SAMPLE_ON','SACCADE_START'};
GenerateAllFiguresPipeline('Classifier','ALL', ResponseControl.Area,0,100,{'Decoding'},ResponseControl.Conditions,ResponseControl.StrTime,ExtraConds{:});
% Collect the figures generated by this analysis 
GenerateAllFiguresPipeline('Classifier','ALL', ResponseControl.Area,7,100,[{'Decoding'}],ResponseControl.Conditions,...
    ResponseControl.StrTime,'Page2SaveNum',{[1 4],[1 4],[1 4]});


%% Run Fine Tune analysis
ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ThisisSinaPC',1,'SweepClassifierConds',0,'ClassifierFunctiononClust','PopulationAnalysisTemp'};

GenerateAllFiguresPipeline('Classifier','ALL', [1],0,100,{'FineTune'},{'3D_FineTuneColorCategoryClassifier'},...
    {'SAMPLE_ON'},ExtraConds{:});


GenerateAllFiguresPipeline('Classifier','ALL', [4],2,100,{'FineTune'},{'3D_ColorCategoryFineTuneSig'},...
    {'SAMPLE_ON'},ExtraConds{:});
GenerateAllFiguresPipeline('Classifier','ALL', [1:5],0,100,{'FineTune'},{'3D_Color_Response_Xgen'},...
    {'SACCADE_START'},ExtraConds{:});
GenerateAllFiguresPipeline('Classifier','ALL', [1],0,100,{'FineTune'},{'3D_Color_Response_XgenBalCong','3D_Color_Response_XgenBalInCong'},...
    {'SAMPLE_ON','SACCADE_START'},ExtraConds{:});

GenerateAllFiguresPipeline('Classifier','Silas', [1],0,100,{'FineTune'},{'3D_Color_Response_Xgen'},...
    {'SAMPLE_ON','SACCADE_START'},ExtraConds{:});
GenerateAllFiguresPipeline('Classifier','Chico', [1],0,100,{'FineTune'},{'3D_Color_Cat_Xgen'},...
    {'SAMPLE_ON','SACCADE_START'},ExtraConds{:});

GenerateAllFiguresPipeline('Classifier','ALL', [1],0,100,{'FineTune'},{'ResponseLocation_BalCorr'},...
    {'SAMPLE_ON','SACCADE_START'},ExtraConds{:});

GenerateAllFiguresPipeline('Classifier','ALL', [1],0,100,{'FineTune'},{'3D_Color_Response_XgenR1R2R3BalInCong'},...
    {'SAMPLE_ON','SACCADE_START'},ExtraConds{:});
%% sweep classifier conditions
ExtraConds={'Classifier_Nrep',500,'DividSpockClassifier',1,'RunDummyFile',0,'SweepClassifierConds',1};
GenerateAllFiguresPipeline('Classifier','ALL', [1:5],0,100,{'Learning'},{'3D_Shared_Color_Response_Entropy'},{'SAMPLE_ON'},ExtraConds{:});

%% %%%%%%%%%%%%%%%%%%%%%%%%%<<<<<<<<<<<subspace geomerty>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate subspace geomerty results
GenerateAllFiguresPipeline('Subspace','ALL', [1:5],3,100,[{'SubspaceLearning'}],[],{'SAMPLE_ON'});
GenerateAllFiguresPipeline('Subspace','ALL', [1:5],3,100,[{'SubspaceSensMotTransform'}],[],{'SAMPLE_ON'});



%% extra PC runs

Entropy.ExtraConds={'Classifier_Nrep',2000,'DividSpockClassifier',3,'RunDummyFile',0,'SweepClassifierConds',0,'ClassifierFunctiononClust','PopulationAnalysisTemp'};
Entropy.Conditions={'3D_Shared_Color_Response_EntropyR3','3D_Shared_Color_Response_EntropyR2','3D_Shared_Color_Response_EntropyR3Bal','3D_Shared_Color_Response_EntropyR2Bal'};
Entropy.Area=1;Entropy.StrTime={'SACCADE_START'};
NTrlRngTrainLearningAreaSet=[-50]; % -75 800 number of trials to take from end of the block for training per area
NTrlRngTestLearningAreaSet=[-50]; % -75 800 number of trials shift for Test
ntrlPerCondAreaSet=[3]; % nu


PopulationAnalysisTemp(0, 0, [], [1], '3D_Shared_Color_Response_EntropyR2Bal', 'ProcessingStep',8,'PopulationAna.PSTHbin',100,'DividSpockClassifier_TrlRng',...
    [16 1 1 1  1 1 1 1],'DividSpockClassifier_Cond',1,'CalShuffleClassifier',0,'DividSpockClassifier',3,'CalShuffTrlOrderClassifier',0,'SpkCntStartFieldName',...
    'SACCADE_START','NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',3)



for i=1:4

PopulationAnalysisTemp(0, 0, [], [1], Entropy.Conditions{i}, 'ProcessingStep',4,'PopulationAna.PSTHbin',100,'DividSpockClassifier_TrlRng',...
    [16 1 1 1  1 1 1 1],'DividSpockClassifier_Cond',1,'CalShuffleClassifier',0,'DividSpockClassifier',3,'CalShuffTrlOrderClassifier',0,'SpkCntStartFieldName',...
    'SACCADE_START','NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',3)
end

