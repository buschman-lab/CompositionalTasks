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

%% plotting single cell results
% plot all GLM result summery for each area(PS16)
GenerateAllFiguresPipeline('GLM','ALL', [1:5],7,100,[],[],{'SAMPLE_ON'})

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%<<<<<<<<<<<Section 2:Compositionality of Reperesentations>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%<<<<<<<<<<<Section 2-1:Color & Shape Information>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot Main figure (figure 2)
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


%% Comparision of encoding of color across areas
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


%% plot entropy results 
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Entropy.Area,2,100,{'Decoding'},Entropy.Conditions,Entropy.StrTime,...
    Entropy.ExtraConds{:},'NTrlRngTrainLearningArea',NTrlRngTrainLearningAreaSet(x),...
    'NTrlRngTestLearningArea',NTrlRngTestLearningAreaSet(x),'ntrlPerCondArea',ntrlPerCondAreaSet(x)),1:length(ntrlPerCondAreaSet));

PopulationAnalysisTemp(0, 0, [], [1], [11],'ProcessingStep',4,'SpkCntStartFieldName','SACCADE_START','CalShuffTrlOrderClassifier',0,'RunCrossTemporalClassifer',0,'CalShuffleClassifier',0,...
'PopulationAna.PSTHbin',100,'SpkCntStartFieldName','SACCADE_START','DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',16*ones(1,7),'DividSpockClassifier',3,'NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',4)

% plot entropy for congruency results '3D_Shared_Color_Response_EntropyR3Bal'
PopulationAnalysisTempCongEntropy(0, 0, [], [1], [11],'ProcessingStep',4,'SpkCntStartFieldName','SACCADE_START','CalShuffTrlOrderClassifier',0,'RunCrossTemporalClassifer',0,'CalShuffleClassifier',0,...
'PopulationAna.PSTHbin',100,'SpkCntStartFieldName','SACCADE_START','DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',16*ones(1,7),'DividSpockClassifier',3,'NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',4)


%% plot learning analysis
%% run main conditions for main figures
SCR_ALT='Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';
SCR_SAME='Learning3D_Shape_Color_Rule_Xgen_SameRule_RB_NoColDup';
CRR_ALT='Learning3D_Color_Response_Rule_Xgen_AltRule';
CRR_SAME='Learning3D_Color_Response_Rule_Xgen_SameRule';

Classifier_NrepShufperFold=250;
Classifier_NrepShuf=50;
ShuffMethod={            'RND','RND','RND','RND','RND' };
NAllConds=length(ShuffMethod);
Num=                      [1    2      3     4     5   ];
ClassifierType=           [1    1      2     3     4   ];
NTrlRngTestLearningArea = [45   40     40    40    40  ];
ntrlPerCondTestArea     =[ones(1,NAllConds)];
IncludeAllClassifierInfo=[ 1    1      1    1       1];%[ones(1,NAllConds)]; %

MaxTrlTestLearningArea  =[120  110    110    110   110 ];
Area                    =[1     1      1      1     1  ];
ClassifierNames=   [ {SCR_ALT}  SCR_ALT SCR_SAME CRR_ALT CRR_SAME ];
Classifier_TrlShuff_GenClassifierSpecs=ones(1,NAllConds);
LimitFromSwitchPerfArea=[nan*ones(1,NAllConds)];
CalShuffTrlOrderClassifier=[nan*ones(1,NAllConds)];
ClassifierNameExtra=[cell(1,NAllConds)];
EqualizeTrialsXConds=[1 1 1 0 0 ];

CodeTypeNum=input('Which code are we using Temp(1)/Final(2)');
CodeType={'PopulationAnalysisTemp','PopulationAnalysis'};
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust',CodeType{CodeTypeNum}};%,'Enforce_nCondsRun_Shuff',{[228,245]}};

Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;
DividSpockClassifier_Cond=1;
Learning.Conditions=ClassifierNames;

SpockTaskCond=2; % 2 for plotting 10 for TDR and processing trial shuffle and 9 for running trial shuffle calculations
if SpockTaskCond==10;Learning.ExtraConds=[Learning.ExtraConds 'EnforceSpockTime',240];end

if SpockTaskCond==0 | SpockTaskCond==3 | SpockTaskCond==2
    CalShuffTrlOrderClassifier=zeros(1,100);
    RunMode2=1;
elseif SpockTaskCond==8
    RunMode2=input('What is the type of data(1:Main,2:Shuffle)');
    CalShuffTrlOrderClassifier=(RunMode2-1)*ones(1,100);
else
    CalShuffTrlOrderClassifier=ones(1,100);
    RunMode2=2;
end
PCProcessingStep=4; % to plot

for CNum=[1:5]
    fprintf(2,'\n Running %s ****ShuffMethod=%s NTrlRngTestLearningArea=%i ntrlPerCondTestArea=%i \n IncludeAllClassifierInfo=%i MaxTrlTestLearningArea=%i LimitFromSwitchPerfArea=%i \n CalShuffTrlOrderClassifier=%i ClassifierNameExtra=%s',...
        ClassifierNames{CNum},ShuffMethod{CNum},NTrlRngTestLearningArea(CNum),ntrlPerCondTestArea(CNum),IncludeAllClassifierInfo(CNum),MaxTrlTestLearningArea(CNum),LimitFromSwitchPerfArea(CNum),CalShuffTrlOrderClassifier(CNum),ClassifierNameExtra{CNum});

    GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,SpockTaskCond,100,{'Learning'},Learning.Conditions(CNum),Learning.StrTime,...
        Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',CalShuffTrlOrderClassifier(CNum),'RunMode2',RunMode2,...
        'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',Classifier_TrlShuff_GenClassifierSpecs(CNum),...
        'EqualizeTrialsXConds',EqualizeTrialsXConds(CNum),'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
        'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',Classifier_NrepShufperFold,'Classifier_NrepShuf',Classifier_NrepShuf,...
        'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
        'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
        'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',DividSpockClassifier_Cond,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum),...
        'ClassifierNameExtra',ClassifierNameExtra{CNum});
end

% concatinate shuffle trl order
PopulationAnalysisTemp(0, 0, [], Learning.Area, ClassifierNames{CNum}, ...
    'IncludeAllClassifierInfo',1,'ProcessingStep',PCProcessingStep,'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',CalShuffTrlOrderClassifier(CNum),'RunMode2',2,...
    'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',EqualizeTrialsXConds(CNum),'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',150,'Classifier_NrepShuf',50,...
    'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
    'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',DividSpockClassifier_Cond,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum),...
    'ClassifierNameExtra',ClassifierNameExtra{CNum});


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%<<<<<Section 4:Compression analysis>>>>>>>>>>%%%%%%%%%%%%%%%%%%%%%%%%%
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysisTemp'};
Learning.Conditions={'Learning3D_Shape_Color_Color_Compression_RB'};
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'ntrlPerCondArea',5);

% plot in pc
PopulationAnalysisTemp(0, 0, [], Learning.Area, Learning.Conditions{1},'ProcessingStep',4,'SpkCntStartFieldName',Learning.StrTime{1},'CalShuffTrlOrderClassifier',0,'RunCrossTemporalClassifer',0,'CalShuffleClassifier',0,...
'PopulationAna.PSTHbin',100,'DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',16*ones(1,7),'DividSpockClassifier',1,'NTrlRngTrainLearningArea',nan,'NTrlRngTestLearningArea',nan,'ntrlPerCondArea',5)

% look at other areas
Learning.Area=2:5;
arrayfun(@(x) GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,2,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'ntrlPerCondArea',x),[5],'UniformOutput',0); % trials could be [12 9 5]


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

PopulationAnalysisTemp(0, 0, [], [ColorShapeEncoding.Area], ColorShapeEncoding.Conditions{1}, 'ProcessingStep',4,'PopulationAna.PSTHbin',100,'DividSpockClassifier_TrlRng',[1],'DividSpockClassifier_Cond',1,'CalShuffleClassifier',0,'CalShuffTrlOrderClassifier',0,...
    'SpkCntStartFieldName',ColorShapeEncoding.StrTime{1},'NTrlRngTrainLearningArea',NTrlRngTrainLearningAreaSet,'NTrlRngTestLearningArea',NTrlRngTestLearningAreaSet,'ntrlPerCondArea',ntrlPerCondAreaSet)

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

%% Axis decoding analysis
%% look at resutls for Axis decocidng 
% we have two conditions CalShuffTrlOrderClassifier_RunAllTrlRng=5 which is
% only for 1-10 and 7:16 and CalShuffTrlOrderClassifier_RunAllTrlRng=0 which is for 1-10 and 64:73
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysisTemp'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

Classifier_NrepShufperFold=250;
Classifier_NrepShuf=50;
Learning.StrTime={'SAMPLE_ON'};
Learning.Conditions={'3D_AxisDecodingSupp'};%'Learning3D_AxisDecoding'};%''};3D_AxisDecodingSupp
Learning.Area=1;

GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,2,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
    'Classifier_NrepShufperFold',Classifier_NrepShufperFold,'Classifier_NrepShuf',Classifier_NrepShuf,...
    'Classifier_TrlShuff_ShuffMethod','RND','Classifier_TrlShuff_TrialSet','FL',...
    'Classifier_TrlShuff_GenClassifierSpecs',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',5,...
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
%% axis suppression analysis
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysis'};%,'Enforce_nCondsRun',{2}};

ColorShapeEncoding.Conditions={'3D_AxisDecodingSupp'};
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1; 

GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,2,100,{'Decoding'},ColorShapeEncoding.Conditions,ColorShapeEncoding.StrTime,...
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
