% scripts to run stat test processing
% run Arima  Important numbers 1 2 8 9 10 11 13 14 16 17
ShuffMethod={           'FEA','RND','FEA','FEA','FEA','FEA','RND','FEA','FEA','FEA','FEA','FEA','RND','RND','FEA','RND','RND',  'FEA','FEA','FEA'};
Num=                     [1    2      3     4     5     6     7     8     9     10    11    12   13     14    15   16     17      18    19    20  ];
NTrlRngTestLearningArea =[nan nan    40    50     60    60    60    40    40    35    45    35   40     40    40   45     35,     40    40    40  ];
ntrlPerCondTestArea     =[nan nan     1     2      2     3     3     1    1     1     1     1    1      1     1    1      1,      1     1     1   ];
IncludeAllClassifierInfo=[1   1       0     0      0     0     1     1    1     1     1     1    1      1     1    1      1,      0     0     0   ];
MaxTrlTestLearningArea  =[nan*ones(1,7)                             110   115  115   115   115   115    110   115  120   120,     115   115   115 ];
LimitFromSwitchPerfArea=[nan*ones(1,11)                                                    0.55  nan    nan   nan  nan   nan      nan   nan   nan ];
CalShuffTrlOrderClassifier=ones(1,length(NTrlRngTestLearningArea));
ClassifierNameExtra=[cell(1,length(NTrlRngTestLearningArea))];%'NQ'
EqualizeTrialsXConds=[ones(1,17) zeros(1,3)];
SCR_ALT='Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';
SCR_SAME='Learning3D_Shape_Color_Rule_Xgen_SameRule_RB';
CRR_ALT='Learning3D_Color_Response_Rule_Xgen_AltRule';
CRR_SAME='Learning3D_Color_Response_Rule_Xgen_SameRule';
ClassifierNames=[repmat({SCR_ALT},1,17) SCR_SAME CRR_ALT CRR_SAME];
 
NConds=1; % which conditions you want to check    
ProcessingStep=17; % what processing step we want to run 
for CNum=19 % what classifier do you want to work with
    fprintf(2,'\n Running ****ShuffMethod=%s NTrlRngTestLearningArea=%i ntrlPerCondTestArea=%i \n IncludeAllClassifierInfo=%i MaxTrlTestLearningArea=%i LimitFromSwitchPerfArea=%i \n CalShuffTrlOrderClassifier=%i ClassifierNameExtra=%s',...
        ShuffMethod{CNum},NTrlRngTestLearningArea(CNum),ntrlPerCondTestArea(CNum),IncludeAllClassifierInfo(CNum),MaxTrlTestLearningArea(CNum),LimitFromSwitchPerfArea(CNum),CalShuffTrlOrderClassifier(CNum),ClassifierNameExtra{CNum});
    % concatinate shuffle trl order
    PopulationAnalysisTemp(0, 0, [], [1], ClassifierNames{CNum}, ...
        'IncludeAllClassifierInfo',1,'ProcessingStep',ProcessingStep,'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
        'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',EqualizeTrialsXConds(CNum),'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
        'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',150,'Classifier_NrepShuf',50,...
        'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
        'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
        'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',NConds,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum),...
        'ClassifierNameExtra',ClassifierNameExtra{CNum});
end
%% Run Arima Stat Test on Data
ARIMA.ExtraConds={'EnforceSpockTime',120,'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysisTemp','DividSpockClassifier',3,'DividSpockGLM_Cond',1};
ARIMA.Conditions={''};
ARIMA.StrTime={'SAMPLE_ON'}; 
ARIMA.Area=[1];
ARIMA.NConds=[1];
ARIMA.NConds_PC=1; % which conditions are we running on PC

%ARIMA.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','Learning3D_Color_Response_Rule_Xgen_AltRule','Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};
%'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB',...
%'Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Color_Response_Rule_Xgen_AltRule'}; 

for CNum=[1]
     fprintf(2,'\n Running ****ShuffMethod=%s NTrlRngTestLearningArea=%i ntrlPerCondTestArea=%i \n IncludeAllClassifierInfo=%i MaxTrlTestLearningArea=%i LimitFromSwitchPerfArea=%i \n CalShuffTrlOrderClassifier=%i ClassifierNameExtra=%s',...
        ShuffMethod{CNum},NTrlRngTestLearningArea(CNum),ntrlPerCondTestArea(CNum),IncludeAllClassifierInfo(CNum),MaxTrlTestLearningArea(CNum),LimitFromSwitchPerfArea(CNum),CalShuffTrlOrderClassifier(CNum),ClassifierNameExtra{CNum});

    GenerateAllFiguresPipeline('Classifier','ALL',ARIMA.Area,10,100,{'Learning'},ClassifierNames(CNum),ARIMA.StrTime,...
        ARIMA.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
        'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',EqualizeTrialsXConds(CNum),'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
        'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',150,'Classifier_NrepShuf',50,...
        'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
        'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
        'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',ARIMA.NConds,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum),...
        'ClassifierNameExtra',ClassifierNameExtra{CNum});
end

%  if you are running on PC
PopulationAnalysisTemp(0, 0, [], [1], ARIMA.Conditions{1}, ...
    'ProcessingStep',17, 'SpkCntStartFieldName','SAMPLE_ON','PopulationAna.PSTHbin',100,...
    ARIMA.ExtraConds{:},'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',150,'Classifier_NrepShuf',50,...
    'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
    'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',ARIMA.NConds_PC);

%% run shuffle order classifier results or normal
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysis'};%,'Enforce_nCondsRun_Shuff',{[1001]}};
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;

Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB'};

for CNum=17
    % N trial of test range and n trls of test
    GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
        Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',CalShuffTrlOrderClassifier(CNum),'RunMode2',2,...
        'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',1,'EqualizeTrialsXConds',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
        'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',100,'Classifier_NrepShuf',50,...
        'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
        'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
        'NTrlStpLearningArea',nan,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum));
end

% to run same shapecolor and response conds
%% run shuffle order classifier results or normal
 
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','Learning3D_Color_Response_Rule_Xgen_AltRule',...
    'Learning3D_Color_Response_Rule_Xgen_SameRule',...
    'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};% 'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB_EQ',

%% run shuffle
% N trial of test range and n trls of test  
GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,9,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',CalShuffTrlOrderClassifier(CNum),'RunMode2',2,...
    'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',1,'EqualizeTrialsXConds',0,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',100,'Classifier_NrepShuf',50,...
    'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
    'NTrlStpLearningArea',nan,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum));

%% run MeanSTD
Learning.Conditions={'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','Learning3D_Color_Response_Rule_Xgen_AltRule',...
    'Learning3D_Color_Response_Rule_Xgen_SameRule','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB'};% 'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB_EQ',

CalShuffTrlOrderClassifier=zeros(1,length(LimitFromSwitchPerfArea));
IncludeAllClassifierInfo=zeros(1,12);

GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,0,100,{'Learning'},Learning.Conditions,Learning.StrTime,...
    Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',0,'RunMode2',2,...
    'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',1,'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',100,'Classifier_NrepShuf',50,...
    'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
    'NTrlStpLearningArea',nan,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum));


%% run main conditions 

Classifier_NrepShufperFold=250;
Classifier_NrepShuf=50;
ShuffMethod={            'RND','RND','RND','RND','RND','RND','RND','RND','FEA'};
NAllConds=length(ShuffMethod);
Num=                      [1    2      3     4     5    6     7      8     9];
ClassifierType=           [1    1      1     1     2    3     4      3     1];
NTrlRngTestLearningArea = [45   50     45    50    45   45    45     45    45];
ntrlPerCondTestArea     =[ones(1,NAllConds)];
IncludeAllClassifierInfo=[ones(1,sum(ClassifierType==1)) 0 0 0 1 1];
MaxTrlTestLearningArea  =[115  115    120    125   115  115   115 115 115];
Classifier_TrlShuff_GenClassifierSpecs=ones(1,NAllConds);
LimitFromSwitchPerfArea=[nan*ones(1,NAllConds)];
CalShuffTrlOrderClassifier=[nan*ones(1,NAllConds)];
ClassifierNameExtra=[cell(1,NAllConds-1) 'EX'];
SCR_ALT='Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';
SCR_SAME='Learning3D_Shape_Color_Rule_Xgen_SameRule_RB';
CRR_ALT='Learning3D_Color_Response_Rule_Xgen_AltRule';
CRR_SAME='Learning3D_Color_Response_Rule_Xgen_SameRule';
ClassifierNames=[repmat({SCR_ALT},1,4) SCR_SAME CRR_ALT CRR_SAME CRR_ALT SCR_ALT];
EqualizeTrialsXConds=[ones(1,4) 0 0 0 0 1];
%% extra conditions
Classifier_NrepShufperFold=250;
Classifier_NrepShuf=50;
ShuffMethod={            'FEA','RND','RND','RND','RND'};
NAllConds=length(ShuffMethod);
Num=                      [1    2      3    4   5];
ClassifierType=           [1    1      1    1   1];
NTrlRngTestLearningArea = [45   45     40   45  50];
ntrlPerCondTestArea     =[ones(1,NAllConds)];
IncludeAllClassifierInfo=[ones(1,sum(ClassifierType==1)) ];
MaxTrlTestLearningArea  =[110  110    110   120  nan];
Classifier_TrlShuff_GenClassifierSpecs=[1 1 0 0 0];
LimitFromSwitchPerfArea=[nan*ones(1,NAllConds)];
CalShuffTrlOrderClassifier=[nan*ones(1,NAllConds)];
ClassifierNameExtra=[cell(1,NAllConds) ];
SCR_ALT='Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';
SCR_SAME='Learning3D_Shape_Color_Rule_Xgen_SameRule_RB';
CRR_ALT='Learning3D_Color_Response_Rule_Xgen_AltRule';
CRR_SAME='Learning3D_Color_Response_Rule_Xgen_SameRule';
ClassifierNames=[repmat({SCR_ALT},1,NAllConds) ];
EqualizeTrialsXConds=[ones(1,NAllConds)];
Enforce_nCondsRun_Shuff={[],[],[101:251],[101:251],[165:251],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]};

Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysisTemp'};
%,'RunWithDepenency',1,'job_id_dep','36410229'};%,'Enforce_nCondsRun_Shuff',{[1001]}};
Learning.StrTime={'SAMPLE_ON'};
Learning.Area=1;
DividSpockClassifier_Cond=1;
Learning.Conditions=ClassifierNames;
PCProcessingStep=17;
SpockTaskCond=10; % 10 for TDR etc and 9 for trial shuffle
if SpockTaskCond==10;Learning.ExtraConds=[Learning.ExtraConds 'EnforceSpockTime',60];end

if SpockTaskCond==0 | SpockTaskCond==3 | SpockTaskCond==2
    CalShuffTrlOrderClassifier=zeros(1,100);
    RunMode2=1;
else
    CalShuffTrlOrderClassifier=ones(1,100);
end
 
for CNum=[1:3]
     fprintf(2,'\n Running %s ****ShuffMethod=%s NTrlRngTestLearningArea=%i ntrlPerCondTestArea=%i \n IncludeAllClassifierInfo=%i MaxTrlTestLearningArea=%i LimitFromSwitchPerfArea=%i \n CalShuffTrlOrderClassifier=%i ClassifierNameExtra=%s',...
        ClassifierNames{CNum},ShuffMethod{CNum},NTrlRngTestLearningArea(CNum),ntrlPerCondTestArea(CNum),IncludeAllClassifierInfo(CNum),MaxTrlTestLearningArea(CNum),LimitFromSwitchPerfArea(CNum),CalShuffTrlOrderClassifier(CNum),ClassifierNameExtra{CNum});

    GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,SpockTaskCond,100,{'Learning'},Learning.Conditions(CNum),Learning.StrTime,...
        Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',1,'RunMode2',2,...
        'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',Classifier_TrlShuff_GenClassifierSpecs(CNum),...
        'EqualizeTrialsXConds',EqualizeTrialsXConds(CNum),'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
        'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',Classifier_NrepShufperFold,'Classifier_NrepShuf',Classifier_NrepShuf,...
        'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
        'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
        'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',DividSpockClassifier_Cond,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum),...
        'ClassifierNameExtra',ClassifierNameExtra{CNum},'Enforce_nCondsRun_Shuff',Enforce_nCondsRun_Shuff(CNum));
end


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
IncludeAllClassifierInfo=[ 1    1      1    1       0];%[ones(1,NAllConds)]; %

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

SpockTaskCond=2; % 10 for TDR etc and 9 for trial shuffle
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

for CNum=[3]
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
for CNum=[1]
PopulationAnalysisTemp(0, 0, [], Learning.Area, ClassifierNames{CNum}, ...
    'IncludeAllClassifierInfo',1,'ProcessingStep',PCProcessingStep,'SweepClassifierConds',0,'CalShuffTrlOrderClassifier',CalShuffTrlOrderClassifier(CNum),'RunMode2',2,...
    'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',0,'EqualizeTrialsXConds',EqualizeTrialsXConds(CNum),'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
    'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',150,'Classifier_NrepShuf',50,...
    'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
    'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
    'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',DividSpockClassifier_Cond,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum),...
    'ClassifierNameExtra',ClassifierNameExtra{CNum});
end

%% run Rb shuffle for learning for 250 permutations 
SpockTaskCond=9; % 10 for TDR etc and 9 for trial shuffle
CalShuffTrlOrderClassifier=zeros(1,100);
CalShuffleClassifier=ones(1,100);
RunMode2=2;
PCProcessingStep=3; % to plot

for CNum=[1:3]
     fprintf(2,'\n Running %s ****ShuffMethod=%s NTrlRngTestLearningArea=%i ntrlPerCondTestArea=%i \n IncludeAllClassifierInfo=%i MaxTrlTestLearningArea=%i LimitFromSwitchPerfArea=%i \n CalShuffleClassifier=%i CalShuffTrlOrderClassifier=%i ClassifierNameExtra=%s',...
        ClassifierNames{CNum},ShuffMethod{CNum},NTrlRngTestLearningArea(CNum),ntrlPerCondTestArea(CNum),IncludeAllClassifierInfo(CNum),MaxTrlTestLearningArea(CNum),LimitFromSwitchPerfArea(CNum),CalShuffleClassifier(CNum),CalShuffTrlOrderClassifier(CNum),ClassifierNameExtra{CNum});

    GenerateAllFiguresPipeline('Classifier','ALL',Learning.Area,SpockTaskCond,100,{'Learning'},Learning.Conditions(CNum),Learning.StrTime,...
        Learning.ExtraConds{:},'SweepClassifierConds',0,'CalShuffleClassifier',CalShuffleClassifier(CNum),'CalShuffTrlOrderClassifier',CalShuffTrlOrderClassifier(CNum),'RunMode2',RunMode2,...
        'Classifier_TrlShuff_TrialSet','FL','Classifier_TrlShuff_GenClassifierSpecs',Classifier_TrlShuff_GenClassifierSpecs(CNum),...
        'EqualizeTrialsXConds',EqualizeTrialsXConds(CNum),'CalShuffTrlOrderClassifier_RunAllTrlRng',1,...
        'Subtract4TrialRange',0,'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum),'Classifier_NrepShufperFold',Classifier_NrepShufperFold,'Classifier_NrepShuf',Classifier_NrepShuf,...
        'ntrlPerCondTestArea',ntrlPerCondTestArea(CNum),'Classifier_TrlShuff_ShuffMethod',ShuffMethod{CNum},'ntrlPerCondArea',nan,'LimitFromSwitchPerfArea',nan,'NTrlRngTrainLearningArea',nan,...
        'NTrlRngTestLearningArea',NTrlRngTestLearningArea(CNum),'MaxTrlTestLearningArea',MaxTrlTestLearningArea(CNum),...
        'NTrlStpLearningArea',nan,'DividSpockClassifier_Cond',DividSpockClassifier_Cond,'LimitFromSwitchPerfArea',LimitFromSwitchPerfArea(CNum),...
        'ClassifierNameExtra',ClassifierNameExtra{CNum});
end


%% BlkSeq Analysis for Rule 2
ColorShapeEncoding.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'ClassifierFunctiononClust','PopulationAnalysis',...
    'Classifier_NrepShuf',100,'Classifier_NrepShufperFold',1000,'CalShuffTrlOrderClassifier',0,...
    'EnforceSpockTime',400};%,'Enforce_nCondsRun',{2}};

% Learning3D_RuleR2SeqAnalysisLearningCtrlR2onR2
ColorShapeEncoding.Conditions={'Learning3D_RuleR2SeqAnalysisLearningCtrlR2','Learning3D_RuleR2SeqAnalysisLearningCtrl'};
IncludeAllClassifierInfo=[1 1 0];
ColorShapeEncoding.StrTime={'SAMPLE_ON'};
ColorShapeEncoding.Area=1; 

% Run on spock
arrayfun(@(CNum) GenerateAllFiguresPipeline('Classifier','ALL',ColorShapeEncoding.Area,3,100,{'Decoding'},ColorShapeEncoding.Conditions(CNum),ColorShapeEncoding.StrTime,...
    ColorShapeEncoding.ExtraConds{:},'IncludeAllClassifierInfo',IncludeAllClassifierInfo(CNum)),1:3);

% Run on PC
for i=1:2
    PopulationAnalysisTemp(0, 0, [], 1, ColorShapeEncoding.Conditions{i}, 'ProcessingStep',4, ...
        'IncludeAllClassifierInfo',IncludeAllClassifierInfo(i), ...
        'MaxTrlTestLearning',nan,'NTrlStpLearningArea',nan,'CalShuffleClassifier',0,'CalShuffTrlOrderClassifier',0, ...
        'Classifier_TrlShuff_ShuffMethod','FEA','ntrlPerCondTestArea',nan,'CalShuffTrlOrderClassifier_RunAllTrlRng',1, ...
        'Subtract4TrialRange',0,'ClassifierNameExtra','','Classifier_TrlShuff_TrialSet','FL','NeuResamplingSubPopSize',110, ...
        'ntrlPerCondArea',nan,'NTrlRngTrainLearningArea',nan,'NTrlRngTestLearningArea',nan,'PopulationAna.PSTHbin',100, ...
        'DividSpockClassifier_TrlRng',[1 1 1 1  1 1 1 1],'DividSpockClassifier_Cond',1,'SpkCntStartFieldName','SAMPLE_ON')
end

%% do test for problem with saving ?
SpockTaskCond=10; % to plot
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysisTemp','EnforceSpockTime',5,'NcoresSpock',10};

CNum=2;
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


%% look at resutls for Axis decocidng 
% we have two conditions CalShuffTrlOrderClassifier_RunAllTrlRng=5 which is
% only for 1-10 and 7:16 and CalShuffTrlOrderClassifier_RunAllTrlRng=0 which is for 1-10 and 64:73
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust','PopulationAnalysisTemp'};%,'Enforce_nCondsRun_Shuff',{[1001]}};

Classifier_NrepShufperFold=250;
Classifier_NrepShuf=50;
Learning.StrTime={'SAMPLE_ON'};
Learning.Conditions={'3D_AxisDecodingSupp'};%'Learning3D_AxisDecoding'};
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

%% shuffle trial order 

Col_Same_ALT='Learning3D_Color_Same_Alt_Rule_Stat_RB';
Classifier_NrepShufperFold=100;
Classifier_NrepShuf=50;
ShuffMethod={'FEA'};
NAllConds=length(ShuffMethod);
Num=                      [1   ];
ClassifierType=           [1   ];
NTrlRngTestLearningArea = [50   ];
ntrlPerCondTestArea     =[ones(1,NAllConds)];
IncludeAllClassifierInfo=[ 0  ];
MaxTrlTestLearningArea  =[125   ];
Area                    =[4    ];
ClassifierNames=   [ {Col_Same_ALT} ];
Classifier_TrlShuff_GenClassifierSpecs=zeros(1,NAllConds);
LimitFromSwitchPerfArea=[nan*ones(1,NAllConds)];
CalShuffTrlOrderClassifier=[ones(1,NAllConds)];
ClassifierNameExtra=[cell(1,NAllConds)];
EqualizeTrialsXConds=[0 ];

CodeTypeNum=input('Which code are we using Temp(1)/Final(2)');
CodeType={'PopulationAnalysisTemp','PopulationAnalysis'};
Learning.ExtraConds={'Classifier_Nrep',250,'DividSpockClassifier',3,'RunDummyFile',0,...
    'ReWriteClassifierData',1,'ThisIsSinaPC',1,'pvalClassifierAnalysis_ClusterCorrect',0.001,...
    'ClassifierFunctiononClust',CodeType{CodeTypeNum}};%,'Enforce_nCondsRun_Shuff',{[228,245]}};

Learning.StrTime={'SAMPLE_ON'};
Learning.Area=4;
DividSpockClassifier_Cond=1;
Learning.Conditions=ClassifierNames;

SpockTaskCond=9; % 10 for TDR etc and 9 for trial shuffle
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

for CNum=[1]
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