% Generates figures for behavioral analysis for 
% Tafazoli et al, Building compositional tasks with shared neural subspaces
% https://www.biorxiv.org/content/10.1101/2024.01.31.578263v1.article-metrics


%% Figure 1 
% behavior
GenerateBhvPaperFigs % generate Fig. 1e-f and Fig. 4a

% GLM  % Change 'ProcessingStep', 11 if you want to run GLM fits
Generatefigures(0,[0],[],1,[],'ProcessingStep',16,'SpkCntStartFieldName','SAMPLE_ON');% generate Fig. 1h-m

%% Figure 2 % Change 'ProcessingStep', 3 if you want to run classifier fits
Generatefigures(0,[0],[],1,'3D_Color_Response_XgenBalInCongV5','ProcessingStep',4,'SpkCntStartFieldName','SAMPLE_ON'); % generate Fig. for reponse direction and color category (Fig. 2a,c,f,i) 
Generatefigures(0,[0],[],1,'3D_Shape_Cat_Xgen','ProcessingStep',4,'SpkCntStartFieldName','SAMPLE_ON');% generate Fig. for shape category(Fig. 2b)
Generatefigures(0,[0],[],1,'3D_Cat_Color_XgenCol_Compare','ProcessingStep',5,'SpkCntStartFieldName','SAMPLE_ON');% generate Fig. 2g

%% Figure 3 
Generatefigures(0,[0],[],1,'3D_Cat_Color_Response_Compare','ProcessingStep',5,'SpkCntStartFieldName','SAMPLE_ON');% generate Fig. 3a

% generate figures Fig. 3c-e   
Generatefigures(0, 0, [], [1], '3D_Shared_Color_Response_EntropyR3Bal','ProcessingStep',4,'SpkCntStartFieldName','SACCADE_START','CalShuffTrlOrderClassifier',0,'RunCrossTemporalClassifer',0,'CalShuffleClassifier',0,...
'PopulationAna.PSTHbin',100,'SpkCntStartFieldName','SACCADE_START','DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',16*ones(1,7),'DividSpockClassifier',3,'NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',4)

%% Figure 4   % Change 'ProcessingStep', 3 if you want to run the classifier fits
Generatefigures(0,[0],[],1,'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','ProcessingStep',4,'SpkCntStartFieldName','SAMPLE_ON');  % generate Fig. 4c,d,e,f,g,h,i 
Generatefigures(0,[0],[],1,'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB','ProcessingStep',4,'SpkCntStartFieldName','SAMPLE_ON'); % generate Fig. 4d,g,j
Generatefigures(0,[0],[],1,'Learning3D_Color_Response_Rule_Xgen_AltRule','ProcessingStep',4,'SpkCntStartFieldName','SAMPLE_ON');  % generate Fig. 4k,l
Generatefigures(0,[0],[],1,'Learning3D_Color_Response_Rule_Xgen_SameRule','ProcessingStep',4,'SpkCntStartFieldName','SAMPLE_ON'); % generate Fig. 4l

%% Figure 5   % Change 'ProcessingStep', 3 if you want to run the classifier fits
Generatefigures(0,[0],[],1,'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','ProcessingStep',4,'SpkCntStartFieldName','SAMPLE_ON');  % generate Fig. 5a,b,d,e,f
Generatefigures(0,[0],[],1,'Learning3D_Shape_Color_Color_Compression_RB','ProcessingStep',4,'SpkCntStartFieldName','SAMPLE_ON');  % generate Fig. 5c


