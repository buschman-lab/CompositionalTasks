% Generating Figure 1
% Guide: In all code Rule1=S1, Rule2=C2 and Rule3=C1
function SuppFigures

%% Prepare Vaiables (Run this section before everything)
global AnalysisOpts 
AnalysisOpts.Project='Rule Representation';
AnalysisOpts.BhvAna.Kidname='ALL';
AnalysisOpts.AnalysisType='Analysis Pipeline';
KickoffMyfunc(0,0);
SetAnalysisOptions_RuleRepresentation(); % set all of the parameters
 
% Adjust the Path to locate the folder with all data *****
DataPath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Submission Code\Final Code for Nature Paper 2025 Tafazoli et al\Figure Data\';
FigSavePath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Submission Code\Final Code for Nature Paper 2025 Tafazoli et al\Figure Results\';

% Setup the class to plot the data
FigParams=fig_params;
bhvAna=BhvAnalysisFuncs;

%set(gcf, 'Renderer', 'opengl');   % good for static plots

%%  %%%%%%%%%%%%%%%%%% Run for individual figures %%%%%%%%%%%%%%%%%
%% Figure S1a-b
FigSaveFileName='FigS1ab';
% load data 
load([DataPath 'BhvData.mat']);
% generate figure
varargout=cell(3);
[varargout{1:3}]=bhvAna.PlotAvgPSM(AllPSMPerf_MonkeyC); % results for Monkey Ch, Fig. S1a
close([varargout{2} varargout{3}]);varargoutAll=varargout([1]);
[varargout{1:3}]=bhvAna.PlotAvgPSM(AllPSMPerf_MonkeyS); % results for Monkey Si, Fig. S1a
close([varargout{2} varargout{3}]);varargoutAll=[varargoutAll varargout([1])];

%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargoutAll,'SaveEachFrame',1,'enforce_daspect',1)

%% Figure S2a
FigSaveFileName='FigS2a';
% load data 
load([DataPath 'BhvData.mat']);
% generate figure
varargout=cell(1,3);
[varargout{1:3}]=bhvAna.PlotAvgPSMCongInCong(BhvData.IndSamp,BhvData.AllPSMPerf,0);
close([varargout{1} varargout{2}]);varargout=varargout([3]);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

%% Figure S2b
FigSaveFileName='FigS2b';
% load data 
load([DataPath 'BhvData.mat']);
% generate figure
varargout=cell(1,3);
[varargout{1:3}]=bhvAna.PlotAvgPSMCongInCong(BhvData.IndSamp,BhvData.AllPSMPerf,1);
close([varargout{1} varargout{2}]);varargout=varargout([3]);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

%% Figure S2e-f 
FigSaveFileName='FigS2ef';
% load data 
load([DataPath 'BhvData.mat']);
% generate figure
varargout=cell(3);
[varargout{1:3}]=bhvAna.PlotTrlPerf(AllTrlPerf,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist);
close([varargout{1} varargout{2}]);varargout=varargout(3);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

%% Figure S2c,d,e,g,h
FigSaveFileName='FigS2cdgh';
% load data 
load([DataPath 'BhvData.mat']);
% generate figure
varargout=cell(1,5);
[varargout{1:5}]=bhvAna.PlotResponseLocDist(BhvData.IndSamp,BhvData.AllSeqHist,'PaperFigName','FigS2cdgh');
close([varargout{3}]);varargout=varargout([1 2 4 5]);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)


%% Figure S4
FigSaveFileName='FigS4';
h1=PlotGLMmixedselectivityResults(FigSaveFileName,DataPath); % plots results  for Fig.1h-m and Fig. S3
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,h1,'SaveEachFrame',1,'enforce_daspect',1)

% run this if you want to see the whole process of analysis
PopulationAnalysisTemp(0, 0, [], [1], 'SensoryCat','ProcessingStep',16)

%% Figure S5c,d,e
FigSaveFileName={'FigS5c','FigS5d','FigS5e'};
% generate figure
varargout=cellfun(@(x) PlotClassifierResults(x,DataPath),FigSaveFileName,'uniformoutput',0);
%  save the figures we have generated
arrayfun(@(x) FigParams.SaveFigSeries(FigSaveFileName{x}, FigSavePath,varargout{x},'SaveEachFrame',1,'enforce_daspect',1),1:length(varargout))


%% Figure S7abc
FigSaveFileName='FigS7abc';

% generate figure
varargout=cell(1);
[varargout{1}]=PlotDynamicTransformationResults('FigS7',DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

% run this if you want to see the whole process of analysis
% PopulationAnalysisTemp(0, 0, [], [1], '3D_Shared_Color_Response_EntropyR2','ProcessingStep',4,'SpkCntStartFieldName','SACCADE_START','CalShuffTrlOrderClassifier',0,'RunCrossTemporalClassifer',0,'CalShuffleClassifier',0,...
% 'PopulationAna.PSTHbin',100,'SpkCntStartFieldName','SACCADE_START','DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',16*ones(1,7),'DividSpockClassifier',3,'NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',3)


%% Figure S8ab
FigSaveFileName='FigS8ab';
% load data 
load([DataPath 'BhvData.mat']);
% generate figure
varargout=cell(3);
[varargout{1:3}]=bhvAna.PlotTrlPerf(AllTrlPerf,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist);
close([varargout{1} varargout{3}]);varargout=varargout(2);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)


%% figure S8f
FigSaveFileName='FigS8f';
h1=PlotBlockInfoTransferResults('FigS8f1',DataPath);
h2=PlotBlockInfoTransferResults('FigS8f2',DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,[h1 h2],'SaveEachFrame',1,'enforce_daspect',1)

% if you want to run the full analysis pipeline
%PopulationAnalysisTemp(0, 0, [], 1, 'Learning3D_RuleR2SeqAnalysisLearningCtrlR2','ProcessingStep',4,'DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',[1 1],'CalShuffTrlOrderClassifier',0)
%PopulationAnalysisTemp(0, 0, [], 1, 'Learning3D_RuleR2SeqAnalysisLearningCtrl','ProcessingStep',4,'DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',[1 1],'CalShuffTrlOrderClassifier',0)

%% Figure S9ab
FigSaveFileName='FigS9ab';

% generate figure
varargout=cell(1);
[varargout{1}]=PlotClassifierAngleData(FigSaveFileName,DataPath);

%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)


%% Figure S9c,d
FigSaveFileName='FigS9cd';
varargout=cell(1);
[varargout{1}]=PlotAngleAnaTDRresults(DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

%% Figure S9g 
FigSaveFileName='FigS9g';

% generate figure
varargout=cell(3);
[varargout{1}]=PlotCompressionIndexData('FigS9g','FEF',DataPath);
[varargout{2}]=PlotCompressionIndexData('FigS9g','LIP',DataPath);
[varargout{3}]=PlotCompressionIndexData('FigS9g','IT',DataPath);

%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

end
function FS=KickoffMyfunc(RunonCluster,DateNum)
% sets up all of the initial options and path
global AnalysisOpts  

%% ok now first setup everything for cluster first because it needs path
AnalysisOpts.RunOnCluster=RunonCluster;

if AnalysisOpts.RunOnCluster
    RootPath='/jukebox/buschman/';
    FS='/';
else
    if ismac
        RootPath='/Volumes/buschman/';
        FS='/';
    elseif ispc
        RootPath='Z:\';
        FS='\';
    end
end

% add core function path first so we can sturt up everything
addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Rule Representation Project',...
    FS 'Submission Code' FS 'Final Code for Nature Paper 2025 Tafazoli et al' FS 'Analysis Pipeline' FS 'Core functions' FS]));
SetupAllVars(DateNum)  %% set up the path and initialize vars

end

