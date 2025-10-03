% Generating Figure 1
% Guide: In all code Rule1=S1, Rule2=C2 and Rule3=C1
function Figure1

%% Prepare Vaiables
global AnalysisOpts 
AnalysisOpts.Project='Rule Representation';
AnalysisOpts.BhvAna.Kidname='ALL';
AnalysisOpts.AnalysisType='Analysis Pipeline';
KickoffMyfunc(0,0);
SetAnalysisOptions_RuleRepresentation(); % set all of the parameters
 
% Adjust the Path to locate the folder with all data
DataPath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Submission Code\Final Code for Nature Paper 2025 Tafazoli et al\Figure Data\';
FigSavePath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Submission Code\Final Code for Nature Paper 2025 Tafazoli et al\Figure Results\';

% Setup the class to plot the data
FigParams=fig_params;
bhvAna=BhvAnalysisFuncs;
opengl software
set(gcf, 'Renderer', 'opengl');   % good for static plots


%% Figure 4a 
FigSaveFileName='Fig4a';
% load data 
load([DataPath 'DataFig1f.mat']);
% generate figure
varargout=cell(3);
[varargout{1:3}]=bhvAna.PlotTrlPerf(AllTrlPerf,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist);
close([varargout{2} varargout{3}]);varargout=varargout(1);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)


%% figure 4b
FigSaveFileName='Fig4b';
h1=PlotBlockInfoTransferResults('Fig4b1',DataPath);
h2=PlotBlockInfoTransferResults('Fig4b2',DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,[h1 h2],'SaveEachFrame',1,'enforce_daspect',1)

% if you want to run the full analysis pipeline
%PopulationAnalysisTemp(0, 0, [], 1, 'Learning3D_RuleR2SeqAnalysisLearningCtrlR2','ProcessingStep',4,'DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',[1 1],'CalShuffTrlOrderClassifier',0)
%PopulationAnalysisTemp(0, 0, [], 1, 'Learning3D_RuleR2SeqAnalysisLearningCtrl','ProcessingStep',4,'DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',[1 1],'CalShuffTrlOrderClassifier',0)

%% figure 4c-l
h1=PlotClassifierLearningResults('Fig4c',DataPath);
h2=PlotClassifierLearningResults('Fig4d',DataPath);
h3=PlotClassifierLearningResults('Fig4f',DataPath);
h4=PlotClassifierLearningResults('Fig4g',DataPath);
h6=PlotClassifierLearningResults('Fig4i',DataPath);
h7=PlotClassifierLearningResults('Fig4j',DataPath);
h8=PlotClassifierLearningResults('Fig4l',DataPath);

%% Figure 4h 
FigSaveFileName='Fig5f';

% generate figure
varargout=cell(1);
[varargout{1}]=PlotCorrelationAnalysis('Fig4h',DataPath);
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
    FS 'Submission Code' FS 'Final Code for Nature Paper 2025 Tafazoli et al' FS 'Core functions' FS]));
SetupAllVars(DateNum)  %% set up the path and initialize vars

end

