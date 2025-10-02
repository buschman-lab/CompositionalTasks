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

%% Figure 1e 
FigSaveFileName='Fig1e';
% load data 
load([DataPath 'BhvData.mat']);
% generate figure
varargout=cell(3);
[varargout{1:3}]=bhvAna.PlotAvgPSM(AllPSMPerf);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout(1),'SaveEachFrame',1,'enforce_daspect',1)


%% Figure 1e 
FigSaveFileName='Fig1f';
% load data 
load([DataPath 'BhvData.mat']);
% generate figure
varargout=cell(3);
[varargout{1:3}]=bhvAna.PlotTrlPerf(AllTrlPerf,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist);
close([varargout{2} varargout{3}]);varargout=varargout(1);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

%% Figure 1h-m  
FigSaveFileName='Fig1hm';

% generate figure
varargout=cell(2);
[varargout{1:2}]=PlotGLMtaskvarEncoding(FigSaveFileName,DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

% run this if you want to see the whole process of analysis
%PopulationAnalysisTemp(0, 0, [], [1], 'SensoryCat','ProcessingStep',16)

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

