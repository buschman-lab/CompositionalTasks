% Generating Figure 1
% Guide: In all code Rule1=S1, Rule2=C2 and Rule3=C1
function Figure3

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

%% Figure 3cde
FigSaveFileName='Fig3cde';

% generate figure
varargout=cell(1);
[varargout{1}]=PlotDynamicTransformationResults('Fig3',DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

% run this if you want to see the whole process of analysis
% PopulationAnalysisTemp(0, 0, [], [1], '3D_Shared_Color_Response_EntropyR3Bal','ProcessingStep',4,'SpkCntStartFieldName','SACCADE_START','CalShuffTrlOrderClassifier',0,'RunCrossTemporalClassifer',0,'CalShuffleClassifier',0,...
% 'PopulationAna.PSTHbin',100,'SpkCntStartFieldName','SACCADE_START','DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',16*ones(1,7),'DividSpockClassifier',3,'NTrlRngTrainLearningArea',-50,'NTrlRngTestLearningArea',-50,'ntrlPerCondArea',4)





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

