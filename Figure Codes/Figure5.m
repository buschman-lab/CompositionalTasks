% Generating Figure 5
% Guide: In all code Rule1=S1, Rule2=C2 and Rule3=C1
function Figure5(DataPath,FigSavePath)

%% Prepare Vaiables
global AnalysisOpts 
AnalysisOpts.Project='Rule Representation';
AnalysisOpts.BhvAna.Kidname='ALL';
AnalysisOpts.AnalysisType='Analysis Pipeline';
KickoffMyfunc(0,0);
SetAnalysisOptions_RuleRepresentation(); % set all of the parameters
 
% Adjust the Path to locate the folder with all data
%DataPath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Submission Code\Final Code for Nature Paper 2025 Tafazoli et al\Figure Data\';
%FigSavePath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Submission Code\Final Code for Nature Paper 2025 Tafazoli et al\Figure Results\';

% Setup the class to plot the data
FigParams=fig_params;
bhvAna=BhvAnalysisFuncs;
%opengl software
%set(gcf, 'Renderer', 'opengl');   % good for static plots


%% Figure 5a,b,d 
FigSaveFileName='Fig5abd';

h1=PlotClassifierEncodingAxisCompressionIndex('Fig5a',DataPath);
h2=PlotClassifierEncodingAxisCompressionIndex('Fig5b',DataPath);
h3=PlotClassifierEncodingAxisCompressionIndex('Fig5d',DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,[h1 h2 h3],'SaveEachFrame',1,'enforce_daspect',1)

%% Figure 5c 
FigSaveFileName='Fig5c';

% generate figure
varargout=cell(1);
[varargout{1}]=PlotCompressionIndexData('Fig5c','PFC',DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)


%% Figure 5e 
FigSaveFileName='Fig5e';

% generate figure
varargout=cell(1);
[varargout{1}]=PlotCorrelationAnalysis('Fig5e',DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

%% Figure 5f 
FigSaveFileName='Fig5f';

% generate figure
varargout=cell(1);
[varargout{1}]=PlotCorrelationAnalysis('Fig5f',DataPath);
%  save the figures we have generated
FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)

%% Figure 5h 
FigSaveFileName='Fig5h';

% generate figure
varargout=cell(1);
[varargout{1}]=PlotClassifierLearningResults('Fig5h',DataPath);
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

AnalysisOpts.FS=FS;

% add core function path first so we can sturt up everything
% addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Rule Representation Project',...
%     FS 'Submission Code' FS 'Final Code for Nature Paper 2025 Tafazoli et al' FS 'Core functions' FS]));
% SetupAllVars(DateNum)  %% set up the path and initialize vars

end

