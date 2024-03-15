function GenerateBhvPaperFigs(varargin)
% Generates figures for behavioral analysis for 
% Tafazoli et al, Building compositional tasks with shared neural subspaces
% https://www.biorxiv.org/content/10.1101/2024.01.31.578263v1.article-metrics
 
global AnalysisOpts
Animal='ALL';
AnalysisOpts.Project='Rule Representation';
AnalysisOpts.BhvAna.Kidname=Animal;
AnalysisOpts.AnalysisType='Analysis Pipeline';
KickoffMyfunc(0,0);
SetAnalysisOptions_RuleRepresentation(varargin) % set up the parameters used in the analysis

%% Load Behavioral Data
BhvData=BhvAnalysisPiplines('LoadBhvData',Animal,[]);

%% Reprocess data if we need to
BhvData=BhvAnalysisPiplines('ReprocessData',Animal,BhvData);

%% generate figures for Psychometic curve
BhvAnalysisPiplines('Psychometric_Curve',Animal,BhvData);

%% generate figures for TrialPerformance
BhvAnalysisPiplines('TrialPerformance',Animal,BhvData);
end

function varargout=BhvAnalysisPiplines(AnalysisType,Animal,BhvData)
global AnalysisOpts

bhvAna=BhvAnalysisFuncs;
FigParams=fig_params;

switch AnalysisType
    case 'PreprocessData'
        varargout{1}=PreprocessBhvData(Animal);
    case 'ReprocessData' % reprocess behavioral data
        [BhvData.AllPSMPerf,BhvData.AllTrlPerf,BhvData.IndSamp,BhvData.AllTrlCount,BhvData.AllTrlCountDay,...
            BhvData.NBlocksDay,BhvData.RewardPulse,BhvData.NCorrectTrl,BhvData.NumRewards,BhvData.AllSeqHist] =...
            bhvAna.CancatinateInfoDays(BhvData.Perf);
        varargout{1}=BhvData;
    case 'LoadBhvData'
        varargout{1}=LoadBhvData(Animal); % load behavioral data
    case 'Psychometric_Curve' % plot Psychometric curves for each animal
        varargout=cell(1,4);
        [varargout{1:4}]=bhvAna.PlotAvgPSM(BhvData.AllPSMPerf);
    case 'TrialPerformance'
        varargout=cell(1,2);
        [varargout{1:2}]=bhvAna.PlotTrlPerf(BhvData.AllTrlPerf,BhvData.AllTrlCount,BhvData.AllTrlCountDay,...
            BhvData.NBlocksDay,BhvData.RewardPulse,BhvData.NCorrectTrl,BhvData.NumRewards,BhvData.AllSeqHist);
end
% %%  save the figures we have generated in case
% FigSaveFileName=sprintf(['%s_%s'],AnalysisType,Animal);
% FigSavePath=[AnalysisOpts.ResultsSavePath AnalysisOpts.FS 'Behavior' AnalysisOpts.FS];
% if ~exist(FigSavePath,'file');mkdir(FigSavePath);end
% FigParams.SaveFigSeries(FigSaveFileName, FigSavePath,varargout,'SaveEachFrame',1,'enforce_daspect',1)
end

function FS=KickoffMyfunc(RunonCluster,DateNum)
% sets up all of the initial options and path
global AnalysisOpts

%% ok now first setup everything for cluster first because it needs path
AnalysisOpts.RunOnCluster=RunonCluster;

if ismac
    RootPath='/Volumes/';
    FS='/';
elseif ispc
    RootPath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Submission Code\Submission to Nature March 2024\Tafazoli et al 2024 Code\';
    FS='\';
end

% add core function path first so we can sturt up everything
addpath(genpath([RootPath 'Core functions' FS]));
SetupAllVars(DateNum)  %% set up the path and initialize vars
end

function BhvData=LoadBhvData(Animal)
global AnalysisOpts
if strcmp(Animal,'ALL')
    BhvData=load([AnalysisOpts.BhvMdlPath Animal '_BhvData.mat'],...
        'AllPSMPerf','Perf','AllTrlPerf','IndSamp','AllTrlCount','AllTrlCountDay','NBlocksDay','RewardPulse','NCorrectTrl','NumRewards','AllSeqHist');
else   
    EndDate=AnalysisOpts.BhvAna.([Animal '_EndDate']);
    InitialDate=AnalysisOpts.BhvAna.([Animal '_InitialDate']);
    
    BhvData=load([AnalysisOpts.BhvMdlPath Animal '_BhvData_' char(datetime(InitialDate)) '_till_' char(datetime(EndDate)) '.mat'],...
        'AllPSMPerf','Perf','AllTrlPerf','IndSamp','AllTrlCount','AllTrlCountDay','NBlocksDay','RewardPulse','NCorrectTrl','NumRewards','AllSeqHist');
end
end

