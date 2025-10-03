function GenerateBhvPaperFigs(AnalysisType,Animal,varargin)
%GENERATEBHVPAPERFIGS generates figures for behavioral analysis for the
%paper 
% Sina Tafazoli Januray 2023
% This function uses functions from
% Z:\Projects\Rule_Representation\Behavior AnalysisCode\Analysis_Code and
% they need to be added to the path along with Core functions

% @FigureType: depending on the figures and subfigures run scripts 
% can be "Psychometric_Curve"
global AnalysisOpts 
AnalysisOpts.Project='Rule Representation';
AnalysisOpts.BhvAna.Kidname=Animal;
AnalysisOpts.AnalysisType='Analysis Pipeline';
KickoffMyfunc(0,0);
SetAnalysisOptions_RuleRepresentation(varargin)
 
%% generate plots is a sequence now 
if ~isempty(AnalysisType) && ~strcmp(AnalysisType,'ALL')
    BhvAnalysisPiplines(AnalysisType,Animal,[]);
    return;
end

%% Preprocess data 
%BhvAnalysisPiplines('PreprocessData',Animal,[]);

%% Load Behavioral Data
BhvData=BhvAnalysisPiplines('LoadBhvData',Animal,[]);

%% Reprocess data if we need to 
BhvData=BhvAnalysisPiplines('ReprocessData',Animal,BhvData);

%% generate figures for Psychometic curve
BhvAnalysisPiplines('Psychometric_Curve',Animal,BhvData);

%% generate figures for TrialPerformance 
BhvAnalysisPiplines('TrialPerformance',Animal,BhvData);

%% generate figures for SamplePerfInfo 
BhvAnalysisPiplines('SamplePerfInfo',Animal,BhvData);

%% generate figures for PSMCongInCong from the end of the block
BhvAnalysisPiplines('PSMCongInCongEndBlk',Animal,BhvData);

%% generate figures for PSMCongInCong from the all of the trials
BhvAnalysisPiplines('PSMCongInCongAllTrls',Animal,BhvData);

%% generate figures for Response Location 
BhvAnalysisPiplines('ResponseLocDist',Animal,BhvData);

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
        varargout=cell(1,3);
        [varargout{1:3}]=bhvAna.PlotAvgPSM(BhvData.AllPSMPerf);
    case 'TrialPerformance'
        varargout=cell(1,3);
        [varargout{1:3}]=bhvAna.PlotTrlPerf(BhvData.AllTrlPerf,BhvData.AllTrlCount,BhvData.AllTrlCountDay,...
            BhvData.NBlocksDay,BhvData.RewardPulse,BhvData.NCorrectTrl,BhvData.NumRewards,BhvData.AllSeqHist);
    case 'SamplePerfInfo'
        varargout=cell(1,1);
        [varargout{1}]=bhvAna.PlotSampPerfInfo(BhvData.IndSamp);
    case 'PSMCongInCongEndBlk'
        varargout=cell(1,4);
        [varargout{1:4}]=bhvAna.PlotAvgPSMCongInCong(BhvData.IndSamp,BhvData.AllPSMPerf,0);
    case 'PSMCongInCongAllTrls'
        varargout=cell(1,4);
        [varargout{1:4}]=bhvAna.PlotAvgPSMCongInCong(BhvData.IndSamp,BhvData.AllPSMPerf,1);
    case 'ResponseLocDist' % plots distribution of response locations 
        varargout=cell(1,5);
        [varargout{1:5}]=bhvAna.PlotResponseLocDist(BhvData.IndSamp,BhvData.AllSeqHist);

end
%%  save the figures we have generated in case
FigSaveFileName=sprintf(['%s_%s'],AnalysisType,Animal);
FigSavePath=[AnalysisOpts.ResultsSavePath 'Behavior' AnalysisOpts.FS];
if ~exist(FigSavePath,'file');mkdir(FigSavePath);end
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

