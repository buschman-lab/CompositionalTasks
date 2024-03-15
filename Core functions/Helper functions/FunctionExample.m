 function FunctionExample(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
% exmaple function that will be same to initialize  
global AnalysisOpts AnalysisData

%% define primary vars before kick off
AnalysisOpts.RunonCluster=RunonCluster;
AnalysisOpts.Project='Learning attentional templates';  % can be 'Learning attentional templates' or 'Rule Representaion'
AnalysisOpts.Animal='Scooter';
AnalysisOpts.AnalysisType='Data Preprocessing'; % can be 'Oscillations Breakdown ToolBox';
AnalysisOpts.AnalysisFocus1  ='';  % 'ExplorParams','PlotMotifs'
AnalysisOpts.SubAnalysisType1='';  % define type of sub analsysis  can be also 'LFP_Trial'
AnalysisOpts.AnalysisPathName='';  % this is the name of the folder used in input-output data

%% kick off the function
FS=KickoffMyfunc(RunonCluster,DateNum);

%% define analysis options in a gobal varibales
AnalysisOpts.PairNum=PairNum; % what pair we are looking at 
AnalysisOpts.Ch_2look=ChNum;% how many channels do you want to look at, leave empty if non
AnalysisOpts.AreaNum=AreaNum;
ParseParams(varargin) % add all of the additional parameters we have added to the function

%% define classes necessary for this analysis 
TlBoxHelp=ToolboxHelpers;
[FigPrms,ManData,TrialFunc,MotifFunc,TimeFreq]=TlBoxHelp.LoadAnalysisClasses; % load all of the classes we want 
[Fs,FsLFP,FsWave,FsWaveTarg]=TlBoxHelp.SetupGeneralFuncs; % set up all general functions we want
% get block and trial information we need
[TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet]=TrialFunc.InitializeTrialFuncs; % run common trial functions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Starting the analysis of this function %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 end
%% put the functions to kick off here
function FS=KickoffMyfunc(RunonCluster,DateNum)
% sets up all of the initial options and path
global AnalysisOpts AnalysisFuncs AnalysisData

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
addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Core functions' FS]));
SetupAllVars(DateNum)  %% set up the path and initialize vars

end
 
