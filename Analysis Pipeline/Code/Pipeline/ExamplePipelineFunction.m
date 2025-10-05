function ExamplePipelineFunction(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
%% Explanation of the function goes here
% define some global functions to carry options and data
% input RunonCluster: runs on cluster
% varargin(1) number of recording date
% Sina Tafazoli timeless
global AnalysisOpts 

%% define primary vars before kick off
AnalysisOpts.RunonCluster=RunonCluster;
AnalysisOpts.Project='Rule Representation';  
AnalysisOpts.Animal='';
AnalysisOpts.AnalysisType='Analysis Pipeline'; 
AnalysisOpts.AnalysisFocus1  ='';  % define focus of this analysis
AnalysisOpts.SubAnalysisType1='';  % define type of sub analysis
AnalysisOpts.AnalysisPathName='';  % this is the name of the folder used in input-output data

%% kick off the function
FS=KickoffMyfunc(RunonCluster,DateNum);

%% define analysis options in a gobal varibales
AnalysisOpts.PairNum=PairNum; % what pair we are looking at
AnalysisOpts.Ch_2look=ChNum;% how many channels do you want to look at, leave empty if non
AnalysisOpts.AreaNum=AreaNum;
ParseParams(varargin) % add all of the additional parameters we have added to the function

%% define classes necessary for this analysis
HelperFunc=HelperFuncs;
[TimeFreq,FigPrms,ManData,TrialFunc,FilterFunc,NeuAnaFunc,CrosFreqCop,AnalysisFunc]=HelperFunc.LoadAnalysisClasses; % load all of the classes we want
[Fs,FsLFP,FsWave,FsWaveTarg]=HelperFunc.SetupGeneralFuncs; % set up all general functions we want

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Starting the analysis of this function %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% details of the analysis will go here


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
    RootPath='Z:\';
    FS='\';
end
% add core function path first so we can sturt up everything
addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Core functions' FS]));
SetupAllVars(DateNum)  %% set up the path and initialize vars

end


