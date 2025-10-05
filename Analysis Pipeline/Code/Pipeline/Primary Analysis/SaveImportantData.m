function [TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet,TimeEpoch]=SaveImportantData(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
%% Explanation of the function goes here
% Pipline to save data for specific conditions so we read them later 
% define some global functions to carry options and data
% input RunonCluster: runs on cluster
% varargin(1) number of recording date
% Sina Tafazoli timeless
global AnalysisOpts AnalysisData
tic
%% define primary vars before kick off
AnalysisOpts.RunonCluster=RunonCluster;
AnalysisOpts.Project='Rule Representation';  
AnalysisOpts.Animal='';
AnalysisOpts.AnalysisType='Analysis Pipeline';
AnalysisOpts.AnalysisFocus1  ='';  % define focus of this analysis
AnalysisOpts.SubAnalysisType1='';  % define type of sub analysis
AnalysisOpts.AnalysisPathName='';  % this is the name of the folder used in input-output data

%% kick off the function
FS=KickoffMyfunc(RunonCluster,DateNum(1));

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
% get block and trial information we need
[TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet,BlockSpec,StimMap]=TrialFunc.InitializeTrialFuncs; % run common trial functions
fprintf('\nGetting recording info for: %s',AnalysisOpts.RecDate)

%% get the number of time epochs based on how long the recording would be
AnalysisOpts.SpkParams.PeriodLength=1.5;
AnalysisOpts.SpkParams.BaselineDelay=0.6; %%% look at x ms before the start
AnalysisOpts.Bhv.MaxNTrial=50;
AnalysisOpts.Bhv.NTrial2Swtch=50;
AnalysisOpts.Area2look=AnalysisOpts.AreaNames(AreaNum);
if isempty(AreaNum);AnalysisOpts.Area2look=AnalysisOpts.AreaNames;end
if isempty(AnalysisOpts.Rule_2look);AnalysisOpts.Rule_2look=[1 2 3];end

% report details of the analysis 
ManData.ReportAnalysisDetails;

%% Define what channels you want
if DateNum==0;DateNum=1:34;end
Recs=DateNum; 
[Chns,RecSpecs]=ManData.RunFuncOnRec(TrialFunc,'FindWantedSpkChSet',Recs,{[AnalysisOpts.SpikeQuality2Look],AnalysisOpts.Area2look},4);
if ~isempty(ChNum);Chns=cellfun(@(x) x(ChNum),Chns,'uniformoutput',0); end % is we want 1 channel
TrialFunc.ParseChannelInformation(Chns,RecSpecs);

PSTHRaster=[];TrialSpikeData=cell(1,2);
%% Grab all of the spiking data from channels
[SpikeData,RecSpecs]=ManData.RunFuncOnRec(HelperFunc,'getSpikeData',Recs,Chns(:,1:3),1);

%% Define what type of trials you want
[WantedBlockSpec]=ManData.RunFuncOnRec(TrialFunc,'FindWantedTrialTypeBlockSpec',Recs,{[],AnalysisOpts.Rule_2look,{'ALL','ALL'},{'ALL','ALL'},'ALL'},1);

%% Take the data from the trials we care about
[TrialSpikeData,RecSpecs]=ManData.RunFuncOnRec(TrialFunc,'GrabWantedSpksWantedTrial',Recs,[WantedBlockSpec,SpikeData],2);
TrialSpikeData=ManData.ConcatinateStructs(TrialSpikeData,1); % concatinate data across recordings

%% start analysis
[PSTHRaster_150ms]=NeuAnaFunc.ProcessSpikeCountingRec(TrialSpikeData,AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.TrlSpkTimeFieldName,'PSTH_Bin',0.150);
[PSTHRaster_100ms]=NeuAnaFunc.ProcessSpikeCountingRec(TrialSpikeData,AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.TrlSpkTimeFieldName,'PSTH_Bin',0.100);
[PSTHRaster_50ms]=NeuAnaFunc.ProcessSpikeCountingRec(TrialSpikeData,AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.TrlSpkTimeFieldName,'PSTH_Bin',0.050);


%% save data for this condition  
SaveFileName=[AnalysisOpts.Animal '_' AnalysisOpts.Area2look{1} '_' AnalysisOpts.TrlSpkTimeFieldName '_' AnalysisOpts.SpkCntStartFieldName]; 
save([AnalysisOpts.DataSavePath 'Core Data' filesep 'SpikingData' filesep SaveFileName],...
    'Chns','RecSpecs','SpikeData','WantedBlockSpec','TrialSpikeData','PSTHRaster_150ms',...
    'PSTHRaster_100ms','PSTHRaster_50ms','-v7.3');

ElapsedTime=toc;
fprintf('\nElapsed time to run this pipeline:%0.1f(secs)',ElapsedTime)
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

