function [TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet,TimeEpoch]=GetRecordingInfo(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
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
% get block and trial information we need
[TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet,BlockSpec,StimMap]=TrialFunc.InitializeTrialFuncs; % run common trial functions
fprintf('\nGetting recording info for: %s',AnalysisOpts.RecDate)

%% get the number of time epochs based on how long the recording would be
%SpikeDataq=HelperFunc.getSpikeData([1],ones(1,3),ChannelArea) ;
%[ChNumAll,ClustNumAll,AreaIndAll]=TrialFunc.FindWantedSpkChSet([1 2],{'PFC','LIP'});
%[NewBlockSpec]    =TrialFunc.FindWantedTrialTypeBlockSpec(BlockSpec,[1],{'ALL','ALL'},{'ALL','ALL'},'ALL');
%[TrialSpikeTime,CriticalTimes]=GrabWantedSpksWantedTrial(obj,WantedBlockSpec,SpikeData,varargin) 
%AnalysisOpts.TrlSpkTimeFieldName='AllTrials';
AnalysisOpts.TrlSpkTimeFieldName='FromSwitch';
%AnalysisOpts.SpkCntStartFieldName='FIXATE_ACQUIRED';%'SAMPLE_ON';
AnalysisOpts.SpkParams.PeriodLength=2;
AnalysisOpts.SpkParams.BaselineDelay=0.5; %%% look at 500ms before the start
AnalysisOpts.Bhv.MaxNTrial=50;
AnalysisOpts.Bhv.NTrial2Swtch=50;
AnalysisOpts.Area2look=AnalysisOpts.AreaNames{AreaNum};
if isempty(AnalysisOpts.Rule_2look);AnalysisOpts.Rule_2look=[1 2 3];end

Recs=1;%AnalysisOpts.ChicoRecNums;
%% Define what channels you want 
[Chns,RecSpecs]     =ManData.RunFuncOnRec(TrialFunc,'FindWantedSpkChSet',Recs,{[1],{AnalysisOpts.Area2look}},3);
if ~isempty(ChNum);Chns=cellfun(@(x) x(ChNum),Chns,'uniformoutput',0); end % is we want 1 channel 
%% Grab all of the spiking data from channels 
[SpikeData,RecSpecs]=ManData.RunFuncOnRec(HelperFunc,'getSpikeData',Recs,Chns,1);
%% Define what type of trials you want 

[WantedBlockSpec]   =ManData.RunFuncOnRec(TrialFunc,'FindWantedTrialTypeBlockSpec',Recs,{[],AnalysisOpts.Rule_2look,{'ALL','ALL'},{'ALL','ALL'},'ALL'},1);
%% Take the data from the trials we care about 

[TrialSpikeData,RecSpecs]=ManData.RunFuncOnRec(TrialFunc,'GrabWantedSpksWantedTrial',Recs,[WantedBlockSpec,SpikeData],2);
TrialSpikeData=ManData.ConcatinateStructs(TrialSpikeData,1); % concatinate data across recordings 

%% start analysis 
[PSTHRaster]=NeuAnaFunc.ProcessSpikeCountingRec(TrialSpikeData,AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.TrlSpkTimeFieldName);

%% do single cell analysis only in time 
SingCellAnaFigs=cell(1,7);
AnalysisOpts.RuleTxt=strrep(num2str(AnalysisOpts.Rule_2look),' ','');
[SingCellAnaFigs{:}]=NeuAnaFunc.SetupSingleCellAnalysis(PSTHRaster,TrialSpikeData{2});
% save figure for single cells 
[~,~,SingCellAnaFigsFileName]=ManData.GetFileName([],['SingCellAna_' AnalysisOpts.Area2look '_' AnalysisOpts.SpkCntStartFieldName 'R' AnalysisOpts.RuleTxt],'SaveInResults',1,'WantedDate','ALL');
FigPrms.SaveFigSeries([],SingCellAnaFigsFileName,[SingCellAnaFigs],'SaveEachFrame',1)

%% do sinlge cell analysis in time and trials
SingCellAnaFigs_TrialTime=cell(1);
SingCellAnaFigs_TrialTime{1}=NeuAnaFunc.SetupSingleCellAnalysis_TrialTime(PSTHRaster,TrialSpikeData{2});

% save figure for single cells 
[~,~,SingCellAnaFigsFileName]=ManData.GetFileName([],['SingCellAna_TrialTime' AnalysisOpts.Area2look '_' AnalysisOpts.SpkCntStartFieldName 'R' AnalysisOpts.RuleTxt],'SaveInResults',1,'WantedDate','ALL');
FigPrms.SaveFigSeries([],SingCellAnaFigsFileName,[SingCellAnaFigs_TrialTime],'SaveEachFrame',1)

% 
% %% cluster spike shapes 
% SpkShp=cellfun(@(y) cell2mat(cellfun(@(x) mean(x.SpikeShape,1)',y,'UniformOutput',0))',Funcout,'UniformOutput',0);
% 
% %% start by plotting the number of channels in each area
% [Count]=hist(ChannelArea(ChannelArea~=0),1:5);
% FigPrms.BarPlot(1:5,Count,'b','Area Number','Number of Chs',['Ch count for Rec' AnalysisOpts.RecDate ]); % plots whatever
% xticklabels(AnalysisOpts.AreaNames)

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

