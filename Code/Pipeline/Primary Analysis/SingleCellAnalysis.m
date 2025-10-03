function [TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet,TimeEpoch]=SingleCellAnalysis(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
%% Explanation of the function goes here
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
%SpikeDataq=HelperFunc.getSpikeData([1],ones(1,3),ChannelArea) ;
%[ChNumAll,ClustNumAll,AreaIndAll]=TrialFunc.FindWantedSpkChSet([1 2],{'PFC','LIP'});
%[NewBlockSpec]    =TrialFunc.FindWantedTrialTypeBlockSpec(BlockSpec,[1],{'ALL','ALL'},{'ALL','ALL'},'ALL');
%[TrialSpikeTime,CriticalTimes]=GrabWantedSpksWantedTrial(obj,WantedBlockSpec,SpikeData,varargin) 
%AnalysisOpts.SpkCntStartFieldName='SAMPLE_ON';
AnalysisOpts.SpkParams.PeriodLength=1.5;
AnalysisOpts.SpkParams.BaselineDelay=0.6; %%% look at x ms before the start
AnalysisOpts.Bhv.MaxNTrial=20;
AnalysisOpts.Bhv.NTrial2Swtch=20;
AnalysisOpts.Area2look=AnalysisOpts.AreaNames(AreaNum);
if isempty(AreaNum);AnalysisOpts.Area2look=AnalysisOpts.AreaNames;end
if isempty(AnalysisOpts.Rule_2look);AnalysisOpts.Rule_2look=[1 2 3];end

% report details of the analysis 
ManData.ReportAnalysisDetails;

Recs=DateNum;%13%AnalysisOpts.ChicoRecNums;
%% Define what channels you want 
[Chns,RecSpecs]     =ManData.RunFuncOnRec(TrialFunc,'FindWantedSpkChSet',Recs,{[AnalysisOpts.SpikeQuality2Look],AnalysisOpts.Area2look},4);
if ~isempty(ChNum);Chns=cellfun(@(x) x(ChNum),Chns,'uniformoutput',0); end % is we want 1 channel 
TrialFunc.ParseChannelInformation(Chns,RecSpecs);

PSTHRaster=[];TrialSpikeData=cell(1,2);
if sum(AnalysisOpts.ProcessingStep==[1 2])
    %% Grab all of the spiking data from channels
    [SpikeData,RecSpecs]=ManData.RunFuncOnRec(HelperFunc,'getSpikeData',Recs,Chns(:,1:3),1);

    %% Define what type of trials you want
    [WantedBlockSpec]=ManData.RunFuncOnRec(TrialFunc,'FindWantedTrialTypeBlockSpec',Recs,{[],AnalysisOpts.Rule_2look,{'ALL','ALL'},{'ALL','ALL'},'ALL'},1);

    %% Take the data from the trials we care about
    [TrialSpikeData,RecSpecs]=ManData.RunFuncOnRec(TrialFunc,'GrabWantedSpksWantedTrial',Recs,[WantedBlockSpec,SpikeData],2);
    TrialSpikeData=ManData.ConcatinateStructs(TrialSpikeData,1); % concatinate data across recordings

    %% start analysis
    [PSTHRaster]=NeuAnaFunc.ProcessSpikeCountingRec(TrialSpikeData,AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.TrlSpkTimeFieldName);
end

%% do single cell analysis only in time 
ChAreas=unique(cell2mat(Chns(:,3)')); % find the areas of the channels we have 
nar=length(ChAreas);
AnalysisOpts.RuleTxt=strrep(num2str(AnalysisOpts.Rule_2look),' ','');
AnalysisOpts.NfigProcesStep=[2 1 3*3 3*2*nar+1 7+3*2*nar];

SingCellAnaFigs=cell(1,AnalysisOpts.NfigProcesStep(AnalysisOpts.ProcessingStep));
TrialFunc.RevertCh_2look;
[SingCellAnaFigs{:}]=NeuAnaFunc.SetupSingleCellAnalysis(PSTHRaster,TrialSpikeData{2},'ZscoreFlag',0,...
    'GLMnMdlCompRuns',AnalysisOpts.PairNum,'UseFakeNeurons',AnalysisOpts.UseFakeNeurons,'GLMLambda',[0], ...
    'NormalizebyMax',0,'n_movavg',10,'TargetArea',ChAreas,'PlotResults',0);

% save figure for single cells 
if isempty(ChNum);ChTxt='';WantedDate='ALL';else;ChTxt=['_Ch' num2str(ChNum)];WantedDate=AnalysisOpts.RecDate; end
if nar>1; AreaTxt='ALL';else;AreaTxt=AnalysisOpts.Area2look{ChAreas};end
[~,~,SingCellAnaFigsFileName]=ManData.GetFileName([],['SingCellAna' AnalysisOpts.ExtraStr '_' AreaTxt ChTxt '_' AnalysisOpts.SpkCntStartFieldName 'Prstp' num2str(AnalysisOpts.ProcessingStep)],'SaveInResults',1,'WantedDate',WantedDate);
FigPrms.SaveFigSeries([],SingCellAnaFigsFileName,[SingCellAnaFigs],'SaveEachFrame',0);

%% do sinlge cell analysis in time and trials
if ~strcmp(AnalysisOpts.TrlSpkTimeFieldName,'AllTrials')
    SingCellAnaFigs_TrialTime=cell(1);
    SingCellAnaFigs_TrialTime{1}=NeuAnaFunc.SetupSingleCellAnalysis_TrialTime(PSTHRaster,TrialSpikeData{2},'ZscoreFlag',0);
    % save figure for single cells
    [~,~,SingCellAnaFigsFileName]=ManData.GetFileName([],['SingCellAna_TrialTime' AnalysisOpts.Area2look{1} '_' AnalysisOpts.SpkCntStartFieldName 'R' AnalysisOpts.RuleTxt],'SaveInResults',1,'WantedDate','ALL');
    FigPrms.SaveFigSeries([],SingCellAnaFigsFileName,[SingCellAnaFigs_TrialTime],'SaveEachFrame',1)
end
ElapsedTime=toc;
fprintf('\Elapsed time to run this pipeline:%0.2f(secs)',ElapsedTime)
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

