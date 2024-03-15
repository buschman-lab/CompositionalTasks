function [TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet]=Generatefigures(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
%% Explanation of the function goes here
% Pipline to run neural population analysis codes
% define some global functions to carry options and data
% input RunonCluster: runs on cluster
% varargin(1) number of recording date
global AnalysisOpts AnalysisData
tic
close all
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
AnalysisOpts.DateNum=DateNum;
ParseParams(varargin) % add all of the additional parameters we have added to the function

%% define classes necessary for this analysis
HelperFunc=HelperFuncs;
[TimeFreq,FigPrms,ManData,TrialFunc,FilterFunc,NeuAnaFunc]=HelperFunc.LoadAnalysisClasses; % load all of the classes we want
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
AnalysisOpts.Bhv.MaxNTrial=25;
AnalysisOpts.Bhv.NTrial2Swtch=25;
AnalysisOpts.Area2look=AnalysisOpts.AreaNames(AreaNum);
if isempty(AreaNum);AnalysisOpts.Area2look=AnalysisOpts.AreaNames;end
if isempty(AnalysisOpts.Rule_2look);AnalysisOpts.Rule_2look=[1 2 3];end

%% Detemine recording dates and report details of the analysis
DateNum=ManData.DetermineDateNum(DateNum); % what is the date and animal we are using
ManData.ReportAnalysisDetails;

%% load information about channels
Recs=DateNum;ExtraPSTHRaster=[];
DataFileName=[AnalysisOpts.DataSavePath AnalysisOpts.Animal '_' AnalysisOpts.Area2look{1} '_' AnalysisOpts.TrlSpkTimeFieldName '_' AnalysisOpts.SpkCntStartFieldName '.mat'];
if sum(AnalysisOpts.ProcessingStep==[1 3 10])
    % if we are looking at specific channels then get those channels
    [Chns,PSTHRaster,TrialSpikeData,ChAreas,RecSpecs]=TrialFunc.GetSpecificChData(DataFileName);
    if AnalysisOpts.ProcessingStep==10 % then make sure to load SACCADE_START information for the same neurons as well
        [~,ExtraPSTHRaster]=TrialFunc.GetSpecificChData([AnalysisOpts.DataSavePath 'Core Data' filesep 'SpikingData' filesep AnalysisOpts.Animal '_' AnalysisOpts.Area2look{1} '_' AnalysisOpts.TrlSpkTimeFieldName '_SACCADE_START.mat']);
    end
    fprintf('\nNumber of channels in area:%s =%i',AnalysisOpts.AreaNames{AreaNum},length(PSTHRaster))
    TrialFunc.ParseChannelInformationTrialSpikeData(TrialSpikeData);
    %TrialFunc.ParseChannelInformation(Chns,RecSpecs);  
    SpkCountTim=PSTHRaster(1).SpkCountTim{1};RasterTim=PSTHRaster(1).RasterTim{1};
    % now delete the fields we are not going to use 
    PSTHRaster=rmfield(PSTHRaster,{'RasterTim','SpkCountTim'});
    TrialSpikeData{1}=[]; % we are nt gonna use it
elseif sum(AnalysisOpts.ProcessingStep==[11 12]) %% Define what channels you want( saves memory but loads each channel)
    load(DataFileName,'Chns','RecSpecs');
    TrialFunc.ParseChannelInformation(Chns,RecSpecs);
    Recs=ManData.GetDateNum(AnalysisOpts.Ch_2look_RecDate);
    [Chns,RecSpecs]=ManData.RunFuncOnRec(TrialFunc,'FindWantedSpkChSet',Recs,{[AnalysisOpts.SpikeQuality2Look],AnalysisOpts.Area2look},4);
    if ~isempty(ChNum);RecChNum=AnalysisOpts.Ch_2look_RecChNum;Chns=cellfun(@(x) x(RecChNum),Chns,'uniformoutput',0); end % is we want 1 channel    
    PSTHRaster=[];TrialSpikeData=cell(1,2);
    %% Grab all of the spiking data from channels
    [SpikeData,RecSpecs]=ManData.RunFuncOnRec(HelperFunc,'getSpikeData',Recs,Chns(:,1:3),1);
    
    %% Define what type of trials you want
    [WantedBlockSpec]=ManData.RunFuncOnRec(TrialFunc,'FindWantedTrialTypeBlockSpec',Recs,{[],AnalysisOpts.Rule_2look,{'ALL','ALL'},{'ALL','ALL'},'ALL'},1);
    
    %% Take the data from the trials we care about
    [TrialSpikeData,RecSpecs]=ManData.RunFuncOnRec(TrialFunc,'GrabWantedSpksWantedTrial',Recs,[WantedBlockSpec,SpikeData],2);
    TrialSpikeData=ManData.ConcatinateStructs(TrialSpikeData,1); % concatinate data across recordings
    
    %% start analysis
    [PSTHRaster]=NeuAnaFunc.ProcessSpikeCountingRec(TrialSpikeData,AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.TrlSpkTimeFieldName,'PSTH_Bin',0.001*AnalysisOpts.PopulationAna.PSTHbin);
    ChAreas=unique(cell2mat(Chns(:,3)')); % find the areas of the channels we have
    nar=length(ChAreas);
    fprintf('\nNumber of channels in area:%s =%i',AnalysisOpts.AreaNames{AreaNum},length(PSTHRaster))
    SpkCountTim=PSTHRaster(1).SpkCountTim{1};RasterTim=PSTHRaster(1).RasterTim{1};
else % inistialize variables
    load(DataFileName,'Chns','RecSpecs');
    PSTHRaster=[];TrialSpikeData=cell(2);ChAreas=unique(cell2mat(Chns(:,3)'));
    TrialFunc.ParseChannelInformation(Chns,RecSpecs);
    SpkCountTim=[];RasterTim=[];
end

%% set up extra varibales for neural analysis
AnalysisOpts.DividSpockClassifier_Cond=AnalysisOpts.DividSpockClassifier_Cond+ AnalysisOpts.ChunkMaxArray*AnalysisOpts.ChunkArrayBiasNum;
AnalysisOpts.NfigProcesStep=[2*ones(1,20)];

if ischar(AnalysisOpts.PairNum) & (sum(AnalysisOpts.ProcessingStep==[3 4 8]))
     SetProcessingConds(AnalysisOpts.ProcessingStep,ManData.GetInd(AnalysisOpts.PopulationAna.Classifier_TaskNameSet,AnalysisOpts.PairNum));
elseif ischar(AnalysisOpts.PairNum) & (sum(AnalysisOpts.ProcessingStep==[5]))
    SetProcessingConds(AnalysisOpts.ProcessingStep,ManData.GetInd(AnalysisOpts.PopulationAna.Classifier_ComparisionNameSet,AnalysisOpts.PairNum));
end

if AnalysisOpts.ProcessingStep==4 % if we are plotting classifiers then preload them
    [ClassifierResults,ClassifierOpts]=NeuAnaFunc.LoadClassiferResults(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'CalShuff',0);
  
    if ~(strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'Modified_MannKendall') & contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'Learning3D'))
        [ClassifierResults_Observed,ClassifierOpts_Shuff,~,ClassifierResults_Shuff]=NeuAnaFunc.LoadClassiferResults(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'CalShuff',1);
    else
        ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
    end

elseif AnalysisOpts.ProcessingStep==3 & ...
    contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'Learning3D')

    AnalysisOpts.PaperSpec.StrTime_SAMPLE_ON=-0.6; % start time before sample on for plotting
    AnalysisOpts.PaperSpec.EndTime_SAMPLE_ON=0.61; % end time before sample on for plotting
    ClassifierResults=[];ClassifierOpts=[];ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
    AnalysisOpts.UseRep4Cluster=1; % use repetition for cluster runs only for shuffle

elseif AnalysisOpts.ProcessingStep==3 & ...
      contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'BalInCongV')
 
    ClassifierResults=[];ClassifierOpts=[];ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
    AnalysisOpts.UseRep4Cluster=1; % use repetition for cluster runs only for shuffle      

elseif AnalysisOpts.ProcessingStep==5 % plot classifier comparision
    Task2CompareName=AnalysisOpts.PopulationAna.Classifier_ComparisionNameSet{AnalysisOpts.PairNum};
    Tasks2Compare=NeuAnaFunc.ClassifierComparisionOpts(Task2CompareName);
    % plot classifier analysis we have and compare them
    for TskCompId=1:Tasks2Compare.NComparisons
        Ind=strcmp(AnalysisOpts.PopulationAna.Classifier_TaskNameSet,Tasks2Compare.TaskInd{TskCompId});
        ClassifierOptsTemp=NeuAnaFunc.DefineClassifierTestOptions(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{Ind});
        FileNameSyntax=['_' ClassifierOptsTemp.Name '_' Tasks2Compare.Area{TskCompId} '_' Tasks2Compare.SpkCntStartFieldName{TskCompId} '_' Tasks2Compare.TrlSpkTimeFieldName{TskCompId} '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
        ClassifierResults{TskCompId}=ManData.LoadVar('Classifier','ClassifierResults',FileNameSyntax,0,'WantedDate','ALL');
        ClassifierOpts{TskCompId}=ManData.LoadVar('Classifier','ClassifierOpts',FileNameSyntax,0,'WantedDate','ALL');
    end
    ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];

else 
    ClassifierResults=[];ClassifierOpts=[];ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
end

% define the options we want to run this analysis
NeuralAnalysisOpts={'ZscoreFlag',0,'ZscoreFlag',0,'GLMLambda',0,...
    'UseFakeNeurons',AnalysisOpts.UseFakeNeurons,'n_movavg',10,'TargetArea',ChAreas,'PSTHTimRef',AnalysisOpts.SpkParams.PSTHTimRef,...
    'PlotResults',0,'ExtraPSTHRaster',ExtraPSTHRaster,'GLMnMdlCompRuns',AnalysisOpts.DividSpockGLM_Cond,...
    'Classifier_TaskName',AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},...
    'SubspaceAna_TaskName',AnalysisOpts.PopulationAna.SubspaceAna_TaskNameSet{AnalysisOpts.PairNum},...
    'ZscoreFactorData',AnalysisOpts.ZscoreFactorData,'DetrendFactorData',AnalysisOpts.DetrendFactorData,...
    'PSTH_Bin',0.001*AnalysisOpts.PopulationAna.PSTHbin,'MaxMatchTrialConds',AnalysisOpts.PopulationAna.MaxMatchTrialConds,...
    'RunCrossTemporalClassifer',AnalysisOpts.RunCrossTemporalClassifer,'CalShuff',AnalysisOpts.CalShuffleClassifier,...
    'CalShuffTrlOrder',AnalysisOpts.CalShuffTrlOrderClassifier,'ClassifierResults_Loaded',ClassifierResults,...
    'ClassifierOpts_Loaded',ClassifierOpts,'ClassifierOptsShuff_Loaded',ClassifierOpts_Shuff,....
    'ClassifierResults_Shuff',ClassifierResults_Shuff,'ClassifierResults_Observed',ClassifierResults_Observed,...
    'MaxMatchTrialConds',AnalysisOpts.PopulationAna.MaxMatchTrialConds};

TrialFunc.RevertCh_2look;
NeuAnaFunc.SetupNeuralAnalysis(PSTHRaster,TrialSpikeData{2},SpkCountTim,RasterTim,NeuralAnalysisOpts{:});


ElapsedTime=toc;
fprintf('\nElapsed time to run this pipeline:%0.4f(secs)',ElapsedTime);
end
%% put the functions to kick off here
function FS=KickoffMyfunc(RunonCluster,DateNum)
% sets up all of the initial options and path
global AnalysisOpts AnalysisFuncs AnalysisData

%% ok now first setup everything for cluster first because it needs path
AnalysisOpts.RunOnCluster=RunonCluster;

if isfield(AnalysisOpts,'KickoffMyfuncRunned')
    if AnalysisOpts.KickoffMyfuncRunned
        FS=filesep;
        return;
    end % we have aleady runned this 
end
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
addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Rule Representation Project' FS 'Submission Code' FS 'Submission to Nature March 2024' FS 'Tafazoli et al 2024 Code' FS 'Core functions' FS]));
SetupAllVars(DateNum)  %% set up the path and initialize vars
if ispc;AnalysisOpts.KickoffMyfuncRunned=1;end
end

function SetProcessingConds(ProcessingStep,PairNum)
global AnalysisOpts

AnalysisOpts.ProcessingStep=ProcessingStep;
AnalysisOpts.PairNum=PairNum;

end

