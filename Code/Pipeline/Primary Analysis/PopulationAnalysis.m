function [TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet]=PopulationAnalysis(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
%% Explanation of the function goes here
% Pipline to run neural population analysis codes
% define some global functions to carry options and data
% input RunonCluster: runs on cluster
% varargin(1) number of recording date
% Sina Tafazoli timeless
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
DataFileName=[AnalysisOpts.DataSavePath 'Core Data' filesep 'SpikingData' filesep AnalysisOpts.Animal '_' AnalysisOpts.Area2look{1} '_' AnalysisOpts.TrlSpkTimeFieldName '_' AnalysisOpts.SpkCntStartFieldName '.mat'];
if sum(AnalysisOpts.ProcessingStep==[1 3 10]) & ~AnalysisOpts.GenerateCPIinfo
    % if we are looking at specific channels then get those channels
    [Chns,PSTHRaster,TrialSpikeData,ChAreas,RecSpecs]=TrialFunc.GetSpecificChData(DataFileName);
   % eval(['PSTHRaster_' num2str(AnalysisOpts.PopulationAna.PSTHbin) 'ms=PSTHRaster;']);
    if AnalysisOpts.ProcessingStep==10 % then make sure to load SACCADE_START information for the same neurons as well
        [~,ExtraPSTHRaster]=TrialFunc.GetSpecificChData([AnalysisOpts.DataSavePath 'Core Data' filesep 'SpikingData' filesep AnalysisOpts.Animal '_' AnalysisOpts.Area2look{1} '_' AnalysisOpts.TrlSpkTimeFieldName '_SACCADE_START.mat']);
    end
   % clear PSTHRaster
  %  eval(sprintf('PSTHRaster=PSTHRaster_%ims;',AnalysisOpts.PopulationAna.PSTHbin));
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
  %  TrialFunc.ParseChannelInformation(Chns,RecSpecs);
    
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
% if RunonCluster==0 
%     AnalysisOpts.DividSpockSubspace=0;AnalysisOpts.DividSpockClassifier=0;
%     AnalysisOpts.CalShuffleClassifier=0;AnalysisOpts.DividSpockClassifier_Cond=1;
% end
%if RunonCluster==1;warning off;end
% add bias to spock cond
AnalysisOpts.DividSpockClassifier_Cond=AnalysisOpts.DividSpockClassifier_Cond+ AnalysisOpts.ChunkMaxArray*AnalysisOpts.ChunkArrayBiasNum;
AnalysisOpts.NfigProcesStep=2*ones(1,20);
%PopulationAnaFigs=cell(1,AnalysisOpts.NfigProcesStep(AnalysisOpts.ProcessingStep));
ThisTestManual='';%'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';%'Learning3D_AxisDecoding';%'3D_Color_Shape_AxisV2';%'3D_Color_Shape_Axis';%'3D_Color_Shape_Response_BalRespDir';%;%'3D_Color_Shape_Axis';%'3D_FineTuneColorCategoryClassifier';%'3D_AxisDecodingSupp';%'';%'3D_AxisDecodingSupp';%'3D_AxisDecodingSupp';3D_Color_Response_XgenBalInCongV5_AreaComp%'3D_AxisDecodingSupp';%'Learning3D_AxisDecoding';%'3D_Response_Xgen_BalInCongAxisSupp';%'Learning3D_AxisDecoding';%'Learning3D_AxisDecoding';%'Learning3D_RuleR2SeqAnalysisLearningCtrlLearning';%'';%'Learning3D_Color_Color_Shape_Xgen_Angle_BlkEnd_Final';%'';%'';%'Learning3D_Color_ColorRB_Response_Xgen_Angle_BlkEnd';%'Learning3D_Color_Color_Response_Xgen_Angle';%'Learning3D_Shape_Color_Color_Xgen_Angle';%'' '3D_Color_Response_XgenBalInCongV5_AreaComp';%'Learning3D_Shape_Color_Response_Xgen_Orhto';%'3D_Color_Response_XgenBalInCongV5';%Learning3D_Shape_Color_Color_Compression_RB';%'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';%'Learning3D_Shape_Color_Rule_Xgen_AltRuleR1_RB';%'Learning3D_Shape_Color_Color_Compression_RB';%'3D_Color_Cat_Xgen';%'';%'';%'3D_Shared_Color_Response_EntropyR3';%'3D_Shared_Color_Response_EntropyR3';%Learning3D_Shape_Color_Rule_Xgen_AltRuleR1_RB';%;'3D_Color_Response_XgenR2R3Entropy';%'';%'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';%'Learning3D_Color_Response_Rule_Xgen_AltRule';%'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';%'3D_Color_Response_XgenBalInCongV2Test';%'3D_Color_Cat_Xgen';%'3D_Color_Response_XgenBalInCongV2Test';%'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';%'3D_Color_Response_XgenBalInCongV2';%'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB';%'Learning3D_Color_ShapeCorrInCorr_Rule_Xgen_AltRule';%'3D_Shared_Color_Response_Entropy';%'Learning3D_Color_Shape_Rule_Xgen_AltRule';%'3D_Color_Response_XgenBalInCong';%'Learning3D_Color_Color_Rule_Xgen_AltRuleR1';%'3D_ColorCategoryFineTuneSig';%'ResponseLocation_BalCorr';%'3D_ColorCategoryFineTuneSig';%'3D_Color_Response_Xgen';%'3D_Shape_Cat_Xgen_BalRespDir';%'Learning3D_Color_Shape_Rule_Xgen_AltRule';%'3D_Color_Cat_Xgen';%;%'Learning3D_Color_Shape_Rule_Xgen_AltRule';%'3D_FineTuneColorCategoryClassifier';%'3D_Color_Response_XgenR2R3Entropy';

if ~isempty(ThisTestManual) & RunonCluster==0 & sum(AnalysisOpts.ProcessingStep==[3 4 8])
    SetProcessingConds(AnalysisOpts.ProcessingStep,ManData.GetInd(AnalysisOpts.PopulationAna.Classifier_TaskNameSet,ThisTestManual));
    if  contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'Learning3D') & RunonCluster==0 & isempty(AnalysisOpts.DividSpockClassifier_TrlRng)
        AnalysisOpts.DividSpockClassifier_TrlRng=16*ones(1,7);
    else
        AnalysisOpts.DividSpockClassifier_TrlRng=ones(1,7);
    end
end

if ischar(AnalysisOpts.PairNum) & (sum(AnalysisOpts.ProcessingStep==[3 4 8 17]))
     SetProcessingConds(AnalysisOpts.ProcessingStep,ManData.GetInd(AnalysisOpts.PopulationAna.Classifier_TaskNameSet,AnalysisOpts.PairNum));
elseif ischar(AnalysisOpts.PairNum) & (sum(AnalysisOpts.ProcessingStep==[5]))
    SetProcessingConds(AnalysisOpts.ProcessingStep,ManData.GetInd(AnalysisOpts.PopulationAna.Classifier_ComparisionNameSet,AnalysisOpts.PairNum));
elseif ischar(AnalysisOpts.PairNum) & (sum(AnalysisOpts.ProcessingStep==[11 14 16 17]))
    SetProcessingConds(AnalysisOpts.ProcessingStep,ManData.GetInd(AnalysisOpts.SingCellAna.GLMMdlNameSet,AnalysisOpts.PairNum));
end

if AnalysisOpts.ProcessingStep==4 % if we are plotting classifiers then preload them
    [ClassifierResults,ClassifierOpts]=NeuAnaFunc.LoadClassiferResults(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'CalShuff',0,'CalShuffTrlOrder',0,'IncludeAllReps4ShuffFold',0);
     
    if contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'AxisDecoding')  %AnalysisOpts.CalShuffTrlOrderClassifier & 
        [ClassifierResults_ObservedTrlOrder,ClassifierOpts_ShuffTrlOrder_Loaded,~,ClassifierResults_ShuffTrlOrder]=NeuAnaFunc.LoadClassiferResults(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'CalShuff',1,...
            'MinimizeFileName',AnalysisOpts.MinimizeFileName4CalShuffTrlOrder,'CalShuffTrlOrder',1,'IncludeAllReps4ShuffFold',AnalysisOpts.IncludeAllReps4ShuffFold);
    else
            ClassifierResults_ObservedTrlOrder=[];ClassifierOpts_ShuffTrlOrder_Loaded=[];ClassifierResults_ShuffTrlOrder=[];
    end

    if (~strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'Modified_MannKendall') | ...
            (contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'Learning3D') & ...
            contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'Angle'))) & ...
            ~AnalysisOpts.SkipLoadingShuffleFile & ~contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'AxisDecoding')
        [ClassifierResults_Observed,ClassifierOpts_Shuff,~,ClassifierResults_Shuff]=NeuAnaFunc.LoadClassiferResults(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'CalShuff',1);
    else
        ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
     end

elseif AnalysisOpts.ProcessingStep==3 & ...
    contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'Learning3D')

    AnalysisOpts.PaperSpec.StrTime_SAMPLE_ON=-0.6; % start time before sample on for plotting
    AnalysisOpts.PaperSpec.EndTime_SAMPLE_ON=0.61; % end time before sample on for plotting
    ClassifierResults=[];ClassifierOpts=[];ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
    ClassifierResults_ObservedTrlOrder=[];ClassifierOpts_ShuffTrlOrder_Loaded=[];ClassifierResults_ShuffTrlOrder=[];
    AnalysisOpts.UseRep4Cluster=1; % use repetition for cluster runs only for shuffle

elseif AnalysisOpts.ProcessingStep==1 & AnalysisOpts.SingCellAna.SkipShuffleNMdlComp==1 % change timing if we doing some test only
     AnalysisOpts.PaperSpec.StrTime_SAMPLE_ON=-0.6; % start time before sample on for plotting
    AnalysisOpts.PaperSpec.EndTime_SAMPLE_ON=0.61; % end time before sample on for plotting   
  
    AnalysisOpts.PaperSpec.StrTime_SACCADE_START=-0.6; % start time before sample on for plotting
    AnalysisOpts.PaperSpec.EndTime_SACCADE_START=0.61; % end time before sample on for plotting

elseif AnalysisOpts.ProcessingStep==3 & ...
      contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'BalInCongV')
 
    ClassifierResults=[];ClassifierOpts=[];ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
    ClassifierResults_ObservedTrlOrder=[];ClassifierOpts_ShuffTrlOrder_Loaded=[];ClassifierResults_ShuffTrlOrder=[];
    AnalysisOpts.UseRep4Cluster=1; % use repetition for cluster runs only for shuffle      

elseif AnalysisOpts.ProcessingStep==5 % plot classifier comparision
    Task2CompareName=AnalysisOpts.PopulationAna.Classifier_ComparisionNameSet{AnalysisOpts.PairNum};
    Tasks2Compare=NeuAnaFunc.ClassifierComparisionOpts(Task2CompareName);
    CompAreaFlag=~isempty(Tasks2Compare.CompRef);
    
    % plot classifier analysis we have and compare them
    for TskCompId=1:Tasks2Compare.NComparisons
        Ind=strcmp(AnalysisOpts.PopulationAna.Classifier_TaskNameSet,Tasks2Compare.TaskInd{TskCompId});
        if CompAreaFlag & contains(Tasks2Compare.CompRef.Area,Tasks2Compare.Area{TskCompId}) & ~isempty(Tasks2Compare.CompRef.AreaRepNum)
            ClassifierOptsTemp=NeuAnaFunc.DefineClassifierTestOptions([AnalysisOpts.PopulationAna.Classifier_TaskNameSet{Ind} '_AreaComp']);
            FileNameSyntax=['_NeuRS' num2str(Tasks2Compare.CompRef.AreaNeuNum) '_C' num2str(Tasks2Compare.CompRef.AreaRepNum) '_' ClassifierOptsTemp.Name '_' Tasks2Compare.Area{TskCompId} '_' Tasks2Compare.SpkCntStartFieldName{TskCompId} '_' Tasks2Compare.TrlSpkTimeFieldName{TskCompId} '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
        elseif CompAreaFlag & contains(Tasks2Compare.CompRef.Area,Tasks2Compare.Area{TskCompId}) & isempty(Tasks2Compare.CompRef.AreaRepNum)
            ClassifierOptsTemp=NeuAnaFunc.DefineClassifierTestOptions([AnalysisOpts.PopulationAna.Classifier_TaskNameSet{Ind} '_AreaComp']);
            FileNameSyntax=['_NeuRS' num2str(Tasks2Compare.CompRef.AreaNeuNum) '_' ClassifierOptsTemp.Name '_' Tasks2Compare.Area{TskCompId} '_' Tasks2Compare.SpkCntStartFieldName{TskCompId} '_' Tasks2Compare.TrlSpkTimeFieldName{TskCompId} '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
        else
            ClassifierOptsTemp=NeuAnaFunc.DefineClassifierTestOptions(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{Ind});
            FileNameSyntax=['_' ClassifierOptsTemp.Name '_' Tasks2Compare.Area{TskCompId} '_' Tasks2Compare.SpkCntStartFieldName{TskCompId} '_' Tasks2Compare.TrlSpkTimeFieldName{TskCompId} '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
        end
        fprintf(2,'\nLoading clasifier comparision data %s',FileNameSyntax);
        ClassifierResults{TskCompId}=ManData.LoadVar('Classifier','ClassifierResults',FileNameSyntax,0,'WantedDate','ALL');
        ClassifierOpts{TskCompId}=ManData.LoadVar('Classifier','ClassifierOpts',FileNameSyntax,0,'WantedDate','ALL');
        ClassifierOpts{TskCompId}.FileNameSyntax=FileNameSyntax;
    end
    ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
    ClassifierResults_ObservedTrlOrder=[];ClassifierOpts_ShuffTrlOrder_Loaded=[];ClassifierResults_ShuffTrlOrder=[];
else 
    ClassifierResults=[];ClassifierOpts=[];ClassifierResults_Shuff=[];ClassifierResults_Observed=[];ClassifierOpts_Shuff=[];
    ClassifierResults_ObservedTrlOrder=[];ClassifierOpts_ShuffTrlOrder_Loaded=[];ClassifierResults_ShuffTrlOrder=[];
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
    'ClassifierOpts_ShuffTrlOrder_Loaded',ClassifierOpts_ShuffTrlOrder_Loaded,....
    'ClassifierResults_ShuffTrlOrder',ClassifierResults_ShuffTrlOrder,'ClassifierResults_ObservedTrlOrder',ClassifierResults_ObservedTrlOrder,...
    'MaxMatchTrialConds',AnalysisOpts.PopulationAna.MaxMatchTrialConds,'TakeMeanHyperplane',1,'IncludeAllReps4ShuffFold',AnalysisOpts.IncludeAllReps4ShuffFold};

%SetProcessingConds(3,43) % set processing conditions
TrialFunc.RevertCh_2look;
if AnalysisOpts.ProcessingStep~=17
    NeuAnaFunc.SetupNeuralAnalysis(PSTHRaster,TrialSpikeData{2},SpkCountTim,RasterTim,NeuralAnalysisOpts{:});
else

%% Run Angle Analysis with TDR
%     AngleAnalysisWithTDR(AnalysisOpts.PairNum,'SAMPLE_ON', 'PFC', [4], 1, 90);

%% Run TDR Analysis
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},4,1,60);
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},4,1,70);
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},4,1,80);
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},4,1,90);
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},4,1,'all');
%     
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},[],1,'all');
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},[],1,90);
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},[],1,80);
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},[],1,70);
    
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},4,1,'all');
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},5,1,'all');
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},3,1,'all');
%     TargetedDimReductionAnalysis(AnalysisOpts.SpkCntStartFieldName,AnalysisOpts.Area2look{1},[],1,'all');

%% Run Fitting Arima Model Analysis
% load data for shuffle
if 1
    if AnalysisOpts.CalShuffTrlOrderClassifier
        %% define parameters we want to check
        DimSet    =[0 0 1 2 3 1 2 1 2 1 2 0];
        LagSet    =[2 2 2 2 2 2 2 2 2 2 2 2];
        NStgSet   =[0 0 0 0 0 0 0 0 0 0 0 14];
        AvgDim2Set=[1 5 1 1 1 4 4 5 5 6 6 1];
        Ind=AnalysisOpts.DividSpockClassifier_Cond;
        Lag=LagSet(Ind);Dim=DimSet(Ind);AvgDim2=AvgDim2Set(Ind);WantedNStg=NStgSet(Ind);

        %% load and process data
        [ClassifierOpts_Shuff]=NeuAnaFunc.DefineClassifierTestOptions(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'CalShuff',1,'CalShuffTrlOrder',1);

        % define path options and check if we have preprocessed stat data
        [~,TargetFactorDim]=arrayfun(@(dim) NeuAnaFunc.getClassifierDimInfo(ClassifierOpts_Shuff,dim),[0 1 3],'UniformOutput',0);
        FileNameSpecs=@(x) [ClassifierOpts_Shuff.Name '_' TargetFactorDim{x}];
        AnalysisOpts.CodeTestDataPath=[AnalysisOpts.DataSavePath 'Code Test Data' filesep];
        if Dim==0;Dim2Look=1:3;else;Dim2Look=Dim;end
        Hsum=[];TargetSps_sum=[];RawDataFigTot=[];
        CorrMethodSet=[];CorrMethodSet={'TFPW_Y_p','TFPW_Y_slop'};
        for dim=Dim2Look
            SaveFileName=[AnalysisOpts.CodeTestDataPath FileNameSpecs(dim) '_D' num2str(dim) 'ShData'];
            DataFileExists=exist([SaveFileName '.mat'],'file');

            if ~DataFileExists
                if isempty(ClassifierResults_Observed)
                    % load shuffle trial order file first
                    [ClassifierResults_Observed,ClassifierOpts_Shuff,~,ClassifierResults_Shuff]=NeuAnaFunc.LoadClassiferResults(AnalysisOpts.PopulationAna.Classifier_TaskNameSet{AnalysisOpts.PairNum},'CalShuff',1,...
                        'MinimizeFileName',AnalysisOpts.MinimizeFileName4CalShuffTrlOrder,'CalShuffTrlOrder',1,'IncludeAllReps4ShuffFold',AnalysisOpts.IncludeAllReps4ShuffFold);
                end

                % process data
                ProcessTrlShuffStatTests(ClassifierResults_Observed,ClassifierResults_Shuff,ClassifierOpts_Shuff,Lag,AvgDim2,dim,FileNameSpecs,SaveFileName)
            elseif dim==Dim2Look(1)
                FileNameSyntax=['_Shuf_'  ClassifierOpts_Shuff.Name '_' AnalysisOpts.Area2look{1} '_' AnalysisOpts.SpkCntStartFieldName '_' AnalysisOpts.TrlSpkTimeFieldName '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
                ClassifierOpts_Shuff=ManData.LoadVar('Classifier','ClassifierOpts',FileNameSyntax,0,'WantedDate','ALL');
            end

            %% run correlation stat tests as well
           
            % plot results
            [Hsum,TargetSps_sum,NStgTxt,Ext,RawDataFig]=RunAllStattestsonData(ClassifierOpts_Shuff,'Accuracy',AvgDim2,dim,WantedNStg,CorrMethodSet,Hsum,TargetSps_sum);
            RawDataFigTot=[RawDataFigTot RawDataFig];

            % do CPI analysis if we are processing shape color rule
            if strcmp(ClassifierOpts_Shuff.TestName,'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB') & ...
                    ClassifierOpts_Shuff.IncludeAllClassifierInfo & dim==1
                CPIfeature={'CPI','ShapeAvgDist','ColorAvgDist'};
                for i=1:3
                    [Hsum,TargetSps_sum,~,~,RawDataFig]=RunAllStattestsonData(ClassifierOpts_Shuff,CPIfeature{i},AvgDim2,1,WantedNStg,CorrMethodSet,Hsum,TargetSps_sum);
                    RawDataFigTot=[RawDataFigTot RawDataFig];
                end
            end
        end
        if (strcmp(ClassifierOpts_Shuff.TestName,'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB') || ...
                strcmp(ClassifierOpts_Shuff.TestName,'Learning3D_Color_Response_Rule_Xgen_AltRule') || ...
                strcmp(ClassifierOpts_Shuff.TestName,'Learning3D_Shape_Color_Rule_Xgen_SameRule_RB_NoColDup')) 

            ProcessCorrelationStatTests(ClassifierOpts_Shuff,AvgDim2)
        end

       FigPrms.SaveFigSeries(['StatsSummary_' NStgTxt Ext '_' ClassifierOpts_Shuff.Name ],[AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep],[Hsum arrayfun(@(x) RawDataFigTot(x),1:length(RawDataFigTot),'uniformoutput',0)]);
      
       % save all figures so we don't miss anthing 
       mkdir([AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep ClassifierOpts_Shuff.Name 'AllStatTests'] )
       FigPrms.SaveCurrentFigs([ClassifierOpts_Shuff.Name 'CorrStatTests'],[AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep  ClassifierOpts_Shuff.Name 'AllStatTests' filesep],'SaveEachFrame',1)

    else
        error('This Analysis is only for trend order statisitcal test')
    end
end

%% generate test data so we can see if we can save some data
if 0
    % Desired size in GB
    desiredSizeGB = 36;

    % Number of bytes in 1 GB
    bytesPerGB = 1024^3;

    % Total size in bytes
    totalSizeBytes = desiredSizeGB * bytesPerGB;

    % Number of elements required for a double precision variable (8 bytes per element)
    numElements = totalSizeBytes / 8;

    % Create a large array of double precision numbers (fill with NaN or any value)
    largeVar = NaN(1, floor(numElements));  % NaN initializes the array with NaN values

    % Check the size of the variable
    disp(['The variable "largeVar" has been created with a size of ', num2str(totalSizeBytes / bytesPerGB), ' GB.']);
    disp(['Size of largeVar: ', num2str(numel(largeVar) * 8 / bytesPerGB), ' GB']);

    ManData.SaveVar('Classifier',largeVar,'largeVar', ...
        ['SizeTest' num2str(desiredSizeGB) 'GB'],'WantedDate','ALL','SaveFileVer','-v7.3');

end
% save figure for population analysis
% if isempty(ChNum);ChTxt='';WantedDate='ALL';else;ChTxt=['_Ch' num2str(ChNum)];WantedDate=AnalysisOpts.RecDate; end
% if nar>1; AreaTxt='ALL';else;AreaTxt=AnalysisOpts.Area2look{ChAreas};end
% [~,~,PopulationAnaFigsFileName]=ManData.GetFileName([],['SingCellAna' AnalysisOpts.ExtraStr '_' AreaTxt ChTxt '_' AnalysisOpts.SpkCntStartFieldName 'Prstp' num2str(AnalysisOpts.ProcessingStep)],'SaveInResults',1,'WantedDate',WantedDate);
% FigPrms.SaveFigSeries([],PopulationAnaFigsFileName,[PopulationAnaFigs],'SaveEachFrame',0);

% if we want to save factorize date
% FactorizedData=NeuAnaFunc.GetFactorizedData4ThisArea(PSTHRaster,TrialSpikeData{2},'UseFakeNeurons',AnalysisOpts.UseFakeNeurons,'n_movavg',10,'TargetArea',ChAreas);
% save PFC_ALL_SACCADE_START.mat FactorizedData AnalysisOpts -v7.3
end

%% end the function 
ElapsedTime=toc;
fprintf('\nElapsed time to run this pipeline:%0.4f(secs)',ElapsedTime);
end


%% put the functions to kick off here
function FS=KickoffMyfunc(RunonCluster,DateNum)
% sets up all of the initial options and path
global AnalysisOpts AnalysisFuncs AnalysisData

%% ok now first setup everything for cluster first because it needs path
AnalysisOpts.RunOnCluster=RunonCluster;

% if isfield(AnalysisOpts,'KickoffMyfuncRunned')
%     if AnalysisOpts.KickoffMyfuncRunned
%         FS=filesep;
%         return;
%     end % we have aleady runned this 
% end
if AnalysisOpts.RunOnCluster
    RootPath='/jukebox/buschman/';
    FS='/';
    warning on
else
    if ismac
        RootPath='/Volumes/buschman/';
        FS='/';
    elseif ispc
        RootPath='Z:\';
        FS='\';
    end
    warning on
end
% add core function path first so we can sturt up everything
addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Core functions' FS]));
SetupAllVars(DateNum)  %% set up the path and initialize vars
if ispc;AnalysisOpts.KickoffMyfuncRunned=1;end
end

function SetProcessingConds(ProcessingStep,PairNum)
global AnalysisOpts

AnalysisOpts.ProcessingStep=ProcessingStep;
AnalysisOpts.PairNum=PairNum;

end

