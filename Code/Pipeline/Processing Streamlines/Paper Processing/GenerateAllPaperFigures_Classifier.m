function GenerateAllPaperFigures_Classifier(Animal,AreaNum,SpkCntStartFieldNames,PSTHbin,AnalysisType,RunMode,varargin) % generates all of the data and figures of the paper (YOOOOOHOOOOO)
%AnalysisType is the analysis type for this run can be
%'ALL','Decoding','CrossDecoding','Learning','Comparision'
%RunMode:0 dont run shuffle, run main and plot; 1= run missing files only; 2=plot only 3=Run main but don't plot
global AnalysisOpts

AnalysisOpts.Project='Rule Representation';
AnalysisOpts.Animal='';
AnalysisOpts.RecDate='ALL';
AnalysisOpts.AnalysisType='Analysis Pipeline';

if strcmp(AnalysisType,'ALL');AnalysisType={'Decoding','CrossDecoding','Learning','Comparision','FineTune'};end

% paper conditions to look at per area Refer to
% SetAnalysisOptions_RuleRepresentations.m for full set of conditions
%% Color shape and response encoding
AnalysisOpts.PaperFigs.Decoding.SpkCntStartFieldName=SpkCntStartFieldNames;
AnalysisOpts.PaperFigs.Decoding.PSTHbin=PSTHbin;
%% crossdecoding conditions
AnalysisOpts.PaperFigs.CrossDecoding.SpkCntStartFieldName=SpkCntStartFieldNames;
AnalysisOpts.PaperFigs.CrossDecoding.PSTHbin=PSTHbin;
%% Fine Tune conditions
AnalysisOpts.PaperFigs.FineTune.SpkCntStartFieldName=SpkCntStartFieldNames;
AnalysisOpts.PaperFigs.FineTune.PSTHbin=PSTHbin;
%% Learning conditions
AnalysisOpts.PaperFigs.Learning.SpkCntStartFieldName=SpkCntStartFieldNames;
AnalysisOpts.PaperFigs.Learning.PSTHbin=PSTHbin;
%% Comparision Conditions
AnalysisOpts.PaperFigs.Comparision.SpkCntStartFieldName=SpkCntStartFieldNames;
AnalysisOpts.PaperFigs.Comparision.PSTHbin=PSTHbin;

% if we want to add extra conditions then process them right now
ParseParams(varargin);

if RunMode==6;return;end % we have updated what we want
%% now go through each stage and generate data for each animal
% loop through analysistype
for nAnaType=1:length(AnalysisType)
    AnaType=AnalysisType{nAnaType};
    fprintf(2,'\nGenerating data and figures for Analysis Type:%s',AnaType)
    AnalysisOpts.PaperFigs.(AnaType).JobID=cellfun(@(y) cellfun(@(x) ProcessData4Paper_Classifier(Animal,AreaNum,x,RunMode,AnaType,'SpkCntStartFieldName',y,'PopulationAna.PSTHbin',AnalysisOpts.PaperFigs.(AnaType).PSTHbin,varargin{:}),...
        AnalysisOpts.PaperFigs.(AnaType).TestName,'UniformOutput',0),AnalysisOpts.PaperFigs.(AnaType).SpkCntStartFieldName,'UniformOutput',0);
end
end
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
addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Core functions' FS]));
SetupAllVarsLessPath(DateNum)  %% set up the path and initialize vars
if ispc;AnalysisOpts.KickoffMyfuncRunned=1;end
end
