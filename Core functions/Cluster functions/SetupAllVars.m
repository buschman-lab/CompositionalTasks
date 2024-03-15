function SetupAllVars(DateNum)
%%% sets up the path for all purpuses takes care of cluster analysis as
%%% well for various projects
global AnalysisOpts  
ManData=ManipulateData;

[AnalysisOpts.RootPath,AnalysisOpts.FS]=FileSepIdentifier(AnalysisOpts.RunOnCluster); % identify filesep

%% load all of the default options first then you can always overwrite after 
SetAnalysisOptions_RuleRepresentation;
AnalysisOpts.PreDateTxt='';
DateNum=ManData.DetermineDateNum(DateNum); % what is the date and animal we are using
 
if ~ischar(DateNum) & DateNum>0
   AnalysisOpts.RecDate=[AnalysisOpts.PreDateTxt AnalysisOpts.DateSet{DateNum} ];
elseif ischar(DateNum) & str2double(DateNum)<40 & str2double(DateNum)>0
    AnalysisOpts.RecDate=[AnalysisOpts.PreDateTxt AnalysisOpts.DateSet{str2double(DateNum)} ];
elseif DateNum==0 | isempty(setdiff(AnalysisOpts.DateSet_2look,AnalysisOpts.DateNum))
    AnalysisOpts.RecDate='ALL';
else
    AnalysisOpts.RecDate=DateNum;
end
AnalysisOpts.RecDate=upper(AnalysisOpts.RecDate);

if ~strcmp(AnalysisOpts.RecDate,'ALL')
    AnalysisOpts.Animal=AnalysisOpts.AnimalSet{strcmp(AnalysisOpts.RecDate,AnalysisOpts.DateSet)};
else
    AnalysisOpts.Animal='ALL';
end
%% set random number generator seed
rng('shuffle');

%% generate path 
GenerateClusterPath(AnalysisOpts.AnalysisType) % genereate the path we need for this analysis 
addpath(AnalysisOpts.ProjectPath) % add the path forThis Project
addpath(genpath(AnalysisOpts.CodePath)) % add the path forThis Project
addpath(genpath(AnalysisOpts.ToolboxPath)) % add the path for toolboxes
 