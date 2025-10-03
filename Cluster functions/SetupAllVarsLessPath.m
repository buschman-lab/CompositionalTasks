function SetupAllVarsLessPath(DateNum)
%%% sets up the path for all purpuses takes care of cluster analysis as
%%% well for various projects
global AnalysisOpts  
ManData=ManipulateData;

[AnalysisOpts.RootPath,AnalysisOpts.FS]=FileSepIdentifier(AnalysisOpts.RunOnCluster); % identify filesep

%% load all of the default options first then you can always overwrite after 
if strcmpi(AnalysisOpts.Project,'Rule Representation')
    SetAnalysisOptions_RuleRepresentation;   
    AnalysisOpts.PreDateTxt='';
elseif strcmpi(AnalysisOpts.Project,'Learning attentional templates')
    SetAnalysisOptions;
    AnalysisOpts.PreDateTxt='18';
elseif strcmpi(AnalysisOpts.Project,'Cortical Dynamics')
    SetAnalysisOptions_CorticalDynamics;   
    AnalysisOpts.PreDateTxt='';
end

DateNum=ManData.DetermineDateNum(DateNum); % what is the date and animal we are using
%DateNum=DateNum(1);

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
    if strcmpi(AnalysisOpts.Project,'Rule Representation')
            AnalysisOpts.Animal=AnalysisOpts.AnimalSet{strcmp(AnalysisOpts.RecDate,AnalysisOpts.DateSet)};
    elseif strcmpi(AnalysisOpts.Project,'Learning attentional templates')
            AnalysisOpts.Animal=AnalysisOpts.AnimalSet{strcmp(AnalysisOpts.RecDate(3:end),AnalysisOpts.DateSet)};
    end
else
    AnalysisOpts.Animal='ALL';
end
%% set random number generator seed
rng('shuffle');

%% generate path 
GenerateClusterPath(AnalysisOpts.AnalysisType) % genereate the path we need for this analysis 
%addpath(genpath(AnalysisOpts.ClusterPath)) % add the path for all coding

% addpath(AnalysisOpts.ProjectPath) % add the path forThis Project
% addpath(genpath(AnalysisOpts.CodePath)) % add the path forThis Project
addpath(genpath([AnalysisOpts.ToolboxPath 'ssh2_v2_m1_r7' AnalysisOpts.FS])) % add the path for toolboxes

%rmpath(genpath([AnalysisOpts.ToolboxPath AnalysisOpts.FS 'Fieldtrip']))

%% set up monitors so that we show figures where we want 
% close all
% f=figure('Units','normalized')
% %place it wherever you want in the 1st monitor
% f
% %check units
% %retreive location info
% pos1 = f.Position
% close all
% f=figure('Units','normalized')
% %place it wherever you want in the 1st monitor
% pos2 = f.Position
% close all
% %Evaluate the following code or add it in the startup.m for permanent settings 
% %in case you want to add it permanently, replace the pos1,2 with the values (hard-coded)
load([AnalysisOpts.DataSavePath filesep 'Core Data' filesep 'MonitorSettings.mat']);
MP = get(0, 'MonitorPositions');
pos=[ 0.1 0.1 0.8 0.8];
if size(MP, 1) == 1  % 1 monitor  
    set(0, 'defaultFigureUnits', 'normalized','defaultFigurePosition', pos)
else %2nd monitor
    set(0, 'defaultFigureUnits', 'normalized','defaultFigurePosition', pos) 
end