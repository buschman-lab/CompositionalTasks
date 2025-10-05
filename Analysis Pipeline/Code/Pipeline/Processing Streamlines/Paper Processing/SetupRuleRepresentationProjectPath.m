function FS=SetupRuleRepresentationProjectPath
% sets up all of the initial options and path
global AnalysisOpts AnalysisFuncs AnalysisData

RunonCluster=0;DateNum=0;
%% ok now first setup everything for cluster first because it needs path
AnalysisOpts.RunOnCluster=RunonCluster;
AnalysisOpts.Project='Rule Representation';
AnalysisOpts.Animal='';
AnalysisOpts.AnalysisType='Analysis Pipeline';
AnalysisOpts.AnalysisFocus1  ='';  % define focus of this analysis
AnalysisOpts.SubAnalysisType1='';  % define type of sub analysis
AnalysisOpts.AnalysisPathName='';

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
addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Rule Representation Project' FS 'Analysis Pipeline' FS 'Code' FS]))
SetupAllVarsLessPath(DateNum)  %% set up the path and initialize vars
if ispc;AnalysisOpts.KickoffMyfuncRunned=1;end
end
