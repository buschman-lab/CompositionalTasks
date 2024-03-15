function GenerateClusterPath(varargin)
%generates all of the path directories 
global AnalysisOpts

RootPath=AnalysisOpts.RootPath;FS=AnalysisOpts.FS;

AnalysisOpts.ClusterPath=RootPath;
AnalysisOpts.ProjectPath=RootPath;
AnalysisOpts.ToolboxPath=[RootPath 'Toolboxes' FS ];
AnalysisOpts.RecDatePath=['Z:\Projects\Rule_Representation\Submission data\' ];%[RootPath 'Data' FS];
AnalysisOpts.DataPath=['Z:\Projects\Rule_Representation\Submission data\' ];%[RootPath  'Data' FS];
AnalysisOpts.SpikeDataPath=[];
AnalysisOpts.DataSavePath='Z:\Projects\Rule_Representation\Submission data\';%[RootPath 'Input Output Data' FS];
AnalysisOpts.ResultsSavePath=[RootPath 'Results' FS];
AnalysisOpts.CodePath=[RootPath  FS 'Code' FS];
AnalysisOpts.CuesPath=[RootPath 'Projects' FS 'Rule_Representation' FS 'Monkeys' FS 'Chico' FS 'Final Recording Code' FS 'ColorShapeMorphlineCorrectedFinal' FS];
AnalysisOpts.BhvMdlPath=AnalysisOpts.DataPath;
