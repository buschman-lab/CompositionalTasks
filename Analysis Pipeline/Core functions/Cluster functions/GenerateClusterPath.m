function GenerateClusterPath(AnalysisType)
% generates all of the path directories we need for various projects
global AnalysisOpts

RootPath=AnalysisOpts.RootPath;FS=AnalysisOpts.FS;


AnalysisOpts.ClusterPath=[RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Rule Representation Project' FS 'Submission Code' FS 'Final Code for Nature Paper 2025 Tafazoli et al' FS];
AnalysisOpts.ProjectPath=[RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Rule Representation Project' FS 'Submission Code' FS 'Final Code for Nature Paper 2025 Tafazoli et al' FS];
AnalysisOpts.ToolboxPath=[RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Toolboxes' FS ];
AnalysisOpts.RecDatePath=[RootPath 'Projects' FS 'Rule_Representation' FS 'Data' FS AnalysisOpts.Animal '_Recording' FS 'Neural_Recording' FS  AnalysisOpts.RecDate FS 'ANALYSIS' FS];
AnalysisOpts.DataPath=[RootPath 'Projects' FS 'Rule_Representation' FS 'Data' FS];
AnalysisOpts.SpikeDataPath=[];
AnalysisOpts.DataSavePath=[RootPath 'Projects'  FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Rule Representation Project' FS AnalysisType FS 'Input Output Data' FS];
AnalysisOpts.ResultsSavePath=[RootPath 'Projects'  FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Rule Representation Project' FS 'Submission Code' FS 'Final Code for Nature Paper 2025 Tafazoli et al' FS AnalysisType FS 'Results' FS];
AnalysisOpts.CodePath=[AnalysisOpts.ProjectPath 'Analysis Pipeline' FS 'Code' FS];
AnalysisOpts.CuesPath=[RootPath 'Projects' FS 'Rule_Representation' FS 'Monkeys' FS 'Chico' FS 'Final Recording Code' FS 'ColorShapeMorphlineCorrectedFinal' FS];
AnalysisOpts.BhvMdlPath=[AnalysisOpts.DataSavePath 'Behavioral Model' FS];
AnalysisOpts.CodeTestDataPath=[AnalysisOpts.DataSavePath 'Code Test Data' FS];
