% runs a sequence of processing steps to process data for figure
function JobID=ProcessData4Paper_Classifier(Animal,AreaNum,TestName,RunMode,AnaType,varargin)
% @RunMode:0 dont run shuffle, run main and plot; 1= run missing files main and shuffle and then plot
% 2=plot only 3=Run main, no shuffle no plot 4=run missing main only no shuffle or plot (in subspace it is run shuffle only)
% 5: plot comparision only 6: collect comparision results
% 7: GLM: plot summery for each area 8: concatinate remaining files 9: run shuffle only
% 10: runs Targeted Dim reduction analysis
%TestName is the name of classifier test
%AnaType type of analysis 'ALL','Decoding','CrossDecoding','Learning','Comparision'

global AnalysisOpts

if contains(AnalysisOpts.ClassifierFunctiononClust,'Temp')
    NeuAnaFunc=NeuralAnalysisFuncsTemp;
else
    NeuAnaFunc=NeuralAnalysisFuncs;
end

ManData=ManipulateData;
CL=ClusterFuncs;
ParseParams(varargin)
%initialize variables
nCondsRun=cell(1,length(AreaNum));
nCondsRun_Shuff=cell(1,length(AreaNum));
nXTrlPnt=[]; JobID=[];
AnalysisOpts.AreaNum=AreaNum;

if RunMode==10 % run Tagrgeted Dimendionality analysis
    JobID.TDRana=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',17,varargin{:},...
        'RunWithDepenency',0),AreaNum,'UniformOutput',0);
    return
elseif ~contains(AnaType,'comparision','IgnoreCase',1) % if we are not doing comparisiom
    % Process this condition
    [ClassifierOpts]=NeuAnaFunc.DefineClassifierTestOptions(TestName); % gte the options for this test
    ParseParams(varargin)
    if ischar(TestName);TestName=ManData.GetInd(AnalysisOpts.PopulationAna.Classifier_TaskNameSet,TestName);end
    % get the conditions we are running on spock
    [nConds,nCondsShuffle,nXTrlPnt]=NeuAnaFunc.GetClassifierSpockRunConds(ClassifierOpts);

    % check what are the files that are missing and check for those
    for Ar=1:length(AreaNum)
        if (RunMode==1 | RunMode==4 | RunMode==9) & ~AnalysisOpts.ReWriteClassifierData% find missing files
            AnalysisOpts.Area2look=AnalysisOpts.AreaNames(AreaNum(Ar));
            if isempty(AnalysisOpts.Enforce_nCondsRun)
                [~,nCondsRun{Ar}]=NeuAnaFunc.ExistClassifierFiles(TestName,'CalShuff',0);
            end
            if isempty(AnalysisOpts.Enforce_nCondsRun_Shuff)
                [~,nCondsRun_Shuff{Ar}]=NeuAnaFunc.ExistClassifierFiles(TestName,'CalShuff',1,'CalShuffTrlOrder',AnalysisOpts.CalShuffTrlOrderClassifier);
            end
        else
            nCondsRun{Ar}=nConds;nCondsRun_Shuff{Ar}=nCondsShuffle;
        end
    end
    PlotProcessingStep=4;
else  % if we are plotting classifier comparision
    Tasks2Compare=NeuAnaFunc.ClassifierComparisionOpts(TestName);
    [ClassifierOpts]=NeuAnaFunc.DefineClassifierTestOptions(Tasks2Compare.TaskInd{1});
    ParseParams(varargin)
    if ischar(TestName);TestName=ManData.GetInd(AnalysisOpts.PopulationAna.Classifier_ComparisionNameSet,TestName);end
    PlotProcessingStep=5;
end

if AnalysisOpts.SweepClassifierConds % if we are sweeping conditions
    nCondsRun=repmat({1},length(AreaNum));
    nCondsRun_Shuff=cell(1,length(AreaNum)); 
elseif RunMode==8  
    nCondsRun=cell(1,length(AreaNum));
    nCondsRun_Shuff=cell(1,length(AreaNum));
elseif RunMode==0 || RunMode==3 || RunMode==4   % don't run shuffle
    nCondsRun_Shuff=cell(1,length(AreaNum));
elseif RunMode==9 % run shuffle only
    nCondsRun=cell(1,length(AreaNum));
elseif RunMode==2 % just plot or concatinate
    nCondsRun=cell(1,length(AreaNum));
    nCondsRun_Shuff=cell(1,length(AreaNum));
end

if ~isempty(AnalysisOpts.Enforce_nCondsRun)
    nCondsRun=AnalysisOpts.Enforce_nCondsRun;
end
if ~isempty(AnalysisOpts.Enforce_nCondsRun_Shuff)
    nCondsRun_Shuff=AnalysisOpts.Enforce_nCondsRun_Shuff;
end

%% run main classifier analysis
varargin=[varargin 'ExchangeableCalShuffClassifier',AnalysisOpts.ExchangeableCalShuffClassifier,...
    'DividSpockClassifier',AnalysisOpts.DividSpockClassifier,'GetOnlyShuffLabelsClassifier',0];

% run this condition without shuffling
AnalysisOpts.CalShuffleClassifier=0;
JobID.Run=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',3,varargin{:},...
    'DividSpockClassifier_Cond',nCondsRun{x==AreaNum},'DividSpockClassifier_TrlRng',cell2mat(nXTrlPnt),'RunWithDepenency',0),AreaNum,'UniformOutput',0);
if AnalysisOpts.SweepClassifierConds;return;end

% Concatinate files together
JobID.Concatinate=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',8,'CalShuffleClassifier',0,varargin{:},...
    'RunWithDepenency',1,'job_id_dep',JobID.Run{x==AreaNum}),AreaNum,'UniformOutput',0);

%% shuffle
% if we are running shuffle, then generate the shuffle distribution first and save it
% JobID.RunShuffleLabel=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',3,varargin{:},...
%     'CalShuffleClassifier',0,'GetOnlyShuffLabelsClassifier',1,'RunWithDepenency',0,'DividSpockClassifier_Cond',CheckContent(nCondsRun_Shuff{x==AreaNum}),'DividSpockClassifier_TrlRng',cell2mat(nXTrlPnt)),AreaNum,'UniformOutput',0);

% run this condition with shuffle now
% JobID.RunShuffle=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',3,varargin{:},...
%     'RunWithDepenency',1,'job_id_dep',JobID.RunShuffleLabel{x==AreaNum},'CalShuffleClassifier',1,'DividSpockClassifier_Cond',nCondsRun_Shuff{x==AreaNum},'DividSpockClassifier_TrlRng',cell2mat(nXTrlPnt)),AreaNum,'UniformOutput',0);

if AnalysisOpts.Classifier_TrlShuff_GenClassifierSpecs & RunMode==9 & AnalysisOpts.CalShuffTrlOrderClassifier==1
    varargin=[varargin 'Classifier_TrlShuff_GenClassifierSpecs',1];

    % then we are generating the file for classifier specs
    JobID.GenClassSpec=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',3,varargin{:},...
        'RunWithDepenency',0,'CalShuffleClassifier',1,'DividSpockClassifier_Cond',1,'DividSpockClassifier_TrlRng',cell2mat(nXTrlPnt)),AreaNum,'UniformOutput',0);
   
    % now run the classifiers
    varargin=[varargin 'Classifier_TrlShuff_GenClassifierSpecs',0];
    JobID.RunShuffle=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',3,varargin{:},...
        'RunWithDepenency',1,'job_id_dep',CL.ConcatinateJobIdDep(JobID.GenClassSpec{x==AreaNum}),'CalShuffleClassifier',1,'DividSpockClassifier_Cond',nCondsRun_Shuff{x==AreaNum},'DividSpockClassifier_TrlRng',cell2mat(nXTrlPnt)),AreaNum,'UniformOutput',0);
else
    varargin=[varargin 'Classifier_TrlShuff_GenClassifierSpecs',0];
    JobID.RunShuffle=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',3,varargin{:},...
        'RunWithDepenency',0,'CalShuffleClassifier',1,'DividSpockClassifier_Cond',nCondsRun_Shuff{x==AreaNum},'DividSpockClassifier_TrlRng',cell2mat(nXTrlPnt)),AreaNum,'UniformOutput',0);
end
%if AnalysisOpts.ExchangeableCalShuffClassifier==1;JobID.RunShuffle=JobID.Run;end % we are running shuffle with main job
% Concatinate shuffle files together
JobID.ConcatinateShuffle=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',8,varargin{:},...
    'CalShuffleClassifier',1,'RunWithDepenency',1,'job_id_dep',CL.ConcatinateJobIdDep(JobID.RunShuffle{x==AreaNum})),AreaNum,'UniformOutput',0);

%% just concatinate the results
% RunMode2 defines which detail operation we are performing in each case
% ConcatinateClassifier Files 1:only main 2: only shuffle 3: both shuffle and main
if RunMode==8
    if AnalysisOpts.RunMode2==1 % run main
        JobID.Concatinate=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',8,'CalShuffleClassifier',0,varargin{:},...
            'RunWithDepenency',0,'job_id_dep',JobID.Run{x==AreaNum}),AreaNum,'UniformOutput',0);
    elseif AnalysisOpts.RunMode2==2  % run shuffle
        JobID.ConcatinateShuffle=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',8,varargin{:},...
            'CalShuffleClassifier',1,'RunWithDepenency',0,'job_id_dep',JobID.RunShuffle{x==AreaNum}),AreaNum,'UniformOutput',0);
    elseif AnalysisOpts.RunMode2==3  % run main and shuffle
        JobID.Concatinate=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',8,'CalShuffleClassifier',0,varargin{:},...
            'RunWithDepenency',0,'job_id_dep',JobID.Run{x==AreaNum}),AreaNum,'UniformOutput',0);
        JobID.ConcatinateShuffle=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',8,varargin{:},...
            'CalShuffleClassifier',1,'RunWithDepenency',0,'job_id_dep',JobID.RunShuffle{x==AreaNum}),AreaNum,'UniformOutput',0);
    end
end

%% plot results
% Plot results with significant tests for this analysis
if RunMode==3 | RunMode==4 | RunMode==9 | RunMode==8 ;JobID.PlotResults=[];return;end
% run after everything has run
JobID.PlotResults=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',PlotProcessingStep,varargin{:},...
    'RunWithDepenency',1,'job_id_dep',CL.ConcatinateJobIdDep([JobID.ConcatinateShuffle(x==AreaNum) JobID.Concatinate(x==AreaNum)])),AreaNum,'UniformOutput',0);

% JobID.PlotResults=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',4,varargin{:}),...
%   AreaNum,'UniformOutput',0);
end
function out=CheckContent(data)% checks if there is content inside the matrix
if isempty(data);out=[];else;out=1;end
end