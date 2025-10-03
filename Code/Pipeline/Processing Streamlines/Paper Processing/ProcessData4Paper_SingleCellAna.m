% runs a sequence of processing steps to process data for figure
function JobID=ProcessData4Paper_SingleCellAna(Animal,AreaNum,TestName,RunMode,AnaType,varargin)
% GLM analysis pipleline 
% Runmode 0: fit each neuron without shuffle (PS11)
% Runmode 1: fit each neuron with it's shuffle (PS11), concatinate (PS15) it and create a summery(PS14) for each neuron(No Plots) 
% Runmode 2: Plot compact singlecell cahractristics with GLM(PS10)
% Runmode 5: run model comparision(PS12)
% Runmode 6 plot model comparision results(PS13)
% Runmode 7 plot all GLM result summery for each area(PS16)
% Runmode 8 concatinate GLM result summery for each area(PS14)
%Runmode 9 concatinate GLM result summery for each single Neurons(PS15)

%AnaType type of analysis 'ALL' 'GLM'

global AnalysisOpts
NeuAnaFunc=NeuralAnalysisFuncs; 
ManData=ManipulateData;
TrialFunc=TrialFuncs;
ParseParams(varargin)
%initialize variables
nCondsRun=cell(1,length(AreaNum));
nCondsRun_Shuff=cell(1,length(AreaNum));
nXTrlPnt=[]; JobID=[];JobID.PlotResults=[];

% Process this condition
if ischar(TestName);TestName=ManData.GetInd(AnalysisOpts.SingCellAna.GLMMdlNameSet,TestName);end

%% check what are the files that are missing and check for those
AnalysisOpts.DateNum=Animal;
for Ar=1:length(AreaNum)
    [ThisDateNum,RecChns,ThisRecChstxt,NChs]=TrialFunc.GetRecordingChannels(Animal,AreaNum(Ar),[],1);
    if RunMode==1  % find missing files
        [nCondsRun{Ar},nCondsRun_ShuffNeu{Ar},nCondsRun_Shuff{Ar}]=NeuAnaFunc.ExistGLMFiles(RecChns,TestName,0);
         nCondsRun{Ar}=find(~nCondsRun{Ar});nCondsRun_ShuffNeu{Ar}=find(~nCondsRun_ShuffNeu{Ar});
         nCondsRun_Shuff{Ar}=cellfun(@(x) find(~x),nCondsRun_Shuff{Ar},'UniformOutput',0);
    else
        nCondsRun{Ar}=RecChns;nCondsRun_ShuffNeu{Ar}=RecChns;
        nCondsRun_Shuff{Ar}=arrayfun(@(x) 1:AnalysisOpts.SingCellAna.GLMShuffleRuns,1:NChs,'UniformOutput',0);
    end
end

if RunMode==0   % don't run shuffle 
    nCondsRun_Shuff=cell(1,length(AreaNum));
elseif RunMode==2 % just plot
    nCondsRun=cell(1,length(AreaNum));
    nCondsRun_ShuffNeu=cell(1,length(AreaNum));
    nCondsRun_Shuff=cell(1,length(AreaNum));%arrayfun(@(x) cell(1),1:length(AreaNum),'uniformoutput',0);
end

%% Runmode 0: fit each neuron without shuffle or number of shuffles less or equal to 100 (PS11)
if  RunMode==0  
    % run this condition without shuffling
    AnalysisOpts.CalShuffleGLM=0;
    JobID.Run=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, nCondsRun{x==AreaNum}, x, TestName,'ProcessingStep',11,varargin{:},...
        'DividSpockGLM',1,'DividSpockGLM_Cond',nCondsRun{x==AreaNum},'DividSpockGLM_TrlRng',1),AreaNum,'UniformOutput',0);
    % % Concatinate files together
    % JobID.Concatinate=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',15,'CalShuffleGLM',0,varargin{:},...
    %     'DividSpockGLM',1,'RunWithDepenency',1,'job_id_dep',JobID.Run{x==AreaNum}),AreaNum,'UniformOutput',0);
end
%% Runmode 1: fit each neuron with it's shuffle (PS11), concatinate (PS15) it and create a summery(PS14) for each neuron(No Plots) 
if sum(RunMode==1)
    % run this condition with shuffle now
    [JobID.RunShuffle,RecChns]=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, nCondsRun_ShuffNeu{Ar}, x, TestName,'ProcessingStep',11,varargin{:},...
        'CalShuffleGLM',1,'CalShuffleGLM_Cond',nCondsRun_Shuff{x==AreaNum}(nCondsRun_ShuffNeu{Ar}),'CalShuffleGLM_TrlRng',1),AreaNum,'UniformOutput',0);
    
    % Concatinate shuffle files together and create a summery for each neuron
    [JobID.ConcatinateShuffle]=arrayfun(@(x) cellfun(@(Jobid) ProcessingStreamLinePopulationAnalysis(0, Animal, RecChns{x==AreaNum}(strcmp(Jobid,JobID.RunShuffle{x==AreaNum})), x, TestName,'ProcessingStep',15,varargin{:},...
        'CalShuffleGLM',1,'RunWithDepenency',1,'job_id_dep',Jobid),JobID.RunShuffle{x==AreaNum},'UniformOutput',0),AreaNum,'UniformOutput',0);
    
    % run this when you want to run concatinate independantly
    %[JobID.ConcatinateShuffle]=ProcessingStreamLinePopulationAnalysis(0, Animal, setdiff(RecChns, nCondsRun_ShuffNeu{Ar}), 3, TestName,'ProcessingStep',15,varargin{:},'CalShuffleGLM',1,'RunWithDepenency',0);
        
    % create a summery file for each neuron after we are done with
    % shuffleconcatinate and main data
%     [JobID.RunShuffle,RecChns]=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,...
%         'ProcessingStep',14,'RunWithDepenency',1,varargin{:}),AreaNum,'UniformOutput',0);
end

%% Runmode 2: Plot compact singlecell cahractristics with GLM(PS10)
if sum(RunMode==2)
    JobID.PlotResults=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',10,varargin{:},...
        'RunWithDepenency',0),AreaNum,'UniformOutput',0);
    %     JobID.PlotResults=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',4,varargin{:}),...
    %         AreaNum,'UniformOutput',0);
end

%% Runmode 5: run model comparision(PS12)
if sum(RunMode==5)
    AnalysisOpts.CalShuffleGLM=0;
    JobID.Run=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',12,varargin{:},...
        'DividSpockGLM',1,'DividSpockGLM_Cond','','DividSpockGLM_TrlRng',1),AreaNum,'UniformOutput',0);
end

%% Runmode 6 plot model comparision results(PS13)
if RunMode==6
    JobID.Run=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',13,varargin{:},...
        'DividSpockGLM',1,'DividSpockGLM_Cond','','DividSpockGLM_TrlRng',1),AreaNum,'UniformOutput',0);
end

%% Runmode 7 plot all GLM result summery for each area(PS16)
if sum(RunMode==7)
    JobID.PlotResults=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',16,varargin{:},...
        'RunWithDepenency',0),AreaNum,'UniformOutput',0);
end

%% Runmode 8 concatinate GLM result summery for each area(PS14)
if sum(RunMode==8)
    JobID.PlotResults=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',14,varargin{:},...
        'RunWithDepenency',0),AreaNum,'UniformOutput',0);
end

%% Runmode 9 concatinate GLM result summery for each single Neurons(PS15)
if sum(RunMode==9)
    AnalysisOpts.CalShuffleGLM=0;
    % concatinate for single cells
    JobID.ConcateSingle=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, nCondsRun{x==AreaNum}, x, TestName,'ProcessingStep',15,varargin{:},...
        'DividSpockGLM',1,'DividSpockGLM_Cond',nCondsRun{x==AreaNum},'DividSpockGLM_TrlRng',1),AreaNum,'UniformOutput',0);
  % concatinate for area
 %     JobID.PlotResults=arrayfun(@(x) ProcessingStreamLinePopulationAnalysis(0, Animal, [], x, TestName,'ProcessingStep',14,varargin{:},...
 %       'RunWithDepenency',1,'job_id_dep',JobID.ConcateSingle{x==AreaNum}),AreaNum,'UniformOutput',0);
end
