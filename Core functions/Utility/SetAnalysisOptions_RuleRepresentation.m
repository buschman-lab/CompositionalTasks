function SetAnalysisOptions_RuleRepresentation(varargin)

% Sinopsys: Lookat_% indicate the variable we want to look at in the set

global AnalysisOpts

%% define classes we might need 
FigParams=fig_params;

%% define project 
AnalysisOpts.Project='Rule Representation';
AnalysisOpts.PreDateTxt='';

%% Path params
AnalysisOpts.DateSet={'061321','062321','062521','070421','070921','071121','071321','072521','073021','080121','080421','080621','081021','081721','081921',...
    '060818','061118','061418','061518','061818','062018','062618','062718','063018','070418','070518','070818','071018','071118','071218','071718','071818','072018','072118'}; % Chico and Silas dataset
AnalysisOpts.AnimalSet={'Chico','Chico','Chico','Chico','Chico','Chico','Chico','Chico','Chico','Chico','Chico','Chico','Chico','Chico','Chico',...
    'Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas','Silas'};
AnalysisOpts.SilasRecNums=find(strcmp(AnalysisOpts.AnimalSet,'Silas'));
AnalysisOpts.ChicoRecNums=find(strcmp(AnalysisOpts.AnimalSet,'Chico'));
AnalysisOpts.DateSet_2look=[1:length(AnalysisOpts.DateSet)];
AnalysisOpts.ReRefChannels=0; % No rereferencing
%AnalysisOpts.ReRefChannels=1; % rereference channels to the mean of each well with only with spike channels
%AnalysisOpts.ReRefChannels=2; % rereference channels to the mean of each well with all of the channels
%AnalysisOpts.ReRefChannels=3; % rereference channels to one of the channels that doesn't have a neuron
AnalysisOpts.ReRefChannelsStr={'','ReRef','ReRefAll','ReRefSing'}; % names corresponding to each of the above conditions of ReRefing
AnalysisOpts.ReWriteData=0; % should we rewrite data even if it exists
AnalysisOpts.ReWriteMotifData=0; % should we rewrite motif data even if it exists
AnalysisOpts.ReWriteWaveData=0; % should we rewrite wavelet data even of ot exists
AnalysisOpts.UseFakeNeurons=0;% are we using fake data 
AnalysisOpts.ReWriteClassifierData=1; % are we rewriting classifier data
% initialize current channel data 
AnalysisOpts.CurrentCh_RecDate=[];
AnalysisOpts.CurrentCh_Animal=[];
AnalysisOpts.CurrentCh =[];
AnalysisOpts.CurrentChClust=[];
AnalysisOpts.DateNum=[];
AnalysisOpts.ExtaStrAnalysis=''; % extra string to add to the name of the file 
if ~isfield(AnalysisOpts,'KickoffMyfuncRunned')
   AnalysisOpts.KickoffMyfuncRunned=0; % have we already runned kick off my func
end
AnalysisOpts.Neu2Send2Spock=[]; % occacianlly if you want to only send specific neurons to spock instead of all
AnalysisOpts.UseIndividualChName=[]; % when we are saving or loading are we using indivdual channel names 
AnalysisOpts.RunDummyFile=0; % if we are running a dummy file to just test the data 

%% General Parameters
AnalysisOpts.UseDataPointer=0; % use data pointer (matfile) instead of loading the file
AnalysisOpts.ProcessingStep=1; % what is our processing step. Depends on the function 
AnalysisOpts.SaveEachFrame=1; % are we saving each frame 
AnalysisOpts.SaveEachSubplots=1; % are we saving each subplot

%% Cluster computing
AnalysisOpts.local_bucket='/jukebox/buschman/Users/Sina/Rule\ Representation\ Project/Server\ Analysis\ forlders/Rule\ Representation\ Cluster\ Analysis\ Functions/';
AnalysisOpts.script_path='Z:\Users\Sina\Rule Representation Project\Server Analysis forlders\Rule Representation Cluster Analysis Functions\';
AnalysisOpts.function_path='/jukebox/buschman/Projects/Rule_Representation/ElecPhys_Analysis/Rule Representation Project/Analysis Pipeline/Code/Pipeline/';
AnalysisOpts.CompResource='spock'; % can be della or spock
AnalysisOpts.RunWithDepenency=0;
AnalysisOpts.SpockPass=''; % spock password
AnalysisOpts.username=''; %spock username
AnalysisOpts.job_id_dep='';
AnalysisOpts.Spock_TaskName=''; % what is the current task name we are sending
AnalysisOpts.Spock_ExtraInfo=''; % what is the current extra information we are providing
AnalysisOpts.UseRep4Cluster=0; % are we using repetition for cluster. This is ONLY used with cluster condition 3 and shuffle 
AnalysisOpts.ClassifierFunctiononClust='PopulationAnalysis'; % what function are we running on cluster for classifier
AnalysisOpts.Enforce_nCondsRun=[]; % enforce these runs for main condition 
AnalysisOpts.Enforce_nCondsRun_Shuff=[]; % enforce these runs for shuffle
AnalysisOpts.ThisIsSinaPC=0; % Is this task being run on my PC? then don't as for crednetials

%% General Parameters
AnalysisOpts.UseDataPointer=0; % use data pointer (matfile) instead of loading the file

%% Area and Trial  paramters 
AnalysisOpts.SpikeQuality2Look=1; % 1 single cell 2% multi 3% usable 4% low quality or time cut
AnalysisOpts.RuleSet=[1 2 3]; % which sets of Rules do you want to look at
AnalysisOpts.ChannelArea.PFC=''; % cahnnels belongin gto each area
AnalysisOpts.ChannelArea.Str='';
AnalysisOpts.ChannelArea.IT='';
AnalysisOpts.ChannelArea.FEF='';
AnalysisOpts.ChannelArea.LIP='';
AnalysisOpts.PairNum=''; % this is a free var used to sweep any parameter
AnalysisOpts.CurrentCh=''; % which channel is being used right now
AnalysisOpts.CurrentChClust=NaN; % what is the cluster of the current channel
AnalysisOpts.Ch='';
AnalysisOpts.ChArea='';
AnalysisOpts.RecDate='';
AnalysisOpts.Animal='';
AnalysisOpts.Rule_2look=''; % which rule to look at
AnalysisOpts.File_2look=''; % file to open and look into instead of redoing everything
AnalysisOpts.Ch_2look='';  % which Ch do you want to look at in this analysis
AnalysisOpts.Area_2look={'PFC','Striatum','IT','FEF','LIP' }; %(DO NOT CHANGE THIS)
AnalysisOpts.AreaNames={'PFC','Striatum','IT','FEF','LIP' };%(DO NOT CHANGE THIS)
AnalysisOpts.Trial.TrialType='ALL';   % which type of trials do we want to look at(DO NOT CHANGE THIS)
AnalysisOpts.Trial.TrialsToAnalyz='ALL';  %% put 'ALL' if you want all of the trials (DO NOT CHANGE THIS)
AnalysisOpts.Trial.MaxTrialsToAnalyz='ALL';   %% how many trials to load (DO NOT CHANGE THIS)
AnalysisOpts.Trial.Rule=1;  %Look at rule 1
AnalysisOpts.AreaWells=[1 1 2 1 2]; % 1 is frontal well and 2 is parietal well
AnalysisOpts.WellNames={'FrontalWell','ParietalWell'}; % name of the wells

%% plotting 
AnalysisOpts.Plotting.NRunningAverage=20; %% how many trials to use for running average
AnalysisOpts.Plotting.SubplotStruct=[4 4]; % how many windos for subplot
AnalysisOpts.Plotting.font_size=15;  % plotting font size
AnalysisOpts.SavePlot=1;   % save the plots
AnalysisOpts.SaveData=1;   % save the data we are generating
AnalysisOpts.ShowPlot=1;   % show the plots we are generating
AnalysisOpts.ShowStatPvalinPlot=1; % are we showing the statistical test p-val in the plots

%% Neural Analysis
AnalysisOpts.NeuralAnalysis.StdTH2look=4;     %%%% what is the TH of spike detection we are looking at
AnalysisOpts.NeuralAnalysis.SamplingFrequency=40000;
AnalysisOpts.NeuralAnalysis.SpikeQualityTH=2; 
AnalysisOpts.NeuralAnalysis.UseWorkingPeriod=1;  %%% are we only using the period of working to do the spike sorting 
 AnalysisOpts.TrialTimesFields=[{'START_TRIAL'}, {'END_TRIAL'},'FIXATE_ACQUIRED','SAMPLE_ON','SACCADE_START','SACCADE_STOP',...
     'CORRECT_TRIAL','REWARD_START','REWARD_END','RULE_SWITCH_ON','RULE_SWITCH_OFF','RT_EyeTracker','CONDITION_NUM_OFFSET',...
     'SAMPLE_NUM_OFFSET','RESLOC_NUM_OFFSET'];
 
%% SpikeCount and spk detection Parameters
AnalysisOpts.SpkParams.PSTH_Bin=0.1;%% bin size
AnalysisOpts.SpkParams.PSTH_BinShift=0.01;  %% shif of bins
AnalysisOpts.SpkParams.PeriodLength=1.5;  %% % NEVER CHANGE THIS
AnalysisOpts.SpkParams.BaselineDelay=0.6; %%% look at 500ms before the start % NEVER CHANGE THIS
AnalysisOpts.SpkParams.RuleSelectvityPeriod=0.25;
AnalysisOpts.SpkParams.StartFieldName='SAMPLE_ON';
AnalysisOpts.SpkParams.StopFiledName='CORRECT_TRIAL'; % correct trial and reward are almost at the same time 'Reward was added on 063018 for Silas
AnalysisOpts.SpkParams.PSTHTimRef='leading'; % which timing are we locking our PSTH? 'leading','trailing','centered'
AnalysisOpts.ExtractSpikes=1; 
AnalysisOpts.ChracPlot.detect_fmin=300;
AnalysisOpts.ChracPlot.detect_fmax=3000;
%[AnalysisOpts.b,AnalysisOpts.a]=ellip(2,0.1,40,[AnalysisOpts.ChracPlot.detect_fmin AnalysisOpts.ChracPlot.detect_fmax]*2/AnalysisOpts.NeuralAnalysis.SamplingFrequency); 

%% Trial Timing Analysis
AnalysisOpts.TrialTiming.StartFieldName='SAMPLE_ON';
AnalysisOpts.TrialTiming.StopFiledName='CORRECT_TRIAL'; % correct trial and reward are almost at the same time 'Reward was added on 063018 for Silas
AnalysisOpts.TrialTiming.PeriodLength=1.5;  %%% Whole period duration
AnalysisOpts.TrialTiming.BaselineDelay=1; %%% look at x ms before the start
AnalysisOpts.TrialTiming.RuleSelectvityPeriod=0.25;

%% LFP Analysis Filter specs
AnalysisOpts.LFPParams.FilterOpts.PassBand=[4 10];  % Passband Frequency
AnalysisOpts.LFPParams.FilterOpts.StopBand=200;  % Stopband Frequency
AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency=1000;  
AnalysisOpts.LFPParams.FilterOpts.FilterOrder=4;   %% filter order
AnalysisOpts.LFPParams.FilterOpts.FilterType = 'IIR'; %Can be FIR or IIR
AnalysisOpts.LFPParams.FilterOpts.PassBandRipple = 0.1; % Passband Ripple (dB)
AnalysisOpts.LFPParams.FilterOpts.StopBandAttenuation = 40; % Stopband Attenuation (dB)
AnalysisOpts.LFPParams.FilterOpts.FilterDesign='bandpass';
AnalysisOpts.LFPParams.FreqToPlot = [10.^[log10(200):0.02:log10(3000)]];
%% Data Preprocessing
AnalysisOpts.PreProcess.MaxNTrial=[]; % N trials how long we want the maximum motif time to be 
AnalysisOpts.PreProcess.MaxTime  ='ALL'; % sec how long we want the maximum motif time to be. 'ALL' for entire data
AnalysisOpts.PreProcess.StrTime  =0;%'ALL'; % sec when to start looking at the data

%% Motif and Time Freq Analysis
AnalysisOpts.MotifAnalysis.FsWave=100; %  wavelet frequency sampling
AnalysisOpts.MotifAnalysis.FsWaveTarg=50; % target wavelet frequency
AnalysisOpts.MotifAnalysis.FsLFP=AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency;
AnalysisOpts.MotifAnalysis.MaxNTrial=''; % sec how long we want the maximum motif time to be 
% motif discovery
AnalysisOpts.MotifAnalysis.PowerMethod='NormPower'; % cwt power calc method can be 'Power','NormPower' also
AnalysisOpts.MotifAnalysis.CwtFreqScale='LinScale';   % how are we placing frequencies together 'LinScale' or 'LogScale'
AnalysisOpts.MotifAnalysis.Navg=1;
AnalysisOpts.MotifAnalysis.FreqStep=1;   % number for freq steps for wavelet 
AnalysisOpts.MotifAnalysis.flimit=[0 100]; % limits of Freq we want to look at
AnalysisOpts.MotifAnalysis.WaveletWidth=310; % wavelength width that we have saved the file with
AnalysisOpts.MotifAnalysis.WaveletMethod='VariableWidth';
AnalysisOpts.MotifAnalysis.K=15; % K
AnalysisOpts.MotifAnalysis.L_ms=1; % L is in sec
AnalysisOpts.MotifAnalysis.maxiter=300; % max iterations
AnalysisOpts.MotifAnalysis.maxiterRefit=300; % max iterations for refit
AnalysisOpts.MotifAnalysis.lambda=0.0002; %lambda
AnalysisOpts.MotifAnalysis.lambdaL1H=1;    %L1 sparsity parameter; Increase to make H's more sparse
AnalysisOpts.MotifAnalysis.lambdaL1W=0;    %L1 sparsity parameter; Increase to make W's more sparse
AnalysisOpts.MotifAnalysis.lambdaOrthoH=1; %Encourages events-based factorizations
AnalysisOpts.MotifAnalysis.lambdaOrthoW=0; % ||Wflat^TWflat||_1,i~=j; ; Encourages parts-based factorizations
AnalysisOpts.MotifAnalysis.MaxTime=60;     % chunk size of motif discovery 
AnalysisOpts.MotifAnalysis.Shift=0;        %Shift factors to center; Helps avoid local minima
AnalysisOpts.MotifAnalysis.ChunkData=1;    %if we have to chunk data intop MaxTime chunks
AnalysisOpts.MotifAnalysis.TestGeneralize=0; % do we want to test if we can generalize to test
AnalysisOpts.MotifAnalysis.ValidateNrep=1; % how many repetition for validation 
% motif Clustring
AnalysisOpts.MotifAnalysis.Phenograph_K=5;   % K of phenograph
AnalysisOpts.MotifAnalysis.MaxClust=30;    % Max number of clusters for heirarchical clustring 
AnalysisOpts.MotifAnalysis.XCorrTensorMethod ='';% 'DTW','' for 'Xcorr', how to calculate correlation between motifs 
AnalysisOpts.MotifAnalysis.ClustringMethod ='Hierarchical';%'EncoderPhenograph';%'PhenoCluster';  % how are we doing the clustring 
AnalysisOpts.MotifAnalysis.TemplateMotif ='Best';%'random' how do we choose template motif;"best" min dist to rest;  
AnalysisOpts.MotifAnalysis.AlignMotifs ='xcorr';%'xcorr' how do we align motifs to create core motifs "dtw' or 'max xcorr'
AnalysisOpts.MotifAnalysis.AlignMotifs ='xcorr';%'xcorr' how do we align motifs to create core motifs "dtw' or 'max xcorr'
AnalysisOpts.MotifAnalysis.CoreMotifsPerc=0.1; % percentage of motifs in the cluster that are avgeraged to find core motifs
AnalysisOpts.MotifAnalysis.NchunksXCorrTensor=1000; % number of chnuks for XCorr on cluster
AnalysisOpts.MotifAnalysis.NSampledMotifs=30000;   % how many motifs are we sampling instead of all of them(has to be EVEN)
AnalysisOpts.MotifAnalysis.NSampledMotifsALL=20000;   % how many motifs are we storing from all of the recordings (Has to be EVEN)
AnalysisOpts.MotifAnalysis.NSampledMotifsSubSet=10000;   % how many motifs are we sampling for clustering(has to be EVEN)
AnalysisOpts.MotifAnalysis.MotifSamplingMethod='WellEqual'; % how do we sample the motifs for clustring 'WellEqual' 
%is equal number of samples per well, 'Random' take motifs randomly from everywhere
AnalysisOpts.MotifAnalysis.ProcessStep=1; % which step of processing are we
AnalysisOpts.MotifAnalysis.XcorrTimInt=0.5; % what is the time interval we care about H xcorr
AnalysisOpts.MotifAnalysis.ShiftWH2ZeroTH=0.3; % threshold to find the starting point of a motif
% param sweep 
AnalysisOpts.MotifAnalysis.ParamSweepFile=0; % should we use the parameters in param sweep file?

%% Coherence Analysis
AnalysisOpts.CohAnalysis.NAvgCoh=5; % number of trials to average for coehrence analysis 

%% behavior with respect to neural analysis 
AnalysisOpts.Bhv.UseModel=1; % are we using behavioral model 
AnalysisOpts.Bhv.MaxNTrial=20;
AnalysisOpts.Bhv.NTrial2Swtch=20;
AnalysisOpts.Bhv.BhvModelTypes={'Hybrid','InferenceAxisFeature'}; % type of behavioral models Flora has generated

%% Single cell Analysis 
AnalysisOpts.SingCellAna.GLMnMdlCompRuns=20; % how many runs of GLM model for model comparision 
AnalysisOpts.SingCellAna.GLMShuffleRuns=1000; % how many runs of GLM for shuffling
AnalysisOpts.SingCellAna.GLM_UseSingFactorShuff=1; % are we using a single factor to shuffle the data for GLM
AnalysisOpts.SingCellAna.nModels2Test=14; % how many GLM models in total are we testing 
AnalysisOpts.SingCellAna.GLMLamdaVals=[1 25 50 75 99]; % list of lambda indexes we want to look at
AnalysisOpts.SingCellAna.GLMLamdaSet=[0.000151616578650318,0.000166399007897004,0.000182622705746232,0.000200428194107491,0.000219969695602974,0.000241416469370147,0.000264954277103429,0.000290786992033118,0.000319138364777787,0.000350253961365093,0.000384403290207242,0.000421882136453918,0.000463015123941747,0.000508158526930682,0.000557703355981194,0.000612078744700694,0.000671755666693551,0.000737251014909138,0.000809132078721281,0.000888021457517436,0.000974602453356731,0.00106962498940552,0.00117391210541308,0.00128836708648821,0.00141398128692269,0.00155184271682794,0.00170314546595805,0.00186920004634406,0.00205144474332217,0.00225145807327330,0.00247097245597682,0.00271188922000190,0.00297629507110602,0.00326648016628192,0.00358495795000196,0.00393448692447167,0.00431809454245623,0.00473910342962801,0.00520116016356065,0.00570826685863947,0.00626481583046175,0.00687562763997312,0.00754599284686044,0.00828171783384876,0.00908917509881084,0.00997535845029513,0.0109479435845496,0.0120153545687308,0.0131868368061435,0.0144725371155008,0.0158835916178146,0.0174322221921474,0.0191318423356824,0.0209971733450192,0.0230443718250023,0.0252911696295030,0.0277570274462550,0.0304633033560217,0.0334334378260743,0.0366931567403054,0.0402706942245284,0.0441970371969683,0.0485061937621234,0.0532354877726974,0.0584258821109561,0.0641223334896171,0.0703741818453888,0.0772355776978981,0.0847659511755862,0.0930305267710504,0.102100888284397,0.112055598847870,0.122980881402122,0.134971365518077,0.148130907033002,0.162573487600086,0.178424201942995,0.195820341366521,0.214912582908185,0.235866294430657,0.258862966959068,0.284101786669447,0.311801360144090,0.342201607837194,0.375565842151060,0.412183048122148,0.452370386471222,0.496475940697783,0.544881732012821,0.598007028223789,0.656311975231779,0.720301582595187,0.790530097680913,0.867605806289907,0.952196301337780,1.04503426522527,1.14692381598267,1.25874747215557,1.38147379675824,1.51616578650318];
AnalysisOpts.SingCellAna.GLMLamda=exp(-8:0.5:-4);%[1e-4 1e-3 1e-2];%[0 1e-6 1e-5 1e-4 1e-3 1e-2]; % list of lambda values 
AnalysisOpts.SingCellAna.SaveGLMmdls=0; % are we saving the generated GLM models?
AnalysisOpts.SingCellAna.GLMMdlNameSet={'SensoryMLSing','SensoryMLCatSing','SensoryMLCatObjSing','SensoryMotorInteractMdl','HybridbhvMdlFull','InferAFbhvMdlFull'}; % what are the names of GLM models we are testing
AnalysisOpts.SingCellAna.GLMMdlName2Test={'SensoryCat'};%,'SensoryCat'};%{'SensoryMotorInteractMdl'};
AnalysisOpts.SingCellAna.GLMMdlComp2Test={'CompareSensoryCatModel'};%CompareSensoryModels
AnalysisOpts.CalShuffleGLM=0; % are we running shuffle test on GLM 
AnalysisOpts.CalShuffleGLM_Cond=[]; % what is the condition we are running the shuffle GLM now
AnalysisOpts.CalShuffleGLM_TrlRng=[]; % what is the trial range we are runningn shuffle GLM now
AnalysisOpts.DividSpockGLM=1; % are we running  GLM on spock
AnalysisOpts.DividSpockGLM_Cond=[]; % what is the condition we are running the GLM now
AnalysisOpts.DividSpockGLM_TrlRng=[]; % what is the trial range we are runningn GLM now
AnalysisOpts.ShuffleStr=''; % string for shuffled data
AnalysisOpts.ExtraGLMStr=''; % extraStr added to anyfile for GLM
AnalysisOpts.ReWriteGLMData=1; % are we rewriting GLM data?
AnalysisOpts.GLMmdlCompStr=''; % what is the current repetition of model comparision for GLM
AnalysisOpts.GLMSpockFitTime=500; % fit time for GLM Spock per model 
AnalysisOpts.GLM_dependent_samples=0; % usually we have independent samples for shuffles
AnalysisOpts.GLM_p_threshold=0.05; % threshold 
AnalysisOpts.GLM_p_threshold_pop=0.001; % threshold for population analysis
AnalysisOpts.GLM_Time4SigStimInfo=0.6; % we look from 0 to this time to determine if a neuron has significant information about shape or color 
AnalysisOpts.GLM_UseR2forBestLambda=1; % are we using R2 to get the best lambda for each neuron(1) or we are using the max number of significant factors(0)
AnalysisOpts.GLM_two_sided=1;  % is this a two sided test
AnalysisOpts.GLM_num_clusters =inf; % max number of clusters
AnalysisOpts.GLM_npoints4indNeuSig=1; % number of points that we require for each neuron to have significant value 
AnalysisOpts.GLM_num_permutations=1000; % num of permutations for shuffle 
AnalysisOpts.GLM_SkipSingleCellCharctristics=0; % are we skipping generating single cell charactristics
AnalysisOpts.GLM_SkipPopulationCharctristics=0; % are we skipping generating population  charactristics
AnalysisOpts.GLM_pval_binomial_AreaCell=0.05; % what is the pval for number of significant neurons in an area
%% Population Analysis 
AnalysisOpts.RuleMainFeature={'ShapeCat','ColorCat','ColorCat'};
AnalysisOpts.RuleCorrectCatRepLoc={[1 1;2 2],[1 4;2 3],[1 2;2 1]};% what is the correct reponse for each category for each rule 
AnalysisOpts.RuleInCorrectCatRepLoc={[1 2;2 1],[1 3;2 4],[1 1;2 2]};% what is the incorrect reponse for each category for each rule 

AnalysisOpts.PopulationAna.ProcessingStepNames={'SaveImportantData','RunSubSpace','PlotSubSpace',... %0 1 2 
    'RunClassifer','PlotClassifier','ClassifierComparision','','CharactrizeSubSpace',... %3 4 5 6 7
    'ConcatinateClassifier','ConcatinateSubSpace',... %8, 9
    'GenSummeryPlotSingleCell','FitGLM','GLMMdlComparision','PlotGLMMdlComparison',... % 10, 11, 12 ,13 
    'CreatSummeryFile','ConcatinateGLMFiles','GenSummeryPlotArea'}; %Processing Step Names % 14 15 16

AnalysisOpts.ZscoreFactorData=0; % are we z-scoring factor data ?
AnalysisOpts.DetrendFactorData=0; % are we detrending factor data?
AnalysisOpts.RunCrossTemporalClassifer=0; % are we running a cross temporal classifier
AnalysisOpts.CalShuffleClassifier=0; % are we running shuffle test on classifier 
AnalysisOpts.CalShuffTrlOrderClassifier=0; % are we running trial order shuffle on the classifier
AnalysisOpts.CalShuffleClassifier_Cond=[]; % what is the condition we are running the shuffle classifier now
AnalysisOpts.CalShuffleClassifier_TrlRng=[]; % what is the trial range we are runningn shuffle classifier now
AnalysisOpts.DividSpockClassifier=3; % are we running  classifier on spock 1% to divide based on Condition and TrialRng 2: to divide baed on repetitions 3: divide based on rep / cond and trial range
AnalysisOpts.DividSpockClassifier_Cond=[]; % what is the condition we are running the  classifier now
AnalysisOpts.DividSpockClassifier_TrlRng=[]; % what is the trial range we are runningn  classifier now
AnalysisOpts.ExchangeableCalShuffClassifier=1; %<<KEEP it 1 always>> we are using the same set of trials for for calculating shuffle and observed in classifer 
AnalysisOpts.CurrentClassifierOpts=[]; % opetions for the current classifier 
AnalysisOpts.Classifier_FileDateTimeTh=datetime(2023,12,21);%is the input to datetime  be YYYY,MM,DD for mat or anything else( keep one week from now) 
AnalysisOpts.Classifier_Nrep=nan; % are we imposing number of reps for classifier analysis 
AnalysisOpts.Classifier_FileCond=''; % the condition of the currecnt file 
AnalysisOpts.PopulationAna.MaxMatchTrialConds=1; % are we maximally matching classifier trial types
if AnalysisOpts.PopulationAna.MaxMatchTrialConds==2;warning('Max Match condiiotn is 2 change to 1');end
AnalysisOpts.SweepClassifierConds=0; % are we looking into different combinations of conditions for classifer to find the best number of neurons?
AnalysisOpts.IncludeAllClassifierInfo=0; % are we including AUC bias and beta
AnalysisOpts.NRep2Use4StatTest=[]; % index of repetitions are we using to calculate the stat test(empty=All reps)
AnalysisOpts.ClassifierLearningTrlSweepDir=1; % -1: move back from switch trial. +1 move forward from switch trial
AnalysisOpts.GetOnlyShuffLabelsClassifier=0; % are we only generating the shuffle label for classifier and saing it?
AnalysisOpts.ChunkMaxArray=2500; % how many max arry in the chunk 
AnalysisOpts.ChunkArrayBiasNum=0; % how much should be add to the so that we have the correct array number 
AnalysisOpts.DeleteClassifierShuffleFiles=0; % are we deleting classifier shuffle files after we Classifier_TaskConcateSpockTimenate?
AnalysisOpts.SkipClassifierLearningPerfStatTest=0; % are we skipping plot the significance for performance of the classifier?
AnalysisOpts.PopulationAna.UseAvg_BeleifRespScoreCorr=1; % are we using trial by trail basis 0 or average blocks 1for corrleating response encoing and Belief
AnalysisOpts.PopulationAna.UseAvg_BeleifStimScoreCorr=1;% are we using trial by trail basis 0 or average blocks1 for corrleating response encoing and Belief
%% classifier tasks 
AnalysisOpts.PopulationAna.Classifier_TaskNameSet={'3D_Color_Cat_Xgen','3D_Shape_Cat_Xgen',... % Decoding Xdecoding conds
    '3D_Color_Cat_Xgen_BalRespDir','3D_Shape_Cat_Xgen_BalRespDir',...
    '3D_Response_Xgen_BalCorr','3D_Response_Xgen_BalCong','3D_Color_Response_Xgen',...
    'Learning3D_Color_Shape_Rule_Xgen_AltRule','Learning3D_Color_Shape_Rule_Xgen_AltRule_HiTh',... % learning conds
    'Learning3D_Color_Shape_Rule_Xgen_SameRule','Learning3D_Color_Color_Rule_Xgen_AltRuleR1',...
    'Learning3D_Color_Shape_Shape_Xgen','Learning3D_Color_Shape_Response_Xgen',...
    '3D_Color_Response_XgenTest','3D_Color_Response_XgenR2',...
    '3D_Color_Response_XgenR2R3','3D_Color_Response_XgenR1R2R3','3D_Color_Response_XgenR1R2R3BalCong',...
    '3D_Color_Response_XgenR1R2R3BalInCong','3D_Response_Xgen_BalInCong','3D_Color_Xgen_Bhv',...
    '3D_Color_Response_XgenR2R3Entropy','3D_Shared_Color_Response_EntropyOld',...
    '3D_Color_Response_Xgen_R2R1','3D_Color_Response_Xgen_R2R1_V2',...
    'Learning_3D_Color_Response_Xgen_R2R1_V2','Learning_3D_Color_Response_Xgen_R2R1_V2_Rev',...
    '3D_Color_Response_XgenR1R2R3BalInCongV2','3D_Color_Response_Xgen_R2R1_V3','3D_FineTuneColorCategoryClassifier',...
    'XTemp_AllObjects','3D_ColorCategoryFineTuneSig', ...
    '3D_Color_Cat_Xgen_BalRespDirNoTrR2','3D_Shape_Cat_Xgen_BalRespDirNoTrR2','3D_Shape_Response_Xgen',...
    '3D_Color_Response_XgenNMS','3D_Color_Response_XgenBalCong','ResponseLocation_BalCorr','3D_Color_Response_XgenBalInCong',...
    '3D_Color_Response_XgenBalInCongPrototype',...
    '3D_Shared_Color_Response_EntropyR3','3D_Shared_Color_Response_EntropyR2','3D_Shared_Color_Response_EntropyInCongR3','3D_Shared_Color_Response_EntropyInCongR2',...
    '3D_Color_Response_XgenBalInCongV2','3D_Color_Response_XgenBalInCongV2Test','3D_Color_Response_XgenBalInCongV2BalSh','Learning3D_Color_Shape_Rule_Xgen_AltRule_SameRule',...
    '3D_Color_Response_XgenBalInCongV3','Learning3D_Color_ShapeCorrInCorr_Rule_Xgen_AltRule',...
    'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB','Learning3D_Shape_Color_Rule_Xgen_SameRule_RB',...
    'Learning3D_Shape_Color_Rule_Xgen_AltRule_HiTh_RB',...
    'Learning3D_Shape_Color_Rule_Xgen_AltRuleR1_RB','Learning3D_Shape_Color_Rule_Xgen_SameRuleR1_RB',...
    'Learning3D_Color_Response_Rule_Xgen_AltRule','Learning3D_Color_Response_Rule_Xgen_SameRule',...
    'Learning3D_Shape_Color_Color_Compression_RB','3D_Shared_Color_Response_EntropyR3Bal','3D_Shared_Color_Response_EntropyR2Bal',...
    'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB_CutShuff','3D_Color_Response_XgenBalInCongV4',...
    '3D_Color_Response_XgenBalInCongV5','3D_Color_Response_XgenBalInCongV6','Learning3D_Shape_Color_Rule_Xgen_AltRule_RB_Step',...
    'Learning3D_Shape_Color_Response_Xgen_Orhto','3D_Shared_Color_Response_EntropyR3Bal_BalCongTest',...
    '3D_Color_Color_BalRespDir','3D_Response_Xgen_BalInCongRBC',...
    'Learning3D_Color_Response_Rule_Xgen_AltRuleRBC','Learning3D_Color_Response_Rule_Xgen_SameRuleRBC'};     % conditions for infomration evolution plots

if AnalysisOpts.DividSpockClassifier==1;SpockBaseTime=120;LearningCoef=1;  % if we are using conditions and trialrng
elseif AnalysisOpts.DividSpockClassifier==2;SpockBaseTime=30;LearningCoef=16; % if we are using repetitions
elseif AnalysisOpts.DividSpockClassifier==3;SpockBaseTime=30;LearningCoef=16;end % if we are using repetitions
% average run time for 3D classifier is 50secs with 5 CV folds and 3
% conditions and shuffling the time needed is 30 mins. Multiply this by the
% number of learning trial range to obtain learning time estimate
AnalysisOpts.PopulationAna.Classifier_TaskSpockBaseTime=SpockBaseTime;
BT=AnalysisOpts.PopulationAna.Classifier_TaskSpockBaseTime;
nClassifierTasks=length(AnalysisOpts.PopulationAna.Classifier_TaskNameSet);
AnalysisOpts.PopulationAna.Classifier_TaskSpockTime=[BT*ones(1,nClassifierTasks)]; % this is time for tasks that don't have learning
AnalysisOpts.PopulationAna.Classifier_TaskSpockTime(contains(AnalysisOpts.PopulationAna.Classifier_TaskNameSet,'Learning'))=LearningCoef*BT; % for learning tasks multiply by 16 which the average number of Trl range
AnalysisOpts.PopulationAna.Classifier_TaskConcateSpockTime=600*ones(1,nClassifierTasks);

AnalysisOpts.Xtemp3dPlotFunctions={'PerfTripleInd','PerfTripleSuperimposed','Scores','CompareTiming','ProjScores2D'};

AnalysisOpts.PopulationAna.Classifier_ComparisionNameSet={'3D_Cat_Color_Area_Summery','3D_Cat_Color_Response_Compare','3D_Cat_Color_XgenCol_Compare','3D_Response_XgenResponse_Compare',...
    'Cat_Color_Area_Summery','Cat_Color_Xgen_Area_Summery','Cat_Abstract_Color_Xgen_Area_Summery',...
    'Cat_Color_BalRespDir_Area_Summery','Cat_Color_Xgen_BalRespDir_Area_Summery',...
    'Cat_Shape_Area_Summery','Cat_Shape_Xgen_Area_Summery',...
    'Stim_Color_Area_Summery','Stim_Color_Xgen_Area_Summery',...
    'Stim_Shape_Area_Summery','Stim_Shape_Xgen_Area_Summery',...
    'AllObjs_Area_Summery',...
    'Resp_Xgen_Area_Summery',...
    'Col_Resp_Xgen_BalCong_Area_Summery','Shape_Resp_Xgen_BalCong_Area_Summery',...
    'Col_Resp_Xgen_BalCorr_Area_Summery','Shape_Resp_Xgen_BalCorr_Area_Summery',...
    'Learning_Cat_Color_Resp_Xgen_Area_Summery',...
    'FineTune_Cat_Color_Xgen_Area_Summery','FineTune_Stim_Color_Area_Summery','FineTune_Stim_Color_Kernel_Area_Summery',...
    'Xgen_Col_Resp_BalCong','Xgen_Col_Resp_BalCorr','Xgen_ColBalRespDir_Resp_BalCong','Xgen_Index_Color','Xgen_Index_Resp',...
    'Xgen_Col_Learning',...
    '2D_Cat_Color_Shape','2D_Cat_Color_Resp_Xgen','Xgen_Cat_Color_Area','Xgen_Resp_Area',...
    'Cat_Color_Area','ResponseLocation_Area'};
nClassifier_ComparisionTasks=length(AnalysisOpts.PopulationAna.Classifier_ComparisionNameSet);
AnalysisOpts.PopulationAna.ClassifierComparision_SpockTime=[60*ones(nClassifier_ComparisionTasks)];
%% Paper conditions to test for classifier analysis
AnalysisOpts.PaperFigs.Decoding.TestName={'3D_Color_Cat_Xgen','3D_Shape_Cat_Xgen',...
    '3D_Color_Cat_Xgen_BalRespDir','3D_Shape_Cat_Xgen_BalRespDir','3D_Color_Cat_Xgen_BalRespDirNoTrR2','3D_Shape_Cat_Xgen_BalRespDirNoTrR2',...
    '3D_Response_Xgen_BalCorr','3D_Response_Xgen_BalCong','3D_Response_Xgen_BalInCong',...
    '3D_Color_Response_Xgen','3D_Shape_Response_Xgen'};

AnalysisOpts.PaperFigs.CrossDecoding.TestName={};
AnalysisOpts.PaperFigs.FineTune.TestName={'3D_Color_Response_XgenR1R2R3BalInCong','3D_Color_Response_XgenR1R2R3BalInCongV2',...
    '3D_Color_Response_Xgen_R2R1','3D_Color_Response_Xgen_R2R1_V3'};%{'XTemp_AllObjects','FineTune_Cat_Color_Xgen','FineTune_Stim_Color','FineTune_Stim_Color_Kernel'};

AnalysisOpts.PaperFigs.Learning.TestName={'Learning3D_Color_Shape_Rule_Xgen_AltRule',...
    'Learning3D_Color_Shape_Rule_Xgen_AltRule_HiTh','Learning3D_Color_Shape_Rule_Xgen_SameRule',...
    'Learning3D_Color_Color_Rule_Xgen_AltRuleR1'};

AnalysisOpts.PaperFigs.Comparision.TestName={'Cat_Color_Area_Summery','Cat_Color_Xgen_Area_Summery',...
    'Cat_Color_BalRespDir_Area_Summery','Cat_Color_Xgen_BalRespDir_Area_Summery',...
    'Cat_Shape_Area_Summery','Cat_Shape_Xgen_Area_Summery',...
    'Resp_Xgen_Area_Summery',...
    'Col_Resp_Xgen_BalCong_Area_Summery','Shape_Resp_Xgen_BalCong_Area_Summery',...
    'Col_Resp_Xgen_BalCorr_Area_Summery','Shape_Resp_Xgen_BalCorr_Area_Summery'};
AnalysisOpts.RunMode=[]; % current run mode for processing 
AnalysisOpts.RunMode2=[]; % current run mode for processing more details
% RunMode2 defines which detail operation we are performing in each case
% ConcatinateClassifier Files 1:only main 2: only shuffle 3: both shuffle and main
AnalysisOpts.NTrlRngTrainLearningArea=nan; % number of trials to take from end of the block for training per area
AnalysisOpts.NTrlRngTestLearningArea=nan; % number of trials shift for Test
AnalysisOpts.ntrlPerCondArea=nan; % number of training trails per area
AnalysisOpts.LimitFromSwitchPerfArea=nan; % if we are limiting the perfromance
AnalysisOpts.NTrlStpLearningArea=nan; % steps for learning
AnalysisOpts.ClassifierLearningParamsFields={'ntrlPerCondArea','NTrlRngTrainLearningArea',...
    'NTrlRngTestLearningArea','LimitFromSwitchPerfArea',...
    'CalShuffTrlOrderClassifier','GetOnlyShuffLabelsClassifier','pvalClassifierAnalysis_ClusterCorrect','NTrlStpLearningArea'};
%% determine which model we are testing and what run we are doing
AnalysisOpts.RunLocator=cell2mat(arrayfun(@(x) [x*ones(1,AnalysisOpts.SingCellAna.nModels2Test) ; 1:AnalysisOpts.SingCellAna.nModels2Test],1:AnalysisOpts.SingCellAna.GLMnMdlCompRuns,'UniformOutput',0))';
AnalysisOpts.RunLocatorDims={'RunNum','ModelNum'}; % what are the columns of the the runlocator 
AnalysisOpts.factornames={'ColorML','ShapeML','ColorCat','ShapeCat','Rule','ResponseLoc','Reward','RT','Time',...
                'Hybrid_Q1','Hybrid_Q2','Hybrid_Q3','Hybrid_Q4','Hybrid_RPE','Hybrid_W_Color1','Hybrid_W_Color2',...
                'Hybrid_W_Shape1','Hybrid_W_Shape2','Hybrid_Baxes',...
                'InferAF_Q1','InferAF_Q2','InferAF_Q3','InferAF_Q4','InferAF_RPE','InferAF_Baxes','InferAF_Bfeature',...
                'InferRule_Q1','InferRule_Q2','InferRule_Q3','InferRule_Q4','InferRule_RPE','InferRule_Brule','Hybrid_W_Diff1','Hybrid_W_Diff2',...
                'Congruency','Feature','Axis','SeqHist','TrialNum','TrialNumReverse','BlkOrder','FromSwitchBhvPerf',...
                'BhvPerfShape50','BhvPerfColor50','BhvPerfShape10','BhvPerfColor10'};
AnalysisOpts.FactorInds2Keep={'ColorML','ShapeML','Reward','ResponseLoc','RuleInfo','BhvPerfShape50','BhvPerfColor50','BhvPerfShape10','BhvPerfColor10'}; % factors we are keeping for further analysis 

AnalysisOpts.NonMatchingRecs=[]; % recording that bhv and bhv model don't match in the number of trials 
AnalysisOpts.CheckAxis1Labels=@(x) sum(x==1)+sum(x==2); % check if these labels belong to Axis 1
AnalysisOpts.CheckAxis2Labels=@(x) sum(x==3)+sum(x==4); % check if these labels belong to Axis 2
AnalysisOpts.CheckBothAxisLabels=@(x) AnalysisOpts.CheckAxis1Labels(x) & AnalysisOpts.CheckAxis2Labels(x);
AnalysisOpts.PlotEntropySig=1; % are we superimposing significance for entropy analysis
%% Statistical tests and comparisions
AnalysisOpts.pvalEntropyAnalysis=[0.05 0.001]; % pval to detect correlations
AnalysisOpts.pvalClassifierAnalysis_ClusterCorrect=0.001; % cluster correction pvalue
AnalysisOpts.pvalClassifierAnalysisTrlShuff_ClusterCorrect=0.1; % cluster correction pvalue for trial shuff
AnalysisOpts.baseline_percentage_RiseTime=0.2; % what is the baseline percantage to calculate rise time
AnalysisOpts.peak_percentage_RiseTime=0.8;% what is the peak percantage to calculate rise time
AnalysisOpts.NPnts_SubtractBaseLine='auto'; % number of point to subtract base line if 'auto' then uses specified timiing AnalysisOpts.PaperSpec.NPnts_SubtractBaseLine
AnalysisOpts.pval_TimeTrendTest=0.05; %pvalue when we do trend test over time and use bonferroni
AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod='Modified_MannKendall'; % can be regression or Kendall(corretion) Modified_MannKendall
AnalysisOpts.Modified_MannKendall_significance_value_tau = 0.05; % significance level above which kendall tau values are statistically significant. One-tailed test.
AnalysisOpts.Modified_MannKendall_significance_value_ac = [0.05]; %predetermined level for selecting only those autocorrelation lag values that are statistically significant. Two-tailed test
AnalysisOpts.Classifier_TrlShuff_TrendCorrMethodCode=1; % which code we are using(KEEP IT 1)
AnalysisOpts.Classifier_TrlShuff_UseBonferroni=2; % are we using bonferroni correction for modidied MannKendall Trend test 1: bonferroni 2: is Benjamini-Hotchberg correction 3: holm bonferroni
AnalysisOpts.Classifier_TrlShuff_BenjaminiHochberg_FalseDiscoveryRate=0.05;
%% subspace Analysis
AnalysisOpts.PopulationAna.SubspaceAna_TaskNameSet=[{'Learning_Cat_Shape_Color_Xgen_AltRule','Learning_Cat_Color_Color_Xgen_AltRule','Sensory_Motor_Transformation'},cell(1,150)];
AnalysisOpts.PopulationAna.PSTHbin=100; % 100ms PSTHbin
AnalysisOpts.CalShuffleSubspace=0; % are we running shuffle test on subspace 
AnalysisOpts.CalShuffleSubspace_Cond=[]; % what is the condition we are running the shuffle subspace now
AnalysisOpts.CalShuffleSubspace_TrlRng=[]; % what is the trial range we are runningn shuffle subspace now
AnalysisOpts.DividSpockSubspace=1; % are we deviding processing on spock for subspace 
AnalysisOpts.DividSpockSubspace_Cond=[]; % what is the condition we are running  subspace now
AnalysisOpts.DividSpockSubspace_TrlRng=[]; % what is the trial range we are runningn  subspace now
AnalysisOpts.CalShuffleAxB=0; % are we calculating shuffle for AxB analysis 
AnalysisOpts.GetFullAxBdata=0; % are we getting the full data for AxB( to make Movies)
%% Trial data structure 
AnalysisOpts.MaxTrialsToAnalyz='ALL';
AnalysisOpts.UseWorkingPeriod=0;  
AnalysisOpts.UseManualSpkSorting=1;  % are using manual spike sorting to create data structures?
AnalysisOpts.UseMedianReferenced=0; % are we using Median Referenced signal 
AnalysisOpts.WorkingPeriodTxt=''; % text to include in file names
AnalysisOpts.ManualSpkSortingTxt=''; % text to include in file names
AnalysisOpts.MedianReferencedTxt=''; % test to include in file names 
AnalysisOpts.BlockSpecsFeilds={'FromSwitch','ToSwitch','AroundSwitch','AllTrials'};
AnalysisOpts.TrlSpkTimeFieldName='AllTrials';
AnalysisOpts.SpkCntStartFieldNames={'FIXATE_ACQUIRED','SAMPLE_ON','SACCADE_START'}; %DO NOT CHANGE THE ORDER
AnalysisOpts.SpkCntStartFieldName='SAMPLE_ON';%'SAMPLE_ON';
AnalysisOpts.StimulusMorphLevels=[0 30 50 70 100 130 150 170];
AnalysisOpts.StimulusMorphLevelsNo50=[0 30 70 100 130 170];
AnalysisOpts.GetIndStrFieldName=(@(x) find(strcmp(AnalysisOpts.SpkCntStartFieldNames,x))); % function to find which episode we are locking the results

%% colors 
HEX2RGB=@(x) transpose(cell2mat(cellfun(@(y) sscanf(y(2:end),'%2x%2x%2x',[1 3])'/255,x,'uniformoutput',0)));
% double checked images used in the experiment with the https://imagecolorpicker.com/en
AnalysisOpts.MorphlevelsColRGB= HEX2RGB({'#d7536b','#c67442','#84962b','#339c64','#709784','#bf7686'});   %[247 19 19;247 171 19;153 153 0;0 153 76;0 153 153;255 153 204]/255; % morph level colors without 50%
AnalysisOpts.MorphlevelsColRGBInc50=HEX2RGB({'#d7536b','#c67442','#aa8a13','#84962b','#339c64','#709784','#9e8b8b','#bf7686'});    %[247 19 19;247 171 19;245 209 66;153 153 0;0 153 76;0 153 153;245 209 66;255 153 204]/255; % includes 50%
AnalysisOpts.MorphlevelsShpRGBInc50=copper(8);%[0 0 0;0.3125 0.1953 0.1244; 0.6250 0.3906 0.2487;0.9375 0.5859,0.3731;1.0000 0.7812 0.4975;...
  %  0.9375 0.5859,0.3731;0.6250 0.3906 0.2487;0.3125 0.1953 0.1244]; % shape color map
AnalysisOpts.ProtoypeImgRGB=[245, 66, 66; 78, 245, 66;66, 96, 245; 239, 66, 245];


%% define my color paletts 
AnalysisOpts.ColorPalett1=HEX2RGB({'#264653', '#2a9d8f', '#8ab17d', '#e9c46a', '#f4a261', '#e76f51'}); % charcoal persian green olivine saffron sandy brown burnt siennahttps://coolors.co/264653-2a9d8f-8ab17d-e9c46a-f4a261-e76f51
AnalysisOpts.ColorPalett2=HEX2RGB({'#3d5a80', '#98c1d9', '#e0fbfc', '#e7b4a5', '#ee6c4d', '#293241'}); %yimin blue, powder blue, light cyan melon burnt sienna Gunmetal URL https://coolors.co/3d5a80-98c1d9-e0fbfc-e7b4a5-ee6c4d-293241
AnalysisOpts.ColorPalett3=HEX2RGB({'#f72585', '#b5179e', '#7209b7', '#3a0ca3', '#4361ee', '#4cc9f0'}); %same as AnalysisOpts.ColorPalett5 URL https://coolors.co/f72585-b5179e-7209b7-3a0ca3-4361ee-4cc9f0
AnalysisOpts.ColorPalett4=HEX2RGB({'#ffcdb2', '#ffb4a2', '#e5989b', '#b5838d', '#6d6875'});
AnalysisOpts.ColorPalett5=HEX2RGB({'#f72585', '#7209b7', '#3a0ca3', '#4361ee', '#4cc9f0'});% Rose, Grape, Zaffre, Noon blue, sky blue https://coolors.co/f72585-7209b7-3a0ca3-4361ee-4cc9f0
AnalysisOpts.ColorPalett6=HEX2RGB({'#ffbe0b', '#fb5607', '#ff006e', '#8338ec', '#3a86ff'});% Amber, Orange, Rose, Blue, Aure https://coolors.co/ffbe0b-fb5607-ff006e-8338ec-3a86ff
AnalysisOpts.ColorPalett7=HEX2RGB({'#03045e', '#0077b6', '#00b4d8', '#90e0ef', '#caf0f8'});% Extention of blue colors https://coolors.co/03045e-0077b6-00b4d8-90e0ef-caf0f8
AnalysisOpts.ColorPalettCamden= [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];[0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1 ]; [1 .3 .7]; [.2 .8 .2];[.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]];
AnalysisOpts.CurrColorPalett=AnalysisOpts.ColorPalett1;
AnalysisOpts.ColorPalett1COSYNE23=HEX2RGB({'#ff006e', '#ff006e', '#ff006e', '#ff006e', '#ff006e', '#ff006e'}); % https://coolors.co/264653-2a9d8f-8ab17d-e9c46a-f4a261-e76f51
AnalysisOpts.ColorPalett2COSYNE23=HEX2RGB({'#293241', '#293241', '#293241', '#293241', '#293241', '#293241'}); %URL https://coolors.co/3d5a80-98c1d9-e0fbfc-e7b4a5-ee6c4d-293241
AnalysisOpts.ColorPalett3COSYNE23=HEX2RGB({'#5D576B', '#5D576B', '#5D576B', '#5D576B', '#5D576B', '#5D576B'}); %URL https://coolors.co/f72585-b5179e-7209b7-3a0ca3-4361ee-4cc9f0
%% define factor colro Sets 
AnalysisOpts.ColorCatColors_ColorName='YlGnBu7'; % color names that we use with othercolor function
AnalysisOpts.ShapeCatColors_ColorName='YlGn6';% color names that we use with othercolor function
AnalysisOpts.ResponseLocColors_ColorName='YlOrRd6';% color names that we use with othercolor function
AnalysisOpts.RuleColors=HEX2RGB({'#f99b38','#0d71b3','#558235'});%HEX2RGB({'#EB6B3A','#7865A9','#E64279'});%[3 15 252;252 157 3;3 252 236]/255; % color of the rules
AnalysisOpts.ColorMLColors=AnalysisOpts.MorphlevelsColRGBInc50; % color of the colorMLs
AnalysisOpts.ShapeMLColors=AnalysisOpts.MorphlevelsShpRGBInc50; % color of the shapeMls
AnalysisOpts.ColorCatColors=AnalysisOpts.RuleColors;%FigParams.getOthercolormap(AnalysisOpts.ColorCatColors_ColorName,5,1);%Purple AnalysisOpts.ColorMLColors([3 1 6],:);
AnalysisOpts.ShapeCatColors=AnalysisOpts.RuleColors;%FigParams.getOthercolormap(AnalysisOpts.ShapeCatColors_ColorName,5,1);%Green  AnalysisOpts.ShapeMLColors([3 1 6],:);
AnalysisOpts.ResponseLocColors=AnalysisOpts.RuleColors;%FigParams.getOthercolormap(AnalysisOpts.ResponseLocColors_ColorName,5,1);%AnalysisOpts.ColorPalett1([1 2 5 6],:);
AnalysisOpts.RewardColors=AnalysisOpts.RuleColors;%AnalysisOpts.ColorPalett3([1 6],:);
AnalysisOpts.RTColors=AnalysisOpts.ColorPalett3;
AnalysisOpts.TimeColors=AnalysisOpts.ColorPalett3;
AnalysisOpts.AreaColors=HEX2RGB({'#100404','#dca67f','#9e7b42','#78472b','#ffd770'});%AnalysisOpts.ColorPalettCamden;
AnalysisOpts.ClassifierAccuracyColormap=parula(256);%FigParams.getOthercolormap('RdYlBu5',256,1);%caxisLimitsPRGn10
AnalysisOpts.DimTxt=[0 2 3];% text for classifier dimensions
AnalysisOpts.QuadrantColors=HEX2RGB({'#E22121','#218717','#C1F7DC','#F26419'}); %[red bunny, green tee, green bunny, red tee]
% define factor line markers 
AnalysisOpts.FactorMarkers={'none','d','none','x','s'};%'Color','Shape','Response','Reward','Rule'
AnalysisOpts.FactorLineStyle={'-','--',':','--','-'};
%% define specifics for the paper 
AnalysisOpts.PaperSpec.LimitTimeAxis=1; % are we limiting time axis for paper plots
AnalysisOpts.PaperSpec.StrTime_SAMPLE_ON=-0.2; % start time before sample on for plotting
AnalysisOpts.PaperSpec.EndTime_SAMPLE_ON=0.61; % end time before sample on for plotting
AnalysisOpts.PaperSpec.StrTime_SACCADE_START=-0.4; % start time before saccade on for plotting
AnalysisOpts.PaperSpec.EndTime_SACCADE_START=0.41; % end time before saccade on for plotting
AnalysisOpts.PaperSpec.NPnts_SubtractBaseLine=[-0.2 0]; % Time points to subtract the base line
AnalysisOpts.PaperSpec.TimeAxisSteps=0.2; % what is the steps for time axis
AnalysisOpts.ThisTimeAxisStart=[]; %optional start time
AnalysisOpts.ThisTimeAxisEnd=[]; % optional end time
AnalysisOpts.CurrentAxisLimits=@(AnalysisOpts) [AnalysisOpts.PaperSpec.(['StrTime_' AnalysisOpts.SpkCntStartFieldName]),...
   AnalysisOpts.PaperSpec.(['EndTime_' AnalysisOpts.SpkCntStartFieldName]) ];

%% Charactristic Plot Neuron 
AnalysisOpts.ClusterNeuronFeature='wav'; % what feature to plot in the charactristic plot for the neuron
AnalysisOpts.max_spikes_plot=1000; % maximum number of spikes to plot to look at 

%% behavioral analysis (pure behavior 
AnalysisOpts.BhvAna.Chico_InitialDate=[2021 06 12];
AnalysisOpts.BhvAna.Chico_EndDate=[2021 08 20];
AnalysisOpts.BhvAna.Silas_InitialDate=[2018 06 07];
AnalysisOpts.BhvAna.Silas_EndDate=[2018 07 22];
AnalysisOpts.BhvAna.IgnoreLastBlk=0; % are we ignoring the last block
AnalysisOpts.BhvAna.NTrlBlkPerf=15;
%% preset some of the valiables 
%% manual spike sorting 
if AnalysisOpts.UseWorkingPeriod
AnalysisOpts.WorkingPeriodTxt='WorkPeriod';
end
if AnalysisOpts.UseManualSpkSorting
    AnalysisOpts.ManualSpkSortingTxt='_ManuSpkSrt';
end
if AnalysisOpts.UseMedianReferenced
    AnalysisOpts.MedianReferencedTxt='_Median';
end

%% replace values given in varargin
ParseParams(varargin) % add all of the additional parameters we have added to the function


