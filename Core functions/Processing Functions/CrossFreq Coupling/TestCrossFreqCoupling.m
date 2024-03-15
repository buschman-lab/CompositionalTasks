function TestCrossFreqCoupling(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
% Discovers frequency temporal motifs in electrophysiology data; then look
% at the cross frequency couple charactristics of compound motifs
global AnalysisOpts AnalysisData

%% define primary vars before kick off
AnalysisOpts.RunonCluster=RunonCluster;
AnalysisOpts.Project='Learning attentional templates';  % can be 'Learning attentional templates' or 'Rule Representaion'
AnalysisOpts.Animal='Scooter';
AnalysisOpts.AnalysisType='Data Preprocessing'; % can be 'Oscillations Breakdown ToolBox';
AnalysisOpts.AnalysisFocus1  ='PlotMotifs';  % 'ExplorParams','PlotMotifs'
AnalysisOpts.SubAnalysisType1='Wavelet_Trial';  % define type of sub analsysis  can be also 'LFP_Trial'
AnalysisOpts.AnalysisPathName='Motif_CFreqCouple';  % this is the name of the folder used in input-output data

%% kick off the function
FS=KickoffMyfunc(RunonCluster,DateNum);

%% define analysis options in a gobal varibales
AnalysisOpts.PairNum=PairNum; % what pair we are looking at
AnalysisOpts.Ch_2look=ChNum;% how many channels do you want to look at, leave empty if non
AnalysisOpts.AreaNum=AreaNum;
ParseParams(varargin) % add all of the additional parameters we have added to the function

%% define classes necessary for this analysis
TlBoxHelp=ToolboxHelpers;
[FigPrms,ManData,TrialFunc,MotifFunc,TimeFreq]=TlBoxHelp.LoadAnalysisClasses; % load all of the classes we want
[Fs,FsLFP,FsWave,FsWaveTarg]=TlBoxHelp.SetupGeneralFuncs; % set up all general functions we want
% get block and trial information we need
[TrialTimes,RuleBlockTrials,ChannelInfo,ChannelArea,ChsSet]=TrialFunc.InitializeTrialFuncs; % run common trial functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Starting the analysis of this function %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update some parameters
AnalysisOpts.MotifAnalysis.L=AnalysisOpts.MotifAnalysis.L_ms.*FsWaveTarg;

%% load parameters that match what we have here based on parameter sweep
[~,~,ParamSweepFileName]=GenerateFileName(FS,AnalysisOpts.DataSavePath,'Parameter_Sweep',[],[],[],[],...
    'SelfName',1,'SelfNameTxt',['ParamSweepIntersect_Ttot' num2str(AnalysisOpts.MotifAnalysis.MaxTime) '_KMax' num2str(AnalysisOpts.MotifAnalysis.K)]);
ParamSweepData=load(ParamSweepFileName,'Params','LamInterSect');
ThisLind=ParamSweepData.Params.L_ms==AnalysisOpts.MotifAnalysis.L_ms; %find ind of desired L
AnalysisOpts.MotifAnalysis.lambda=ParamSweepData.LamInterSect(ThisLind); % replace lambda
%

% get the Xcorr file name from ALL of the recordings
[~,~,AnalysisData.XCorrFileName]=ManData.GetFileName('Motif_Cluster',AnalysisOpts.ExtraParamTxtWrite,'WantedDate','ALL');
% load Core Motifs, plot them and save it off
[AnalysisData.CoreMotifs,AnalysisData.SizeW]=MotifFunc.ClassifyMotifs_Area(AnalysisData.XCorrFileName);

%% Loop through channels and find motifs
extraStr=['_L' num2str(AnalysisOpts.MotifAnalysis.L_ms) 'T' num2str(AnalysisOpts.MotifAnalysis.MaxTime)];
for Ch=ChsSet
    AnalysisOpts.CurrentCh=Ch;
    if ManData.IsVarExistinFileAnalysis(AnalysisOpts.AnalysisPathName,'AnalysisData',extraStr) & ~AnalysisOpts.ReWriteData
        fprintf('\nChannel %i motifs already exist; skippking this channel ...',Ch)
    else
        
        fprintf('\n processing Motifs channel %s \n',num2str(Ch));
    %    [AnalysisData.RawData]=TlBoxHelp.getRawEphysData(Ch,AnalysisOpts.UseDataPointer); AnalysisData.RawData=AnalysisData.RawData{1};
    %    [AnalysisData.cwt,AnalysisData.cwt_f]=TlBoxHelp.getWaveletData(Ch,AnalysisOpts.UseDataPointer);AnalysisData.cwt=AnalysisData.cwt{1};
        if 0
            %% Get the CWT data
            AnalysisData.cwt=transpose(downsample(AnalysisData.cwt',FsWave/FsWaveTarg));
            AnalysisData.Time_org=0:1/FsWaveTarg:(size(AnalysisData.cwt,2)-1)/FsWaveTarg;
            AnalysisData.TimeRaw_org=0:1/FsLFP:(length(AnalysisData.RawData)-1)/FsLFP;
            AnalysisData.RawData_org=AnalysisData.RawData;
            
            flimitind=AnalysisData.cwt_f>=AnalysisOpts.MotifAnalysis.flimit(1) & AnalysisData.cwt_f<=AnalysisOpts.MotifAnalysis.flimit(2);
            AnalysisData.cwt_f=AnalysisData.cwt_f(flimitind);
            AnalysisData.cwt_org=AnalysisData.cwt(flimitind,:);
            
            
            % load refit data for this channel
            AnalysisOpts.ExtraParamTxtRead=['_L' num2str(AnalysisOpts.MotifAnalysis.L_ms) 'T' num2str(AnalysisOpts.MotifAnalysis.MaxTime) ];
            AnalysisOpts.ExtraParamTxtRefit=[AnalysisOpts.ExtraParamTxtRead '_MotifDistArea'  num2str(AreaNum) '_rf' num2str(Ch) ];  % MC indicates it is a multi-channel
            AnalysisOpts.CurrentCh='';
            [~,~,RefitFileName]=ManData.GetFileName('Motif_Cluster',AnalysisOpts.ExtraParamTxtRefit);
            load(RefitFileName);
            AnalysisData.MotifSpecs_Refit{1}=MotifSpecs_Refit;
            
            % get H and W, shift Ws to zero and shift Hs
            %    H=arrayfun(@(x) AnalysisData.MotifSpecs_Refit{x}.H,1:length(AnalysisData.Ch),'UniformOutput',0); % retrive H and W
            W=AnalysisData.MotifSpecs_Refit{1}.W;
            H=AnalysisData.MotifSpecs_Refit{1}.H;
            % [W,H,AnalysisData.Lags]=MotifFunc.ShiftWH2Zero(W,H,AnalysisOpts.MotifAnalysis.ShiftWH2ZeroTH); % shift H values so that everyting is centered on start of Motifs
            
            TlBoxHelp.getSpikeData(ChsSet,ChannelArea);  % get spike data we want
            XfreqCoup=cell(1,4);
            [XfreqCoup{:}]=MotifFunc.CrosFreqAnalysis(AnalysisData.cwt_org,AnalysisData.RawData_org,W,...
                H,AnalysisData.MotifSpecs_Refit{1}.loadings,AnalysisData.cwt_f,AnalysisData.Time_org,...
                AnalysisData.SpkData,0.4,'LogScale',0);
            
            [FigSaveFileName,Path]=GenerateFileName(FS,AnalysisOpts.ResultsSavePath,AnalysisOpts.AnalysisPathName,AnalysisOpts.Animal,AnalysisOpts.RecDate,Ch,'ExampleXfeqCoupFit');
            FigPrms.SaveFigSeries(FigSaveFileName,Path,[XfreqCoup  ])
            
        end
        % try loading motifs and reftting them                
        if 1
            MaxTime=3000; % fit to 60 seconds of data
            load('RecDateInfo.mat','ChsSet')
            ChSet=ChsSet{9}(end);
            load Exmp2CoreMotifsPhenoClust
            [MotifSpecs_Refit ]=MotifFunc.RefitMotifsGeneral(Ch,9,CoreMotifs,MaxTime,AnalysisData.SizeW);
            W= MotifSpecs_Refit.W;
            H= MotifSpecs_Refit.H;
            % [W,H,AnalysisData.Lags]=MotifFunc.ShiftWH2Zero(W,H,AnalysisOpts.MotifAnalysis.ShiftWH2ZeroTH); % shift H values so that everyting is centered on start of Motifs
             
         %   save Exmp2CoreMotifsPhenoClust MotifSpecs_Refit -append
             XfreqCoup1=cell(1,3);
            [XfreqCoup1{:}]=MotifFunc.CrosFreqAnalysis(AnalysisData.cwt,AnalysisData.RawData,W,...
                H, MotifSpecs_Refit.loadings,AnalysisData.cwt_f,AnalysisData.Time_cwt,...
                [7],0.998,'LogScale',0);
            
             XfreqCoup2=cell(1,3);
            [XfreqCoup2{:}]=MotifFunc.CrosFreqAnalysis(AnalysisData.cwt,AnalysisData.RawData,W,...
                H, MotifSpecs_Refit.loadings,AnalysisData.cwt_f,AnalysisData.Time_cwt,...
                [8],0.99,'LogScale',0);
            
            
        end
        if 0
        
        %% Get the CWT data
        AnalysisData.cwt=transpose(downsample(AnalysisData.cwt',FsWave/FsWaveTarg));
        AnalysisData.Time_org=0:1/FsWaveTarg:(size(AnalysisData.cwt,2)-1)/FsWaveTarg;
        AnalysisData.TimeRaw_org=0:1/FsLFP:(length(AnalysisData.RawData)-1)/FsLFP;
        AnalysisData.RawData_org=AnalysisData.RawData;
        
        %  AnalysisData.cwt=movmean(AnalysisData.cwt,AnalysisOpts.MotifAnalysis.Navg,2); % smooth the cwt
        flimitind=AnalysisData.cwt_f>=AnalysisOpts.MotifAnalysis.flimit(1) & AnalysisData.cwt_f<=AnalysisOpts.MotifAnalysis.flimit(2);
        AnalysisData.cwt_f=AnalysisData.cwt_f(flimitind);
        AnalysisData.cwt_org=AnalysisData.cwt(flimitind,:);
        
        if ~AnalysisOpts.MotifAnalysis.TestGeneralize % if we are testing the generalization of data as well
            Timetot=AnalysisData.Time_org(end)-1/FsWaveTarg; % how many secs of data
            TimeSteps=0:AnalysisOpts.MotifAnalysis.MaxTime:Timetot; %in sec
            StrTim=TimeSteps(1:end-1)+1/FsWaveTarg;
            StpTim=TimeSteps(2:end)+1/FsWaveTarg; % in sec
            StrTimDisc=StrTim(1:end);StpTimDisc=StpTim(1:end); % discovery times
            % log in our timings
            AnalysisData.StrTim=StrTim;        AnalysisData.StpTim=StpTim;
            AnalysisData.StrTimDisc=StrTimDisc;AnalysisData.StpTimDisc=StpTimDisc;
            
        else % split the data into discovery and test epochs
            
            Timetot=AnalysisData.Time_org(end)-1/FsWaveTarg; % how many secs of data
            TimeSteps=0:AnalysisOpts.MotifAnalysis.MaxTime:(Timetot); %in sec
            StrTim=TimeSteps(1:end-1)+1/FsWaveTarg;
            StpTim=TimeSteps(2:end)+1/FsWaveTarg; % in se
            StrTimDisc=StrTim(1:2:end);StpTimDisc=StpTim(1:2:end); % discovery times
            StrTimTest=StrTim(2:2:end);StpTimTest=StpTim(2:2:end); % test Times
            % log in our timings
            AnalysisData.StrTim=StrTim;        AnalysisData.StpTim=StpTim;
            AnalysisData.StrTimDisc=StrTimDisc;AnalysisData.StpTimDisc=StpTimDisc;
            AnalysisData.StrTimTest=StrTimTest;AnalysisData.StpTimTest=StpTimTest;
        end
        
        for step=1:length(StrTimDisc)
            %% prepare the data for this chunk
            % load discover data
            AnalysisData.cwt_Disc    =AnalysisData.cwt_org(:,(StrTimDisc(step)*FsWaveTarg):(StpTimDisc(step)*FsWaveTarg-1));
            AnalysisData.RawData_Disc=AnalysisData.RawData_org((StrTimDisc(step)*FsLFP):(StpTimDisc(step)*FsLFP-1));
            AnalysisData.Time_Disc   =AnalysisData.Time_org((StrTimDisc(step)*FsWaveTarg):(StpTimDisc(step)*FsWaveTarg-1));
            
            
            if AnalysisOpts.MotifAnalysis.TestGeneralize
                % load test data
                AnalysisData.cwt_Test    =AnalysisData.cwt_org(:,(StrTimTest(step)*FsWaveTarg):(StpTimTest(step)*FsWaveTarg-1));
                AnalysisData.RawData_Test=AnalysisData.RawData_org((StrTimTest(step)*FsLFP):(StpTimTest(step)*FsLFP-1));
                AnalysisData.Time_Test   =AnalysisData.Time_org((StrTimTest(step)*FsWaveTarg):(StpTimTest(step)*FsWaveTarg-1));
            end
            f_lin=[ManData.BinData(0.25,[],0.5,1) ManData.BinData(1,[],1,20) ManData.BinData(2,[],30,20)];
            f1=0.5;f2=80;nf=49;
             f_lin=f1:(f2-f1)/nf:f2;

            f_log=2.^(log2(f1):(log2(f2)-log2(f1))/nf:log2(f2));
            WAV=ApplyWaveletTransform(AnalysisData.RawData_Disc',FsLFP,'Freqs',f_lin,...
                                'WaveletWidth',7);
            cwt_s=WAV{1}.amp;cwt_f=WAV{1}.f;
            AnalysisData.cwt_Disc=cwt_s;AnalysisData.cwt_f=cwt_f;
            AnalysisData.cwt_Disc=transpose(downsample(AnalysisData.cwt_Disc',FsLFP/FsWaveTarg));
            %% Step 2: Discover and Plot motifs
            fprintf('\n processing Motifs rec %s channel %i step %i ',AnalysisOpts.RecDate,Ch,step);
            MotifFunc.ChNum=Ch;
            MotifFigs=cell(1,AnalysisOpts.MotifAnalysis.K+1);
            %% Discover motifs
            fprintf('\n discovering motifs')
            [AnalysisData.MotifSpecs_Disc{step},MotifFigs{:}]=MotifFunc.PlotMotifs(AnalysisData.cwt_Disc,AnalysisData.cwt_f,[],...
                'SavePlot',AnalysisOpts.SavePlot,'SaveData',AnalysisOpts.SaveData,'ShowPlot',step==1,...
                'DownSampleFactor',1,'ExtraString','Sim_','LogScale',0,...
                'PowerMethod' ,AnalysisOpts.MotifAnalysis.PowerMethod,'maxiter',100,'K',AnalysisOpts.MotifAnalysis.K,'W_fixed',0,...
                'L',AnalysisOpts.MotifAnalysis.L,'lambda',AnalysisOpts.MotifAnalysis.lambda,'lambdaOrthoH',AnalysisOpts.MotifAnalysis.lambdaOrthoH,...
                'lambdaOrthoW',AnalysisOpts.MotifAnalysis.lambdaOrthoW,'lambdaL1W',AnalysisOpts.MotifAnalysis.lambdaL1W,'lambdaL1H',AnalysisOpts.MotifAnalysis.lambdaL1H);
            
            %% test significance with test period
            %       fprintf('\nTesting significance of factors on held-out data')
            %       [AnalysisData.MotifSpecs_Disc{ChInd,step}.pvals,AnalysisData.MotifSpecs_Disc{ChInd,step}.is_significant] = ...
            %       test_significance(AnalysisData.cwt_Test,AnalysisData.MotifSpecs_Disc{ChInd,step}.W,0.05);
            %% remove motifs that don't have any loadings
            AnalysisData.MotifSpecs_Disc{step}=MotifFunc.RemoveNoloadingMotifs(AnalysisData.MotifSpecs_Disc{step});
            
            %% test generalization of the motifs
            if AnalysisOpts.MotifAnalysis.TestGeneralize
                fprintf('\nRefitting to held-out data')
                [AnalysisData.MotifSpecs_Test{step}]=MotifFunc.PlotMotifs(AnalysisData.cwt_Test,AnalysisData.cwt_f,[],...
                    'SavePlot',AnalysisOpts.SavePlot,'SaveData',AnalysisOpts.SaveData,'ShowPlot',step==1,...
                    'DownSampleFactor',1,'ExtraString','Sim_','LogScale',0,'PowerMethod' ,AnalysisOpts.PowerMethod,...
                    'maxiter',AnalysisOpts.MotifAnalysis.maxiter,'K',AnalysisOpts.MotifAnalysis.K,...
                    'L',AnalysisOpts.MotifAnalysis.L,'lambda',0,'lambdaOrthoH',AnalysisOpts.MotifAnalysis.lambdaOrthoH,...
                    'lambdaOrthoW',AnalysisOpts.MotifAnalysis.lambdaOrthoW,'lambdaL1W',AnalysisOpts.MotifAnalysis.lambdaL1W,'lambdaL1H',0,...
                    'W_init',AnalysisData.MotifSpecs_Disc{step}.W,'W_fixed',1,'useWupdate',0);
                %% remove motifs that don't have any loadings
                AnalysisData.MotifSpecs_Test{step}=MotifFunc.RemoveNoloadingMotifs(AnalysisData.MotifSpecs_Test{step});
            end
        end
        %%  Show some example Raw data so we are sure what we are looking at
        if AnalysisOpts.ShowPlot
            RawFigs=cell(1,2);
            AnalysisData.TrialTimes=[];
            [RawFigs{:}]=MotifFunc.PlotRawData(AnalysisData.cwt_Disc,AnalysisData.cwt_f,AnalysisData.Time_Disc,...
                AnalysisData.RawData_Disc,AnalysisData.TrialTimes,'PowerMethod',AnalysisOpts.MotifAnalysis.PowerMethod,'Fs',FsLFP,...
                'LogScale',0,'ChNum',Ch,'DownSampleFactor',1);
            %MotifFunc.ShowRawDataTF(AnalysisData.RawData_Disc,AnalysisData.cwt_Disc,AnalysisData.cwt_f,AnalysisData.TrialTimes,Fs,FsWaveTarg,FsWaveTarg,1)
        end
        %% Plot Collective PEV
        if AnalysisOpts.ShowPlot
            PEVfig=cell(1);
            FigPrms.RenderFigure(1,[])
            PEVfig=MotifFunc.PlotMotifPEV({AnalysisData.cwt_Disc},AnalysisData.MotifSpecs_Disc(end),'PEVmethod','var','freq',AnalysisData.cwt_f,'PowerMethod' ,AnalysisOpts.MotifAnalysis.PowerMethod,'DownSampleFactor',1);
        end
        end
        %% Plot Example Hs and H triggered analysis
        %        if AnalysisOpts.ShowPlot
        ExampleHFig=cell(1,5);
        HTriggAnaFins=cell(1,4);
        TlBoxHelp.getSpikeData(ChsSet,ChannelArea);  % get spike data we want
        [HTriggAnaFins{:}]=MotifFunc.CrosFreqAnalysis(AnalysisData.cwt_Disc,AnalysisData.RawData_Disc,AnalysisData.MotifSpecs_Disc{step}.W,...
            AnalysisData.MotifSpecs_Disc{step}.H,AnalysisData.MotifSpecs_Disc{step}.loadings,AnalysisData.cwt_f,AnalysisData.Time_Disc,...
            4,0.8,'LogScale',0);
        
        %           ExampleHFig=arrayfun(@(x) MotifFunc.PlotExampleHs(AnalysisData.cwt_Disc,AnalysisData.RawData_Test,AnalysisData.MotifSpecs_Disc{step}.W,...
        %           AnalysisData.MotifSpecs_Test{step}.H,AnalysisData.cwt_f,AnalysisData.Time_Test,'LogScale',0,'PowerMethod',AnalysisOpts.PowerMethod),1:5);
        %        end
        %% Save data and plots
        if AnalysisOpts.SavePlot
            [FigSaveFileName,Path]=GenerateFileName(FS,AnalysisOpts.ResultsSavePath,AnalysisOpts.AnalysisPathName,AnalysisOpts.Animal,AnalysisOpts.RecDate,Ch,extraStr);
            FigPrms.SaveFigSeries(FigSaveFileName,Path,[RawFigs MotifFigs {PEVfig} ])
            FigPrms.CloseFigs([RawFigs MotifFigs {PEVfig}])
        end
        %%
        
        if AnalysisOpts.SaveData
            % remove the fields we don't need so we save space
            AnalysisData=ManData.rmfieldExept(AnalysisData,{'MotifSpecs_Disc','MotifSpecs_Test','cwt_f','StrTim','StpTim','StrTimDisc','StpTimDisc','StrTimTest','StpTimTest'});
            ManData.SaveVar(AnalysisOpts.AnalysisPathName,AnalysisOpts,'AnalysisOpts',extraStr)
            ManData.SaveVar(AnalysisOpts.AnalysisPathName,AnalysisData,'AnalysisData',extraStr)
        end
    end
end


end

%% put the functions to kick off here
function FS=KickoffMyfunc(RunonCluster,DateNum)
% sets up all of the initial options and path
global AnalysisOpts

%% ok now first setup everything for cluster first because it needs path
AnalysisOpts.RunOnCluster=RunonCluster;

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
SetupAllVars(DateNum)  %% set up the path and initialize vars

end

