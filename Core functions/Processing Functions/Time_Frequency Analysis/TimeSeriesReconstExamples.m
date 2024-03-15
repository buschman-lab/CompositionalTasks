function TimeSeriesReconstExamples(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
% loads example data, runs motifs discovery and then tries to reconstrcut
% the orginal time series of data
global AnalysisOpts AnalysisData

%% define primary vars before kick off
AnalysisOpts.RunonCluster=RunonCluster;
AnalysisOpts.Project='Learning attentional templates';  % can be 'Learning attentional templates' or 'Rule Representaion'
AnalysisOpts.Animal='Scooter';
AnalysisOpts.AnalysisType='Data Preprocessing'; % can be 'Oscillations Breakdown ToolBox';
AnalysisOpts.AnalysisFocus1  ='';  % 'ExplorParams','PlotMotifs'
AnalysisOpts.SubAnalysisType1='';  % define type of sub analsysis  can be also 'LFP_Trial'
AnalysisOpts.AnalysisPathName='';  % this is the name of the folder used in input-output data

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
[~,~,ParamSweepFileName]=GenerateFileName(FS,AnalysisOpts.DataSavePath,'Parameter_Sweep','ALL','ALL',[],[],...
'SelfName',1,'SelfNameTxt','Parameter_Sweep_ALL_ALL__Params.mat');

if exist(ParamSweepFileName,'file')
    ParamSweepData=load(ParamSweepFileName,'Params','LamInterSect');
    ThisLind=ParamSweepData.Params.L_ms==AnalysisOpts.MotifAnalysis.L_ms; %find ind of desired L
    ThisMaxTim=ParamSweepData.Params.MaxTimeSet==AnalysisOpts.MotifAnalysis.MaxTime;
    AnalysisOpts.MotifAnalysis.lambda=ParamSweepData.LamInterSect(ThisMaxTim,ThisLind); % replace lambda
end

load('ExampleData.mat','data');
opts.WaveletMethod='matlab';
opts.q1=3;opts.q2=14;
opts.FrqStp=[1];
opts.frstp=1;
NRaw=length(opts.FrqStp)+1;Ncol=2;

N=AnalysisOpts.MotifAnalysis.MaxTime*FsLFP; % number of samples to look at
f1=1;f2=100;
%FigParams.RenderFigure(1,[]);
f_lin=sort(unique([0.25 0.5 0.75 ManData.BinData(opts.frstp,[],f1,f2)])); % define frequency axis
q=opts.q1:(opts.q2-opts.q1)/(length(f_lin)-1):opts.q2;
% cal cwt
[AnalysisData.cwt_data,AnalysisData.cwt_f,~] = TimeFreq.ComputeTimeFreq(data,'Fs',FsLFP,...
    'WaveletMethod',opts.WaveletMethod,'FreqToPlot',f_lin,'WaveletWidth',q,'VoicesPerOctave',8);

%[AnalysisData.cwt_f,AnalysisData.cwt_f_ind]=sort(AnalysisData.cwt_f);
%AnalysisData.cwt_data=AnalysisData.cwt_data(AnalysisData.cwt_f_ind,1:N);
AnalysisData.cwt_data=AnalysisData.cwt_data(:,1:N);
% cal power of data 
%AnalysisData.cwtPower=ManData.NormPower(AnalysisData.cwt_data(:,1:N),AnalysisData.cwt_f);
%      CWTPower=smoothdata(CWTPower,2);

AnalysisData.cwt_data=transpose(downsample(AnalysisData.cwt_data',FsLFP/FsWaveTarg)); % downsample to 50hz

% discover motifs
[AnalysisData.MotifSpecs,MotifFigs]=MotifFunc.PlotMotifs(AnalysisData.cwt_data,AnalysisData.cwt_f,[],...
    'SavePlot',AnalysisOpts.SavePlot,'SaveData',AnalysisOpts.SaveData,'ShowPlot',1,...
    'DownSampleFactor',1,'ExtraString','','LogScale',0,...
    'PowerMethod' ,AnalysisOpts.MotifAnalysis.PowerMethod,'maxiter',AnalysisOpts.MotifAnalysis.maxiter,'K',AnalysisOpts.MotifAnalysis.K,'W_fixed',0,...
    'L',AnalysisOpts.MotifAnalysis.L,'lambda',AnalysisOpts.MotifAnalysis.lambda,'lambdaOrthoH',AnalysisOpts.MotifAnalysis.lambdaOrthoH,...
    'lambdaOrthoW',AnalysisOpts.MotifAnalysis.lambdaOrthoW,'lambdaL1W',AnalysisOpts.MotifAnalysis.lambdaL1W,'lambdaL1H',AnalysisOpts.MotifAnalysis.lambdaL1H);
 AnalysisData.MotifSpecs =MotifFunc.RemoveNoloadingMotifs(AnalysisData.MotifSpecs);
 
% reconstruct the data using these motifs 
cwtPower_hat=MotifFunc.reconstruct(AnalysisData.MotifSpecs.W,AnalysisData.MotifSpecs.H);
AnalysisData.cwtPower=ManData.NormPower(AnalysisData.cwt_data,AnalysisData.cwt_f);

% devide this data by f 
cwtPower_hat=cwtPower_hat./AnalysisData.cwt_f;

Ampcwt_hat=sqrt(cwtPower_hat);
anglecwt_hat= (angle(AnalysisData.cwt_data));

cwt_data_hat=Ampcwt_hat.*(cos(anglecwt_hat)+1i*sin(anglecwt_hat));
figure
data_hat=icwt(cwt_data_hat);
data_ds=downsample(data(1:N),FsLFP/FsWaveTarg);
subplot(121)
hold on
Time=0:1/FsWaveTarg:(length(data_ds)-1)/FsWaveTarg;
FigPrms.Plot(Time,data_ds,2,'Time','Volts','X and Xhat');
FigPrms.Plot(Time,data_hat,1,'Time','Volts','X and Xhat');
legend({'X','Xhat'})
title(sprintf('X and Xhat. PEV:%s',num2str(MotifFunc.PEV(data_ds,data_hat'),2)))
subplot(122)
hold on
Time=0:1/FsWaveTarg:(length(data_ds(1:200))-1)/FsWaveTarg;
FigPrms.Plot(Time,data_ds(1:200),2,'Time','Volts','X and Xhat');
FigPrms.Plot(Time,data_hat(1:200),1,'Time','Volts','X and Xhat');
legend({'X','Xhat'})

% plot each motif with an example of phase taken from the real data 
H=AnalysisData.MotifSpecs.H;W=AnalysisData.MotifSpecs.W;
NMotifs=size(H,1);NTim=size(W,3);
TimeMtf=0:1/FsWaveTarg:(NTim-1)/FsWaveTarg;
stdH=std(H,0,2);
figure
for m=1:NMotifs
    ThisMotif=squeeze(W(:,m,:));
    Hind=find(H(1,:)>0.2,1,'last');
    
    Phase=angle(AnalysisData.cwt_data(:,Hind:Hind+NTim-1));
    ThisMotif=ThisMotif./AnalysisData.cwt_f;
    ThisMotif=sqrt(ThisMotif);
    ThisMotif_hat=ThisMotif.*(cos(Phase)+1i*sin(Phase));
    ThisMotif_hat=icwt(ThisMotif_hat);
    subplot(2,NMotifs,m)
    plot(ThisMotif_hat);
    FigPrms.Plot(TimeMtf,ThisMotif_hat,1,'Time','Volts',sprintf('Motis %i TS',m));

    subplot(2,NMotifs,m+NMotifs)
    MotifFunc.PlotCoreMotifs(W,m,[],'LogScale',0)
end

% 
% subplot(NRaw,Ncol,ww+(find(frstp==opts.FrqStp)-1)*Ncol)
% helperCWTTimeFreqPlot(CWTPower,(0:N-1)/FsWave,f_Linear,'justplot1',['q Wave:' opts.wwsTitle{ww}, 'Freq Stp:' num2str(frstp)],'Time(s)','f',0)
% axis square
% % end
% subplot(NRaw,Ncol,ww+(find(frstp==opts.FrqStp))*Ncol)
% FigParams.Plot(f_Linear,FWHM,ww,'Freq(Hz)','FWHM(Hz)',opts.wwsTitle{ww});
% axis square

end
%% put the functions to kick off here
function FS=KickoffMyfunc(RunonCluster,DateNum)
% sets up all of the initial options and path
global AnalysisOpts AnalysisFuncs AnalysisData

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






