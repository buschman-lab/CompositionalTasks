function PrepareEphysData(RunonCluster,DateNum,ChNum,AreaNum,PairNum,varargin)
%% takes each channel and computes LFP data and applies wavelet transform
% define some global functions to carry options and data
% input RunonCluster: runs on cluster
% varargin(1) number of recording date
% Sina Tafazoli timeless
global AnalysisOpts 

%% define primary vars before kick off
AnalysisOpts.RunonCluster=RunonCluster;
AnalysisOpts.Project='Rule Representation';  
AnalysisOpts.Animal='';
AnalysisOpts.AnalysisType='Analysis Pipeline'; 
AnalysisOpts.AnalysisFocus1  ='';  % define focus of this analysis
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
HelperFunc=HelperFuncs;
[TimeFreq,FigPrms,ManData,TrialFunc,FilterFunc,NeuAnaFunc,CrosFreqCop,AnalysisFunc]=HelperFunc.LoadAnalysisClasses; % load all of the classes we want
[Fs,FsLFP,FsWave,FsWaveTarg]=HelperFunc.SetupGeneralFuncs; % set up all general functions we want

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Starting the analysis of this function %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define filtering functions first
f1=0.5;f2=150;nf=49; % range and binning of frequencies to look at 
% f_log=2.^(log2(f1):(log2(f2)-log2(f1))/nf:log2(f2)); % log frequency
% scale
 f_lin=f1:(f2-f1)/nf:f2; % linear frequency scale 
%f_lin=[0.25 0.5 0.75 ManData.BinData(AnalysisOpts.MotifAnalysis.FreqStep,[],f1,f2)]; % define frequency axis


%%% loop thourgh Channels
extraStr=['Pr' num2str(AnalysisOpts.ProbeNum) 'Width' num2str(AnalysisOpts.MotifAnalysis.WaveletWidth)]; % add this to end of everyfile we write

for Ch= ChNum %%% loop thourgh Useable Channels
    % if this channel already exists and we don't want to rewrite just
    % skip it
    AnalysisOpts.CurrentCh=Ch;
    if ManData.IsVarExistinFileAnalysis('Wavelet_Data','Wavelet_Linear',extraStr) & ~AnalysisOpts.ReWriteData
        fprintf('\nChannel %i  already exists; skippking this channel ...',Ch)
    else
        %%% now load the raw data of the channel
        tic
        %% get extra data about this channel
        AnalysisOpts.CurrentCh=Ch;
        % AreaNum=TrialFunc.FindChArea(Ch,ChannelInfo); % for future to
        % retrive the area name
        [ThisChData,ChDepth,Opts]=HelperFunc.getLFPEphysData(ChNum);
        ThisChData=ThisChData{1};
        Fs=Opts.lfp_freq;
        AnalysisOpts.CurrentChDepth=ChDepth{1};

        fprintf('\n ***************************************************');
        fprintf('\n Started processing channel %s in area %s, recording %s  \n',num2str(Ch),AnalysisOpts.Area_2look{AreaNum},AnalysisOpts.RecDate);
        %% do we need to rereference the data 
%         if AnalysisOpts.ReRefChannels 
%             ThisChWell=AnalysisOpts.AreaWells(AreaNum); % find the well of this channel
%             ThisChData=ThisChData-MeanWellRawData{ThisChWell};
%         end
        %% preprocess for Waveletdata
        data= detrend(transpose(ThisChData));
    %   data = ApplyNotchFilter(data,Fs); % notch filter the data
        data = ApplyLFPLowPass(data,Fs,'PassBand',f2,'StopBand',f2+10);  % apply low pass filter to only look at lfp
        data = downsample(data,Fs/FsLFP);        % downsample the data to a reasonable range
        RawData=data';
        
        %% do wavelet transform with linear freq scale
        [Wavelet_Linear,f_Linear,~] = TimeFreq.ComputeTimeFreq(data,'Fs',FsLFP,...
            'WaveletMethod','Code','FreqToPlot',f_lin,'WaveletWidth',AnalysisOpts.MotifAnalysis.WaveletWidth);
         % plot an example snippet of this channel
        RawDataFig=cell(1);
        RawDataFig{:}=AnalysisFunc.PlotRawData(Wavelet_Linear,f_Linear,0:1/FsLFP:2,data,[],'PowerMethod','NormPower');
%         figure
%         [Sxx,faxis,df,fNQ]=TimeFreq.CalPSD(data',FsLFP,1); 
%         TimeFreq.PlotPSD(Sxx,faxis,[0 150],'r')
        % now downsample the data and save important variables in the file 
        Wavelet_Linear=transpose(downsample(Wavelet_Linear',FsLFP/FsWave)); % downsample to FsWave
        ParamSpace.Fs=Fs;
        ParamSpace.FsWave=FsWave;
        ParamSpace.FsLFP=FsLFP;
        ParamSpace.Channel=Ch;
        ParamSpace.AreaNum=AreaNum;
        
        %             %% do wavelet transform with lof freq scale
        %             [DataSpace.Wavelet_Log,DataSpace.f_Log,~] = TimeFreq.ComputeTimeFreq(data,'Fs',FsLFP,...
        %                 'WaveletMethod','Code','FreqToPlot',f_log);
        %             DataSpace.Wavelet_Log=transpose(downsample(DataSpace.Wavelet_Log',FsLFP/FsWave)); % downsample to 100hz
        %

         % some motif analysis 
        [motif_onset,st_norm] = ReadMotifs(DateNum,1);
        TrigFigs=cell(1);
        TrigFigs{:}=AnalysisFunc.HTriggeredAnalysis(Wavelet_Linear,f_Linear,RawData,motif_onset,st_norm,0.5);
        [~,~,AnalysisData.MotTrigWav]=ManData.GetFileName([],['MotTrigWavCh' num2str(Ch) 'Pr' num2str(AnalysisOpts.ProbeNum)],'SaveInResults',1);
        FigPrms.SaveFigSeries([],AnalysisData.MotTrigWav,[TrigFigs])

        %saveas(gcf,['H_trig_Ch' num2str(Ch)],'jpg')
        %% save data
        ManData.SaveVar('Wavelet_Data',Wavelet_Linear,'Wavelet_Linear',extraStr);
        ManData.SaveVar('Wavelet_Data',f_Linear,'f_Linear',extraStr);
        ManData.SaveVar('Wavelet_Data',RawData,'RawData',extraStr);
        ManData.SaveVar('Wavelet_Data',ParamSpace,'ParamSpace',extraStr);
        ManData.SaveVar('Wavelet_Data',AnalysisOpts,'AnalysisOpts',extraStr);
        
        clear data ParamSpace RawData
        ElapsedTime=num2str(toc);
        fprintf('\n Finished this Channel. Elapsed time: %s \n',ElapsedTime)
       
    end
end

%% after this is done produce a recording summery
%ProduceRecordingSummery(RunonCluster,DateNum,ChsSet,[],PairNum,varargin)

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
    RootPath='Z:\';
    FS='\';
end
% add core function path first so we can sturt up everything
addpath(genpath([RootPath 'Projects' FS 'Rule_Representation' FS 'ElecPhys_Analysis' FS 'Core functions' FS]));
SetupAllVars(DateNum)  %% set up the path and initialize vars

end


