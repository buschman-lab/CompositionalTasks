function CheckRawData(data,Fs,varargin)
%CHECKRAWDATA Summary of this function goes here
% checks rawdata, plots some traces and simple freqency analysis 
%   Detailed explanation goes here
if size(data,2)>size(data,1)  % samples should be on the raw
    data=data';
end

opts.verbose=1;
opts.TotalRecTime=size(data,1)/Fs;
opts.Nsamp=size(data,1);
opts.ExmpL=1;  % length of example data points in sec
opts.ExmpN=5;  % number of examples 
opts.PreProcessData=0;  % should we apply notch filter and detrend to this data?

opts = ParseOptionalInputs(opts,varargin);
% setup relevant functions 
fp=fig_params;
ManData=ManipulateData;
TimeFreq=TimeFreqAnalysis;  % time fre analysis

% grab some initial data

data=double(data);
% detrend the data first
if opts.PreProcessData
    data = detrend(data);
    data = ApplyNotchFilter(data,Fs);
    data = ApplyLFPLowPass(data,Fs,'PassBand',120,'StopBand',120+10);
    data=data';
end
Time=0:1/Fs:opts.TotalRecTime-1/Fs;
% plot different chunks of data to double check 
% take three sample data points along the way 
ExmpTimPint=1:floor(opts.Nsamp/opts.ExmpN):opts.Nsamp;
Figs=fp.RenderFigure(1,[]);

    for i=1:opts.ExmpN
        subplot(3,opts.ExmpN,i)
        TimeInd=[ExmpTimPint(i):ExmpTimPint(i)+floor(opts.ExmpL*Fs)];
        ThisTime=Time(TimeInd);ThisData=data(TimeInd);
        plot(ThisTime,ThisData);axis tight
        
        % plot Time_Freq
        subplot(3,opts.ExmpN,i+opts.ExmpN)
       [cwt_s,cwt_f] = cwt(ThisData,'amor',Fs,'VoicesPerOctave',4,'FrequencyLimits',[1 120]);
       % [cwt_s,cwt_f]=TimeFreq.ComputeTimeFreq(ThisData,'Fs',Fs,'WaveletMethod','code','WaveletWidth',7);
        TF_power=ManData.Power(cwt_s);
     %   TF_power=ManData.NormPower(cwt_s,cwt_f);
        %TF_power=transpose(downsample(TF_power',Fs/50));
        helperCWTTimeFreqPlot(TF_power,ThisTime,cwt_f,'justplot1','PSD','Time','Freq',0)
        
        % plot Psd
        subplot(3,opts.ExmpN,i+2*opts.ExmpN)
        TimeFreq.PlotPSD(ThisData,'Fs',Fs)   
    end
end

