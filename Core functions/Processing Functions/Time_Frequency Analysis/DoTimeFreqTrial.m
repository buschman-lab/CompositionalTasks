function [cwt_all,cwt_f,TrueTime]=DoTimeFreqTrial(data,varargin)

global AnalysisOpts
%% Parse optional inputs
%Default save filename is based off of the first file name passed
opts.FreqToPlot = 0.5:0.5:120;%AnalysisOpts.LFPParams.FreqToPlot; %frequencies at which to calculate the PSD
opts.PSDMethod = 'FFT'; %What algorithm to use to estimate PSD.  Can be 'FFT' or 'pwelch'
opts.TimeFreqMethod='Wavelet';  % which method we are using wavelet or shortterm FFT 'STFFT';
opts.Fs= AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency;    
opts.navg=10; % how many trials to average
opts.TrialTimesFields=AnalysisOpts.TrialTimesFields;

%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin),
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end
     
%% Do all of the calculations of power and FFT here 
L=size(data,1); Ntrials=size(data,2);
%%% spectrogram parameters all in ms
Step=1;  % 1ms of step
Window=500;noverlap=(Window-Step)*opts.Fs/1000;
PeriodLength=AnalysisOpts.TrialTiming.PeriodLength*1000;  %%% Whole period duration
StrTime=AnalysisOpts.TrialTiming.BaselineDelay*1000; %%% look at x ms before the start
StpTime=PeriodLength-StrTime;
Time=-1*StrTime:Step:(StpTime-Window); %this is used for windowed short term FFT

TrueTime=-1*StrTime:1000/opts.Fs:StpTime-1000/opts.Fs; %% what was atually in the experiment

for i=1:Ntrials  % loop on differernt trials or channels   
    
         %% Compute power spectrum 
        if strcmp(opts.PSDMethod,'pwelch')
         [temp_amp_psd,temp_amp_f]=pwelch(data(:,i), [],[],L,opts.Fs,'power'); 
        elseif strcmp(opts.PSDMethod,'FFT')   
         [temp_amp_psd,temp_amp_f]=periodogram(data(:,i), [],L,opts.Fs,'power'); 
        end
        
  %Interpolate to specific frequency values
        if ~isempty(opts.FreqToPlot)
           temp_amp_psd = interp1(temp_amp_f, temp_amp_psd, opts.FreqToPlot, 'linear', NaN);
           temp_amp_f = opts.FreqToPlot;
        end
        amp_psd(:,i) = temp_amp_psd.*temp_amp_f; % make sure this multiplication is correct
        amp_f = temp_amp_f;
        
        %% compute Spectectorgram 
        if strcmp(opts.TimeFreqMethod,'Wavelet')  % use wavelet 
            
            [cwt_s,cwt_f] = cwt(data(:,i),'amor',opts.Fs,'VoicesPerOctave',4,...
               'FrequencyLimits',[0 opts.FreqToPlot(end)] );
            cwt_all(:,:,i)=cwt_s;
           
        elseif strcmp(TimeFreqMethod,'STFFT')    %use STFFT      
             [s,w,t]=spectrogram(data(:,i),Window,noverlap,[],opts.Fs);
             cwt_all(:,:,i)=s;
        end
       
end