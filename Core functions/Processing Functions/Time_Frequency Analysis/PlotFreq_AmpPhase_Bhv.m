function PlotFreq_AmpPhase_Bhv(data,TrialTimes,varargin) 
%%% plots phase and amplitude of freq spectrum signal and behavioral
%%% performance 
%%% inputs 
%%% data             data time series. Trial*Time*Channel
%%% FreqToPlot       Cut off frequency to look 
%%% PSDMethod        The method to calculte the power spectrum of the data
%%% Fs               Sampling Frequency of data

global AnalysisOpts
%% Parse optional inputs
saveMoviePath=[AnalysisOpts.ResultsSavePath sprintf('Angle_Ch_%s_Area_%s',num2str(AnalysisOpts.CurrentChNum),AnalysisOpts.CurrentChArea)];
%Default save filename is based off of the first file name passed
opts.FreqToPlot = [0.5:0.5:80];%AnalysisOpts.LFPParams.FreqToPlot; %frequencies at which to calculate the PSD
opts.PSDMethod = 'FFT'; %What algorithm to use to estimate PSD.  Can be 'FFT' or 'pwelch'
opts.TimeFreqMethod='Wavelet';  % which method we are using wavelet or shortterm FFT 'STFFT';
opts.Fs= AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency;    
opts.SubPlot=[2 2];
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
L=size(data,1); 
%%% spectrogram parameters all in ms
Step=1;  % 1ms of step
Window=500;noverlap=(Window-Step)*opts.Fs/1000;
PeriodLength=AnalysisOpts.TrialTiming.PeriodLength*1000;  %%% Whole period duration
StrTime=AnalysisOpts.TrialTiming.BaselineDelay*1000; %%% look at x ms before the start
StpTime=PeriodLength-StrTime;
Time=-1*StrTime:Step:(StpTime-Window); %this is used for windowed short term FFT

TrueTime=-1*StrTime:1000/opts.Fs:StpTime-1000/opts.Fs; %% what was atually in the experiment

for i=1:size(data,2)  % loop on differernt trials or channels   
    
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
            s_all(:,:,i)=cwt_s;
           
        elseif strcmp(TimeFreqMethod,'STFFT')    %use STFFT      
             [s,w,t]=spectrogram(data(:,i),Window,noverlap,[],opts.Fs);
              s_all(:,:,i)=s;
        end
       
end
figure;set(gcf,'Units','Normalized','Position',[0 0 0.9 0.9])
%compute frequency spectrum using FFT
Yfft=fft(data); 
         
%Y(1)=0; %% remove DC
%% Plot Amplitude in db
P2=abs(Yfft);
P1=P2(1:L/2+1,:);
P1(2:end-1,:)=2*P1(2:end-1,:);
f = [ opts.Fs*(0:(L/2))/L]';  %% frequency vector
LFP_F=f>0.5 & f<=opts.FreqToPlot(end);

subplot(opts.SubPlot(1),opts.SubPlot(2),1)
FFT_Amp=P1(LFP_F,:);
PlotMeanStd(FFT_Amp,f(LFP_F),'f (Hz)','|P1(f)| db',1)
title('Single-Sided Amplitude Spectrum of X(t)')


%% plot Phase
phaseY =  (angle(Yfft(1:L/2+1,:)));
subplot(opts.SubPlot(1),opts.SubPlot(2),2)
PlotMeanStd(phaseY(LFP_F,:),f(LFP_F),'f (Hz)','radians',2)
title('Phase of X(t)')

%% plot power spectrum
subplot(opts.SubPlot(1),opts.SubPlot(2),3)
hold on
PlotMeanStd(amp_psd,amp_f,'f (Hz)','Power(dbW)',3)
title('Power spectrum of X(t)')

%% Plot Spectectogram
subplot(opts.SubPlot(1),opts.SubPlot(2),4)

helperCWTTimeFreqPlot(mean(s_all,3),TrueTime,cwt_f,'surf','PSD','Time from Stim Onset(ms)','Frequency(Hz)')

% w=w(w<=opts.FreqToPlot(end));   % cut the frequency 
% CutPSD=PSD_all(1:length(w),:,:);
% imagesc((mean(CutPSD,3)));colorbar
% Lt=length(t);Lw=length(w);
% XTik=floor([ 1 find(Time==0) Lt]);
% YTik=floor([ 1 Lw/4 Lw/2 Lw*3/4 Lw]);
% xticks(XTik);yticks(YTik);
% xticklabels(arrayfun(@(x) num2str(Time(x),3),XTik,'Uniformoutput',0))
% yticklabels(arrayfun(@(x) num2str(w(x),3),YTik,'Uniformoutput',0))
% xlabel('Time from Stim Onset (ms)')
% ylabel('PSD')
% title('All Trails Average PSD')

mvFrame(1) = getframe(gcf);
%% now make a movie of each trial's psd during learning for each neuron 
MeanPSD=movmean(s_all,opts.navg,3,'Endpoints','shrink');
MeanLFPsig=movmean(data,opts.navg,3,'Endpoints','shrink');
figure;set(gcf,'Units','Normalized','Position',[0 0 0.9 0.9])

Perf=PlotBehavioralPerf(TrialTimes,opts);  % behavior
%%find correltion of PSD to behavior
for i=1:size(MeanPSD,1)
    for j=1:size(MeanPSD,2)
        A=squeeze(MeanPSD(i,j,:));
        Bhv_PSD_Corr(i,j)=corr(Perf',abs(A).^2);       
    end
end


%%
colormap(jet);PSDaxisLim=[0 5e-5]; CorraxisLim=[0 0.5];

for i=1:size(MeanPSD,3)
    
   %% plot PSD here
   subplot(opts.SubPlot(1),opts.SubPlot(2),1)
   % helperCWTTimeFreqPlot(MeanPSD(:,:,i).*repmat(cwt_f,1,size(MeanPSD,2)),TrueTime,cwt_f,'surf','PSD','Time from Stim Onset(ms)','Frequency(Hz)')
   helperCWTTimeFreqPlot(MeanPSD(:,:,i),TrueTime,cwt_f,'surf','PSD','Time from Stim Onset(ms)','Frequency(Hz)')
   caxis(PSDaxisLim)
   
   %% plot Behavior here
   subplot(opts.SubPlot(1),opts.SubPlot(2),2);hold on
   plot(1:length(Perf),Perf,'k');
   plot(i,Perf(i),'r*')
   xlabel('Trial')
   ylabel('Percent Correct')
   title(' Behavioral Performance')
   
   %% plot average LFP waveform 
   subplot(opts.SubPlot(1),opts.SubPlot(2),3); 
   plot(TrueTime,MeanLFPsig(:,i),'r');
   xlabel('Time from Stim Onset(ms)');
   ylabel('Volatage');
   title(sprintf('Avg LFP; Ch %s Area %s Trial %s',num2str(AnalysisOpts.CurrentChNum),AnalysisOpts.CurrentChArea,num2str(i)));
  
   %% plot correlation of behavior to PSD
   subplot(opts.SubPlot(1),opts.SubPlot(2),4); 
   helperCWTTimeFreqPlot(Bhv_PSD_Corr,TrueTime,cwt_f,'Corrlation','Correlation of PSD to Behavior','Time from Stim Onset(ms)','Frequency(Hz)')
  % caxis(CorraxisLim)   
    pause(0.05) 
    mvFrame(i+1) = getframe(gcf);
end
% now make the movie 
fps=2;
MakeMovieFromFrames(mvFrame,fps,saveMoviePath)

end
%%% code that is using multi taper to estimate PSD over time 
% movingwin=[0.5 0.001];
% params .Fs=1000; params .fpass=[0 opts.FreqToPlot(end)];
% params .tapers=[5 9]; params .trialave=1; params .err=0;
% [S1,t,f]=mtspecgramc (data(:,i) ,movingwin,params);
% plot_matrix(S1,t,f);  
% colorbar;
%     
% 
function PlotMeanStd(amp,amp_f,Xlabel,Ylabel,ColInd)
if nargin<5
    ColInd=3;
end
co = get(gca, 'ColorOrder');
plot(amp_f, nanmean(10*log10(amp), 2), '-', 'Color', co(mod(ColInd-1, size(co, 1))+1, :)); hold on;
plot(amp_f, nanmean(10*log10(amp), 2) + nanstd(10*log10(amp), [], 2), ':', 'Color', co(mod(ColInd-1, size(co, 1))+1, :));
plot(amp_f, nanmean(10*log10(amp), 2) - nanstd(10*log10(amp), [], 2), ':', 'Color', co(mod(ColInd-1, size(co, 1))+1, :));  
xlabel(Xlabel);ylabel(Ylabel)
grid on;
axis tight
end

function Perf=PlotBehavioralPerf(TrialTimes,opts)

CorrectInd=strcmp(opts.TrialTimesFields,'CORRECT_TRIAL');%% we don't use this here but still good to have it

Bhv=~isnan([TrialTimes(:,CorrectInd)]');
sm_kern = ones(1, opts.navg); 
sm_kern = sm_kern./sum(sm_kern);
Perf=convn(Bhv, sm_kern, 'same');

end
