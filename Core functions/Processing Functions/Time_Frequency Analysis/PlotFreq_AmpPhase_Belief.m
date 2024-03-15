function PlotFreq_AmpPhase_Belief(data,TrialTimes,belief,saveMoviePath,varargin) 
%%% plots phase and amplitude of freq spectrum signal and behavioral
%%% performance 
%%% inputs 
%%% data             data time series. Trial*Time*Channel
%%% FreqToPlot       Cut off frequency to look 
%%% PSDMethod        The method to calculte the power spectrum of the data
%%% Fs               Sampling Frequency of data

global AnalysisOpts
%% Parse optional inputs
fp=fig_params;
close all
%Default save filename is based off of the first file name passed
opts.FreqToPlot = [0.5:0.5:80];%AnalysisOpts.LFPParams.FreqToPlot; %frequencies at which to calculate the PSD
opts.PSDMethod = 'FFT'; %What algorithm to use to estimate PSD.  Can be 'FFT' or 'pwelch'
opts.TimeFreqMethod='Wavelet';  % which method we are using wavelet or shortterm FFT 'STFFT';
opts.Fs= AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency;    
opts.SubPlot=[3 3];
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

helperCWTTimeFreqPlot(mean(s_all,3),TrueTime,cwt_f,'power','Average PSD','Time from Stim Onset(ms)','Frequency(Hz)')

mvFrame(1) = getframe(gcf);
%% now make a movie of each trial's psd during learning for each neuron 
MeanPSD=movmean(s_all,opts.navg,3,'Endpoints','discard');
MeanLFPsig=movmean(data,opts.navg,2,'Endpoints','discard');
belief=PlotBelief(belief,opts);
figure;set(gcf,'Units','Normalized','Position',[0 0 0.9 0.9])

Perf=PlotBehavioralPerf(TrialTimes,opts);  % behavior
%%find correltion of PSD to behavior and belief
for i=1:size(MeanPSD,1)
    for j=1:size(MeanPSD,2)
        A=squeeze(MeanPSD(i,j,:));
        [Bhv_PSD_Corr(i,j),p_Bhv_PSD_Corr(i,j)]    =corr(Perf',abs(A).^2);
        [Bhv_Phase_Corr(i,j),p_Bhv_Phase_Corr(i,j)]=corr(Perf',angle(A));  
        
        for b=1:3
            [Belief_PSD_Corr(i,j,b),p_Belief_PSD_Corr(i,j,b)]    =corr(belief(:,b),abs(A).^2) ;          
            [Belief_Phase_Corr(i,j,b),p_Belief_Phase_Corr(i,j,b)]=corr(belief(:,b),angle(A));
        end
    end
end
Belief_PSD_Corr(p_Belief_PSD_Corr>0.01)=0; % just look at significant values
Bhv_PSD_Corr(p_Bhv_PSD_Corr>0.001)=0;
%% Project the PSD data to low dimensions and plot them 
%as the first step take a mean acorss time 
Trial_s_all=transpose(squeeze(mean(MeanPSD,2)));
[score] = tsne(abs(Trial_s_all).^2);


%%
colormap(jet);PSDaxisLim=[0 5e-5]; CorraxisLim=[-0.6 0.6];pvalAxisLim=[0 1];
GrayMap=colormap(gray(Ntrials));
for i=1:size(MeanPSD,3)
    
   %% plot PSD here
   subplot(opts.SubPlot(1),opts.SubPlot(2),1)
   % helperCWTTimeFreqPlot(MeanPSD(:,:,i).*repmat(cwt_f,1,size(MeanPSD,2)),TrueTime,cwt_f,'surf','PSD','Time from Stim Onset(ms)','Frequency(Hz)')
   helperCWTTimeFreqPlot(MeanPSD(:,:,i),TrueTime,cwt_f,'power','PSD','Time from Stim Onset(ms)','Frequency(Hz)')
   caxis(PSDaxisLim)
   
   %% plot Phase here
   subplot(opts.SubPlot(1),opts.SubPlot(2),2)
   helperCWTTimeFreqPlot(MeanPSD(:,:,i),TrueTime,cwt_f,'phase','Phase','Time from Stim Onset(ms)','Frequency(Hz)')

    %% plot low dimensional plot here
   sp3=subplot(opts.SubPlot(1),opts.SubPlot(2),3);cla(sp3);
   hold on;sp3.Colormap=GrayMap;
  % helperCWTTimeFreqPlot(MeanPSD(:,:,i),TrueTime,cwt_f,'powerphase','Phase-Power','Time from Stim Onset(ms)','Frequency(Hz)')
   arrayfun(@(x) plot(score(x,1),score(x,2),'.','Markersize',15,'color',GrayMap(x,:)),1:size(score,1))
   plot(score(i,1),score(i,2),'r*','MarkerSize',10);
   xlabel('Tsne1');
   ylabel('Tsne2');
   title('PSD projection to TSNE')
   %% plot Behavior and belief here
   
   sp4=subplot(opts.SubPlot(1),opts.SubPlot(2),4);cla(sp4)
   hold on   
   plot(1:length(Perf),Perf,'k','linewidth',fp.p_line_width); %plots performance
   plot(1:length(Perf),belief,'linewidth',fp.p_line_width);  
   plot([i i],[0 1],'r')
   legend({'Perf','BRule1','BRule2','BRule3'},'Location','northeastoutside')
   xlabel('Trial')
   ylabel('% Correct/Belief')
   title([' Beh Perf for Rule' num2str(opts.CurrentRule) 'avg Trial ' num2str(opts.navg)])
      
   %% plot average LFP waveform 
   subplot(opts.SubPlot(1),opts.SubPlot(2),5); 
   plot(TrueTime,MeanLFPsig(:,i),'r');
   xlabel('Time from Stim Onset(ms)');
   ylabel('Volatage');
   title(sprintf('Avg LFP; Ch %s Area %s Trial %s',num2str(AnalysisOpts.CurrentChNum),AnalysisOpts.CurrentChArea,num2str(i)));
  
   %% plot correlation of behavior to PSD
   subplot(opts.SubPlot(1),opts.SubPlot(2),6); 
   helperCWTTimeFreqPlot(Bhv_PSD_Corr,TrueTime,cwt_f,'Corrlation','Correlation of PSD to Behavior','Time from Stim Onset(ms)','Frequency(Hz)')
   caxis(CorraxisLim)   
  
  %% plot correlation of PSD to Belief about Rule 1,2,3
  for r=1:3
       subplot(opts.SubPlot(1),opts.SubPlot(2),6+r); 
       helperCWTTimeFreqPlot(Belief_PSD_Corr(:,:,r),TrueTime,cwt_f,'Corrlation',['Corr of PSD to Belief of Rule ',num2str(r)],'Time from Stim Onset(ms)','Frequency(Hz)')
       caxis(CorraxisLim)  
       %% plot pvalues
%        subplot(opts.SubPlot(1),opts.SubPlot(2),9+r); 
%        helperCWTTimeFreqPlot(p_Belief_PSD_Corr(:,:,r),TrueTime,cwt_f,'pval',['pval of PSD to Belief of Rule ',num2str(r)],'Time from Stim Onset(ms)','Frequency(Hz)')
%       caxis(pvalAxisLim)  
  end
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
Perf=convn(Bhv, sm_kern, 'valid');

end


function belief=PlotBelief(Belief,opts)

sm_kern = ones(1, opts.navg); 
sm_kern = sm_kern./sum(sm_kern);
belief=transpose(convn(Belief', sm_kern, 'valid'));

end
