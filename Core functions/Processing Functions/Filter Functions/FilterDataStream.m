function [TrialLFP_Data,Time,opts]=FilterDataStream(RawData,TrialTimes)
%%%filters data in prescribed trials in various frequencies 
%%% inputs: Raw data: matfile ,atrix pointing to the data 
%%% TrialTimes: Time of the trials 

global AnalysisOpts
  
%% filter the data now 

%[Filter,impResp]=DesignLFPfilter(AnalysisOpts.LFPParams.FilterOpts);
%FiltData=ApplyFilter(RawData,Filter,impResp,AnalysisOpts.LFPParams.FilterOpts);


%bandpass(RawData,[100 200],AnalysisOpts.LFPParams.FilterOpts.SamplingFrequency)

SamplingFrequency = 2000;
FreqRatio=SamplingFrequency/AnalysisOpts.LFPParams.FilterOpts.SamplingFrequency;
opts.FilterOrder = 4;StepFreq=20;
opts.BandSet=[1:StepFreq:140;(1:StepFreq:140)+StepFreq]';
for FrBand=1:size(opts.BandSet,1)   %% loop on frequency bands
    opts.PassBand(1) = opts.BandSet(FrBand,1);
    opts.PassBand(2) = opts.BandSet(FrBand,2);
    fprintf('\n Band pass filtering data in band: %s to %s Hz',num2str(opts.PassBand(1)),num2str(opts.PassBand(2)));
    d=fdesign.bandpass('N,F3dB1,F3dB2', opts.FilterOrder, opts.PassBand(1), opts.PassBand(2), SamplingFrequency);
    Hd = design(d, 'butter');
    impResp = impz(Hd, [], SamplingFrequency);
    imp_len = length(impResp); % Get the length of the impulse response -- we'll use this to avoid edge effects
     %Apply filter
    amp = downsample(RawData,1/FreqRatio);
    amp_len = length(amp);
    amp = filtfilt(Hd.SOSMatrix, Hd.ScaleValues, cat(2, zeros(size(amp, 1), imp_len), amp, zeros(size(amp, 1), imp_len)));
    LFPdata(FrBand,:) = amp(:, (imp_len+1):(end-imp_len)); 
end
%Time=0:FreqRatio:(amp_len-1)*FreqRatio;
Time=downsample(0:1/40000:(length(RawData)-1)/40000,1/FreqRatio);

%% this section cuts the trials based on their timing 
TrialLFP_Data=GrabDataTimeTrial(data,TrialTimes,opts.BandSet);