 function out_data = ApplyWaveletTransform(data,Fs, varargin)

% Does a wavelet transformation of the amplifier signal.
%
% Inputs:
%  stream           stream structure being processed.
%  step_id          ID # of step in stream.
%
% Optional Inputs (must be specified in key-value pairs):
%  Channel          which channel to process. Will only transform the specified
%                   channel.
%
%  TimeRange        Do wavelet transformation on a specific time range --
%                   will expand this to ensure there are no edge effects.
%                   This supports parallelization over time.
%
%  Freqs            Center frequencies of wavelets.
%                   Default is 2.^[0:0.25:7].
%
%  WaveletWidth     The width to use for wavelets.  This sets up the
%                   time-frequency trade-off to use. Default is 5.  Can be
%                   a vector to specify different widths for different
%                   frequencies of interest.

global AnalysisOpts
%% Parse optional inputs

%Default save filename is based off of the first file name passed
opts.Channel = '';
opts.TimeRange = [];
opts.Freqs = 2.^[0:0.25:7];
opts.WaveletWidth = 5;
opts.WaveletTimePaddingWidth = 3;
opts.Verbose = 1; %whether to be chatty about what we are doing
opts.SamplingFrequency=Fs;  % Sampling Frequency of data

%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin),
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end
  
 
%% Get channel/trial list
if isempty(opts.Channel),
    chan_list = 1:size(data,2);
else
    chan_list = opts.Channel;
end
%% Loop through channels, applying the notch in turn

%Get wavelets
[wlt, sf, st] = GenerateWavelets(opts.Freqs, opts.WaveletWidth, opts.WaveletTimePaddingWidth, opts.SamplingFrequency);

%Loop through channels
for cur_chan = 1:length(chan_list),
   
    in_data.amp =[data(:,cur_chan)]';
    amp_len = size(in_data.amp,2);
    
     
    %Limit to specified time range?
    if isempty(opts.TimeRange),
        out_amp_ind = [1 amp_len];
    else        
        out_amp_ind = opts.TimeRange(1:2);
    end
    out_amp_ind = min(max(out_amp_ind, 1), amp_len);
    
    %Loop through wavelets
    for cur_freq = 1:length(opts.Freqs),
        %Get the input index
        inp_amp_ind = out_amp_ind + [-1 1]*length(wlt{cur_freq});
        inp_amp_ind = min(max(inp_amp_ind, 1), amp_len);
        
        %Do wavelet transformation
        if length(wlt{cur_freq}) < (diff(inp_amp_ind)+1),
            temp_amp = conv(in_data.amp(1, inp_amp_ind(1):inp_amp_ind(2)), wlt{cur_freq}(:)', 'same');
          %  temp_amp(1:(length(wlt{cur_freq})-1)) = NaN;
          %  temp_amp((end - length(wlt{cur_freq}) + 1):end) = NaN;
        else
            temp_amp = NaN*ones(1, diff(inp_amp_ind)+1);
        end
         
 
%         %Open up output matfile (initializing it if necessary)
%         if ~exist(chan_out_fn, 'file'),
%             out_data = matfile(chan_out_fn, 'Writable', true);
%             out_data.amp = complex(NaN*ones(1, amp_len));
%             %Save options and wavelet parameters to file
%             out_data.opts = opts;
%             out_data.f = opts.Freqs;
%             out_data.wlt = wlt;
%             out_data.sf = sf;
%             out_data.st = st;
%         else
%             out_data = matfile(chan_out_fn, 'Writable', true);
%         end
        
        out_data{cur_chan}.amp(cur_freq, [out_amp_ind(1):out_amp_ind(2)]) = temp_amp([out_amp_ind(1):out_amp_ind(2)] - inp_amp_ind(1) + 1);
        
    end %frequency loop
      % Save options and wavelet parameters to file
        out_data{cur_chan}.opts = opts;
        out_data{cur_chan}.f = opts.Freqs;
        out_data{cur_chan}.wlt = wlt;
        out_data{cur_chan}.sf = sf;
        out_data{cur_chan}.st = st;
end %channel loop

 