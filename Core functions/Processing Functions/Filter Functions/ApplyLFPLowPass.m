function Allamp = ApplyLFPLowPass(data,Fs, varargin)

% Isolates the LFP from the signal by filtering out the lower frequenices.
% It then downsamples the signal to a more reasonable sampling rate.
%
% Inputs:
%  stream           stream structure being processed.
%  step_id          ID # of step in stream.
%
% Optional Inputs (must be specified in key-value pairs):
%  Channel          which channel/trial to process. Will only lowpass the specified
%                   channel.
%
%  TimeRange        limit lowpass to a specific time range.
%
%  FilterType       The type of filter to use. Default is 'IIR' (butterworth)
%                   but you can also specify 'FIR' for an equiripple
%                   low-pass filter.
%
%  StopBand         The frequency at which there is a -3dB attenuation.
%                   Default is 150 Hz.
%
%  FiltOrder        The order of the filter to use. Too high causes edge
%                   effects. An order of 20 (default) seems to be a nice
%                   compromise with attenuation beginning around 120 Hz
%                   (-5*10^-4 dB attenuation)

global AnalysisOpts
%% Parse optional inputs

%Default save filename is based off of the first file name passed
opts.Channel = '';
opts.TimeRange = [];
opts.FilterType = 'IIR'; %Can be FIR or IIR
opts.FilterOrder = 4;
opts.PassBand = 150; % Passband Frequency
opts.StopBand = 200; % Stopband Frequency
opts.PassBandRipple = 0.1; % Passband Ripple (dB)
opts.StopBandAttenuation = 40; % Stopband Attenuation (dB)
opts.Verbose = 1; %whether to be chatty about what we are doing
opts.SaveData =0; % save the data

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

if isempty(opts.FilterOrder),
    if opts.StopBand <= opts.PassBand,
        error('StopBand must be greater than the PassBand');
    end
end


 
%% Get channel/trial list
if isempty(opts.Channel),
    chan_list = 1:size(data,2);
else
    chan_list = opts.Channel;
end


%% Loop through channels, applying the notch in turn

%Loop through channels/Trials
for cur_chan = 1:length(chan_list),
    
    in_data.amp =[data(:,cur_chan)]';
    amp_len = size(in_data.amp,2);
    
    %Limit to specified time range?
    if ~isempty(opts.TimeRange),
        overlap_start_offset = round((opts.TimeRange(2) - opts.TimeRange(1))/3);
        overlap_start_offset = min(opts.TimeRange(1) - 1, overlap_start_offset); %make sure we don't go before beginning
        overlap_end_offset = round((opts.TimeRange(2) - opts.TimeRange(1))/3);
        overlap_end_offset = min(amp_len - opts.TimeRange(2), overlap_end_offset); %make sure we don't go past end
        amp = in_data.amp(1, (opts.TimeRange(1) - overlap_start_offset):(opts.TimeRange(2) + overlap_end_offset));
    else
        amp = in_data.amp;
    end
    clear('in_data');
    
    %Design filter
    if isempty(opts.FilterOrder),
        %Use the full specification
        d=fdesign.lowpass('Fp,Fst,Ap,Ast', opts.PassBand, opts.StopBand, opts.PassBandRipple, opts.StopBandAttenuation, ...
            opts.SamplingFrequency);
    else
%         d=fdesign.lowpass('N,F3dB', opts.FilterOrder, opts.PassBand, opts.SamplingFrequency);
         d=fdesign.lowpass('Fp,Fst,Ap,Ast', opts.PassBand, opts.StopBand, opts.PassBandRipple, opts.StopBandAttenuation, ...
            opts.SamplingFrequency);
    end
    if strcmpi(opts.FilterType, 'FIR'),
        Hd = design(d, 'equiripple');
    elseif strcmpi(opts.FilterType, 'IIR'),
        if isempty(opts.FilterOrder),
            Hd = design(d, 'butter', 'MatchExactly', 'passband');
        else
            Hd = design(d, 'butter');
        end
    end
    
    %Get the length of the impulse response -- we'll use this to avoid edge
    %effects
    impResp = impz(Hd, [], opts.SamplingFrequency);
    imp_len = length(impResp);
    
    %If we have NaNs, we need to interpolate over them and then add them
    %back (with padding)
    fill_non_finite_gap_starts = []; fill_non_finite_gap_ends = [];
    if any(~isfinite(amp)),
        non_finite_ind = find(~isfinite(amp));
        non_finite_gap_starts = non_finite_ind([1 (find(diff(non_finite_ind) ~= 1) + 1)]);
        non_finite_gap_ends = non_finite_ind([find(diff(non_finite_ind) ~= 1) length(non_finite_ind)]);
        
        %We are going to end up filling back with NaNs -- more because of filter carrying over
        fill_non_finite_gap_starts = non_finite_gap_starts - 2*imp_len;
        fill_non_finite_gap_starts = max(1, min(length(amp), fill_non_finite_gap_starts));
        fill_non_finite_gap_ends = non_finite_gap_ends + 2*imp_len;
        fill_non_finite_gap_ends = max(1, min(length(amp), fill_non_finite_gap_ends));
        
        %Let's interpolate over the gaps
        for i = 1:length(non_finite_gap_starts),
            if non_finite_gap_starts(i) == 1,
                %At the beginning, just fill with the first non-NaN value
                amp(non_finite_gap_starts(i):non_finite_gap_ends(i)) = amp(non_finite_gap_ends(i)+1);
            elseif non_finite_gap_ends(i) == length(amp),
                %Goes to the end, just fill with the last non-NaN value
                amp(non_finite_gap_starts(i):non_finite_gap_ends(i)) = amp(non_finite_gap_starts(i)-1);
            else
                %Linearly interpolate between the two
                amp(non_finite_gap_starts(i):non_finite_gap_ends(i)) = interp1([non_finite_gap_starts(i)-1 non_finite_gap_ends(i)+1], amp([non_finite_gap_starts(i)-1 non_finite_gap_ends(i)+1]), [non_finite_gap_starts(i):non_finite_gap_ends(i)]);
            end
        end
    end
    
    %Apply our filter
    if strcmpi(opts.FilterType, 'FIR'),
        %Apply filter
        if opts.Verbose, fprintf('Applying FIR (equirriple) low-pass filter to amplifier signal...\n'); end
        amp = filtfilt(Hd.Numerator, 1, amp);
    elseif strcmpi(opts.FilterType, 'IIR'),
        %Apply filter
        if opts.Verbose, fprintf('Applying IIR (butterworth) low-pass filter to amplifier signal...\n'); end
        %Apply filter (padding to avoid some edge effects of high-order filters) and save
        amp = filtfilt(Hd.SOSMatrix, Hd.ScaleValues, cat(2, zeros(size(amp, 1), imp_len), amp, zeros(size(amp, 1), imp_len)));
        amp = amp(:, (imp_len+1):(end-imp_len));
    else
        error('Unknown filter specified.');
    end
    
    %If needed, we should fill back in the gaps
    if ~isempty(fill_non_finite_gap_starts),
        for i = 1:length(fill_non_finite_gap_starts),
            amp(fill_non_finite_gap_starts(i):fill_non_finite_gap_ends(i)) = NaN;
        end
    end
     
    
      
    %Was I limited to a specified time range?
    if ~isempty(opts.TimeRange),    
        %Remove edges use to avoid edge effects
        amp = amp((1 + overlap_start_offset):(end - overlap_end_offset));
        %Write to file
         sprintf('Applied LFP filter to channel %s, time range %d-%d.', chan_list(cur_chan), opts.TimeRange);
    end
        Allamp(:,cur_chan)=amp; %% save data for all of channels or trials

    
end %channel loop

