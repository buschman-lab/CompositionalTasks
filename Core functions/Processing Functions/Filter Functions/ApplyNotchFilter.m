function [Allamp,opts,cur_NotchFrequencies] = ApplyNotchFilter(data,Fs, varargin)

% Notches specified frequencies.  Can be used to notch noise by first
% detecting whether there is a peak.
%
% Inputs:
%  stream        stream structure being processed.
%  step_id          ID # of step in stream.
%
% Optional Inputs (must be specified in key-value pairs):
%  Channel          which channel to process. Will only notch the specified
%                   channel.
%
%  TimeRange        limit notch to a specific time range (specified in the
%                   # of samples of the signal).
%
%  NotchFrequencies	The list of frequencies to notch.
%
%  NotchType        Default is 'notch' which generates a notch filter
%                   using the built-in Matlab filters.  Alternatively,
%                   'FFTNotch' will apply the notch in FFT space and
%                   'SineFit' will fit a sine wave and then subtract it.
%
%  NotchOnlyPeaks   If true (1) then this will only apply notch filter if
%                   there is a peak at the specified frequency.
%
%  SubSampleTime    If empty (default), uses all data. Otherwise uses a
%                   limited section of the data to calculate FFT (assuming
%                   that noise is constant).  The size of this chunk is
%                   specified by the time here (in seconds).
%
%  NotchOrder       these are parameters for creating the notch filter, if
%  NotchBW          using built-in filters.

global AnalysisOpts


%% Parse optional inputs

%Default save filename is based off of the first file name passed
opts.Channel = '';
opts.TimeRange = [];
opts.NotchFrequencies = [60 120 180]; %frequencies to notch out noise
opts.NotchType = 'notch'; %what type of notching to do. Can be 'notch' (using built-in filters), 'FFTNotch' to use fft/ifft approach, or 'SineFit' to fit sinusoid
opts.NotchOnlyPeaks = 0; %only apply notch if there is a peak here in the PSD
opts.SubSampleTime = []; %amount of time (in seconds) to use to estimate PSD.  Also used to fit sine if using SineFit notching
opts.Verbose = 1; %whether to be chatty about what we are doing
opts.SaveData=0; % whethere to save these data
%Notch filter options
opts.NotchOrder = 12; %Order of the filter
opts.NotchBW = 0.1; %Width of the notch

SamplingFrequency=Fs;  % Sampling Frequency of data

%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin),
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end


 
%% Get channel list
if isempty(opts.Channel),
    chan_list = 1:size(data,2);
else
    chan_list = opts.Channel;
end

chan_out_fn = ['NotchFiltData.mat'];

%% Loop through channels, applying the notch in turn

%Loop through channels
for cur_chan = 1:length(chan_list),
     
    
    %Open up input matfile -- shouldn't have to lock as previous steps
    %should all be done in order for me to do this.
    in_data.amp = [data(:,cur_chan)]';
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
    
    %Do we only notch peaks in the PSD?
    if opts.NotchOnlyPeaks,
        if opts.Verbose, fprintf('Limiting notch to peaks in PSD...\n'); end
        
        %Estimate the PSD
        if isempty(opts.SubSampleTime) || (opts.SubSampleTime >= (amp_len/SamplingFrequency)),
            amp_ind = [1:length(amp)];
        else
            amp_ind = floor((amp_len - (SamplingFrequency*opts.SubSampleTime))/2) + [1:floor(SamplingFrequency*opts.SubSampleTime)];
        end
                
        %If we have NaNs, we need to pick the largest part of our time
        %range that doesn't have NaNs
        if any(~isfinite(amp(amp_ind))),
            finite_ind = find(isfinite(amp(amp_ind)));
            finite_range_starts = finite_ind([1 (find(diff(finite_ind) ~= 1) + 1)]);
            finite_range_ends = finite_ind([find(diff(finite_ind) ~= 1) length(finite_ind)]);
            finite_range_size = finite_range_ends - finite_range_starts;
            
            [~, mi] = max(finite_range_size);
            amp_ind = amp_ind(finite_range_starts(mi(1)):finite_range_ends(mi(1)));
        end
        
        %Calculate PSD
        [amp_psd, amp_f] = periodogram(amp(amp_ind), [], [], SamplingFrequency);
        amp_psd = 10*log10(amp_psd); %dB/Hz now..
        
        %Find the peaks in the PSD -- must be at most 1 Hz wide
        ver_info = ver;
        if str2double(ver_info(strcmpi({ver_info.Name}, 'MATLAB')).Version) > 8.1,
            %Use the advanced version of findpeaks
            [~,locs,~,p] = findpeaks(amp_psd, amp_f, 'MaxPeakWidth', 1);
            %Exclude anything that doesn't have a minimum height of 6 dB (roughly 4x) above
            %neighbors
            locs = locs(p >= 6);
        else
            %Use the older version of findpeaks
            [~,locs] = findpeaks(amp_psd, 'THRESHOLD', 6);
            locs = amp_f(locs);
        end
        if opts.Verbose,
            if length(locs) > 100, fprintf('\tFound %d peaks > 6dB (too many to print out).\n', length(locs)); else
                fprintf('\tFound %d peaks > 6dB at several frequencies: %s\n', length(locs), mat2str(locs)); end
        end
        
        %Notch if within 1 Hz of passed frequency
        do_notch = false(length(opts.NotchFrequencies), 1);
        for i = 1:length(opts.NotchFrequencies),
            do_notch(i) = any(abs(opts.NotchFrequencies(i) - locs) <= 1);
        end
        cur_NotchFrequencies = opts.NotchFrequencies(do_notch);
        if opts.Verbose, fprintf('\tLimited to notching when peak within 1 Hz: %s\n', mat2str(cur_NotchFrequencies)); end
        clear('locs', 'p', 'amp_f', 'amp_psd');
    else
        cur_NotchFrequencies = opts.NotchFrequencies;
    end
    
    %Apply notch
    if strcmpi(opts.NotchType, 'notch'),
        if opts.Verbose, fprintf('Applying Butterworth Notch filter...\n'); end
        %Create notch filters for the specified frequencies
        for notch_ind = 1:length(cur_NotchFrequencies),
            cur_notch_freq = cur_NotchFrequencies(notch_ind);
            d = fdesign.notch('N,F0,BW', opts.NotchOrder, cur_notch_freq, opts.NotchBW, SamplingFrequency);
            notchHd(notch_ind) = design(d); %default is a butterworth filter
            
            %See length of impulse response
            impResp = impz(notchHd(notch_ind), [], SamplingFrequency);
            imp_len(notch_ind) = length(impResp);
        end
        %Pad with the maximum length filter
        imp_len = max(imp_len);
        
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
        
        %Apply notch filters
        for notch_ind = 1:length(cur_NotchFrequencies),
            cur_notch_freq = cur_NotchFrequencies(notch_ind);
            if opts.Verbose, fprintf('\tNotching %4.1f Hz...\n', cur_notch_freq); end
            %Apply filter (padding to avoid some edge effects of high-order filters) and save
            amp = filtfilt(notchHd(notch_ind).SOSMatrix, notchHd(notch_ind).ScaleValues, ...
                cat(2, zeros(size(amp, 1), imp_len), amp, zeros(size(amp, 1), imp_len)));
            amp = amp(:, (imp_len + 1):(end-imp_len));
        end
        
        %If needed, we should fill back in the gaps
        if ~isempty(fill_non_finite_gap_starts),
            for i = 1:length(fill_non_finite_gap_starts),
                amp(fill_non_finite_gap_starts(i):fill_non_finite_gap_ends(i)) = NaN;
            end
        end
        
    elseif strcmpi(opts.NotchType, 'FFTNotch'),
        if opts.Verbose, fprintf('Applying FFT Notch filter...\n'); end
        if any(~isfinite(amp)), error('Currently can''t use FFTNotch approach with NaNs in signal.'); end
        %Do FFT based notching
        L = length(amp); %Length of signal
        notch_temp_fft = fft(amp, L); %FFT
        notch_temp_f = SamplingFrequency/2*linspace(0,1,L/2+1); %frequencies
        %Loop through to-be-notched frequencies
        for notch_ind = 1:length(cur_NotchFrequencies),
            cur_notch_freq = cur_NotchFrequencies(notch_ind);
            if opts.Verbose, fprintf('\tNotching %4.1f Hz...\n', cur_notch_freq); end
            [~, mi] = min(abs(notch_temp_f - cur_notch_freq));
            mi = mi(1);
            notch_temp_fft([mi (L-mi+2)]) = 0;
        end
        amp = ifft(notch_temp_fft, L);
    elseif strcmpi(opts.NotchType, 'SineFit'),
        if opts.Verbose, fprintf('Applying Sine-Fit filter...\n'); end
        %Fit sine wave -- based on Womelsdorf, 2006 -- should be
        %identical (but slower?) to above?
        chunk_size = length(amp);
        if ~isempty(opts.SubSampleTime),
            chunk_size = max(0, min(length(amp), round(opts.SubSampleTime*SamplingFrequency)));
        end
        
        for cur_chunk = 1:chunk_size:length(amp),
            cur_chunk_ind = cur_chunk + [0:(chunk_size-1)];
            cur_chunk_ind = cur_chunk_ind((cur_chunk_ind > 0) & (cur_chunk_ind <= length(amp)));
            
            for notch_ind = 1:length(cur_NotchFrequencies),
                cur_notch_freq = cur_NotchFrequencies(notch_ind);
                %Fit our sine wave
                %options = optimset(optimset('fmincon'), 'display', 'none');
                x = fmincon(@(x) nanmean(sqrt((amp(cur_chunk_ind) - x(1)*sin([0:(length(cur_chunk_ind)-1)]*2*pi*cur_notch_freq/SamplingFrequency + x(2))).^2)), [nanmean(abs(amp(cur_chunk_ind))) pi], [], [], [], [], [0 0], [max(abs(amp(cur_chunk_ind))) 2*pi]);%, options);
                %Subtract our fitted wave from our signal
                amp(cur_chunk_ind) = amp(cur_chunk_ind) - x(1)*sin([0:(length(cur_chunk_ind)-1)]*2*pi*cur_notch_freq/SamplingFrequency + x(2));
            end %notch frequency loop
        end %chunk loop
        clear('x', 'notch_temp_fft', 'notch_temp_f', 'mi');
    end %notch type
     
     
    %Was I limited to a specified time range?
    if isempty(opts.TimeRange),
        %Save data to our output file
        if opts.SaveData;save(chan_out_fn, 'amp', 'opts', 'cur_NotchFrequencies', '-v7.3');end
        %Add note saying we spike-band passed
        cur_note = sprintf('Notched channel %s at %s Hz.', chan_list(cur_chan), mat2str(cur_NotchFrequencies));
    else
        %Open up output matfile (initializing it if necessary)
        if opts.SaveData
            if ~exist(chan_out_fn, 'file') ,
                out_data = matfile(chan_out_fn, 'Writable', true);
                out_data.amp = [];
                out_data.opts = opts;
                out_data.cur_NotchFrequencies = cur_NotchFrequencies;
            else
                out_data = matfile(chan_out_fn, 'Writable', true);
            end
        end
        %Remove edges use to avoid edge effects
        amp = amp((1 + overlap_start_offset):(end - overlap_end_offset));
        %Write to file
        if opts.SaveData;out_data.amp(1, opts.TimeRange(1):opts.TimeRange(2)) = amp;end;
        %Add note saying we spike-band passed
        cur_note = sprintf('Notched channel %s, time range %d-%d, at %s Hz.', chan_list(cur_chan), opts.TimeRange, mat2str(cur_NotchFrequencies));
        clear out_data;
    end
    
    Allamp(:,cur_chan)=amp; %% save data for all of channels or trials
     
end %channel loop

 