
classdef FilterFuncs
    %FILTERFUNCS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FilterOrder=4;
    end
    
    methods
        function obj = FilterFuncs(varargin)
            if nargin~=0 % initialize vars
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
        end
        function obj=ParseParams(obj,InputArgs)
            %Process optional inputs
            if mod(length(InputArgs), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(InputArgs)
                try
                    obj.(InputArgs{i}) = InputArgs{i+1};
                catch
                    error('Couldn''t set option ''%s''.', InputArgs{2*i-1});
                end
            end
        end
        function [output,BandSet]=BandPassFilterOld(obj,data,Fs,BandSet,varargin)
            
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
            
            if size(BandSet,2)~=1;BandSet=BandSet';end
            if size(data,1)~=1;data=data';end
            BandSet=sort(BandSet);
            BandSet=[BandSet(1:end-1) BandSet(2:end)]; %create set of bands
            
            for FrBand=1:size(BandSet,1)   %% loop on frequency bands
                PassBand(1) = BandSet(FrBand,1);
                PassBand(2) = BandSet(FrBand,2);
                fprintf('\n Band pass filtering data in band: %s to %s Hz',num2str(PassBand(1)),num2str(PassBand(2)));
                d=fdesign.bandpass('N,F3dB1,F3dB2',obj.FilterOrder,PassBand(1),PassBand(2), Fs);
                Hd = design(d, 'butter');
                impResp = impz(Hd, [], Fs);
                imp_len = length(impResp); % Get the length of the impulse response -- we'll use this to avoid edge effects
                %Apply filter
                filtdata = filtfilt(Hd.SOSMatrix, Hd.ScaleValues, cat(2, zeros(size(data, 1), imp_len), data, zeros(size(data, 1), imp_len)));
                output(FrBand,:) = filtdata(:, (imp_len+1):(end-imp_len));
            end
            
            
        end
        
        function amp=BandPassFilter(obj,amp,Fs, varargin)
            global AnalysisOpts
            % Isolates the spike band from the signal by filtering out the frqeuencies
            % with action potential waveform energy.
            %
            % Inputs:
            %
            % Optional Inputs (must be specified in key-value pairs):
            %  Channel          which channel to process. Will only filter the specified
            %                   channel.
            %
            %
            %  PassBand         The edges of the pass-band. Default is [300 3000] Hz.
            %
            %  StopBand         The frequencies at which the stop band is in full effect
            %                   (determines roll-off). Default is [250 4000] Hz.
            %
            %  StopBandAtten    The attentuation of the signal at and above StopBand
            %                   (in dB). Default is 80 dB
            %
            %  PassBandRipple   Allowed degree of ripple in passband (in dB). Default
            %                   is 0.01 dB.
            %
            %  FilterDesign     What form of filter to use. Must be valid for built in
            %                   designer from Matlab filter toolbox. Default is
            %                   'butter'.
            
            
            %% Parse optional inputs
            
            %Default save filename is based off of the first file name passed
            opts.Channel = '';
            opts.TimeRange = [];
            opts.FilterOrder = 4; %recommended to use a higher order for elliptic
            opts.PassBand = [1 10];%[300 3000]
            opts.StopBand = [0.5 11];%[250 4000]
            opts.StopBandAtten = 40; %in dB
            opts.PassBandRipple = 0.1; %in dB
            opts.FilterDesign = 'butter'; %what type of filter to use, can be butter, equiripple or elliptic
            opts.Verbose = 1; %whether to be chatty about what we are doing
            opts.ShowFilter=0; % do we want to visualize filter with fvtool?
            
            %Process optional inputs
            if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(varargin),
                try
                    opts.(varargin{i}) = varargin{i+1};
                catch
                    error('Couldn''t set option ''%s''.', varargin{2*i-1});
                end
            end
            
            if opts.PassBand(1) >= opts.PassBand(2),
                error('PassBand must be increasing.');
            end
            
            if isempty(opts.FilterOrder),
                if opts.StopBand(1) >= opts.StopBand(2),
                    error('StopBand must be increasing.');
                end
                if opts.StopBand(1) >= opts.PassBand(1),
                    error('Can''t have the lower stop band frequency higher than the pass band frequency (this is a low-pass filter).');
                end
                if opts.StopBand(2) <= opts.PassBand(2),
                    error('Can''t have the upper stop band frequency below the pass band frequency (this is a low-pass filter).');
                end
            end
            
            if size(amp,1)~=1;amp=amp';end
            %Design filter for this signal
            impResp = [];
            if any(strcmpi(opts.FilterDesign, {'butter', 'equiripple'})),
                if isempty(opts.FilterOrder),
                    %Use the full specification
                    d=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
                        opts.StopBand(1), opts.PassBand(1), opts.PassBand(2), opts.StopBand(2), ...
                        opts.StopBandAtten, opts.PassBandRipple, opts.StopBandAtten, ...
                        Fs);
                else
                    d=fdesign.bandpass('N,F3dB1,F3dB2', opts.FilterOrder, opts.PassBand(1), opts.PassBand(2), Fs);
                end
                Hd = design(d, opts.FilterDesign);
                impResp = impz(Hd, [], Fs);
            elseif strcmpi(opts.FilterDesign, 'elliptic'),
                d = designfilt('bandpassiir', 'FilterOrder', opts.FilterOrder, ...
                    'PassbandFrequency1', opts.PassBand(1), 'PassbandFrequency2', opts.PassBand(2), ...
                    'PassbandRipple', opts.PassBandRipple, 'StopbandAttenuation1', opts.StopBandAtten, 'StopbandAttenuation2', opts.StopBandAtten, ...
                    'SampleRate', Fs);
                impResp = impz(d, Fs);
            else
                error('Unknown filter design.');
            end
            % show filter design 
            if opts.ShowFilter
                fvtool(Hd)
            end
            %Get the length of the impulse response -- we'll use this to avoid edge
            %effects
            if isempty(impResp),
                error('Couldn''t get filter impulse response.');
            end
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
            amp_len = length(amp);
            if strcmpi(opts.FilterDesign, 'butter'),
                %Apply filter
                if opts.Verbose, fprintf('Applying butterworth filter to signal... %iHz-%iHz\n', opts.PassBand(1), opts.PassBand(2)); end
                amp = filtfilt(Hd.SOSMatrix, Hd.ScaleValues, cat(2, zeros(size(amp, 1), imp_len), amp, zeros(size(amp, 1), imp_len)));
                amp = amp(:, (imp_len+1):(end-imp_len));
            elseif strcmpi(opts.FilterDesign, 'equiripple'),
                %Apply filter
                if opts.Verbose, fprintf('Applying equiripple filter to signal... %iHz-%iHz\n', opts.PassBand(1), opts.PassBand(2)); end
                amp = filtfilt(Hd.Numerator, 1, cat(2, zeros(size(amp, 1), imp_len), amp, zeros(size(amp, 1), imp_len)));
                amp = amp(:, (imp_len+1):(end-imp_len));
            elseif strcmpi(opts.FilterDesign, 'elliptic'),
                amp = filtfilt(d, cat(2, zeros(size(amp, 1), imp_len), amp, zeros(size(amp, 1), imp_len)));
                amp = amp(:, (imp_len+1):(end-imp_len));
            else
                error('Unknown filter specified.');
            end
           
        end                     
    end
end

