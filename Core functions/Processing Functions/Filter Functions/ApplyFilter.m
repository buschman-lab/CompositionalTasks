function amp_downSampled=ApplyFilter(amp,Filter,impResp,FilterOpts)
%%% applies specified filter that is desinged by DesignLFPfilter.m
%%% inputs
%%% data: Data strem to be filtered
%%% Filter: filter desgin
%%% output Filtered data(FiltData)
imp_len = length(impResp);

 %Apply our filter
if strcmpi(FilterOpts.FilterType, 'FIR'),
    %Apply filter
    fprintf('Applying FIR (equirriple) %s filter to amplifier signal...\n',FilterOpts.FilterDesign);  
    amp = filtfilt(Filter.Numerator, 1, amp);
elseif strcmpi(FilterOpts.FilterType, 'IIR'),
    %Apply filter
    fprintf('Applying IIR (butterworth) %s filter to amplifier signal...\n',FilterOpts.FilterDesign); 
    %Apply filter (padding to avoid some edge effects of high-order filters) and save
    amp = filtfilt(Filter.SOSMatrix, Filter.ScaleValues, cat(2, zeros(size(amp, 1), imp_len), amp, zeros(size(amp, 1), imp_len)));
    amp = amp(:, (imp_len+1):(end-imp_len));
else
    error('Unknown filter specified.');
end

amp_downSampled = downsample(amp,FilterOpts.SamplingFrequency/FilterOpts.TargetLFPSamplingFrequency);
