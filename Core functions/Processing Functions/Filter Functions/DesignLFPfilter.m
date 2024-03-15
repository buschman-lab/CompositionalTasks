function [Hd,impResp]=DesignLFPfilter( FilterOpts)

%% DESIGN FILTER

%SamplingFrequency = 200;
%FilterOrder = 4;
% PassBand(1) = 40;
% PassBand(2) = 80;

 %Design filter
if isempty(FilterOpts.FilterOrder),
    %Use the full specification
    d=fdesign.lowpass('Fp,Fst,Ap,Ast', FilterOpts.PassBand, FilterOpts.StopBand, FilterOpts.PassBandRipple, FilterOpts.StopBandAttenuation, ...
        FilterOpts.SamplingFrequency);
else
    if strcmp(FilterOpts.FilterDesign,'lowpass')
        d=fdesign.lowpass('N,F3dB', FilterOpts.FilterOrder, FilterOpts.PassBand, FilterOpts.SamplingFrequency);
    elseif strcmp(FilterOpts.FilterDesign,'bandpass')
        d=fdesign.bandpass('N,F3dB1,F3dB2', FilterOpts.FilterOrder, FilterOpts.PassBand(1), FilterOpts.PassBand(2), FilterOpts.SamplingFrequency);
    end
end

if strcmpi(FilterOpts.FilterType, 'FIR'),
    Hd = design(d, 'equiripple');
elseif strcmpi(FilterOpts.FilterType, 'IIR'),
    if isempty(FilterOpts.FilterOrder),
        Hd = design(d, 'butter', 'MatchExactly', 'passband');
    else
        Hd = design(d, 'butter');
    end
end
% d=fdesign.bandpass('N,F3dB1,F3dB2', FilterFilterOpts.FilterOrder, FilterFilterOpts.PassBand(1), FilterFilterOpts.PassBand(2), FilterFilterOpts.SamplingFrequency);
% Hd = design(d, 'butter');
impResp = impz(Hd, [], FilterOpts.SamplingFrequency); %impulse response 
