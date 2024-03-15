function [wlt, sf, st] = GenerateWavelets(freqoi, width, gwidth, fsample)

% This is borrowed from the fieldtrip toolbox.  We just don't do padding, etc.

if numel(width) == 1
  width = ones(1,length(freqoi)) * width;
end

wlt = cell(length(freqoi), 1);
for ifreqoi = 1:length(freqoi)
    dt = 1/fsample; %time step
    sf = freqoi(ifreqoi) / width(ifreqoi); %frequency std
    st = 1/(2*pi*sf); %time std
    toi2 = -gwidth*st:dt:gwidth*st; %Our time vector
    
    %Envelope
    A = 1/sqrt(st*sqrt(pi));
    tap = (A*exp(-toi2.^2/(2*st^2)))';
    
    % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
    ind  = (-(size(tap,1)-1)/2 : (size(tap,1)-1)/2)'   .*  ((2.*pi./fsample) .* freqoi(ifreqoi));
    
    % create wavelet
    prezer = zeros(10, 1); pstzer = zeros(10, 1);
    wavelet = complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer));
    %wltspctrm{ifreqoi} = complex(zeros(1,length(wavelet)));
    %wltspctrm{ifreqoi} = fft(wavelet,[],1)';
    wlt{ifreqoi} = wavelet;
end