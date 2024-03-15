function [SIm,SIpRad,SIpDeg,PowerTS_Upper_FFT, faxis] = CalSynchronizationIndex(UpperSig,LowerSig,Fs)
%Calculates synchronization index between two signals 
% refere to methods of Micheal X Cohen 2008 J neuro methods for deta
% Sig 1 and sig 2 have to be filtered at desired frequencies 
Nsamp=length(UpperSig);
% take Hilbert transform of both signal s


Hilb_UpperSig=hilbert(UpperSig);
Hilb_LowerSig=hilbert(LowerSig);
% calculate powertime series of upper signal
PowerTS_Upper=real(Hilb_UpperSig).^2+imag(Hilb_UpperSig).^2;
% take fft of this power time series 
NFFT=2^nextpow2(Nsamp);
[PowerTS_Upper_FFT, faxis] = pmtm(PowerTS_Upper', [], NFFT, Fs);
% detrend power time  
PowerTS_Upper=detrend(PowerTS_Upper);
% get the phases 
Hilbert_PowerTS_Upper=hilbert(PowerTS_Upper);
Hilbert_PowerTS_Upper_Phase=angle(Hilbert_PowerTS_Upper);
Hilb_LowerSig_Phase=angle(Hilb_LowerSig);
% calculate synchronization index
SI=0;
for k=1:Nsamp
    SI=SI+exp(1i*(Hilb_LowerSig_Phase(k)-Hilbert_PowerTS_Upper_Phase(k)));
end
SI=SI/Nsamp;
%
SIm=abs(SI); % magnitude of SI
[SIpRad]=angle(SI); % preferred phase of sinchronizaiton
SIpDeg=NaN;
end

