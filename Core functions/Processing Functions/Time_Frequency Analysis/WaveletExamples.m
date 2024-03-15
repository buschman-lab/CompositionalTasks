% shows example of wavaelet transform with different parameters
load ExampleData.mat

TimeFreq=TimeFreqAnalysis;
ManData=ManipulateData;
FigParams=fig_params;
opts.WaveletMethod='code';
opts.q1=3;opts.q2=10;
opts.wws={[3],[7],'fixes std time',sprintf('variable width(%i-%i)',opts.q1,opts.q2),'matlab'};
opts.wwsTitle={'3','7','fixes std time',sprintf('variable width(%i-%i)',opts.q1,opts.q2),'matlab'};
opts.FrqStp=[1];
opts.st=0.1; % time std of 50ms 
opts.sf=1/(2*pi*opts.st); % freq std
NRaw=length(opts.FrqStp)+2;Ncol=length(opts.wws);


%wts=[1 2 5 10];
FsWave=100;
N=10*FsWave; % number of samples to look at
f1=1;f2=100;
%f_lin=sort(unique([0.25 0.5 0.75 ManData.BinData(opts.FrqStp,[],f1,f2)])); % define frequency axis

Fig=FigParams.RenderFigure(1,[]);
for frstp=opts.FrqStp
    f_lin=sort(unique([0.25 0.5 0.75 ManData.BinData(frstp,[],f1,f2)])); % define frequency axis
    
    for ww=1:length(opts.wws) %wavelet width
        %  for wt=wts %wavelet time padding
        if strcmpi(opts.wws{ww},'fixes std time')
            q=f_lin/opts.sf;
        elseif strcmpi(opts.wws{ww},sprintf('variable width(%i-%i)',opts.q1,opts.q2))
            q=opts.q1:(opts.q2-opts.q1)/(length(f_lin)-1):opts.q2;
        elseif strcmpi(opts.wws{ww},'matlab')
            opts.WaveletMethod='matlab';
            q=nan;
        else
            q=opts.wws{ww};
        end
        
        [Wavelet_Linear,f_Linear,~] = TimeFreq.ComputeTimeFreq(data,'Fs',FsLFP,...
            'WaveletMethod',opts.WaveletMethod,'FreqToPlot',f_lin,'WaveletWidth',q,'VoicesPerOctave',8);
        
        FWHM = 1*sqrt(2*log(2))*f_Linear./q; %frequency std
        
        
        CWTPower=ManData.NormPower(Wavelet_Linear(:,1:N),f_Linear);
  %      CWTPower=smoothdata(CWTPower,2);
        
        Wavelet_Linear=transpose(downsample(Wavelet_Linear',FsLFP/FsWave)); % downsample to 100hz
        subplot(NRaw,Ncol,ww+(find(frstp==opts.FrqStp)-1)*Ncol)
        helperCWTTimeFreqPlot(CWTPower,(0:N-1)/FsWave,f_Linear,'justplot1',['q Wave:' opts.wwsTitle{ww}, 'Freq Stp:' num2str(frstp)],'Time(s)','f',0)
       % axis square
        
        subplot(NRaw,Ncol,ww+(find(frstp==opts.FrqStp))*Ncol)
        helperCWTTimeFreqPlot(CWTPower,(0:N-1)/FsWave,f_Linear,'justplot1',['q Wave:' opts.wwsTitle{ww}, 'Freq Stp:' num2str(frstp)],'Time(s)','f',0)
     %   axis square
        set(gca,'YScale','log')
         
        subplot(NRaw,Ncol,ww+(find(frstp==opts.FrqStp)+1)*Ncol)
        FigParams.Plot(f_Linear,FWHM,ww,'Freq(Hz)','FWHM(Hz)',opts.wwsTitle{ww});
        axis square

        %     %% show cross correlation of signal itself
        %     WavCrosCor=xcorr2(CWTPower,CWTPower);
        %     subplot(NRaw,Ncol,find(ww==opts.wws)+Ncol)
        %     helperCWTTimeFreqPlot(WavCrosCor,(0:size(WavCrosCor,2)-1)/FsWave,1:size(WavCrosCor,1),'justplot1',['Autocorr'],'Time(ms)','f',0)
        %
    end
end