classdef TimeFreqAnalysis
    %TIMEFREQANALYSIS Summary of this class goes here
    %   Aggregate of functions to do Time frequency analysis
    
    properties
        %Default save filename is based off of the first file name passed
        FreqToPlot = 0.5:7:100;%AnalysisOpts.LFPParams.FreqToPlot; %frequencies at which to calculate the PSD
        PSDMethod = 'MTM'; %What algorithm to use to estimate PSD.  Can be 'FFT' or 'pwelch' or'MTM'
        TimeFreqMethod='Wavelet';  % which method we are using wavelet or shortterm FFT 'STFFT';
        WaveletMethod='code';  % can be using Matlab or our own code ('MATLAB','Code')
        WaveletWidth          % width of wavelet using code 
        WaveletWidth1=3 % in case of variable wavelet lower
        WaveletWidth2=10% in case of variable wavelet higher
        WaveletTimePaddingWidth=3 % timepadding of wavelet 
        VoicesPerOctave=4 % voices per octave
        navg=10; % how many trials to average   
        NormPSD=1; % should we normalize PSD? 
        ShowSTD=0; % should we show STD  
    end
    properties (Access=protected)
        Fs
        TrialTimesFields
        PeriodLength
        BaselineDelay
        FigParams=fig_params;
        FiltFuncs=FilterFuncs;
    end
 
    methods
        function  obj=ParseParams(obj,InputArgs)
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

        function obj = TimeFreqAnalysis(varargin)
            %TIMEFREQANALYSIS Construct an instance of this class
            %   Detailed explanation goes here           
            global AnalysisOpts
            if isfield(AnalysisOpts,'TrialTimesFields') % 
                obj.TrialTimesFields =AnalysisOpts.TrialTimesFields;
                obj.PeriodLength=AnalysisOpts.TrialTiming.PeriodLength*1000;  %%% Whole period duration
                obj.BaselineDelay=AnalysisOpts.TrialTiming.BaselineDelay*1000; %%% look at x ms before the start
                obj.Fs= AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency;     
            end
            %Process additional optional inputs
            if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
            for i = 1:2:length(varargin),
                try
                    obj.(varargin{i}) = varargin{i+1};
                catch
                    error('Couldn''t set option ''%s''.', varargin{2*i-1});
                end
            end
        end  % initialize all of the values we need
        
        function [cwt_all,cwt_f,TrueTime] = ComputeTimeFreq(obj,data,varargin)
            obj=obj.ParseParams(varargin) ; %%Process optional inputs            
            %% Do all of the calculations of power and FFT here 
            L=size(data,1); Ntrials=size(data,2);
            %%% spectrogram parameters all in ms
            Step=1;  % 1ms of step
            StrTime=obj.BaselineDelay;
            Window=500;noverlap=(Window-Step)*obj.Fs/1000;           
            StpTime=obj.PeriodLength-StrTime;
            Time=-1*StrTime:Step:(StpTime-Window); %this is used for windowed short term FFT
            TrueTime=-1*StrTime:1000/obj.Fs:StpTime-1000/obj.Fs; %% what was atually in the experiment

            for i=1:Ntrials  % loop on differernt trials or channels   

                    %% compute Spectectorgram 
                    if strcmpi(obj.TimeFreqMethod,'Wavelet')  % use wavelet 
                        
                        if strcmpi(obj.WaveletMethod,'MATLAB') % log scale
                       
                           [cwt_s,cwt_f] = cwt(data(:,i),'amor',obj.Fs,'VoicesPerOctave',obj.VoicesPerOctave,...
                           'FrequencyLimits',[obj.FreqToPlot(1) obj.FreqToPlot(end)] );
                       
                        elseif strcmpi(obj.WaveletMethod,'Code') %linear scale
                            WAV=ApplyWaveletTransform(data(:,i),obj.Fs,'Freqs',obj.FreqToPlot,...
                                'WaveletWidth',obj.WaveletWidth);
                            cwt_s=WAV{1}.amp;cwt_f=WAV{1}.f;
                        elseif strcmpi(obj.WaveletMethod,'VariableWidth')  % wavelet with variable width
                            obj.WaveletWidth=obj.WaveletWidth1:(obj.WaveletWidth2-obj.WaveletWidth1)/(length(obj.FreqToPlot)-1):obj.WaveletWidth2;
                            WAV=ApplyWaveletTransform(data(:,i),obj.Fs,'Freqs',obj.FreqToPlot,...
                                'WaveletWidth',obj.WaveletWidth);
                            cwt_s=WAV{1}.amp;cwt_f=WAV{1}.f;
                        end
                       
                        cwt_all(:,:,i)=cwt_s;

                    elseif strcmpi(obj.TimeFreqMethod,'STFFT')    %use STFFT      
                         [s,cwt_f,t]=spectrogram(data(:,i),Window,noverlap,[],obj.Fs);
                         cwt_all(:,:,i)=s;                         
                         %% Compute power spectrum 
                        if strcmp(obj.PSDMethod,'pwelch')
                         [temp_amp_psd,temp_amp_f]=pwelch(data(:,i), [],[],L,obj.Fs,'power'); 
                        elseif strcmp(obj.PSDMethod,'FFT')   
                         [temp_amp_psd,temp_amp_f]=periodogram(data(:,i), [],L,obj.Fs,'power'); 
                        end
                     %Interpolate to specific frequency values
                        if ~isempty(obj.FreqToPlot)
                           temp_amp_psd = interp1(temp_amp_f, temp_amp_psd, obj.FreqToPlot, 'linear', NaN);
                           temp_amp_f = obj.FreqToPlot;
                        end
                        amp_psd(:,i) = temp_amp_psd.*temp_amp_f; % make sure this multiplication is correct
                        amp_f = temp_amp_f;

                    end

            end
        end
        function PlotPSD(obj,PSD,faxis,f_lim,Col,varargin)
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs 
          %  if size(data,1)==1;data=data';end
%             Yfft=fft(data); 
          %    L=size(data,1);
%             %Y(1)=0; %% remove DC
%             %% Plot Amplitude in db
%             P2=abs(Yfft);
%             P1=P2(1:L/2+1,:);
%             P1(2:end-1,:)=2*P1(2:end-1,:);
%             f = [ obj.Fs*(0:(L/2))/L]';  %% frequency vector
%             LFP_F=f>obj.FreqToPlot(1) & f<=obj.FreqToPlot(end);
%             FFT_Amp=P1(LFP_F,:);
%             PlotMeanStd(FFT_Amp,f(LFP_F),'f (Hz)','|P1(f)| db',1)
%             title('Single-Sided Amplitude Spectrum of X(t)')

%             [temp_amp_psd,f]=periodogram(data, [],L,obj.Fs,'power');
%              LFP_F=f>obj.FreqToPlot(1) & f<=obj.FreqToPlot(end);
%              PlotFreqMeanStd(temp_amp_psd(LFP_F).*f(LFP_F),f(LFP_F),'f (Hz)','|P1(f)| db',1)
%              title('PSD of X(t)')
             f_ind=faxis>=f_lim(1) & faxis<=f_lim(2);
             PSD=PSD(:,f_ind);
             faxis=faxis(f_ind);


             if obj.NormPSD % should we normalize these values
                 if size(faxis,1)>1;faxis=faxis';end
                 PSD=PSD.*repmat(faxis,[size(PSD,1),1]); 
             end
             obj.FigParams.Plot( (faxis), nanmean(10*log10(PSD), 1), Col,'f(Hz)','PSD(db)',[]); hold on;
             if obj.ShowSTD
                 obj.FigParams.Plot( (faxis), nanmean(10*log10(PSD), 1) + nanstd(10*log10(PSD), [], 1), Col,[],[],[],'p_line_style',':');
                 obj.FigParams.Plot( (faxis), nanmean(10*log10(PSD), 1) - nanstd(10*log10(PSD), [], 1), Col,[],[],[],'p_line_style',':');  
             end
        end
    
        
        function [Sxx,faxis,df,fNQ]=CalPSD(obj,x,Fs,UseHann,varargin)  % calculates Power spectrum  of signal averages across trials 
               obj=obj.ParseParams(varargin) ;      %%Process optional inputs 
               % x is Trial*Time Matrix        
               % Fs Sampling Frequency
               Nsamp=size(x,2);
               dt=1/Fs; T=Nsamp*dt; 
%            NTrl=size(x,1);
%             if UseHann;Taper=hann(Nsamp)';else;Taper=ones(1,Nsamp);end
%             for k=1:NTrl                
%                % xf = fft(Taper.*x(k,:));      %1. Compute the Fourier transform of x with Hanning window
%                % ThisSxx = (2*dt^2/T )* xf.*conj(xf);   %2. Compute the power spectrum.
%                 [Sxx(k,:),faxis]=pmtm(x(k,:),[],[],Fs);
%              %   Sxx(k,:) = ThisSxx(1:Nsamp/2+1); %3. Ignore negative frequencies.            
%             end           
            obj.PSDMethod='MTM';
               NFFT=2^nextpow2(Nsamp);
               %Apply the desired PSD method
               if strcmpi(obj.PSDMethod, 'FFT')
                   [Sxx, faxis] = periodogram(x', [], NFFT, Fs);
               elseif strcmpi(obj.PSDMethod, 'MTM')
                   %nw=4;
                   [Sxx, faxis] = pmtm(x', [], NFFT, Fs);
               elseif strcmpi(obj.PSDMethod, 'Welch')
                   [Sxx, faxis] = pwelch(x', [], [], NFFT, Fs);
               end
                             
               df = 1/max(T);                   %4. Determine the frequency resolution.
               fNQ=1/dt/2;                      %5. Determine the Nyquist frequency.
              % faxis = (0:df:fNQ);              %6. Construct the frequency axis. 
              % mean_Sxx=mean(Sxx,1);
               Sxx=Sxx'; % take a transpose so each line is a trial still         
        end
        % %%%%*******        % %%%%*******

        function  [cohr,cohrPhase,f]=CalCoherence(obj,x,y,Fs,UseHann,ifPlot)  % calculates coherence 
            
            K = size(x,1);            %Define the number of trials.
            N = size(x,2);            %Define the number of indices per trial.
            dt = 1/Fs;                %Define the sampling interval.
            T = N*dt;                 %Define the duration of data.
            Sxx = zeros(K,N);         %Create variables to save the spectra.
            Syy = zeros(K,N);
            Sxy = zeros(K,N);
            
            if UseHann;Taper=hann(N)';else;Taper=ones(1,N);end
            
            for k=1:K
                FFTX=fft(Taper.*x(k,:));
                FFTY=fft(Taper.*y(k,:));
                %Compute the spectra for each trial.
                Sxx(k,:) = (2*dt^2/T) *FFTX.* conj(FFTX);
                Syy(k,:) = (2*dt^2/T) *FFTY.* conj(FFTY);
                Sxy(k,:) = (2*dt^2/T) *FFTX.* conj(FFTY);
            end
            
            Sxx = Sxx(:,1:N/2+1);            %Ignore negative frequencies.
            Syy = Syy(:,1:N/2+1);
            Sxy = Sxy(:,1:N/2+1);
         
            Sxx = mean(Sxx,1);            %Average the spectra across trials.
            Syy = mean(Syy,1);
            Sxy = mean(Sxy,1);
            
            cohr = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));           %Compute the coherence.
            cohrPhase=angle(Sxy);         %Compute the coherence.
          
            df = 1/max(T);                 %Determine the frequency resolution.
            fNQ = 1/ dt / 2;               %Determine the Nyquist frequency.
            f = (0:df:fNQ);            %Construct frequency axis.          
            if ifPlot
               
                subplot(121)
                plot(f, real(cohr));       %Plot the results
                ylim([0 1]) ;                  %Set the axes limits
                xlabel('Frequency [Hz]')        %Label axes.
                ylabel('Coherence [ ]')
                subplot(122)
                plot(f, cohrPhase);            %Plot the results
                xlabel('Frequency [Hz]')            %Label axes.
                ylabel('Coherence Phase')
            end
        end  
        function  [cohr,cohrPhase]=CalWaveCoherence(obj,x,y,Fs,ifPlot)  % calculates wavelet coherence 
            
            F=size(x,1);
            K = size(x,3);            %Define the number of trials.
            N = size(x,2);            %Define the number of indices per trial.
            
            Sxx = zeros(F,N,K);         %Create variables to save the spectra.
            Syy = zeros(F,N,K);
            Sxy = zeros(F,N,K);
            
            for k=1:K
                %Compute the wavelet spectra for each trial.
                Sxx(:,:,k) =  x(:,:,k) .* conj(x(:,:,k));
                Syy(:,:,k) =  y(:,:,k) .* conj(y(:,:,k));
                Sxy(:,:,k) =  x(:,:,k) .* conj(y(:,:,k));
            end
                                
            Sxx = mean(Sxx,3);            %Average the wavlet spectra across trials.
            Syy = mean(Syy,3);
            Sxy = mean(Sxy,3);
            
            cohr = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));           %Compute the coherence.
            cohrPhase=angle(Sxy);         %Compute the coherence.
          
             
           
        end  
        function [PowXf,f]=CalTFRhilbert(obj,Xraw,Fs,FreqRang,FreqInc,FreqBW,varargin) % calculates time frequecny representation with hilbert transform
           % FreqRang, frequecny range 
           %FreqInc, frequency increments
           %FreqBW, frequency band witth
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
            if size(Xraw,1)>1;Xraw=Xraw';end
            % determine frequency bands
            BandSet=cell2mat(arrayfun(@(x) [FreqRang(x,1):FreqInc(x):FreqRang(x,2);(FreqRang(x,1)+FreqBW(x)):FreqInc(x):(FreqRang(x,2)+FreqBW(x))],...
               1:size(FreqRang,1),'UniformOutput',0));
             Nfreq=length(BandSet);
            % filter data now
             Xf=arrayfun(@(x) obj.FiltFuncs.BandPassFilter(Xraw,Fs,'PassBand',[BandSet(1,x) BandSet(2,x)]),1:size(BandSet,2),'UniformOutput',0);
            normXf=arrayfun(@(x) (Xf{x}-mean(Xf{x}))/std(Xf{x}),1:Nfreq,'UniformOutput',0);
           % normXf=arrayfun(@(x) Xf{x}*(BandSet(1,x)+BandSet(2,x))/2,1:Nfreq,'UniformOutput',0);
            
            PowXf=cell2mat(arrayfun(@(x) (abs(hilbert(normXf{x})).^2)',1:Nfreq,'UniformOutput',0));
            PowXf=PowXf';
            f=mean(BandSet,1);
        end


    end
end

