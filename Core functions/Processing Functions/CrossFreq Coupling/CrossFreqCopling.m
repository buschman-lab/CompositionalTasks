classdef CrossFreqCopling
    %CROSSFREQCOPLING Series of function to look at cross frequency
    %coupling between different frequecies
    % Written by Sina Tafazoli Jan 2021
    properties
        % comodulation analysis parameters
        FpFreqRang=[5 10]; % lower frequency range (Hz)
        FpFreqInc=1; % lower frequncy increments (Hz)
        FpFreqBW=2;    % lower frequency band width (Hz) 
        FaFreqRang=[20 29;30 80]; % higher frequecny range (Hz) 
        FaFreqInc=[1;2]; % lower frequncy increments (Hz)
        FaFreqBW =[2;4];   % lower frequency band width (Hz) 
        
    end
    properties(Access=private)
        FiltFuncs=FilterFuncs; % call filter functions
        ManData=ManipulateData; % obvious
    end
    
    methods
        % starter functions 
        function obj = CrossFreqCopling(varargin)
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
        function MI = CalComodulation(obj,Xraw,Fs,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
          
            %CalComodulation Calculates comodulation index for large number
            %of frequency pairs details of the method are described in 
            %Supporting Information Tort et al. 10.1073/pnas.0810524105
            fprintf('\nCalculating Comodulation index')
            %% (1) First, xraw(t) is filtered at the 2 frequency ranges under analysis (fp and fA).
            % We denote the filtered signals as xfp(t) and xfA(t) (Fig. S2 b and c).
            [CWT,F]=cwt(Xraw,'amor',1000,'VoicesPerOctave',48,...
                           'FrequencyLimits',[0 50] );
            FpBandSet=cell2mat(arrayfun(@(x) [obj.FpFreqRang(x,1):obj.FpFreqInc(x):obj.FpFreqRang(x,2);(obj.FpFreqRang(x,1)+obj.FpFreqBW(x)):obj.FpFreqInc(x):(obj.FpFreqRang(x,2)+obj.FpFreqBW(x))],...
               1:size(obj.FpFreqRang,1),'UniformOutput',0));
            FaBandSet=cell2mat(arrayfun(@(x) [obj.FaFreqRang(x,1):obj.FaFreqInc(x):obj.FaFreqRang(x,2);(obj.FaFreqRang(x,1)+obj.FaFreqBW(x)):obj.FaFreqInc(x):(obj.FaFreqRang(x,2)+obj.FaFreqBW(x))],...
               1:size(obj.FaFreqRang,1),'UniformOutput',0));
            
            % filter data now
            fprintf('\nFiltering phase data ...')
            Xfp=arrayfun(@(x) obj.FiltFuncs.BandPassFilter(Xraw,Fs,'PassBand',[FpBandSet(1,x) FpBandSet(2,x)]),1:size(FpBandSet,2),'UniformOutput',0);
            Xfa=arrayfun(@(x) obj.FiltFuncs.BandPassFilter(Xraw,Fs,'PassBand',[FaBandSet(1,x) FaBandSet(2,x)]),1:size(FaBandSet,2),'UniformOutput',0);
           
            %% (2) The time series of the phases of xfp(t) [denoted as 􏰌fp(t)] is obtained from the standard Hilbert transform of xfp(t) (Fig. S2d).
            %The Hilbert transform is also applied to xfA(t) to extract the time series of the amplitude envelope of xfA(t) [denoted as AfA(t); (Fig. S2e)].
            %The composite time series [􏰌fp(t),AfA(t)] is then constructed, which informs the amplitude of the fA oscillation at each 
            %phase of the fp rhythm (Fig. S2f).
            PHfp=cellfun(@(x) obj.ManData.CalAngle(hilbert(x)),Xfp,'UniformOutput',false); % hilbert transform of Xfp
            Afa =cellfun(@(x) abs(hilbert(x))  ,Xfa,'UniformOutput',false); % hilbert transform of Xfa
           
            FreqBins= obj.ManData.BinData(20,[],0,360);
            nfbin=length(FreqBins)-1;
            % loop through frequency bands and calculate what we want 
            Nfp=length(PHfp);Nfa=length(Afa);
            Hmax=log2(nfbin);
            MI=nan*ones(Nfp,Nfa);
            for i=1:Nfp
                for j=1:Nfa
                    CompositTS=[PHfp{i}' Afa{j}'];
                    %% the phases 􏰌fp(t) are binned into eighteen 20o intervals (0o to 360o), and the mean of AfA over
                    %each phase bin is calculated (Fig. S2g). We denote as 􏰂AfA􏰃􏰌fp (j) the mean AfA value at the phase bin j.
                    MeanAfa=[];
                    for f=1:(nfbin)
                        if f~=(nfbin)
                            PhaseInd=CompositTS(:,1)>=FreqBins(f) & CompositTS(:,1)<FreqBins(f+1);
                        else
                            PhaseInd=CompositTS(:,1)>=FreqBins(f) & CompositTS(:,1)<=FreqBins(f+1);
                         end
                        MeanAfa(f)=mean(CompositTS(PhaseInd,2));                   
                    end
                    %% apply the entropy measure H,
                    SumMeanAfa=sum(MeanAfa);
                    H=0;
                    for f=1:(nfbin)
                        H=H+(MeanAfa(f)/SumMeanAfa)*log2(MeanAfa(f)/SumMeanAfa);
                    end
                    H=-1*H;
                        
                    H2=-1*sum(arrayfun(@(x) (MeanAfa(x)/SumMeanAfa)*log2(MeanAfa(x)/SumMeanAfa),1:nfbin));
                    %% The MI is finally obtained by normalizing H by the maximum possible entropy value (Hmax)
                    MI(i,j)=(Hmax-H)/Hmax;
                end
            end
            figure
           imagesc(mean(FpBandSet,1),mean(FaBandSet,1),MI');
            set(gca,'YDir','normal')
            %  CompositTS=arrayfun(
        end
        function MI = CalComodulationCWT(obj,Xraw,Fs,varargin) % calculates comodulation with CWT
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
          
            %CalComodulation Calculates comodulation index for large number
            %of frequency pairs details of the method are described in 
            %Supporting Information Tort et al. 10.1073/pnas.0810524105
            fprintf('\nCalculating Comodulation index')
            %% (1) First, xraw(t) is filtered at the 2 frequency ranges under analysis (fp and fA).
            % We denote the filtered signals as xfp(t) and xfA(t) (Fig. S2 b and c).
             
            % filter data now
            fprintf('\nFiltering phase data ...')
%             [Xfp,Fp]=cwt(Xraw,'amor',1000,'VoicesPerOctave',20,'FrequencyLimits',[3 20]);
%             [Xfa,Fa]= cwt(Xraw,'amor',1000,'VoicesPerOctave',20,'FrequencyLimits',[20 80]);
%             [Fp,indsortFp]=sort(Fp);
%             Xfp=Xfp(indsortFp,:);
%             [Fa,indsortFa]=sort(Fa);
%             Xfa=Xfa(indsortFa,:);
            
            Fp=obj.ManData.BinData(0.1,[],3,20);
            WAV=ApplyWaveletTransform(Xraw',Fs,'Freqs',Fp,...
                'WaveletWidth',7);
            Xfp=WAV{1}.amp;Fp=WAV{1}.f;
            
            Fa=obj.ManData.BinData(1,[],20,80);
            WAV=ApplyWaveletTransform(Xraw',Fs,'Freqs',Fa,...
                'WaveletWidth',7);
            Xfa=WAV{1}.amp;Fa=WAV{1}.f;

            
            %% (2) The time series of the phases of xfp(t) [denoted as 􏰌fp(t)] is obtained from the standard Hilbert transform of xfp(t) (Fig. S2d).
            %The Hilbert transform is also applied to xfA(t) to extract the time series of the amplitude envelope of xfA(t) [denoted as AfA(t); (Fig. S2e)].
            %The composite time series [􏰌fp(t),AfA(t)] is then constructed, which informs the amplitude of the fA oscillation at each 
            %phase of the fp rhythm (Fig. S2f).
            PHfp=  obj.ManData.CalAngle(Xfp) ; % hilbert transform of Xfp
            Afa =  abs(Xfa)   ; % hilbert transform of Xfa
           
            FreqBins= obj.ManData.BinData(20,[],0,360);
            nfbin=length(FreqBins)-1;
            % loop through frequency bands and calculate what we want 
            Nfp=length(Fp);Nfa=length(Fa);
            Hmax=log2(nfbin);
            MI=nan*ones(Nfp,Nfa);
            for i=1:Nfp
                for j=1:Nfa
                    CompositTS=[PHfp(i,:)' Afa(j,:)'];
                    %% the phases 􏰌fp(t) are binned into eighteen 20o intervals (0o to 360o), and the mean of AfA over
                    %each phase bin is calculated (Fig. S2g). We denote as 􏰂AfA􏰃􏰌fp (j) the mean AfA value at the phase bin j.
                    MeanAfa=[];
                    for f=1:(nfbin)
                        if f~=(nfbin)
                            PhaseInd=CompositTS(:,1)>=FreqBins(f) & CompositTS(:,1)<FreqBins(f+1);
                        else
                            PhaseInd=CompositTS(:,1)>=FreqBins(f) & CompositTS(:,1)<=FreqBins(f+1);
                         end
                        MeanAfa(f)=mean(CompositTS(PhaseInd,2));                   
                    end
                    %% apply the entropy measure H,
                    SumMeanAfa=sum(MeanAfa);
                    H=0;
                    for f=1:(nfbin)
                        H=H+(MeanAfa(f)/SumMeanAfa)*log2(MeanAfa(f)/SumMeanAfa);
                    end
                    H=-1*H;
                        
                    H2=-1*sum(arrayfun(@(x) (MeanAfa(x)/SumMeanAfa)*log2(MeanAfa(x)/SumMeanAfa),1:nfbin));
                    %% The MI is finally obtained by normalizing H by the maximum possible entropy value (Hmax)
                    MI(i,j)=(Hmax-H)/Hmax;
                end
            end
            imagesc( Fp,Fa,MI');
            set(gca,'YDir','normal')
            xlabel('Phase Providing Frequency')
            ylabel('Amplitude Providing Frequency')
            %  CompositTS=arrayfun(
        end

    end
end

