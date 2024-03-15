function AX=helperCWTTimeFreqPlot(cfs,time,freq,PlotType,varargin)
%   This function helperCWTTimeFreqPlot is only in support of
%   CWTTimeFrequencyExample and PhysiologicSignalAnalysisExample. 
%   It may change in a future release.
    FigParams=fig_params;
    ClrMap=parula(128);%magma(128);   %parula(128)
    if size(freq,1)==1;freq=freq';end
    Nf=length(freq);
    NNeuorns=size(cfs,1)/Nf; 
    if mod(NNeuorns,1)==0
    freq=repmat(freq,NNeuorns,1);  
    else
        freq=1:size(cfs,1);
    end

    params = parseinputs(varargin{:});
    if params.uselog==1
        fprintf('Using log frequency for X axis\n')
        %%% Normalize Power
        if size(freq,1)==1;freq=freq';end
      %  cfs=cfs.*repmat(freq,1,size(cfs,2));   % normalize power by 1/f
        %%%  take log of  freq
        FreqChunk=sort(freq(1:Nf));
        LogFreq=log2(FreqChunk);
        YTick    =repmat(LogFreq(1:1:length(LogFreq)),NNeuorns,1);
        YTicksLbl=repmat(FreqChunk(1:1:length(FreqChunk)),NNeuorns,1);
        freq=log2(freq);
        %%% change label
        if regexpi(params.ylab,'freq')
            params.ylab=['log2' params.ylab];
        end
    end
    

    if strcmpi(PlotType,'Power')
        args = {time,freq,abs(cfs).^2};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(ClrMap);
        h = colorbar;
        h.Label.String = 'Power';
            if isempty(params.xlab) && isempty(params.ylab)
                xlabel('Time'); ylabel('Hz');
            else
             xlabel(params.xlab); ylabel(params.ylab);
            end
    elseif strcmpi(PlotType,'Corrlation')
        args = {time,freq,cfs};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(ClrMap);
        h = colorbar;
        h.Label.String = 'Corr Coefficient';
            if isempty(params.xlab) && isempty(params.ylab)
                xlabel('Time'); ylabel('Hz');
            else
             xlabel(params.xlab); ylabel(params.ylab);
            end
            
      elseif strcmpi(PlotType,'Phase')
        args = {time,freq,angle(cfs)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(ClrMap);
        h = colorbar;
        h.Label.String = 'Phase(rad)';
            if isempty(params.xlab) && isempty(params.ylab)
                xlabel('Time'); ylabel('Radian');
            else
             xlabel(params.xlab); ylabel(params.ylab);
            end
                  
      elseif strcmpi(PlotType,'pval')
        args = {time,freq,cfs};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(ClrMap);
        h = colorbar;
        h.Label.String = 'Phase(rad)';
            if isempty(params.xlab) && isempty(params.ylab)
                xlabel('Time'); ylabel('Radian');
            else
             xlabel(params.xlab); ylabel(params.ylab);
            end
            
      elseif strcmpi(PlotType,'PowerPhase')
        args = {time,freq,abs(cfs).^2,angle(cfs)};
        surf(args{:},'edgecolor','none');
        view(0,90);
        axis tight;
        shading interp; colormap(ClrMap);
        h = colorbar;
        h.Label.String = 'Phase(rad)';
            if isempty(params.xlab) && isempty(params.ylab)
                xlabel('Time'); ylabel('Radian');
            else
             xlabel(params.xlab); ylabel(params.ylab);
            end
            
       elseif strcmpi(PlotType,'Justplot1')
            args = {time,freq,cfs};
            AX=surf(args{:},'edgecolor','none');
            view(0,90);
            axis tight;
%             axis equal;
          %  shading interp;
          colormap(ClrMap);
           % h = colorbar;
           % h.Label.String = 'Power';
                if isempty(params.xlab) && isempty(params.ylab)
                    xlabel('Time'); ylabel('Hz');
                else
                 xlabel(params.xlab); ylabel(params.ylab);
                end   
      elseif strcmpi(PlotType,'Justplot2')
            args = {time,freq,cfs};
            surf(args{:},'edgecolor','none');
%             view(0,90);
            axis tight;
%             axis equal;
            shading interp; colormap(ClrMap);
           % h = colorbar;
           % h.Label.String = 'Power';
                if isempty(params.xlab) && isempty(params.ylab)
                    xlabel('Time'); ylabel('Hz');
                else
                 xlabel(params.xlab); ylabel(params.ylab);
                end            
    elseif strcmpi(PlotType,'contour')
        contour(time,freq,abs(cfs).^2);
        grid on; colormap(ClrMap); 
        h = colorbar;
        h.Label.String = 'Power';
            if isempty(params.xlab) && isempty(params.ylab)
                xlabel('Time'); ylabel('Hz');
            else
             xlabel(params.xlab); ylabel(params.ylab);
            end
            
    elseif strcmpi(PlotType,'contourf')
        contourf(time,freq,abs(cfs).^2);
        grid on; colormap(ClrMap); 
        h = colorbar;
        h.Label.String = 'Power';
            if isempty(params.xlab) && isempty(params.ylab)
                xlabel('Time'); ylabel('Hz');
            else
             xlabel(params.xlab); ylabel(params.ylab);
            end
            
    elseif strcmpi(PlotType,'image')
        imagesc(time,freq,cfs);
        colormap(ClrMap); 
        AX = gca;
        shading interp;
        set(gca,'YDir','normal')
      %  h = colorbar;
      %  h.Label.String = 'Power';
            if isempty(params.xlab) && isempty(params.ylab)
                xlabel('Time(s)'); ylabel('Freq(Hz)');
            else
                xlabel(params.xlab); ylabel(params.ylab);
            end
    end
    
    if ~isempty(params.PlotTitle)
        title(params.PlotTitle);
    end
    
    if params.uselog==1
        yticks(YTick);
        yticklabels(arrayfun(@(X) num2str(YTicksLbl,3),1:length(YTicksLbl),'UniformOutput',0));
       
    end
     set(gcf,'renderer','painters')
    FigParams.FormatAxes(gca);
    %----------------------------------------------------------------
    function params = parseinputs(varargin)
        
        params.PlotTitle = [];
        params.xlab = [];
        params.ylab = [];
        params.threshold = -Inf;
        
    
        if isempty(varargin)
            return;
        end
        
        Len = length(varargin);
        if (Len==1)
            params.PlotTitle = varargin{1};
        end
    
        if (Len == 3)
            params.PlotTitle = varargin{1};
            params.xlab = varargin{2};
            params.ylab = varargin{3};
        end
           
        if (Len == 4)
            params.PlotTitle = varargin{1};
            params.xlab = varargin{2};
            params.ylab = varargin{3};
            params.uselog=varargin{4};
        end
        
        
  
 
        