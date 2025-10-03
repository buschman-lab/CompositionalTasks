classdef AnalysisFuncs
    %MOTIFFUNCS perfroms different operations of Motif or data to find
    %motifs
    % Detailed explanation goes here
    
    properties
        %% motif functions
        L
        K
        lambdaL1W
        lambdaL1H
        lambda
        lambdaOrthoH
        lambdaOrthoW
        maxiter
        W_fixed=0;
        W_init
        useWupdate=1;
        %% general vars
        ChNum
        DownSampleFactor;
        SavePlot
        SaveData
        ShowPlot
        verbose=1; 
        ExtraString='';  % extra string in plot or data name
        BlockSpec
        TrialTimes
        PowerMethod      % how are we normalizing the power can be Power, NormPower, NormMedianPower
        PEVmethod='var'         % how are we calculating PEV 'var','pls'
        Fs;       % sampling freq of    data
        LogScale=0   % if plot TF in logscale
        ClassTrainFrac=0.70;  % precetange of training for classifier
        Nrep=10;   % classifer repetition
        MotifClassficationMethod='EncoderPhenograph';  % how to classify motifs 'EncoderPhenograph' or 'phenohgraph'
        ClassificaitonAutoEncLayerSiz=[30]; % hidden layer size of autoencoder
        ClassificaitonAutoEncTargLayer=2;       % Target layer index
        CoreMotifsPerc=0.1; % percentage of motifs for calculted for Core motif
        OverWrite=0; % are we overwriting esxting data?
        FigTitle=''; % title of this figure
        XAxisLabel=''; % label of x axis
        YAxisLabel=''; % label of y axis
        NewFig=0; % creates new figure
        AnalysisPathName=''; % path of the analysis we want to save data
        %% shuffling and data managmenet
        Nshuffle=50; % 50 repetition for shuffling
        Navg=10; % how many samples to average
        UseDataPointer=0; % use data pointer for data retrival
        DataPointerVar=''; % variable to use for data pointer retrival
        %% Time Freq Analysis
        UseHann=1; % use hanning window to calculate FFT
        % trial timing
        StartFieldName={'FIX_ON','TARGET_ON','GIVE_REWARD'}; % what is the zero for time alignment 
        %%
    end
    properties (Access=private)
        ManData=ManipulateData;
        TrialFunc=TrialFuncs;
        FilterFunc=FilterFuncs;
        FigParams=fig_params;
        NeuAnaFunc=NeuralAnalysisFuncs;
        TimFreq=TimeFreqAnalysis;
        CrosFreqCop=CrossFreqCopling;
        freq
    end
    
    methods
        function  obj=AnalysisFuncs(varargin)
            if nargin~=0 % initialize vars
                obj=obj.ParseParams(varargin) ; %%Process optional inputs
            end
        end
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
        function  out=ConcatinateWs(obj,W)
            k=1;
            for i=1:length(W)
                w=W{i};
                for j=1:size(w,2)
                    temp=squeeze(w(:,j,:));
                    out(:,k)=temp(:);
                    k=k+1;
                end
            end
        end
        function  PlotW(obj,W,f,T,varargin)
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %Plots W
            SumWs=sum(sum(W,1),3);
            Nws=find(SumWs~=0);
            for i=1:length(Nws)
                subplot(1,length(Nws),i)
                helperCWTTimeFreqPlot(squeeze(W(:,Nws(i),:)) ,T,f,'justplot1',['W' num2str(i)],'Time','Frequency',obj.LogScale);
                axis square
            end
        end      % plots W
        function  varargout=PlotW_H_RealImag(obj,W,H,f,T,NTrials,varargin) % only for real and imag analysis
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % NTrials is the number of Trials for each Rule
            % W and H are W and H mats seqNMF.
            %f and T are requency and Time
            %%Plots W s first
            TsWave=1/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            Nws=size(W,2);
            NCh= size(W,1)/length(f);
            Ncol=length(NTrials)+1;if isempty(NTrials);Ncol=4;end
            TimH=[1:size(H,2)]*TsWave;
            for i=1:Nws
                %subplot(NSubps(1),NSubps(2),i)
                if sum(sum(squeeze(W(:,i,:))))
                    varargout(i)=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]);
                    % plot example period first
                    
                    titleRealImag={'Real','Imag'};
                    for s=1:NCh
                        
                        Ind=length(f)*(s-1)+1:length(f)*(s);
                        subplot(NCh,Ncol,Ncol*(s-1)+1)
                        ExampleReal=AnalysisData.X(Ind,1:100);
                        helperCWTTimeFreqPlot(ExampleReal ,(1:size(ExampleReal,2))*TsWave,f,'justplot1',['Example TFR' titleRealImag{s}],'Time(s)','Frequency',obj.LogScale);
                        axis square
                        axis tight
                        colorbar
                        subplot(NCh,Ncol,2+Ncol*(s-1))
                        helperCWTTimeFreqPlot(squeeze(W(Ind,i,:)) ,(1:size(W,3))*TsWave,f,'justplot1',['SpatioTemp W ' titleRealImag{s}],'Time(s)','Frequency',obj.LogScale);
                        axis square
                        axis tight
                        colorbar
                        if s==1;RealPrt=squeeze(W(Ind,i,:));
                        elseif s==2;ImagPrt=squeeze(W(Ind,i,:));
                        end
                    end
                    ReconstructW=RealPrt+1i*ImagPrt;
                    AbsPhaseW{1}=abs(ReconstructW);
                    AbsPhaseW{2}=angle(ReconstructW);
                    Pow_W=obj.ManData.NormPower(ReconstructW,f);   % calculat the power of the data
                    titleReconst={'Abs','Phase'};
                    for e=1:2
                        subplot(NCh,Ncol,3+Ncol*(e-1))
                        helperCWTTimeFreqPlot(AbsPhaseW{e},(1:size(W,3))*TsWave,f,'justplot1',['SpatioTemp W ' titleReconst{e} ],'Time(s)','Frequency',obj.LogScale);
                        axis square
                        axis tight
                        colorbar
                    end
                    subplot(NCh,Ncol,4)
                    helperCWTTimeFreqPlot(Pow_W,(1:size(W,3))*TsWave,f,'justplot1',['SpatioTemp W norm Power*freq' ],'Time(s)','Frequency',obj.LogScale);
                    axis square
                    axis tight
                    colorbar
                else
                    varargout{i}=[];
                end
            end
            %%Plot corresponding Hs now
            if ~isempty(NTrials)
                H=obj.ManData.ReshapeTrials(H,NTrials);
                for Rule=1:length(NTrials) % loop on rules
                    for nW=1:Nws   % loop on W s
                        if sum(sum(squeeze(W(:,nW,:))))
                            figure(varargout{nW})
                            %subplot(2,2,Rule+1)
                            subplot(NCh,Ncol,Rule+1)
                            h=obj.ManData.SmoothData(H{nW,Rule});
                            %h=H{nW,Rule};
                            helperCWTTimeFreqPlot(h ,T,TrialOrder,...
                                'justplot1',['H for W' num2str(nW) ' Rule:' num2str(Rule)],...
                                'Time(s)','Trial',obj.LogScale)
                            axis square
                        end
                    end
                end
            else
                for nW=1:Nws   % loop on W s
                    if sum(sum(squeeze(W(:,nW,:))))
                        figure(varargout{nW})
                        subplot(NCh,Ncol,8)
                        h=H(nW,:);
                        bar(TimH,h,0.5);
                        xlabel('Time(s)')
                        ylabel('H')
                        axis tight
                        axis square
                    end
                end
            end
            
            %  varargout=arrayfun(@(x) varargout(x),1:length(varargout),'UniformOutput',0)  ;
        end   % plots W and H
        function  varargout=PlotW_H(obj,W,H,f,T,NTrials,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % NTrials is the number of Trials for each Rule
            % W and H are W and H mats seqNMF.
            %f and T are requency and Time
            %%Plots W s first
            TsWave=1/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            Nws=size(W,2);
            NCh= size(W,1)/length(f);
            Ncol=length(NTrials)+1;if isempty(NTrials);Ncol=2;end
            TimH=[1:size(H,2)]*TsWave;
            for i=1:Nws
                %subplot(NSubps(1),NSubps(2),i)
                if sum(sum(squeeze(W(:,i,:))))
                    varargout(i)=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]);
                    for s=1:NCh
                        subplot(NCh,Ncol,1+Ncol*(s-1))
                        Ind=length(f)*(s-1)+1:length(f)*(s);
                        helperCWTTimeFreqPlot(squeeze(W(Ind,i,:)) ,(1:size(W,3))*TsWave,f,'justplot1',['SpatioTemp W' num2str(i) ' CH' num2str(obj.ChNum)],'Time(s)','Frequency',obj.LogScale);
                        axis square
                        axis tight
                    end
                else
                    varargout{i}=[];
                end
            end
            %%Plot corresponding Hs now
            if ~isempty(NTrials)
                H=obj.ManData.ReshapeTrials(H,NTrials);
                for Rule=1:length(NTrials) % loop on rules
                    for nW=1:Nws   % loop on W s
                        if sum(sum(squeeze(W(:,nW,:))))
                            figure(varargout{nW})
                            %subplot(2,2,Rule+1)
                            subplot(NCh,Ncol,Rule+1)
                            h=obj.ManData.SmoothData(H{nW,Rule});
                            %h=H{nW,Rule};
                            helperCWTTimeFreqPlot(h ,T,TrialOrder,...
                                'justplot1',['H for W' num2str(nW) ' Rule:' num2str(Rule)],...
                                'Time(s)','Trial',obj.LogScale)
                            axis square
                        end
                    end
                end
            else
                for nW=1:Nws   % loop on W s
                    if sum(sum(squeeze(W(:,nW,:))))
                        figure(varargout{nW})
                        subplot(NCh,Ncol,2)
                        h=H(nW,:);
                        bar(TimH,h,0.5);
                        xlabel('Time(s)')
                        ylabel('H')
                        axis tight
                        axis square
                    end
                end
                
            end
            
            %  varargout=arrayfun(@(x) varargout(x),1:length(varargout),'UniformOutput',0)  ;
        end   % plots W and H
        
        function  varargout=PlotW_H_Trial(obj,W,H,f,T,NTrials,TrialOrder,ExpBlkNum,varargin)   % plots W_H with trials
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % NTrials is the number of Trials for each Rule
            % W and H are W and H mats seqNMF.
            % H is 3 dimensional matrix and each block's data is in the 3rd
            % dim
            %f and T are requency and Time
            %%Plots W s first
            TsWave=1/AnalysisOpts.NeuralAnalysis.FsWaveTarg;
            Nws=size(W,2);
            NCh= size(H,1);
            Ncol=length(NTrials)+3;if isempty(NTrials);Ncol=2;end
            ChNums=AnalysisData.Ch;
            for i=1:Nws
                %subplot(NSubps(1),NSubps(2),i)
                if sum(sum(squeeze(W(:,i,:))))
                    varargout{i}=figure;
                    set(gcf,'Units','Normalized','Position',[0 0 0.8 0.8])
                    for s=1:NCh
                        subplot(NCh,Ncol,1+Ncol*(s-1))
                        helperCWTTimeFreqPlot(squeeze(W(:,i,:)) ,(1:size(W,3))*TsWave,f,'justplot1',['SpatioTemp W' num2str(i) ' CH' num2str(ChNums(s))],'Time(s)','Frequency',obj.LogScale);
                        %     axis square
                        axis tight
                    end
                else
                    varargout{i}=[];
                end
            end
            %%Plot corresponding Hs now
            for nW=1:Nws   % loop on W s
                for s=1:NCh
                    for Blk=1:length(NTrials) % loop on rules
                        ThisH=squeeze(H{s,Blk}(nW,:,:));
                        
                        if sum(sum(squeeze(W(:,nW,:))))
                            figure(varargout{nW})
                            subplot(NCh,Ncol,Blk+1+(Ncol)*(s-1))
                            %  h=obj.ManData.SmoothData(H{nW,Rule});
                            h=ThisH';
                            helperCWTTimeFreqPlot(h ,T,TrialOrder,...
                                'image',['Block' num2str(ExpBlkNum(Blk)) ],...
                                'Time(s)','Trial',obj.LogScale)
                            axis square
                        end
                        if Blk==length(NTrials)
                            % in the last one just plot the average over time and over
                            % trials
                            subplot(NCh,Ncol,length(NTrials)+2+(Ncol)*(s-1))
                            PlotMeanStd(T, h,[],'Time','Mean H value','b',1,'');
                            axis square
                            
                            subplot(NCh,Ncol,length(NTrials)+3++(Ncol)*(s-1))
                            PlotMeanStd(TrialOrder, h',[],'Trial','Mean H value','r',1,'');
                            axis square
                        end
                    end
                end
            end
            
            
            
            
            
            %  varargout=arrayfun(@(x) varargout(x),1:length(varargout),'UniformOutput',0)  ;
        end   % plots W and H
        function  varargout=PlotW_H_Trial_Area(obj,W,H,f,T,NTrials,TrialOrder,ExpBlkNum,ChArea,varargin)   % plots W_H with trials for each Area
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %% first get the data for each area
            %% HRaw Organization 1xBlk cell array of (motifsxTime*Trial*Channel)
            HRaw=obj.ManData.ReshapeCell2Mat(H,42); % reshape all channels into a big matrix
            H=[];
            Nblks=length(ExpBlkNum);
            AreaInds=[1 4 5];
            
            for i=AreaInds % go over areas and take average for each block
                for b=1:Nblks
                    Ind=find(i==AreaInds);
                    ThisAreaNeu=ChArea==i;
                    H{Ind,b}=mean(HRaw{b}(:,:,:,ThisAreaNeu),4);
                end
            end
            
            
            %              obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %             FieldName='NumRewards';
            %             CondValSet=0:6;
            %
            %             Nblks=length(BlockSpec.ThisBlkTrialTimes);
            %             BlkSet=[6 9 11 13];%[1:Nblks];%[1 2 12];%[6 9 11 13];
            %             Nch=size(H,1);
            %             CondVal=CondValSet;
            %             Ncond=length(CondVal);
            %             CondInd=strcmpi(AnalysisOpts.TrialTimesFields,FieldName);
            %             TrialOrder=BlockSpec.TrialOrder;
            %             CondCol=colormap(jet(Ncond));
            %             % loop on blocks and cut the data for each condition
            %             for Ch=1:Nch
            %                 for cond=1:Ncond
            %                     tempmatAll=[];
            %                     for b=1:length(BlkSet)
            %                         tempmat=nan*ones(size(H{1,1}));
            %                         ThisCondTrls=BlockSpec.ThisBlkTrialTimes{BlkSet(b)}(:,CondInd)==CondVal(cond);
            %                         tempmat(:,:,ThisCondTrls)=H{Ch,BlkSet(b)}(:,:,ThisCondTrls);
            %                         tempmatAll=cat(4,tempmatAll,tempmat);
            %                     end
            %                     H_cond{Ch,cond}=nanmean(tempmatAll,4);
            %                 end
            %             end
            %
            %             %% now get average of each condition for each area
            %             AreaInds=[1 4 5];
            %             H_resh_cond=obj.ManData.ReshapeCell2Mat(H_cond,42); % reshape all channels into a big matrix
            %             H_area_cond=[]; % inistail H for each area
            %             for i=AreaInds % go over areas and take average for each block
            %                 for cond=1:Ncond
            %                     Ind=find(i==AreaInds);
            %                     ThisAreaNeu=ChArea==i;
            %                     H_area_cond{Ind,cond}=mean(H_resh_cond{cond}(:,:,:,ThisAreaNeu),4);
            %                 end
            %             end
            %             H=H_area_cond;
            % NTrials is the number of Trials for each Rule
            % W and H are W and H mats seqNMF.
            % H is 3 dimensional matrix and each block's data is in the 3rd
            % dim
            %f and T are requency and Time
            %%Plots W s first
            TsWave=1/AnalysisOpts.NeuralAnalysis.FsWaveTarg;
            Nws=size(W,2);
            NCh= size(H,1);
            Ncol=length(NTrials)+3;if isempty(NTrials);Ncol=2;end
            ChNums=AnalysisData.Ch;
            for i=1:Nws
                %subplot(NSubps(1),NSubps(2),i)
                if sum(sum(squeeze(W(:,i,:))))
                    varargout{i}=figure;
                    set(gcf,'Units','Normalized','Position',[0 0 0.8 0.8])
                    for s=1:NCh
                        subplot(NCh,Ncol,1+Ncol*(s-1))
                        helperCWTTimeFreqPlot(squeeze(W(:,i,:)) ,(1:size(W,3))*TsWave,f,'justplot1',['SpatioTemp W' num2str(i) ' ' AnalysisOpts.AreaNames{AreaInds(s)}],'Time(s)','Frequency',obj.LogScale);
                        %     axis square
                        axis tight
                    end
                else
                    varargout{i}=[];
                end
            end
            %%Plot corresponding Hs now
            for nW=1:Nws   % loop on W s
                for s=1:NCh
                    for Blk=1:length(NTrials) % loop on rules
                        ThisH=squeeze(H{s,Blk}(nW,:,:));
                        
                        if sum(sum(squeeze(W(:,nW,:))))
                            figure(varargout{nW})
                            subplot(NCh,Ncol,Blk+1+(Ncol)*(s-1))
                            %  h=obj.ManData.SmoothData(H{nW,Rule});
                            h=ThisH';
                            helperCWTTimeFreqPlot(h ,T,TrialOrder,...
                                'image',['Block' num2str(ExpBlkNum(Blk)) ],...
                                'Time(s)','Trial',obj.LogScale)
                            axis square
                        end
                        if Blk==length(NTrials)
                            % in the last one just plot the average over time and over
                            % trials
                            subplot(NCh,Ncol,length(NTrials)+2+(Ncol)*(s-1))
                            PlotMeanStd(T, h,[],'Time','Mean H value','b',1,'');
                            axis square
                            
                            subplot(NCh,Ncol,length(NTrials)+3++(Ncol)*(s-1))
                            PlotMeanStd(TrialOrder, h',[],'Trial','Mean H value','r',1,'');
                            axis square
                        end
                    end
                end
            end
            
            %  varargout=arrayfun(@(x) varargout(x),1:length(varargout),'UniformOutput',0)  ;
        end
        function  varargout=PlotW_H_Trial_Condition(obj,W,H,f,T,BlockSpec,ChArea,FieldName,CondValSet,varargin)   % plots W_H with trials for each defined condition such as reward etc etc
            %  FieldName: which field to look at in trial times
            %  FieldValSet: what are the valaues for this field we are
            %  intrested in
            
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            FieldName='NumRewards';
            CondValSet=0:6;
            
            Nblks=length(BlockSpec.ThisBlkTrialTimes);
            BlkSet=[6 9 11 13];%[1:Nblks];%[1 2 12];%[6 9 11 13];
            Nch=size(H,1);
            CondVal=CondValSet;
            Ncond=length(CondVal);
            CondInd=strcmpi(AnalysisOpts.TrialTimesFields,FieldName);
            TrialOrder=BlockSpec.TrialOrder;
            CondCol=colormap(jet(Ncond));
            % loop on blocks and cut the data for each condition
            for Ch=1:Nch
                for cond=1:Ncond
                    tempmatAll=[];
                    for b=1:length(BlkSet)
                        tempmat=nan*ones(size(H{1,1}));
                        ThisCondTrls=BlockSpec.ThisBlkTrialTimes{BlkSet(b)}(:,CondInd)==CondVal(cond);
                        tempmat(:,:,ThisCondTrls)=H{Ch,BlkSet(b)}(:,:,ThisCondTrls);
                        tempmatAll=cat(4,tempmatAll,tempmat);
                    end
                    H_cond{Ch,cond}=nanmean(tempmatAll,4);
                end
            end
            
            %% now get average of each condition for each area
            AreaInds=[1 4 5];
            H_resh_cond=obj.ManData.ReshapeCell2Mat(H_cond,42); % reshape all channels into a big matrix
            H_area_cond=[]; % inistail H for each area
            for i=AreaInds % go over areas and take average for each block
                for cond=1:Ncond
                    Ind=find(i==AreaInds);
                    ThisAreaNeu=ChArea==i;
                    H_area_cond{Ind,cond}=mean(H_resh_cond{cond}(:,:,:,ThisAreaNeu),4);
                end
            end
            H=H_area_cond;
            % NTrials is the number of Trials for each Rule
            % W and H are W and H mats seqNMF.
            % H is 3 dimensional matrix and each block's data is in the 3rd
            % dim
            %f and T are requency and Time
            %%Plots W s first
            TsWave=1/AnalysisOpts.NeuralAnalysis.FsWaveTarg;
            Nws=size(W,2);
            NCh= size(H,1);
            Ncol=3;if isempty(Ncond);Ncol=2;end
            ChNums=AnalysisData.Ch;
            for i=1:Nws
                %subplot(NSubps(1),NSubps(2),i)
                if sum(sum(squeeze(W(:,i,:))))
                    varargout{i}=figure;
                    set(gcf,'Units','Normalized','Position',[0 0 0.8 0.8])
                    for s=1:NCh
                        subplot(NCh,Ncol,1+Ncol*(s-1))
                        helperCWTTimeFreqPlot(squeeze(W(:,i,:)) ,(1:size(W,3))*TsWave,f,'justplot1',['SpatioTemp W' num2str(i) ' ' AnalysisOpts.AreaNames{AreaInds(s)}],'Time(s)','Frequency',obj.LogScale);
                        %     axis square
                        axis tight
                    end
                else
                    varargout{i}=[];
                end
            end
            %%Plot corresponding Hs now
            for nW=1:Nws   % loop on W s
                for s=1:NCh
                    for Blk=1:Ncond
                        ThisH=squeeze(H{s,Blk}(nW,:,:));
                        
                        if sum(sum(squeeze(W(:,nW,:))))
                            figure(varargout{nW})
                            %     subplot(NCh,Ncol,Blk+1+(Ncol)*(s-1))
                            %  h=obj.ManData.SmoothData(H{nW,Rule});
                            h=ThisH';
                            %                          helperCWTTimeFreqPlot(h ,T,TrialOrder,...
                            %                         'image',['Block' num2str(ExpBlkNum(Blk)) ],...
                            %                         'Time(s)','Trial',obj.LogScale)
                            %                         axis square
                        end
                        
                        %  if Blk==Ncond
                        % in the last one just plot the average over time and over
                        % trials
                        subplot(NCh,Ncol,2+(Ncol)*(s-1));hold on
                        PlotMeanStd(T, h,[],'Time','Mean H value',CondCol(Blk,:),3,'');
                        axis square
                        if Blk==1 || Blk==Ncond
                            subplot(NCh,Ncol,3+(Ncol)*(s-1));hold on
                            PlotMeanStd(TrialOrder, h',[],'Trial','Mean H value',CondCol(Blk,:),3,'');
                            axis square
                        end
                        %  end
                    end
                end
            end
            
            %  varargout=arrayfun(@(x) varargout(x),1:length(varargout),'UniformOutput',0)  ;
        end
        function  varargout=PlotW_H_EachW_Condition(obj,W,H,f,T,TrialTimes,ChArea,FieldName,CondValSet,varargin) % plots H of all channls for each W across conditions there is no trial info available
            %  FieldName: which field to look at in trial times
            %  FieldValSet: what are the valaues for this field we are
            %  intrested in
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            %  FieldName='NumRewards';
            %  CondValSet=0:6;
            
            Nch=size(H,1);
            Nws=size(W,2);
            CondVal=CondValSet;
            Ncond=length(CondVal);
            CondInd=strcmpi(AnalysisOpts.TrialTimesFields,FieldName);
            if contains(FieldName,'Color');CircConSpac=1;else;CircConSpac=0;end
            ColorMap=[];
            % loop on blocks and cut the data for each condition
            for Ch=1:Nch
                for cond=1:Ncond
                    ThisCondTrls=TrialTimes(:,CondInd)==CondVal(cond);
                    H_cond{Ch,cond}=H{Ch,1}(:,:,ThisCondTrls);
                end
            end
            % now plot this for each channel
            ChFigs=cell(1,Nws);
            [ChFigs{:}]=obj.PlotAvgHCondition(W,H_cond,f,T,[],ChArea,ColorMap,CondValSet,CircConSpac);
            
            %% now get average of each condition for each area
            AreaInds=[1 4 5];
            H_resh_cond=obj.ManData.ReshapeCell2Mat(H_cond,42); % reshape all channels into a big matrix
            H_area_cond=[]; % inistail H for each area
            for i=AreaInds % go over areas and take average for each block
                for cond=1:Ncond
                    Ind=find(i==AreaInds);
                    ThisAreaNeu=ChArea==i;
                    H_area_cond{Ind,cond}=mean(H_resh_cond{cond}(:,:,:,ThisAreaNeu),4);
                end
            end
            % plot the areas now
            AreaFigs=cell(1,Nws);
            [AreaFigs{:}]=obj.PlotAvgHCondition(W,H_area_cond,f,T,1:length(AreaInds),AreaInds,ColorMap,CondValSet,CircConSpac);
            varargout=[ChFigs AreaFigs];
        end
        function  varargout=PlotAvgHCondition(obj,W,H,f,T,ChNum,ChArea,ColorMap,CondValSet,CircConSpac,varargin) % plots average H per condition
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
            % NTrials is the number of Trials for each Rule
            % W and H are W and H mats seqNMF.
            % H is 3 dimensional matrix and each block's data is in the 3rd
            % dim
            %f and T are requency and Time
            %%Plots W s first
            % ChName is the name of the channel or the area
            % CircConSpac =1 if the condition space is circular
            TsWave=1/AnalysisOpts.NeuralAnalysis.FsWaveTarg;
            Nws=size(W,2);
            NCh= size(H,1);
            NCond=size(H,2);
            Nsups=2+NCh;
            if Nsups<=5;Nrow=1;else;Nrow=5;end
            Ncol=ceil(Nsups/Nrow);
            
            % adjust the vals now
            if isempty(ColorMap);ColorMap=colormap(jet(NCond));end
            if isempty(CondValSet);CondValSet=1:NCond;end
            if isempty(T);T=AnalysisData.Time;end
            if isempty(ChNum);ChNum=AnalysisData.Ch;end
            % plots the W on the first subplot, plot the color of each
            % condition, average of condition
            % for each channel in the subsequent subplots
            varargout=obj.FigParams.RenderFigure(Nws,[]);
            for i=1:Nws
                figure(varargout{i});
                subplot(Nrow,Ncol,1);
                helperCWTTimeFreqPlot(squeeze(W(:,i,:)) ,(1:size(W,3))*TsWave,f,'justplot1',['Motif ' num2str(i)],'Time(s)','Freq(Hz)',obj.LogScale);
                %    axis square
                axis tight
                % plot colors of values
                subplot(Nrow,Ncol,2);hold on
                if CircConSpac
                    arrayfun(@(x) plot(cos(CondValSet(x)),sin(CondValSet(x)),'.','color',ColorMap(x,:),'MarkerSize',20),1:NCond);
                else
                    arrayfun(@(x) plot(1,CondValSet(x),'.','color',ColorMap(x,:),'MarkerSize',20),1:NCond);
                end
                title('Condition Space')
            end
            %%Plot corresponding Hs now
            for nW=1:Nws   % loop on W s
                figure(varargout{nW})
                for ch=1:NCh  % loop on Chs
                    subplot(Nrow,Ncol,2+ch);hold on
                    for cond=1:NCond
                        ThisH=squeeze(H{ch,cond}(nW,:,:));
                        h=ThisH';
                        PlotMeanStd(T, h,[],'Time','Mean H',ColorMap(cond,:),3,['Ch ' num2str(ChNum(ch)) ' Ar ' AnalysisOpts.AreaNames{ChArea((ch))}]);
                        %axis square
                        axis tight
                    end
                end
            end
        end
        function  varargout=PlotRawData(obj,cwt,f,Time,RawData,TrialTimes,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
            % make sure verything has the same size
            NSamp=length(Time);
            cwt=cwt(:,1:NSamp);RawData=RawData(1:NSamp);
            %calculate power of data
            if strcmpi(obj.PowerMethod,'NormPower')
                eval(['Pow=obj.ManData.' obj.PowerMethod '(cwt,f);']);   % calculat the power of the data
            else
                eval(['Pow=obj.ManData.' obj.PowerMethod '(cwt);']);   % calculat the power of the data
            end
            NRaw=2;Ncol=1;
            varargout(1)=obj.FigParams.RenderFigure(1,[0 0 0.95 0.95]);;
            h{1}=subplot(NRaw,Ncol,1);
            % plot TF plot first
            helperCWTTimeFreqPlot(Pow ,Time,f,'justplot1',['CWT CH' num2str(obj.ChNum)],...
                'Time(s)','Frequency',obj.LogScale);
            % plot the raw data now
            h{2}=subplot(NRaw,Ncol,2);
            plot(Time,RawData,'k','LineWidth',1)
            xlabel('Time(s)');ylabel('volts')
            axis tight
            if 0 % do we want to plot individual frequencies 
                % take the trail times and stamp fix,stim and reward time
                if ~isempty(TrialTimes)
                    cellfun(@(x) obj.PlotTrialTimes(x,TrialTimes,[],Time),h,'UniformOutput',0);
                end
                % show filtered data in a seperate plot now
                % [FiltRaw,BandSet]=obj.FilterFunc.BandPassFilter(RawData,obj.Fs,f); % filter raw data
                if size(f,1)==1;f=f';end
                BandSet=f;
                BandSet=sort(BandSet);
                BandSet=[BandSet(1:end-1) BandSet(2:end)]; %create set of bands
                Nfreq=size(BandSet,1);
                Col=distinguishable_colors(Nfreq);
                varargout(2)=obj.FigParams.RenderFigure(1,[0 0 0.95 0.95]);
                for i=1:Nfreq  % loop of frequencies
                    filtplt{i}=subplot(Nfreq,1,i);
                    plot(Time,abs(cwt(Nfreq-i+1,:)),'LineWidth',1,'Color',Col(Nfreq-i+1,:))
                    if i~=Nfreq;axis off;end
                    ylabel(num2str(BandSet(Nfreq-i+1,1),2))
                    xlabel('Time(sec)')
                    axis tight
                    AxisProp=get(gca);
                    %   set(AxisProp,'xcolor','none')
                    AxisProp.YAxis.Label.Color=[0 0 0];
                    AxisProp.YAxis.Label.Visible='on';
                    AxisProp.YAxis.Label.FontSize=5;
                end
                cellfun(@(x) obj.PlotTrialTimes(x,TrialTimes,[],Time),h,'UniformOutput',0);
            end
        end
        function  PlotTrialTimes(~,FigHndl,TrialTimes,StrTrl,TimPriod) % plots the timing of trials on each axis
            if isempty(TrialTimes);return;end
            global AnalysisOpts
            FixInd=strcmp(AnalysisOpts.TrialTimesFields,'FIXATE_ACQUIRED');
            SampleInd=strcmp(AnalysisOpts.TrialTimesFields,'SAMPLE_ON');
            %% subtract the trial we want to start from
            if ~isempty(StrTrl)
                IndStart= strcmp(AnalysisOpts.TrialTimesFields,'START_TRIAL');
                IndEnd=strcmp(AnalysisOpts.TrialTimesFields,'END_TRIAL');
                TrialTimes(:,1:5)=TrialTimes(:,1:5)-TrialTimes(StrTrl,IndStart); % get the timing only
            end
            %% only take these times
            if ~isempty(TimPriod)
                WantedTrls=find((TrialTimes(:,FixInd)>=TimPriod(1) & TrialTimes(:,FixInd)<=TimPriod(end)) | (TrialTimes(:,SampleInd)>=TimPriod(1) & TrialTimes(:,SampleInd)<=TimPriod(end)));
                TrialTimes=TrialTimes(WantedTrls,:);
            end
            %%
            NTrl=size(TrialTimes,1);
            AxisLimits=axis(FigHndl);
            
            %   RewardInd=strcmp(AnalysisOpts.TrialTimesFields,'REWARD');
            % getaxis properties
            YLim=AxisLimits(3:4);
            % now plot them
            axes(FigHndl);hold on
            arrayfun(@(x) plot([TrialTimes(x,FixInd) TrialTimes(x,FixInd)],YLim,'r','linewidth',1),1:NTrl,'UniformOutput',0);
            arrayfun(@(x) plot([TrialTimes(x,SampleInd) TrialTimes(x,SampleInd)],YLim,'g','linewidth',1),1:NTrl,'UniformOutput',0);
            %  arrayfun(@(x) plot([TrialTimes(x,RewardInd) TrialTimes(x,RewardInd)],YLim,'b','linewidth',2),1:NTrl,'UniformOutput',0);
            axis tight
            if ~isempty(TimPriod);xlim([ TimPriod(1) TimPriod(end)]);end
        end
        function  Output=PrepapreCWT(obj,data,varargin) % prepares CWT for motif analysis
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
            %% concatinate trials into time
            data=obj.ManData.Cell2Mat(data);   % concatinate data in time
            if strcmpi(obj.PowerMethod,'NormPower')
                if isempty(obj.freq);obj.freq=AnalysisData.cwt_f;end
                eval(['Pow_data=obj.ManData.' obj.PowerMethod '(data,obj.freq);']);   % calculat the power of the data
            elseif strcmpi(obj.PowerMethod,'Power') || strcmpi(obj.PowerMethod,'NormMedianPower')
                eval(['Pow_data=obj.ManData.' obj.PowerMethod '(data);']);   % calculat the power of the data
            elseif strcmpi(obj.PowerMethod,'Angle')
                Pow_data=angle(data);
            elseif strcmpi(obj.PowerMethod,'none')
                Pow_data=double(data);
            elseif strcmpi(obj.PowerMethod,'AbsAngle')% concatinate abs angle together
                AbsData=abs(data);AbsData=AbsData/max(AbsData(:));
                AngleData=angle(data)+pi;AngleData=(AngleData)/(2*pi);
                Pow_data=[AbsData;AngleData];
            elseif strcmpi(obj.PowerMethod,'RealImag') % concatinate real and imaginary together
                RealData=real(data);RealData=obj.ManData.NormalizeMinMax(RealData);
                ImagData=imag(data);ImagData=obj.ManData.NormalizeMinMax(ImagData);
                
                Pow_data=[RealData;ImagData];
            end
            Output=transpose(downsample(Pow_data',obj.DownSampleFactor)); % downsample data we don't need 1 ms resolution
            
        end
        function  [Output,varargout]=FindMotifs(obj,data,varargin)
            fprintf('\nFindMotifs: Discovering Motifs with seqNMF ...')
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ;      %%Process optional inputs
            
            X=obj.PrepapreCWT(data); % prepare CWT data for processing
            %% Fit with seqNMF
            if AnalysisOpts.ShowPlot;varargout{1}=figure;else, varargout{1}=[];end
            if ~isempty(obj.W_init)
                [Output.W,Output.H, Output.Cost,Output.loadings,Output.power]= seqNMF(X,'K',obj.K,'L',obj.L,...
                    'lambda', obj.lambda, 'maxiter', obj.maxiter, 'showPlot', obj.ShowPlot,...
                    'tolerance',0,'shift',0,'lambdaL1W',obj.lambdaL1W,...
                    'lambdaL1H',obj.lambdaL1H,'SortFactors',1,'W_fixed',obj.W_fixed,'W_init',obj.W_init,...
                    'lambdaOrthoH',obj.lambdaOrthoH,'lambdaOrthoW',obj.lambdaOrthoW);
            else
                [Output.W,Output.H, Output.Cost,Output.loadings,Output.power]= seqNMF(X,'K',obj.K,'L',obj.L,...
                    'lambda', obj.lambda, 'maxiter', obj.maxiter, 'showPlot', obj.ShowPlot,...
                    'tolerance',0,'shift',0,'lambdaL1W',obj.lambdaL1W,...
                    'lambdaL1H',obj.lambdaL1H,'SortFactors',1,'W_fixed',obj.W_fixed,...
                    'lambdaOrthoH',obj.lambdaOrthoH,'lambdaOrthoW',obj.lambdaOrthoW);
            end
            
            AnalysisData.X=X;
            [Output.PEVind,Output.PEV]=obj.CalMotifPEV(X,Output.W,Output.H,'PEVmethod','var');
            %         [Output.cost,Output.regularization,~] = helper.get_seqNMF_cost(X,Output.W,Output.H);
            %             if obj.SaveData
            %                save([AnalysisOpts.DataSavePath obj.ExtraString 'MotifFit_Ch' num2str(obj.ChNum)  '.mat'])
            %             end
        end  % finds motifs
        function  [out,varargout]=PlotMotifs(obj,cwt_all,cwt_f,NTrials,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % find motifs and plot them
            [out,MotifFindFig{1}]=obj.FindMotifs(cwt_all,'freq',cwt_f);
            
            Time=-1:.01:(0.5-0.01); %% what was atually in the experiment
            % if we are working on the whole data then get H during
            % trials
            if ~isempty(NTrials) % trial based analysis
                H=cellfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(...
                    out.H ,obj.TrialTimes(x,:),[],'WaveletDownSampleFactor',1,...
                    'Fs',AnalysisOpts.LFPParams.FilterOpts.TargetLFPSamplingFrequency/obj.DownSampleFactor,...
                    'AnalysisType','Wavelet_Trial_H'),obj.BlockSpec.ThisBlkTrials,'UniformOutput',0);
                out.H=obj.ManData.Cell2Mat(H);
            end
            
            if obj.ShowPlot & ~isempty(NTrials)
                Wfigs=cell(1,obj.K);
                [Wfigs{:}]=obj.PlotW_H(out.W,out.H,cwt_f,Time,NTrials) ; % plot Ws
                varargout=[MotifFindFig Wfigs];
                % save an example part of the signal too
                figure  % plot first 5 Trials
                helperCWTTimeFreqPlot(AnalysisData.X(:,1:750),1:750,cwt_f,'justplot1','PSD','Time','Frequency(Hz)',obj.LogScale)
            elseif obj.ShowPlot & isempty(NTrials) % then we dont have trial info
                Wfigs=cell(1,obj.K);
                if strcmp(AnalysisOpts.MotifAnalysis.PowerMethod,'RealImag')
                    [Wfigs{:}]=obj.PlotW_H_RealImag(out.W,out.H,cwt_f,Time,[]); % plot Ws
                else
                    [Wfigs{:}]=obj.PlotW_H(out.W,out.H,cwt_f,Time,[]) ; % plot Ws
                end
                varargout=[MotifFindFig Wfigs];
            else, varargout=cell(1,obj.K+1);
            end
            %             if obj.SavePlot
            %                 saveas(gcf,[AnalysisOpts.ResultsSavePath 'ExmpSpectrogram_Ch' num2str(obj.ChNum)  '.png'],'png')
            %                 arrayfun(@(x) saveas(varargout{x},[AnalysisOpts.ResultsSavePath obj.ExtraString 'ConvNMF_W' num2str(x) 'Ch' num2str(obj.ChNum) 'K' num2str(obj.K) ...
            %                  'L' num2str(obj.L) '.png'],'png'),2:size(out.W,2)+1)
            %             end
            
        end
        function  Output=ExploreMotifParams(obj,data,varargin)
            % explores different parameters of seqNMF algorithm
            %% define params
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            nLambdas = 5; % increase if you're patient
            Ks= 10;
            numfits = 3; %number of fits to compare
            lambdas =[0 sort(logspace(-1,-5,nLambdas), 'ascend')];
            Ls=10:10:100;   % how long is the L considering 100 hz sampling rate
            %% setup Data
            data=obj.ManData.Cell2Mat(data);
            Pow_data=obj.PrepapreCWT(data); % prepare CWT data for processing
            X=transpose(downsample(Pow_data',obj.DownSampleFactor)); % downsample data we don't need 1 ms resolution
            
            %% fit with seqNMF
            for l=1:length(Ls)  % loop on L
                LL=Ls(l);
                display(['Testing L ' num2str(l) '/' num2str(length(Ls))])
                for li = 1:length(lambdas)
                    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
                    for k = 1:length(Ks)
                        for ii = 1:numfits
                            fprintf('running seqNMF with K = %i try %i\n',Ks(k),ii)
                            
                            [Output.W{ii,k,li,l}, Output.H{ii,k,li,l}, ~,Output.loadings{ii,k,li,l},...
                                Output.power{ii,k,li,l}]= seqNMF(X,'K',Ks(k),'L',LL,...
                                'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0);
                            
                            [Output.cost(ii,k,li,l),Output.regularization(ii,k,li,l),~] = helper.get_seqNMF_cost(X,Output.W{ii,k,li,l},Output.H{ii,k,li,l});
                            Output.PEV{ii,k,li,l}=obj.CalMotifPEV(X,Output.W{ii,k,li,l}, Output.H{ii,k,li,l});
                        end
                        inds = nchoosek(1:numfits,2);
                        for i = 1:size(inds,1) % consider using parfor for larger numfits
                            Output.Diss(i,k,li,l) = helper.DISSX(Output.H{inds(i,1),k,li,l},Output.W{inds(i,1),k,li,l},Output.H{inds(i,2),k,li,l},Output.W{inds(i,2),k,li,l});
                        end
                    end
                end
            end
            Output.lambdas=lambdas;
            Output.Ks=Ks;
            %  Output.X=X;
            Output.Ls=Ls;
            if AnalysisOpts.SaveData
                save([AnalysisOpts.DataSavePath obj.ExtraString 'MotifParams' num2str(obj.ChNum) '.mat'],'Output','obj','-v7.3')
            end
        end
        function  varargout=PlotDiss(obj,data)
            %% Plot Diss and choose K with the minimum average diss.
            global AnalysisOpts
            Nlambdas=length(data.lambdas);
            Ks=data.Ks;
            
            
            varargout{1}=figure; set(gcf,'Units','Normalized','Position',[0 0 1 1])
            for li=1:Nlambdas
                subplot(1,Nlambdas,li)
                plot(Ks,data.Diss(:,Ks,li),'ko'), hold on
                h1 = plot(Ks,median(data.Diss(:,Ks,li),1),'k-','linewidth',2);
                h2 = plot([3,3],[0,0.5],'r--');
                if li==1;legend([h1 h2], {'median Diss','est K'},'Location','southoutside');end
                xlabel('K')
                ylabel('Diss')
                title(['Lambda' num2str(data.lambdas(li))])
                axis square
            end
            if obj.SavePlot
                saveas(gcf,[AnalysisOpts.ResultsSavePath 'Diss_Ch' num2str(obj.ChNum)  '.png'],'png')
            end
        end
        function  varargout=PlotLambdaCost(obj,data)
            
            %% define some vars
            global AnalysisOpts
            lambdas           =data.lambdas;
            All_regularization=data.regularization;
            All_cost          =data.cost;
            Ks                =data.Ks;
            Ls                =data.Ls;
            %% plot costs as a function of lambda
            varargout{1}=figure; set(gcf,'Units','Normalized','Position',[0 0 0.9 0.9])
            for l=1:length(Ls)
                for k=1:length(Ks)
                    
                    %  subplot(1,length(Ks),k-Ks(1)+1)
                    regularization=transpose(squeeze(mean(All_regularization(:,k,:,l),1)));
                    cost=transpose(squeeze(mean(All_cost(:,k,:,l),1)));
                    windowSize = 2;
                    b = (1/windowSize)*ones(1,windowSize);
                    a = 1;
                    Rs = filtfilt(b,a,regularization);
                    minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
                    Rs = (Rs-minRs)/(maxRs-minRs);
                    R = (regularization-minRs)/(maxRs-minRs);
                    Cs = filtfilt(b,a,cost);
                    minCs =  prctile(cost,10); maxCs =  prctile(cost,90);
                    Cs = (Cs -minCs)/(maxCs-minCs);
                    C = (cost -minCs)/(maxCs-minCs);
                    % plot results
                    subplot(length(Ks),length(Ls),l+(k-1)*length(Ks))
                    hold on
                    plot(lambdas,Rs, 'b')
                    plot(lambdas,Cs,'r')
                    scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
                    scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
                    xlabel('Lambda'); ylabel('Cost (au)')
                    if k==1; set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'off','Location','southoutside');end
                    set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none','xtick',lambdas)
                    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
                    title(['K=' num2str(Ks(k)) ',L=' num2str(Ls(l)*obj.DownSampleFactor) 'ms'])
                    axis square
                end
            end
            if obj.SavePlot
                saveas(gcf,[AnalysisOpts.ResultsSavePath 'LambdaCost_Ch' num2str(obj.ChNum)  '.png'],'png')
            end
        end
        function varargout=PlotParamSweep_Lambda(obj,ParamVals,data,T_tot,K_max,varargin)% Plots param sweep with respect to lambda
                fprintf('\n PlotParamSweep_Lambda:Plots param sweep with respect to lambda ...');
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %         T_tot=30; % total time
            %         K_max=15;
            n_MaxTim  =find(floor(ParamVals.MaxTimeSet)==T_tot);
            n_L       =length(ParamVals.L_ms) ; % L is in s
            n_K       =find(ParamVals.K ==K_max); % K
            n_lambda  =find(ParamVals.lambda);
            n_lambOrth=find(ParamVals.lambdaOrthoH==1);
            
            lambdas=ParamVals.lambda;
            Nch=length(data); % how many channels do we have
            Nraw=ceil(n_L/5);
            Ncol=5;
            UseIntersectLambda=size(data{1}.Output_Disc.cost,4)==1;
            for l=1:n_L
                for i=1:Nch
                    this_regularization=transpose(squeeze(mean(data{i}.Output_Disc.regularization(n_MaxTim,l,n_K,:,n_lambOrth,:),6)));
                    this_cost=transpose(squeeze(mean(data{i}.Output_Disc.cost(n_MaxTim,l,n_K,:,n_lambOrth,:),6)));
                    %  PEV_ind(i,:)=arrayfun(@(x) sum(data{i}.Output.PEV_ind{n_MaxTim,l,n_K,x,n_lambOrth}),1:length(ParamVals.lambda));
                    PEV_all_Disc(i,:)=transpose(squeeze(mean(data{i}.Output_Disc.PEV_all(n_MaxTim,l,n_K,:,n_lambOrth,:),6)));
                    PEV_all_Test(i,:)=transpose(squeeze(mean(data{i}.Output_Test.PEV_all(n_MaxTim,l,n_K,:,n_lambOrth,:),6)));
                    NMotif(i,:)=transpose(squeeze(mean(data{i}.Output_Disc.NMotifs(n_MaxTim,l,n_K,:,n_lambOrth,:),6)))/K_max;
                    regularization(i,:)=mapminmax(this_regularization,0,1);
                    cost(i,:)=mapminmax(this_cost,0,1);
                end
                
                if UseIntersectLambda % if we are using intersect lambda
                    InterSect=ParamVals.lambda(n_MaxTim,:);
                    IntersectPEV.Disc(:,l)=PEV_all_Disc;
                    IntersectPEV.Test(:,l)=PEV_all_Test;
                    IntersectNMotifs(:,l)=NMotif;
                end
                
                
                % if we only have one lambda then don't plot these resutls;
                % it means we have calculated the intersect lambda
                if ~UseIntersectLambda
                    
                    % plot results
                    if l==1; varargout{1}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);end
                    subplot(Nraw,Ncol,l)
                    hold on
                    [rcl,reg_mean]=obj.FigParams.PlotMeanStd(lambdas, regularization,[],'Lambda','Cost (au)','b',0,'');
                    [ccl,cost_mean]=obj.FigParams.PlotMeanStd(lambdas, cost,[],[],[],'r',0,'');
                    pel=obj.FigParams.PlotMeanStd(lambdas, PEV_all_Disc,[],[],[],'k',0,'');
                    nml=obj.FigParams.PlotMeanStd(lambdas, NMotif,[],[],[],'g',0,'');
                    
                    set(gca, 'xscale', 'log', 'ytick', [0 0.25 0.5 0.75 1], 'color', 'none','xtick',lambdas(1:3:end),'fontsize',7);grid off
                    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
                    title(['L=' num2str(ParamVals.L_ms(l)) 'sec-T=' num2str(T_tot) ])
                    axis square
                    %% fit poly nomials to mean and see where they intersect
                    [InterSect(l),LineFunc] = FindLineInterSection(lambdas,cost_mean,reg_mean);
                    inl=plot((InterSect(l)),LineFunc(InterSect(l)),'k*','markersize',12);
                    
                    %% fit polynomials to each PEV value in discovery and test; see what is the PEV value at the intersect
                    %                 MeanPEV=mean(PEV_all,1);
                    %                 [polyPEV,S_polyPEV]=polyfit(log(lambdas),MeanPEV,7);
                    %                 IntersectPEV(l)=polyval(polyPEV,log(InterSect(l)));
                    %                 plot(InterSect(l),IntersectPEV(l),'r*','markersize',12);
                    for i=1:Nch
                        % fit for discovery epoch
                        [polyPEV_Disc{i},S_polyPEV_Disc(i)]=polyfit(log(lambdas),PEV_all_Disc(i,:),7);
                        IntersectPEV.Disc(i,l)=polyval(polyPEV_Disc{i},log(InterSect(l)));
                        % fit for test epoch
                        [polyPEV_Test{i},S_polyPEV_Test(i)]=polyfit(log(lambdas),PEV_all_Test(i,:),7);
                        IntersectPEV.Test(i,l)=polyval(polyPEV_Test{i},log(InterSect(l)));
                        % for for number of motifs
                        [polyNMotifs{i},S_polyNMotifs(i)]=polyfit(log(lambdas),NMotif(i,:),7);
                        IntersectNMotifs(i,l)=polyval(polyNMotifs{i},log(InterSect(l)));
                    end
                    plot( (InterSect(l)),mean(IntersectPEV.Disc(:,l)),'r*','markersize',12);
                    plot( (InterSect(l)),mean(IntersectNMotifs(:,l)),'b*','markersize',8);
                    %
                    %                 [regfit,reg_gop] =fit(lambdas', reg_mean','poly6');
                    %                 [costfit,cost_gop]=fit(lambdas',cost_mean','poly6');
                    %                 cof=@(x) regfit(x)-costfit(x);
                    %                 InterSect(l)=fzero(cof,0.0002);
                    %                 inl=plot(InterSect(l),regfit(InterSect(l)),'k*','markersize',12);
                    % plot legend
                    if l==1; legend([rcl ccl pel nml inl],{'Correlation cost', 'Reconstruction cost','PEV',['# motif/' num2str(K_max)]}, 'Box', 'off','Location','southoutside');end
                else
                    varargout{1}=1;
                end
            end
            varargout{2}=InterSect;  % send out the intersect
            varargout{3}=IntersectPEV;  % intersect of PEV
            varargout{4}=IntersectNMotifs; % intersect of NMotifs
        end
        function varargout=PlotParamSweep_LambdaFreq(obj,ParamVals,data,K_max,L_ms,T_tot,varargin)% Plots param sweep with respect to lambda and freq sampling steps
            fprintf('\n PlotParamSweep_LambdaFreq:Plots param sweep with respect to lambda and freq sampling steps ...');
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %         T_tot=30; % total time
            %         K_max=15;
            n_MaxTim  =find(ParamVals.MaxTimeSet==T_tot);
            n_L       =find(floor(ParamVals.L_ms)==L_ms);  % L is in s
            n_K       =find(ParamVals.K ==K_max); % K
            n_lambda  =find(ParamVals.lambda);
            n_lambOrth=find(ParamVals.lambdaOrthoH==1);
            nFreqStp=length(ParamVals.ParamSwpFrqStp);
            
            FreqStps=ParamVals.ParamSwpFrqStp;
            lambdas=ParamVals.lambda;
            Nch=length(data); % how many channels do we have
            Nraw=ceil(nFreqStp/5);
            Ncol=5;
            UseIntersectLambda=size(data{1}.Output_Disc.cost,4)==1;
            
            for fs=1:nFreqStp
                for i=1:Nch
                    this_regularization=transpose(squeeze(data{i}.Output_Disc.regularization(n_MaxTim,n_L,n_K,:,n_lambOrth,fs)));
                    this_cost=transpose(squeeze(data{i}.Output_Disc.cost(n_MaxTim,n_L,n_K,:,n_lambOrth,fs)));
                    %  PEV_ind(i,:)=arrayfun(@(x) sum(data{i}.Output.PEV_ind{n_MaxTim,l,n_K,x,n_lambOrth}),1:length(ParamVals.lambda));
                    PEV_all_Disc(i,:)=transpose(squeeze(data{i}.Output_Disc.PEV_all(n_MaxTim,n_L,n_K,:,n_lambOrth,fs)));
                    PEV_all_Test(i,:)=transpose(squeeze(data{i}.Output_Test.PEV_all(n_MaxTim,n_L,n_K,:,n_lambOrth,fs)));
                    NMotif(i,:)=transpose(squeeze(data{i}.Output_Disc.NMotifs(n_MaxTim,n_L,n_K,:,n_lambOrth,fs)))/K_max;
                    regularization(i,:)=mapminmax(this_regularization,0,1);
                    cost(i,:)=mapminmax(this_cost,0,1);
                end
                
                if UseIntersectLambda % if we are using intersect lambda
                    InterSect=ParamVals.lambda(n_MaxTim,:);
                    IntersectPEV.Disc(:,fs)=PEV_all_Disc;
                    IntersectPEV.Test(:,fs)=PEV_all_Test;
                    IntersectNMotifs(:,fs)=NMotif;
                end               
                % if we only have one lambda then don't plot these resutls;
                % it means we have calculated the intersect lambda
                if ~UseIntersectLambda
                    
                    % plot results
                    if fs==1; varargout{1}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);end
                    subplot(Nraw,Ncol,fs)
                    hold on
                    [rcl,reg_mean]=obj.FigParams.PlotMeanStd(lambdas, regularization,[],'Lambda','Cost (au)','b',0,'');
                    [ccl,cost_mean]=obj.FigParams.PlotMeanStd(lambdas, cost,[],[],[],'r',0,'');
                    pel=obj.FigParams.PlotMeanStd(lambdas, PEV_all_Disc,[],[],[],'k',0,'');
                    nml=obj.FigParams.PlotMeanStd(lambdas, NMotif,[],[],[],'g',0,'');
                    
                    set(gca, 'xscale', 'log', 'ytick', [0 0.25 0.5 0.75 1], 'color', 'none','xtick',lambdas(1:3:end),'fontsize',7);grid off
                    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
                    title(['Freq Stp=' num2str(ParamVals.ParamSwpFrqStp(fs)) 'sec-T=' num2str(T_tot) ])
                    axis square
                    %% fit poly nomials to mean and see where they intersect
                    [InterSect(fs),LineFunc] = FindLineInterSection(lambdas,cost_mean,reg_mean);
                    inl=plot((InterSect(fs)),LineFunc(InterSect(fs)),'k*','markersize',12);
                    
                    %% fit polynomials to each PEV value in discovery and test; see what is the PEV value at the intersect                  
                    for i=1:Nch
                        % fit for discovery epoch
                        [polyPEV_Disc{i},S_polyPEV_Disc(i)]=polyfit(log(lambdas),PEV_all_Disc(i,:),7);
                        IntersectPEV.Disc(i,fs)=polyval(polyPEV_Disc{i},log(InterSect(fs)));
                        % fit for test epoch
                        [polyPEV_Test{i},S_polyPEV_Test(i)]=polyfit(log(lambdas),PEV_all_Test(i,:),7);
                        IntersectPEV.Test(i,fs)=polyval(polyPEV_Test{i},log(InterSect(fs)));
                        % for for number of motifs
                        [polyNMotifs{i},S_polyNMotifs(i)]=polyfit(log(lambdas),NMotif(i,:),7);
                        IntersectNMotifs(i,fs)=polyval(polyNMotifs{i},log(InterSect(fs)));
                    end
                    plot( (InterSect(fs)),mean(IntersectPEV.Disc(:,fs)),'r*','markersize',12);
                    plot( (InterSect(fs)),mean(IntersectNMotifs(:,fs)),'b*','markersize',8);
                    %
                    %                 [regfit,reg_gop] =fit(lambdas', reg_mean','poly6');
                    %                 [costfit,cost_gop]=fit(lambdas',cost_mean','poly6');
                    %                 cof=@(x) regfit(x)-costfit(x);
                    %                 InterSect(l)=fzero(cof,0.0002);
                    %                 inl=plot(InterSect(l),regfit(InterSect(l)),'k*','markersize',12);
                    % plot legend
                    if fs==1; legend([rcl ccl pel nml inl],{'Correlation cost', 'Reconstruction cost','PEV',['# motif/' num2str(K_max)]}, 'Box', 'off','Location','southoutside');end
                else
                    varargout{1}=1;
                end
            end
             
            varargout{2}=InterSect;  % send out the intersect
            varargout{3}=IntersectPEV;  % intersect of PEV
            varargout{4}=IntersectNMotifs; % intersect of NMotifs
        end
        function varargout=PlotParamSweep_L(obj,Data,ParamVals,data,InterSect,InterSectPEV,InterSectNMotifs,K_Max,MaxTim,varargin) % sweep parameters with respect to L
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %%%
            varargout{1}=obj.FigParams.RenderFigure(1,[]);hold on
            %  n_MaxTim  =length(ParamVals.MaxTimeSet);
            n_MaxTim  =find(ParamVals.MaxTimeSet ==MaxTim);
            n_L       =ParamVals.L_ms ; % L is in s
            n_K       =find(ParamVals.K ==K_Max); % K
            n_lambOrth=find(ParamVals.lambdaOrthoH==1);
            % find the lambda that is closes to the intersection
            n_lambda=arrayfun(@(x) find(ParamVals.lambda>=InterSect(x),1,'first'),1:length(n_L));
            Col=distinguishable_colors(n_MaxTim);Legends=[];
            mt=n_MaxTim;
            %% plot chnage of L and PEV for different Max times
%             for mt=n_MaxTim
%                 NMotif=[];
%                 for i=1:length(data)
%                     %  PEV(i,:)=arrayfun(@(x) sum(data{i}.Output.PEV{mt,x,n_K,n_lambda,n_lambOrth}),1:length(n_L));
%                     for j=1:length(n_lambda)
%                         PEV_all(i,j)=transpose(squeeze(mean(data{i}.Output_Disc.PEV_all(mt,j,n_K,n_lambda(j),n_lambOrth,:),6)));
%                         PEV_all_test(i,j)=transpose(squeeze(mean(data{i}.Output_Test.PEV_all(mt,j,n_K,n_lambda(j),n_lambOrth,:),6)));
%                         NMotif(j,:,i)=transpose(squeeze((data{i}.Output_Disc.NMotifs(mt,j,:,n_lambda(j),n_lambOrth,1))));
%                     end
%                 end
%             end
            %% plot PEV value for intersect for this MaxTime
            X=n_L;
            Xlbl='Motif Length(sec)';Ylbl='PEV';
            subplot(2,3,[1 2]);axis square
            plt=obj.FigParams.PlotMeanStd(X,InterSectPEV.Disc,[],Xlbl,Ylbl,2*mt-1,0,'');
            plttest=obj.FigParams.PlotMeanStd(X,InterSectPEV.Test,[],Xlbl,Ylbl,2*mt,0,'PEV changes with motif length(L)','--');
            Legends=[Legends plt plttest];
            if obj.LogScale
                set(gca, 'xscale', 'log','xtick',n_L(1:4:end))
            else
                set(gca,'xtick',n_L(1:4:end))
            end
            set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
            axis square
            %   legend(Legends,arrayfun(@(x) ['T=' num2str(x)],ParamVals.MaxTimeSet,'UniformOutput',0),'Location','northeastoutside')
            legend(Legends,cellfun(@(x) ['T=' num2str(ParamVals.MaxTimeSet(mt)) ' ' x],{'Disc','Test'},'UniformOutput',0),'Location','northeastoutside')
            %% plot difference of PEVs between disc and test for different L values
            subplot(2,3,3);axis square
            obj.FigParams.PlotMeanStd(X,InterSectPEV.Disc-InterSectPEV.Test,[],Xlbl,'PEV difference',2*mt-1,0,'Difference of PEV between Dics and Test');
            %% plot number of motifs for best lambda for differrent K across changes in Max times
            for mt=n_MaxTim
                subplot(2,3,3+mt);
                %  plot(X,squeeze(mean(NMotif,3)))
                obj.FigParams.PlotMeanStd(X,InterSectNMotifs*ParamVals.K,[],Xlbl,'Mean number of motifs',Col(mt,:),1,'');
                if obj.LogScale
                    set(gca,'xscale','log','xtick',n_L(1:4:end))
                else
                    set(gca,'xtick',n_L(1:4:end))
                end
                set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
                legend(arrayfun(@(x) ['K=' num2str(x)],ParamVals.K,'UniformOutput',0),'Location','bestoutside');
                title(['T=' num2str(ParamVals.MaxTimeSet(mt))])
                axis square
            end
            %% Plot charctristics of the motifs
            % grab data from all of the channels and across different
            % Motif length(L) values
            Fields={'Duration','FracTimActive','FracSimTimActive','CTAhalflife'};
            NRow=2;NCol=ceil(length(Fields)/NRow);
            XAxisLbl={'Motif Length(s)','Motif Length(s)','Motif Length(s)','Motif Length(s)'};
            YAxisLbl={'Duration(s)','Fraction Time Active','Fraction Time Sim Active','Motif Half Life'};
            Title={'Duration of Motifs','Fraction of Time Motif Active','Fraction of Time 2 or more Motifs Active','Motif Cross Temporal Autoccorelation'};
            if 0
                for mt=n_MaxTim
                    varargout{1+mt}=obj.FigParams.RenderFigure(1,[]);hold on
                    MotifChractr=cellfun(@(y) cell2mat(arrayfun(@(x) obj.CharactrizeMotifActivity(y.AnalysisData.MotifSpecs_Disc{mt,x,n_K,n_lambda(x),n_lambOrth,:},...
                        'ShowPlot',0),1:length(n_lambda),'UniformOutput',0)),Data,'UniformOutput',0);
                    MotifChractr=cell2mat(MotifChractr');
                    % now go through and plot each property in a sperate color
                    arrayfun(@(x) obj.FigParams.PlotMeanStd(X,reshape([MotifChractr.(Fields{x})],size(MotifChractr)),[],XAxisLbl{x},YAxisLbl{x},x,0,Title{x},'ThisSubplot',[NRow NCol x]),1:length(Fields))
                end
                % show a few example charactristics of motifs with different Ls
                ShowCharcLs=[0.1 0.9 2 3];ShowCharcLsInd=arrayfun(@(x) find(n_L==x),ShowCharcLs,'UniformOutput',1);
                [~,CharacFigs1,CharacFigs2]=arrayfun(@(x) obj.CharactrizeMotifActivity(Data{2}.AnalysisData.MotifSpecs_Disc{n_MaxTim,x,n_K,n_lambda(x),n_lambOrth,:},...
                    'ShowPlot',1,'FigTitle',['L='  num2str(n_L(x))]),ShowCharcLsInd,'UniformOutput',0);
            end
            CharacFigs1=1;CharacFigs2=2;
            % update varargout
            varargout=[varargout CharacFigs1 CharacFigs2];
        end
        function varargout=PlotParamSweep_Freq(obj,Data,ParamVals,data,InterSect,InterSectPEV,InterSectNMotifs,K_Max,MaxTim,varargin) % plots sweep parameters with respect to freq
            fprintf('\nPlotParamSweep_Freq:plots sweep parameters with respect to freq ...')
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %%%
            varargout{1}=obj.FigParams.RenderFigure(1,[]);hold on
            %  n_MaxTim  =length(ParamVals.MaxTimeSet);
            n_MaxTim  =find(ParamVals.MaxTimeSet ==MaxTim);
            n_L       =ParamVals.L_ms ; % L is in s
            n_K       =find(ParamVals.K ==K_Max); % K
            n_lambOrth=find(ParamVals.lambdaOrthoH==1);
            nFreqStp=length(ParamVals.ParamSwpFrqStp);
            % find the lambda that is closes to the intersection
            n_lambda=arrayfun(@(x) find(ParamVals.lambda>=InterSect(x),1,'first'),1:nFreqStp);
            Col=distinguishable_colors(n_MaxTim);Legends=[];
            mt=n_MaxTim;
%             %% plot chnage of L and PEV for different Max times
%             for mt=n_MaxTim
%                 NMotif=[];
%                 for i=1:length(data)
%                     %  PEV(i,:)=arrayfun(@(x) sum(data{i}.Output.PEV{mt,x,n_K,n_lambda,n_lambOrth}),1:length(n_L));
%                     for j=1:length(n_lambda)
%                         PEV_all(i,j)=transpose(squeeze(mean(data{i}.Output_Disc.PEV_all(mt,j,n_K,n_lambda(j),n_lambOrth,:),6)));
%                         PEV_all_test(i,j)=transpose(squeeze(mean(data{i}.Output_Test.PEV_all(mt,j,n_K,n_lambda(j),n_lambOrth,:),6)));
%                         NMotif(j,:,i)=transpose(squeeze((data{i}.Output_Disc.NMotifs(mt,j,:,n_lambda(j),n_lambOrth,1))));
%                     end
%                 end
%             end
            %% plot PEV value for intersect for this MaxTime
            X=ParamVals.ParamSwpFrqStp;
            Xlbl='Freq Step';Ylbl='PEV';
            subplot(2,3,[1 2]);axis square
            plt=obj.FigParams.PlotMeanStd(X,InterSectPEV.Disc,[],Xlbl,Ylbl,2*mt-1,0,'');
            plttest=obj.FigParams.PlotMeanStd(X,InterSectPEV.Test,[],Xlbl,Ylbl,2*mt,0,'PEV changes with freq step','--');
            Legends=[Legends plt plttest];
            if obj.LogScale
                set(gca, 'xscale', 'log','xtick',ParamVals.ParamSwpFrqStp )
            else
                set(gca,'xtick',ParamVals.ParamSwpFrqStp )
            end
            set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
            axis square
            %   legend(Legends,arrayfun(@(x) ['T=' num2str(x)],ParamVals.MaxTimeSet,'UniformOutput',0),'Location','northeastoutside')
            legend(Legends,cellfun(@(x) ['T=' num2str(ParamVals.MaxTimeSet(mt)) ' ' x],{'Disc','Test'},'UniformOutput',0),'Location','northeastoutside')
            %% plot difference of PEVs between disc and test for different L values
            subplot(2,3,3);axis square
            obj.FigParams.PlotMeanStd(X,InterSectPEV.Disc-InterSectPEV.Test,[],Xlbl,'PEV difference',2*mt-1,0,'Difference of PEV between Dics and Test');
            %% plot number of motifs for best lambda for differrent K across changes in Max times
            for mt=n_MaxTim
                subplot(2,3,3+mt);
                %  plot(X,squeeze(mean(NMotif,3)))
                obj.FigParams.PlotMeanStd(X,InterSectNMotifs*ParamVals.K,[],Xlbl,'Mean number of motifs',Col(mt,:),1,'');
                if obj.LogScale
                    set(gca,'xscale','log','xtick',n_L(1:4:end))
                else
                    set(gca,'xtick',n_L(1:4:end))
                end
                set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
                legend(arrayfun(@(x) ['K=' num2str(x)],ParamVals.K,'UniformOutput',0),'Location','bestoutside');
                title(['T=' num2str(ParamVals.MaxTimeSet(mt))])
                axis square
            end
            %% Plot charctristics of the motifs
            % grab data from all of the channels and across different
            % Motif length(L) values
            Fields={'Duration','FracTimActive','FracSimTimActive','CTAhalflife'};
            NRow=2;NCol=ceil(length(Fields)/NRow);
            XAxisLbl={'Motif Length(s)','Motif Length(s)','Motif Length(s)','Motif Length(s)'};
            YAxisLbl={'Duration(s)','Fraction Time Active','Fraction Time Sim Active','Motif Half Life'};
            Title={'Duration of Motifs','Fraction of Time Motif Active','Fraction of Time 2 or more Motifs Active','Motif Cross Temporal Autoccorelation'};
            
%             for mt=n_MaxTim
%                 varargout{1+mt}=obj.FigParams.RenderFigure(1,[]);hold on
%                 MotifChractr=cellfun(@(y) cell2mat(arrayfun(@(x) obj.CharactrizeMotifActivity(y.AnalysisData.MotifSpecs_Disc{mt,x,n_K,n_lambda(x),n_lambOrth,:},...
%                     'ShowPlot',0),1:length(n_lambda),'UniformOutput',0)),Data,'UniformOutput',0);
%                 MotifChractr=cell2mat(MotifChractr');
%                 % now go through and plot each property in a sperate color
%                 arrayfun(@(x) obj.FigParams.PlotMeanStd(X,reshape([MotifChractr.(Fields{x})],size(MotifChractr)),[],XAxisLbl{x},YAxisLbl{x},x,0,Title{x},'ThisSubplot',[NRow NCol x]),1:length(Fields))
%             end
            % show a few example charactristics of motifs with different Ls
            ShowCharcFreq=ParamVals.ParamSwpFrqStp;
            ShowCharcFreqInd=arrayfun(@(x) find(ParamVals.ParamSwpFrqStp==x),ShowCharcFreq,'UniformOutput',1);
          %  MotifFig=cell(1,length(ShowCharcFreq));
%             [~,CharacFigs1,CharacFigs2]=arrayfun(@(x) obj.CharactrizeMotifActivity(Data{2}.AnalysisData.MotifSpecs_Disc{n_MaxTim,2,n_K,n_lambda(x),n_lambOrth,x},...
%                 'ShowPlot',1,'FigTitle',['L='  num2str(n_L(x))]),ShowCharcFreqInd,'UniformOutput',0);
             MotifFig=arrayfun(@(x) obj.PlotCoreMotifs(Data{2}.AnalysisData.MotifSpecs_Disc{n_MaxTim,2,n_K,n_lambda(x),n_lambOrth,x}.W,[],0,...
                'NewFig',1,'ShowPlot',1,'FigTitle',['FreqStp='  num2str(ParamVals.ParamSwpFrqStp(x))]),ShowCharcFreqInd,'UniformOutput',0);
            % update varargout
            varargout=[varargout MotifFig];
          
        end
        
        function varargout=PlotParamSweep_L_old(obj,ParamVals,data,InterSect,InterSectPEV,K_Max,MaxTim,varargin) % sweep parameters with respect to L
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %%%
            varargout{1}=obj.FigParams.RenderFigure(1,[]);hold on
            %  n_MaxTim  =length(ParamVals.MaxTimeSet);
            n_MaxTim  =find(ParamVals.MaxTimeSet ==MaxTim);
            n_L       =ParamVals.L_ms ; % L is in s
            n_K       =find(ParamVals.K ==K_Max); % K
            % n_lambda  =find(ParamVals.lambda==ParamVals.lambda(4));
            n_lambOrth=find(ParamVals.lambdaOrthoH==1);
            
            %   n_lambda=arrayfun(@(x) find(ParamVals.lambda<=InterSect(x),1,'last'),1:length(n_L));
            n_lambda=arrayfun(@(x) find(ParamVals.lambda>=InterSect(x),1,'first'),1:length(n_L));
            Col=distinguishable_colors(n_MaxTim);Legends=[];
            %% plot chnage of L and PEV for different Max times
            subplot(2,3,[1 2]);axis square
            for mt=n_MaxTim
                NMotif=[];
                for i=1:length(data)
                    %  PEV(i,:)=arrayfun(@(x) sum(data{i}.Output.PEV{mt,x,n_K,n_lambda,n_lambOrth}),1:length(n_L));
                    for j=1:length(n_lambda)
                        PEV_all(i,j)=transpose(squeeze(mean(data{i}.Output_Disc.PEV_all(mt,j,n_K,n_lambda(j),n_lambOrth,:),6)));
                        PEV_all_test(i,j)=transpose(squeeze(mean(data{i}.Output_Test.PEV_all(mt,j,n_K,n_lambda(j),n_lambOrth,:),6)));
                        NMotif(j,:,i)=transpose(squeeze((data{i}.Output_Disc.NMotifs(mt,j,:,n_lambda(j),n_lambOrth,1))));
                    end
                end
                % plot results
                X=n_L;
                Xlbl='Motif Length(sec)';Ylbl='PEV';
                plt=obj.FigParams.PlotMeanStd(X,PEV_all,[],Xlbl,Ylbl,Col(mt,:),1,'');
                plttest=obj.FigParams.PlotMeanStd(X,PEV_all_test,[],Xlbl,Ylbl,1-Col(mt,:),1,'','--');
                Legends=[Legends plt plttest];
                if obj.LogScale
                    set(gca, 'xscale', 'log','xtick',n_L(1:4:end))
                else
                    set(gca,'xtick',n_L(1:4:end))
                end
                set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
                axis square
            end
            legend(Legends,arrayfun(@(x) ['T=' num2str(x)],ParamVals.MaxTimeSet,'UniformOutput',0),'Location','northeastoutside')
            %% plot number of motifs for best lambda for differrent K across changes in Max times
            for mt=n_MaxTim
                subplot(2,3,3+mt);
                %  plot(X,squeeze(mean(NMotif,3)))
                obj.FigParams.PlotMeanStd(X,squeeze(NMotif)',[],Xlbl,'Mean number of motifs',Col(mt,:),1,'');
                if obj.LogScale
                    set(gca, 'xscale', 'log','xtick',n_L(1:4:end))
                else
                    set(gca,'xtick',n_L(1:4:end))
                end
                set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
                legend(arrayfun(@(x) ['K=' num2str(x)],ParamVals.K,'UniformOutput',0),'Location','bestoutside');
                title(['T=' num2str(ParamVals.MaxTimeSet(mt))])
                axis square
            end
        end
        function [MotifChractr,varargout]=CharactrizeMotifActivity(obj,MotifSpecs,varargin)% charactrizes motifs,
            % 1-how long they are active
            % 2-fratction of time they are active
            % 3-degree of compositionality of motifs
            % 4-cross temporal autocorrelation of motifs
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            ShowPlots=obj.ShowPlot;NRow=2;NCol=3;
            % grab all of the necessary info from Motif specs
            MotifSpecs=obj.RemoveNoloadingMotifs(MotifSpecs);
            W=MotifSpecs.W;
            H=MotifSpecs.H;
            PEVind=MotifSpecs.PEVind;
            NW=size(H,1);
            NTim=size(H,2);
            FsWaveTarg=AnalysisOpts.MotifAnalysis.FsWaveTarg;
            %% how long motifs are active
            % find the times each motif is hight than 1 std of the mean
            ActTh=1;% std of mean
            MovMeanSamp=1;
            Npadd=1;
            stdH=std(H,0,2);
            SmoothedHPadded=padarray(movmean(H,MovMeanSamp,2),[0 Npadd],0); % smmoth and padd H with zeros
            ActTimeH=arrayfun(@(x) SmoothedHPadded(x,:)>ActTh*stdH(x),1:NW,'UniformOutput',0); % active time of H
            DiffActTimeH=cellfun(@diff,ActTimeH,'UniformOutput',0); % what is the diff
            StrTim=cellfun(@(x) find(x==1),DiffActTimeH,'UniformOutput',0);
            StpTim=cellfun(@(x) find(x==-1),DiffActTimeH,'UniformOutput',0);
            % remove the StrTims that are unfinished, usually they happen
            % at the end
            %             StrTim=arrayfun(@(x) StrTim{x}(1:length(StpTim{x})),1:NMotifs,'uniformoutput',0);
            MotifDuration=arrayfun(@(x) (StpTim{x}-StrTim{x})/FsWaveTarg,1:NW,'UniformOutput',0);
            MotifChractr.Duration=mean(cell2mat(MotifDuration));
            if ShowPlots
                % plot an example
                varargout(1)=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]); % create figures
                subplot(NRow,NCol,1)
                expnum=1;
                X=(1:length(ActTimeH{expnum}))/FsWaveTarg;
                obj.FigParams.Plot(X,ActTimeH{expnum}*max(SmoothedHPadded(expnum,:),[],2),'b',' ',' ',' ');
                hold on;
                obj.FigParams.Plot(X,SmoothedHPadded(expnum,:),'r',' ',' ',' ');
                v=axis;
                obj.FigParams.Plot([v(1) v(2)],[ActTh*stdH(expnum) ActTh*stdH(expnum)],'g','Time(s)','H Value','Example duration of Motif');
                % show distribution of values for duration of motif
                subplot(NRow,NCol,2)
                obj.FigParams.HistogramPlot(cell2mat(MotifDuration),[0:0.1:1],1,'Motif Duration(s)','Freq of occurance',...
                    'Distribution of Motif durations');
            end
            %% calculate fraction of time Motifs are active
            ActiveTim=cellfun(@sum ,MotifDuration)/(NTim/FsWaveTarg);
            MotifChractr.FracTimActive=mean(ActiveTim);
            if ShowPlots
                subplot(NRow,NCol,3)
                obj.FigParams.HistogramPlot(ActiveTim,0:0.1:1,2,'Fraction Active Time','Freq of occurance','Total active time of Motifs');
            end
            %% calculate degree of compositionality of Motifs
            % calculate the fraction of time that motifs are active
            % toghether by summing up their activity and detecting where
            % more than two motifs are active at a time
            SumActTimeH=sum(obj.ManData.ReshapeCell2Mat(ActTimeH,2),1);
            SumActTimeH=SumActTimeH(Npadd+1:end-Npadd); % remove padding
            if ShowPlots
                % plot results
                subplot(NRow,NCol,4)
                obj.FigParams.HistogramPlot(SumActTimeH,[1:10],1,'Fraction Simltanous Active Time','Freq of occurance','Co-occurance of Motifs');
            end
            MotifChractr.FracSimTimActive= sum(SumActTimeH>1)/NTim;
            % calculate average pairwise correlation between Hs
            HpairNW=nchoosek(1:NW,2);
            MotifChractr.MeanCorrH=mean(arrayfun(@(x) corr(H(HpairNW(x,1))',H(HpairNW(x,2))'),1:size(HpairNW,1)));
            %% calculate cross temporal autocorrelation of motifs
            ConcW=obj.ConcatinateWs({W})';
            [~,~,~,CTA,Rng]=arrayfun(@(x) TemporalCrossCorr(ConcW(x,:),ConcW(x,:),[size(W,1) size(W,3)]),1:NW,'UniformOutput',0);
            CTA=obj.ManData.ReshapeCell2Mat(CTA,2);
            Rng=Rng{1};ZeroLag=find(Rng==0);
            Rng=Rng(ZeroLag:end)/FsWaveTarg;
            CTA=CTA(:,ZeroLag:end);
            %             % fit exponential to each motif CTA
            %             [CTAexpfit,gof]=arrayfun(@(x) FitExp(Rng',CTA(x,:)'),1:NW,'UniformOutput',0);
            %             MotifSpecs.CTAfita=mean(cellfun(@(x) [x.a],CTAexpfit));
            %             MotifSpecs.CTAfitb=mean(cellfun(@(x) [x.b],CTAexpfit));
            %             MotifSpecs.CTAfitc=mean(cellfun(@(x) [x.c],CTAexpfit));
            %%             fit exp to mean CTA
            [CTAexpfit,gof]= FitExp(Rng',mean(CTA,1)');
            MotifChractr.CTAfita=CTAexpfit.a;
            MotifChractr.CTAfitb=CTAexpfit.b;
            MotifChractr.CTAfitc=CTAexpfit.c;
            MotifChractr.CTAhalflife=log(2)/CTAexpfit.c;
            Exp=@(x) MotifChractr.CTAfita-MotifChractr.CTAfitb*exp(-MotifChractr.CTAfitc*x);
            
            %%
            if ShowPlots
                % plot results
                subplot(NRow,NCol,5);hold on
                obj.FigParams.Plot(Rng,Exp(Rng),2,'','',''); % show exp estimate
                obj.FigParams.PlotMeanStd(Rng,CTA,[],'Temporal Offset within Motif','correlation',1,1,'cross temporal autocorrelation')
                obj.FigParams.Text(Rng(end)/2,0.5,['\tau=' num2str(MotifChractr.CTAhalflife,2) ' s'],'k');
                %%  show all of the motifs for this example
                varargout(2)=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]); % create figures
                obj.PlotCoreMotifs(W,[],0)
            end
        end
        function  varargout=PlotParamCost(obj,data,Lambda,L,Title,varargin)
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            varargout{1}=obj.FigParams.RenderFigure(1,[]);
            data=squeeze(data);
            Siz=size(data);
            for i=1:Siz(3)
                for j=1:Siz(2)
                    A=[];
                    for k=1:Siz(1)
                        A(k)=sum(data{k,j,i});
                    end
                    PEV_mean(j,i)=mean(A);
                    PEV_std=std(A);
                end
            end
            helperCWTTimeFreqPlot(PEV_mean,L,Lambda,'justplot1',Title,...
                'L','Lambda',obj.LogScale);
            colorbar
        end
        %%%% function for clustering of motifs in an area with many motifs that need
        %%%% cluster
        function  [CoreMotifs,SizeW,cwt_f]=ClassifyMotifs_Area(obj,XCorrFileName,varargin) % classify motifs for many channels of area
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % if we have it don't recalculate it
            if obj.ManData.IsVarExistinFile(XCorrFileName,'CoreMotifs') & ~obj.OverWrite
                load(XCorrFileName,'CoreMotifs','SizeW');
                return
            end
            
            NRandMotif=AnalysisOpts.MotifAnalysis.NSampledMotifs; % how many random motifs we want to classify with phenoraoph first
            load(XCorrFileName,'Motifs','DistMatDTW','SizeW','SampledMotifsInd');
            NMotifs=size(Motifs,1);
            X=eye(size(DistMatDTW));X=X+triu(DistMatDTW,1)+triu(DistMatDTW,1)'; % X is full distance matrix
            AnalysisData.MotifSimilarityMat=triu(DistMat)+triu(DistMat)';
            AnalysisData.SizeW=SizeW;
            %% choose a random set of motifs to cluster with phenograph at first
            RandMotifDist=X;
            RandMotif=Motifs(SampledMotifsInd,:);
            
            % cluster motifs with phenograph
            [clust_ind_pheno,~, clust_ind_pheno_alllvls] = Phenograph(RandMotifDist,'Level',1);
            AnalysisData.clust_ind=clust_ind_pheno;
            % find template motifs for thse motifs for now
            obj.FindTemplateMotif(AnalysisData.clust_ind,RandMotif)
            
            % generalize this clustering to withheld  motifs
            %       [clust_ind_all]=obj.GeneralizeMotifClustering(Motifs,RandMotifInd,clust_ind_pheno,AnalysisData.MeanCluster,SizeW);
            
            CoreMotifs=AnalysisData.MeanCluster;
            save(XCorrFileName,'clust_ind_pheno_alllvls','clust_ind_pheno','CoreMotifs','-append')
        end
        function  [CoreMotifs,SizeW]=LoadCoreMotifs(obj,XCorrFileName,varargin) % loads core motifs
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % if we have it don't recalculate it
            if obj.ManData.IsVarExistinFile(XCorrFileName,'CoreMotifs')
                load(XCorrFileName,'CoreMotifs','SizeW');
            else
                CoreMotifs=[];SizeW=[];
            end
        end
        function [varargout]=PlotClusterQualityXcorr(obj,XCorrFileName,varargin) % uses Xcorr distance function
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            opts.UseEuDist=0;% if using Euclidean distance to instead of correlation
            fprintf('\nAnalyzing cluster quality...')
            % X has to be distance metric computed by  euclidean or xcorr
            if opts.UseEuDist
                % load this file's DistMat and Motifs
                load(XCorrFileName,'Motifs','EuDistMat','SizeW','SampledMotifsInd');
                X=zeros(size(EuDistMat));X=X+triu(EuDistMat,1)+triu(EuDistMat,1)'; % X is full distance matrix
                dissimMat=X;
                ExtraStr='Eu'; % show we are using Euclidean distance
            else
                % load this file's DistMat and Motifs
                load(XCorrFileName,'Motifs','DistMat','SizeW','SampledMotifsInd');
                X=eye(size(DistMat));X=X+triu(DistMat,1)+triu(DistMat,1)'; % X is full distance matrix
                dissimMat=1-X;
                ExtraStr='';
            end
            AnalysisData.SizeW=SizeW;
            % remove inf and nan values from x and set them to zero
            NaNInfInd=isnan(X(:)) | isinf(X(:));
            fprintf('Number of NaN or Inf entries in Similarity Matrix:%i',num2str(sum(NaNInfInd)));
            X(NaNInfInd)=0;
            
            %% choose a random set of motifs to cluster with phenograph at first
            Motifs=Motifs(SampledMotifsInd,:); % this is now correspoding indexing to X
            NMotifs=size(Motifs,1);
            % if we need to subsample this motif space do so
            if NMotifs>AnalysisOpts.MotifAnalysis.NSampledMotifsSubSet  
                 SubSampledMotifsInd=randsample(1:NMotifs,AnalysisOpts.MotifAnalysis.NSampledMotifsSubSet);                 
                 Motifs=Motifs(SubSampledMotifsInd,:);
                 X=X(SubSampledMotifsInd,SubSampledMotifsInd);
                 dissimMat=dissimMat(SubSampledMotifsInd,SubSampledMotifsInd);
                 NChar=num2str(AnalysisOpts.MotifAnalysis.NSampledMotifsSubSet/10000);
            else
                 NChar='';
            end
            %% take TSNe and PCA of Motifs first and save them off
            %             if ~obj.ManData.IsVarExistinFile(XCorrFileName,'TSNEMotifsDTW') | obj.OverWrite
            %                 % calculate tsne of motifs in 2D and 3D
            %                 % [TSNEMotifsDTW.Val2D, TSNEMotifsDTW.loss2D]= tsne(Motifs,'Algorithm','barneshut','NumPCAComponents',0,'Distance',@TemporalCrossCorrDTWdist);
            %                 [TSNEMotifsDTW.Val3D, TSNEMotifsDTW.loss3D]= tsne(Motifs,'Algorithm','barneshut','NumPCAComponents',0,'Distance',@TemporalCrossCorrDTWdist,'NumDimensions',3);
            %                 TSNEMotifsDTW.Val2D=TSNEMotifsDTW.Val3D(:,1:2);
            %                 % calculate PCA of motifs
            %                 [~,PCAMotifs]=pca(Motifs);
            %                 save(XCorrFileName,'TSNEMotifsDTW','PCAMotifs','-append'); % save off this data (Take a long time to process)
            %             else
            %                 load(XCorrFileName,'TSNEMotifsDTW','PCAMotifs')
            %             end
            %% try different methods of clustring and save them off in different files
            % create a filename for each clustring method based off of the
            % XCorrFileName
            [PATH,TargFile]=fileparts(XCorrFileName);
            if strcmpi(obj.MotifClassficationMethod,'PhenoCluster')
                fprintf('\n Phenograph clustering...')
                % define the filename to save data
                FileName=[PATH AnalysisOpts.FS TargFile '_' obj.MotifClassficationMethod '_k' num2str(AnalysisOpts.MotifAnalysis.Phenograph_K) '_' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '_N' NChar '.mat'];
                FigFileStr=[obj.MotifClassficationMethod '_k' num2str(AnalysisOpts.MotifAnalysis.Phenograph_K) '_' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '_TM' AnalysisOpts.MotifAnalysis.TemplateMotif '_N' NChar];
                
                if ~obj.ManData.IsVarExistinFile(FileName,'clust_ind_PhenoCluster') | obj.OverWrite
                    % perform phenophraph on similarity matrix
                    clust_ind_PhenoCluster=PhenoCluster(X,'k',AnalysisOpts.MotifAnalysis.Phenograph_K); % cluster with phenograph
                    clust_ind=clust_ind_PhenoCluster(:,1); % take the level 1 for now
                    obj.ManData.SaveVar([],clust_ind_PhenoCluster,'clust_ind_PhenoCluster',FileName,'SelfName',1);
                else
                    load(FileName,'clust_ind_PhenoCluster')
                    clust_ind=clust_ind_PhenoCluster(:,1); % take the level 1 for now
                end
            elseif strcmpi(obj.MotifClassficationMethod,'Hierarchical')
                % use heirarchical clustring to clsuter the data
                % define the filename to save data
                fprintf('\n Hierarchical clustering...')
                FileName=[PATH AnalysisOpts.FS TargFile '_' obj.MotifClassficationMethod '_' num2str(AnalysisOpts.MotifAnalysis.MaxClust) '_' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '_N' NChar '.mat'];
                FigFileStr=[obj.MotifClassficationMethod '_MC' num2str(AnalysisOpts.MotifAnalysis.MaxClust) '_' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '_TM' AnalysisOpts.MotifAnalysis.TemplateMotif '_N' NChar ];
           
                if ~obj.ManData.IsVarExistinFile(FileName,'clust_ind') | obj.OverWrite                    
                    % convert similarity matrix to dissimilarity matrix and
                    % then into pdist form
                    MaxClust=AnalysisOpts.MotifAnalysis.MaxClust;
                    yOut = squareform(dissimMat,'tovector');
                    Z = linkage(yOut,'ward');%'method' can be 'ward' or 'average'
                    % leafOrder = optimalleaforder(Z,yOut);
                    clust_ind = cluster(Z,'maxclust',MaxClust );
                    cutoff = median([Z(end-MaxClust+1,3) Z(end-MaxClust+2,3)]);
                    obj.ManData.SaveVar([],clust_ind,'clust_ind',FileName,'SelfName',1);
                    obj.ManData.SaveVar([],cutoff,'cutoff',FileName,'SelfName',1);
                    obj.ManData.SaveVar([],Z,'Z',FileName,'SelfName',1);
                    obj.ManData.SaveVar([],MaxClust,'MaxClust',FileName,'SelfName',1);                                      
                else % load the information
                    load(FileName,'clust_ind','cutoff','Z','MaxClust')
                end
               % tests what is the maximum number of clusters given the
               % core motifs with hierarchicla clustering
               % AnalysisOpts.PairNum has the current recording number
              if strcmpi(AnalysisOpts.AnalysisFocus1,'HierHypoTestMotifs')                    
                MotifSigTest=obj.HierarchicalHypothesisTestingMotifs(Motifs,SizeW,Z,X);
                [HypoTestFile.Fold,HypoTestFile.Name,HypoTestFile.Ext]=fileparts(FileName);
                 HypoTestFile.FullName=sprintf('%s%s%sHTRecs%s',HypoTestFile.Fold,AnalysisOpts.FS,...
                     HypoTestFile.Name,HypoTestFile.Ext);
                 % save each recording by its name 
                obj.ManData.SaveVar(AnalysisOpts.AnalysisPathName,MotifSigTest,...
                    ['MotifSigTestRec' num2str(AnalysisOpts.PairNum)],HypoTestFile.FullName,'SelfName',1);
                return
              elseif strcmpi(AnalysisOpts.AnalysisFocus1,'PlotHierHypoTestMotifs')  % plot the results
                  obj.PlotHierarchicalHypothesisTestingMotifs();
              end
                
            elseif strcmpi(obj.MotifClassficationMethod,'DBSCAN')
            elseif strcmpi(obj.MotifClassficationMethod,'SOM')
                % perform SOM on similarity matrix and then Phenograph on
                % that
                [clust_ind]=SOMCluster(X,'method','matlab','verbose',0,'clustringMethod',AnalysisOpts.MotifAnalysis.SOMClustringMethod);
                
            elseif strcmpi(obj.MotifClassficationMethod,'EncoderPhenograph')
                % build autoencoder, find features and feed them to
                % clusterign algorithm
                % first align motifs to a random motif
                
                if ~obj.ManData.IsVarExistinFile(FileName,'clust_ind_EncoderPhenograph') | obj.OverWrite
                    %                     RandAlightMotifInd=randsample(size(Motifs,1),1);
                    %                     Motifs_AlignedRand=obj.AlignMotifs(Motifs(RandAlightMotifInd,:),Motifs);
                    AutoEncFeature=BuildAutoEncoderMotifs(Motifs,obj.ClassificaitonAutoEncLayerSiz);
                    AutoEncFeaturesDTW=AutoEncFeature{obj.ClassificaitonAutoEncTargLayer+1}';
                    [clust_ind_EncoderPhenograph,ovr_Q] = Phenograph(AutoEncFeatures,'DistanceMetric','euclidean'); % phenograph clustering
                    
                    obj.ManData.SaveVar([],clust_ind_EncoderPhenograph,'clust_ind_EncoderPhenographDTW',FileName,'SelfName',1);
                    obj.ManData.SaveVar([],AutoEncFeatures,'AutoEncFeaturesDTW',FileName,'SelfName',1);
                    
                else
                    load(FileName,'clust_ind_EncoderPhenographDTW')
                end
                clust_ind=clust_ind_EncoderPhenograph(:,1); % take the level 1 for now
            end
            %% find template motif
            [CoreMotifs,TemplateMotif,Motifs_Aligned,clust_ind_Aligned,clust_ind_notAligned,TestCoreMotifs]=obj.FindTemplateMotif(clust_ind,Motifs,X,'CoreMotifsPerc',1);
            if 1
                %% plot Core motifs
                CoreMotifsFig=cell(1,2);
                [CoreMotifsFig(1:2)]=obj.FigParams.RenderFigure(2,[0 0 0.8 0.8]);
                figure(CoreMotifsFig{1});CoreMotifsFig{1}=obj.PlotCoreMotifs(CoreMotifs,[],1,'FigTitle',[num2str(obj.CoreMotifsPerc*100) '%']);
                figure(CoreMotifsFig{2});CoreMotifsFig{2}=obj.PlotCoreMotifs(TestCoreMotifs.P01,[],1,'FigTitle','10%-#'); % plot core motifs with 100 percent
                
                %  [CoreMotifs,TemplateMotif,Motifs_Aligned,clust_ind_Aligned,clust_ind_notAligned,TestCoreMotifs]=obj.FindTemplateMotif(clust_ind,Motifs,'CoreMotifsPerc',1);
                % compute tsne of aligned motifs
                %             if ~obj.ManData.IsVarExistinFile(FileName,'TsneMotifs_Aligned') | obj.OverWrite
                %                 fprintf('\nCalculating Tsne...')
                %                 TsneMotifs_Aligned=tsne(Motifs_Aligned,'Algorithm','exact','NumPCAComponents',0,'Distance','correlation');
                %                 [~,PCAMotifs_Aligned]=pca(Motifs_Aligned);
                %                 obj.ManData.SaveVar([],TsneMotifs_Aligned,'TsneMotifs_Aligned',FileName,'SelfName',1);
                %                 obj.ManData.SaveVar([],PCAMotifs_Aligned,'PCAMotifs_Aligned',FileName,'SelfName',1);
                %
                %              else
                %                 load(FileName,'TsneMotifs_Aligned','PCAMotifs_Aligned')
                %             end
                %% plot dendrogram
                CWC_Figs=cell(1);DendroFig=cell(1);
                if strcmpi(obj.MotifClassficationMethod,'Hierarchical')
                    [DendroFig{1}]=obj.PlotDendrogramwithMotifs(Z,cutoff,clust_ind,CoreMotifs);
                elseif strcmpi(obj.MotifClassficationMethod,'PhenoCluster')
                    %% see if there are different clusters inside some of the clusters we are seeing CWC(cluster within cluster)
                    Nclust=length(CoreMotifs);
                    CWC_MotifNum=[1:Nclust];
                    ovr_q=nan*ones(1,length(CWC_MotifNum));
                    for n_CWC=[CWC_MotifNum]
                        CWC_Fig=cell(1,2);
                        CWCind=clust_ind==n_CWC; % find motifs in this cluster
                        CWC_X=X(CWCind,CWCind);
                        % get some information about the matrix
                        CWC_X_upper=tril(nan*ones(size(CWC_X)))+triu(CWC_X);
                        CWC_X_mean(n_CWC)=nanmean(CWC_X_upper(:));
                        CWC_X_var(n_CWC)=nanvar(CWC_X_upper(:));
                        %cluster inside this cluster
                        CWC_Motifs=Motifs(CWCind,:);
                        [ovr_q(n_CWC),Entropy_CWC_X(n_CWC),WBdist{n_CWC},CoreMotifsCWC{n_CWC},CWC_Fig{:}]=obj.ClusterWithinCluster(CWC_X,CWC_Motifs,n_CWC);
                        CWC_Figs=[CWC_Figs CWC_Fig];
                    end
                    % plot distincations between different methods
                    [CWC_statsfig]=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]);
                    subplot(221); hold on% plot modularity
                    obj.FigParams.Plot(CWC_MotifNum,ovr_q,'b','Cluster #','Modularity','Clust Modularity Index','p_line_style','--');
                    subplot(222) % plot entropy
                    obj.FigParams.Plot(CWC_MotifNum,Entropy_CWC_X,'b','Cluster #','Bits','Sim Matrix Entropy','p_line_style','--');
                    subplot(223) % plot withinbetween distance
                    obj.FigParams.Plot(CWC_MotifNum,1./cellfun(@(x) x.Ratio,WBdist),'b','Cluster #','nu','WithinClust/BetweenClust dist ratio','p_line_style','--');
                    subplot(224) % plot mean and STD of sim matrix
                    obj.FigParams.PlotMeanStd(CWC_MotifNum,CWC_X_mean,CWC_X_var,'Cluster #','Corr','b',1,'Avg similarity')
                    CWC_Figs=[CWC_Figs CWC_statsfig];
                end
                %% plot similarity matrix
                SimMatFig=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]);
                obj.PlotClusterSimilarityMatrix(X,clust_ind);
                %% reorder to original for tsne and PCA;
                TSNEMotifsDTW=[];PCAMotifs=[];
                %             TSNEMotifsDTW.Val2D= TSNEMotifsDTW.Val2D(clust_ind_notAligned,:);
                %             TSNEMotifsDTW.Val3D= TSNEMotifsDTW.Val3D(clust_ind_notAligned,:);
                %             PCAMotifs=PCAMotifs(clust_ind_notAligned,:);
                %% plot examples of clusters
                ClustExmplFig=cell(1,2+length(CoreMotifs));
                [ClustExmplFig{:}]=obj.PlotMotifClusts(Motifs_Aligned,TemplateMotif,CoreMotifs,clust_ind_Aligned,[],[]);
                %% save figures
                [~,~,ClusterQualityFigName]=obj.ManData.GetFileName([],[ExtraStr 'ClustQual_' TargFile FigFileStr],'SaveInResults',1,'WantedDate','ALL');
                obj.FigParams.SaveFigSeries([],ClusterQualityFigName,[DendroFig CoreMotifsFig SimMatFig ClustExmplFig CWC_Figs ])
                varargout=[CoreMotifsFig SimMatFig ClustExmplFig CWC_Figs];
            end
            %% save core motifs and clusting results
            save(XCorrFileName,'clust_ind','CoreMotifs','-append')
        end
        function MotifsSigTest=HierarchicalHypothesisTestingMotifs(obj,Motifs,SizeW,Z,X,varargin)
           global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %% fit these motifs to some channels how they can explain the data
            opts.MaxTime=AnalysisOpts.MotifAnalysis.MaxTime; % fit to 60 seconds of data
            opts.NSampCh=1; % how many example channel from each area
            opts.MaxClsutSet=[2:30];
            RecData=load('RecDateInfo.mat');
            MotifsSigTest=[];
            for MC=opts.MaxClsutSet
                ChCnt=1;
                clust_indSweep = cluster(Z,'maxclust',MC );
                [CoreMotifsSweep]=obj.FindTemplateMotif(clust_indSweep,Motifs,X,'CoreMotifsPerc',1);
                
                for Rec=AnalysisOpts.PairNum
                    % get two channels from each recording
                    FrontalChs =randsample(find(RecData.ChannelArea{Rec}==1 | RecData.ChannelArea{Rec}==4),opts.NSampCh)';
                    ParietalChs=randsample(find(RecData.ChannelArea{Rec}==5),opts.NSampCh)';
                    ChSet=[FrontalChs ParietalChs];
                    
                    for Ch=ChSet
                        tic
                        fprintf('\n Refitting Motifs ... Ch %i',ChCnt)                        
                        [MotifsSigTest.MotifSpecs_Refit{ChCnt,MC==opts.MaxClsutSet}]=obj.RefitMotifsGeneral(Ch,...
                            ['18' AnalysisOpts.DateSet{Rec}],CoreMotifsSweep,opts.MaxTime,SizeW);
                        MotifsSigTest.PEV(ChCnt,MC==opts.MaxClsutSet)=MotifsSigTest.MotifSpecs_Refit{ChCnt,MC==opts.MaxClsutSet}.PEV;
                        MotifsSigTest.PEVind{ChCnt,MC==opts.MaxClsutSet}=MotifsSigTest.MotifSpecs_Refit{ChCnt,MC==opts.MaxClsutSet}.PEVind;
                        ChCnt=ChCnt+1;
                        toc
                    end                 
                end
            end
        end
        function varargout=PlotHierarchicalHypothesisTestingMotifs(obj,PEV)
           global AnalysisOpts  AnalysisData
           NumClusts=size(PEV,2);
           NumChs=size(PEV,1);
           ClustNum=2:30;
           
           ClustPairs=[1:(NumClusts-1);2:NumClusts];
           [sig,pval]=arrayfun(@(x) ttest(PEV(:,ClustPairs(1,x)),PEV(:,ClustPairs(2,end))),1:(NumClusts-1));
            
           FigParams.PlotMeanStd(ClustNum,PEV,[],'Number of clusters','PEV',1,1,'Hierarchical Clustering of Motifs')
           ylim([0 1.1])
           arrayfun(@(x) AddSig(pval(x),[ClustPairs(1,x)+1 ClustPairs(1,x)+1.5 0.9],1),1:(NumClusts-1),'UniformOutput',0)
        end
        function [varargout]=PlotDendrogramwithMotifs(obj,Z,cutoff,clust_ind,CoreMotifs,varargin)
            fprintf('\nPlotDendrogramwithMotifs: Plotting Dendrograms...')
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % Z is calculated  from linkage
            % cutoff is where we have the clusters
            opts.ImgSiz=0.07; % size of images in pixels
            varargout(1)=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]);
            [H,T,outperm]=dendrogram(Z,0,'ColorThreshold',cutoff);
            hold on
            AxisInfo=axis;
            AxesInfo=get(gca,'Position');
            % create a loop where you plot core motif on the cut off
            Nclust=length(unique(clust_ind));
            ClustSeq=clust_ind(outperm); % sequence of plotted clusters on the plot
            
            for cl=1:Nclust
                % find center of this cluster
                StrInd=find(ClustSeq==cl,1,'first');
                StpInd=find(ClustSeq==cl,1,'last');
                MidIndX=mean([StrInd StpInd])*AxesInfo(3)/AxisInfo(2);
                MidIndY=cutoff*AxesInfo(4)/AxisInfo(4);
                % get position and plot core motif in there
                if Nclust>10
                    VerShift=mod(cl,7);
                else
                    VerShift=1;
                end
                Pos=[AxesInfo(1)+MidIndX-opts.ImgSiz/2,AxesInfo(2)+VerShift*opts.ImgSiz+0.05+MidIndY-opts.ImgSiz/2,...
                    opts.ImgSiz,opts.ImgSiz];
                axes('position',Pos)
                obj.PlotCoreMotifs(CoreMotifs,cl,1,'FigTitle',num2str(cl));
                axis off
                axis square
            end
            
        end
        function [varargout]=PlotClusterQualityDTW(obj,XCorrFileName,varargin) % uses DTW distance function
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            % load this file's DistMat and Motifs
            load(XCorrFileName,'Motifs','DistMatDTW','SizeW','SampledMotifsInd');
            % X has to be distance metric computed by dtw(dynamic time
            % warping
            X=eye(size(DistMatDTW));X=X+triu(DistMatDTW,1)+triu(DistMatDTW,1)'; % X is full distance matrix
            % remove inf and nan values from x and set them to zero
            X(isnan(X(:)) | isinf(X(:)))=0;
            AnalysisData.SizeW=SizeW;
            AnalysisData.MotifSimilarityMat=X;
            %% choose a random set of motifs to cluster with phenograph at first
            Motifs=Motifs(SampledMotifsInd,:); % this is now correspoding indexing to X
            %% take TSNe and PCA of Motifs first and save them off
            %             if ~obj.ManData.IsVarExistinFile(XCorrFileName,'TSNEMotifsDTW') | obj.OverWrite
            %                 % calculate tsne of motifs in 2D and 3D
            %                 % [TSNEMotifsDTW.Val2D, TSNEMotifsDTW.loss2D]= tsne(Motifs,'Algorithm','barneshut','NumPCAComponents',0,'Distance',@TemporalCrossCorrDTWdist);
            %                 [TSNEMotifsDTW.Val3D, TSNEMotifsDTW.loss3D]= tsne(Motifs,'Algorithm','barneshut','NumPCAComponents',0,'Distance',@TemporalCrossCorrDTWdist,'NumDimensions',3);
            %                 TSNEMotifsDTW.Val2D=TSNEMotifsDTW.Val3D(:,1:2);
            %                 % calculate PCA of motifs
            %                 [~,PCAMotifs]=pca(Motifs);
            %                 save(XCorrFileName,'TSNEMotifsDTW','PCAMotifs','-append'); % save off this data (Take a long time to process)
            %             else
            %                 load(XCorrFileName,'TSNEMotifsDTW','PCAMotifs')
            %             end
            %% try different methods of clustring and save them off in different files
            % create a filename for each clustring method based off of the
            % XCorrFileName
            [PATH,TargFile]=fileparts(XCorrFileName);
            if strcmpi(obj.MotifClassficationMethod,'PhenoCluster')
                % define the filename to save data
                FileName=[PATH AnalysisOpts.FS TargFile '_' obj.MotifClassficationMethod '_k' num2str(AnalysisOpts.MotifAnalysis.Phenograph_K) '_' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '.mat'];
                FigFileStr=[obj.MotifClassficationMethod '_k' num2str(AnalysisOpts.MotifAnalysis.Phenograph_K) '_' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '_TM' AnalysisOpts.MotifAnalysis.TemplateMotif ];
                
                if ~obj.ManData.IsVarExistinFile(FileName,'clust_ind_PhenoClusterDTW') | obj.OverWrite
                    % perform phenophraph on similarity matrix
                    clust_ind_PhenoClusterDTW=PhenoCluster(X,'k',AnalysisOpts.MotifAnalysis.Phenograph_K); % cluster with phenograph
                    clust_ind=clust_ind_PhenoClusterDTW(:,1); % take the level 1 for now
                    obj.ManData.SaveVar([],clust_ind_PhenoClusterDTW,'clust_ind_PhenoClusterDTW',FileName,'SelfName',1);
                else
                    load(FileName,'clust_ind_PhenoClusterDTW')
                    clust_ind=clust_ind_PhenoClusterDTW(:,1); % take the level 1 for now
                end
            elseif strcmpi(obj.MotifClassficationMethod,'Hierarchical')
                % use heirarchical clustring to clsuter the data
                % define the filename to save data
                FileName=[PATH AnalysisOpts.FS TargFile '_' obj.MotifClassficationMethod '_' num2str(AnalysisOpts.MotifAnalysis.MaxClust) '_' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '.mat'];
                FigFileStr=[obj.MotifClassficationMethod '_MC' num2str(AnalysisOpts.MotifAnalysis.MaxClust) '_' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '_TM' AnalysisOpts.MotifAnalysis.TemplateMotif ];
                % convert similarity matrix to dissimilarity matrix and
                % then into pdist form
                MaxClust=AnalysisOpts.MotifAnalysis.MaxClust;
                dissimMat=1-X;
                yOut = squareform(dissimMat,'tovector');
                Z = linkage(yOut,'ward');%'method' can be 'ward' or 'average'
                clust_ind = cluster(Z,'maxclust',MaxClust );
                DendroFig=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]);
                cutoff = median([Z(end-MaxClust+1,3) Z(end-MaxClust+2,3)]);
                [H,T,outperm]=dendrogram(Z,0,'ColorThreshold',cutoff);
                
            elseif strcmpi(obj.MotifClassficationMethod,'DBSCAN')
            elseif strcmpi(obj.MotifClassficationMethod,'SOM')
                % perform SOM on similarity matrix and then Phenograph on
                % that
                [clust_ind]=SOMCluster(X,'method','matlab','verbose',0,'clustringMethod',AnalysisOpts.MotifAnalysis.SOMClustringMethod);
                
            elseif strcmpi(obj.MotifClassficationMethod,'EncoderPhenograph')
                % build autoencoder, find features and feed them to
                % clusterign algorithm
                % first align motifs to a random motif
                
                if ~obj.ManData.IsVarExistinFile(FileName,'clust_ind_EncoderPhenograph') | obj.OverWrite
                    %                     RandAlightMotifInd=randsample(size(Motifs,1),1);
                    %                     Motifs_AlignedRand=obj.AlignMotifs(Motifs(RandAlightMotifInd,:),Motifs);
                    AutoEncFeature=BuildAutoEncoderMotifs(Motifs,obj.ClassificaitonAutoEncLayerSiz);
                    AutoEncFeaturesDTW=AutoEncFeature{obj.ClassificaitonAutoEncTargLayer+1}';
                    [clust_ind_EncoderPhenographDTW,ovr_Q] = Phenograph(AutoEncFeatures,'DistanceMetric','euclidean'); % phenograph clustering
                    
                    obj.ManData.SaveVar([],clust_ind_EncoderPhenographDTW,'clust_ind_EncoderPhenographDTW',FileName,'SelfName',1);
                    obj.ManData.SaveVar([],AutoEncFeaturesDTW,'AutoEncFeaturesDTW',FileName,'SelfName',1);
                    
                else
                    load(FileName,'clust_ind_EncoderPhenographDTW')
                end
                clust_ind=clust_ind_EncoderPhenographDTW(:,1); % take the level 1 for now
            end
            % find template motif
            [CoreMotifs,TemplateMotif,Motifs_Aligned,clust_ind_Aligned,clust_ind_notAligned,TestCoreMotifs]=obj.FindTemplateMotif(clust_ind,Motifs,'CoreMotifsPerc',1);
            % compute tsne of aligned motifs
            %             if ~obj.ManData.IsVarExistinFile(FileName,'TsneMotifs_Aligned') | obj.OverWrite
            %                 fprintf('\nCalculating Tsne...')
            %                 TsneMotifs_Aligned=tsne(Motifs_Aligned,'Algorithm','exact','NumPCAComponents',0,'Distance','correlation');
            %                 [~,PCAMotifs_Aligned]=pca(Motifs_Aligned);
            %                 obj.ManData.SaveVar([],TsneMotifs_Aligned,'TsneMotifs_Aligned',FileName,'SelfName',1);
            %                 obj.ManData.SaveVar([],PCAMotifs_Aligned,'PCAMotifs_Aligned',FileName,'SelfName',1);
            %
            %              else
            %                 load(FileName,'TsneMotifs_Aligned','PCAMotifs_Aligned')
            %             end
            % plot Core motifs
            CoreMotifsFig=cell(1,3);
            [CoreMotifsFig(1:3)]=obj.FigParams.RenderFigure(3,[0 0 0.8 0.8]);
            figure(CoreMotifsFig{1});CoreMotifsFig{1}=obj.PlotCoreMotifs(CoreMotifs,[],'FigTitle',[num2str(obj.CoreMotifsPerc*100) '%']);
            figure(CoreMotifsFig{2});CoreMotifsFig{2}=obj.PlotCoreMotifs(TestCoreMotifs.P01,[],'FigTitle','10%'); % plot core motifs with 10 percent
            figure(CoreMotifsFig{3});CoreMotifsFig{3}=obj.PlotCoreMotifs(TestCoreMotifs.P1,[],'FigTitle','100%'); % plot core motifs with 100 percent
            % plot similarity matrix
            SimMatFig{1}=obj.PlotClusterSimilarityMatrix(X,clust_ind);
            % plot examples of clusters
            ClustExmplFig=cell(1,2+length(CoreMotifs));
            % reorder to original for tsne and PCA;
            TSNEMotifsDTW=[];PCAMotifs=[];
            %             TSNEMotifsDTW.Val2D= TSNEMotifsDTW.Val2D(clust_ind_notAligned,:);
            %             TSNEMotifsDTW.Val3D= TSNEMotifsDTW.Val3D(clust_ind_notAligned,:);
            %             PCAMotifs=PCAMotifs(clust_ind_notAligned,:);
            %   [ClustExmplFig{:}]=obj.PlotMotifClusts(Motifs_Aligned,TemplateMotif,CoreMotifs,clust_ind_Aligned,TSNEMotifsDTW,PCAMotifs);
            %% now see if there are different clusters inside some of the clusters we are seeing CWC(cluster within cluster)
            CWC_Figs=[];
            CWC_MotifNum=[1:14];
            ovr_q=nan*ones(1,length(CWC_MotifNum));
            for n_CWC=[CWC_MotifNum]
                CWC_Fig=cell(1,3);
                CWCind=clust_ind==n_CWC; % find motifs in this cluster
                CWC_X=X(CWCind,CWCind);
                CWC_Motifs=Motifs(CWCind,:);
                [ovr_q(n_CWC),CWC_Fig{:}]=obj.ClusterWithinCluster(CWC_X,CWC_Motifs,n_CWC);
                CWC_Figs=[CWC_Figs CWC_Fig];
            end
            %%
            [~,~,ClusterQualityFigName]=obj.ManData.GetFileName([],['ClustQual_' TargFile FigFileStr],'SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.SaveFigSeries([],ClusterQualityFigName,[CoreMotifsFig SimMatFig ClustExmplFig CWC_Figs ])
            varargout=[CoreMotifsFig SimMatFig ClustExmplFig CWC_Figs];
        end
        function [ovr_q,Entropy_CWC_X,WBdist,CoreMotifs,varargout]=ClusterWithinCluster(obj,CWC_X,CWC_Motifs,CWC_MotifNum,sp,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            if nargin <5; sp = [1,1,1]; end
            
            [clust_ind_PhenoCluster,~,ovr_q]=PhenoCluster(CWC_X,'k',AnalysisOpts.MotifAnalysis.Phenograph_K); % cluster with phenograph
            %[clust_ind_PhenoCluster,ovr_q]=community_louvain(CWC_X,1,[],'negative_sym');
            clust_ind=clust_ind_PhenoCluster(:,1); % take the level 1 for now
            [CoreMotifs]=obj.FindTemplateMotif(clust_ind,CWC_Motifs,'CoreMotifsPerc',1);
            [varargout(1)]=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]);
            figure(varargout{1})
            obj.PlotCoreMotifs(CoreMotifs,[],'FigTitle',['CWC ' num2str(CWC_MotifNum)]);
            % plot similarity matrix
            obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]);
            subplot(121)
            [varargout{2}]=obj.PlotClusterSimilarityMatrix(CWC_X,clust_ind);
            title(['Sim Matrix for CWC ' num2str(CWC_MotifNum)])
            % plot entropy
            subplot(122)
            imhist(CWC_X);
            axis square
            Entropy_CWC_X=entropy(CWC_X);
            title(sprintf('Entropy %i'),Entropy_CWC_X);
            % calculate between and within cluster distance
            WBdist=obj.CalWithinBetweenClsutDist(CWC_X,clust_ind);
        end
        function WBdist=CalWithinBetweenClsutDist(obj,X,clust_ind,varargin)
            % calculates within cluster distance, between cluster distance
            % and their ratios
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            NClust=length(unique(clust_ind));
            [sortClust,SortClsutInd]=sort(clust_ind);
            Xsort=X(SortClsutInd,SortClsutInd);
            sizXsort=size(Xsort);
            [ClustIndsStr]=arrayfun(@(x) find(sortClust==x,1,'first'),1:NClust); %find first inds
            [ClustIndsStp]=arrayfun(@(x) find(sortClust==x,1,'last'),1:NClust); %find first inds
            % loop in clusters and calculate between and within clust
            Xsort_upper=tril(nan*ones(size(Xsort)))+triu(Xsort);
            WithinDist=[];BetwinDist=[];
            for cl=1:NClust
                ThisWithinDist= Xsort_upper(ClustIndsStr(cl):ClustIndsStp(cl),ClustIndsStr(cl):ClustIndsStp(cl));
                WithinDist=[WithinDist ThisWithinDist(:)'];
                BetwinDistInd=setdiff(1:sizXsort(1),ClustIndsStr(cl):ClustIndsStp(cl));
                ThisBetweenDist=Xsort_upper(ClustIndsStr(cl):ClustIndsStp(cl),BetwinDistInd);
                BetwinDist=[BetwinDist  ThisBetweenDist(:)'];
            end
            WBdist.WithinDistMean=nanmean(WithinDist(:));
            WBdist.BetwinDistMean=nanmean(BetwinDist(:));
            WBdist.Ratio=WBdist.BetwinDistMean/WBdist.WithinDistMean;
        end
        function varargout=PlotClusterSimilarityMatrix(obj,X,clust_ind,varargin)
            fprintf('\nPlotClusterSimilarityMatrix: Plotting Similarity Matrix ...')
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            %% now plot the quality metrics for these motifs
            % plot similarity matrix first
            NClust=length(unique(clust_ind));
            ClustCol=distinguishable_colors(NClust);
            [sortClust,SortClsutInd]=sort(clust_ind);
            Xsort=X(SortClsutInd,SortClsutInd);
            [ClustIndsStr]=arrayfun(@(x) find(sortClust==x,1,'first'),1:NClust); %find first inds
            [ClustIndsStp]=arrayfun(@(x) find(sortClust==x,1,'last'),1:NClust); %find first inds
            ClustIndsMid=mean([ClustIndsStr;ClustIndsStp],1);
            %% plot
            
            hold on
            imagesc(Xsort) % plot the imagesc
            set(gca,'YDir','reverse')
            colorbar
            % plot lines with color horizontal
            arrayfun(@(x) line([ClustIndsStr(x) ClustIndsStp(x)],[-100 -100],'color',ClustCol(x,:),...
                'LineWidth',5),1:NClust)
            arrayfun(@(x) text([ClustIndsMid(x) ClustIndsMid(x)],[-500 -500],num2str(x),...
                'FontSize',10),1:NClust)
            % plot lines with color vertical
            arrayfun(@(x) line([-100 -100],[ClustIndsStr(x) ClustIndsStp(x)],'color',ClustCol(x,:),...
                'LineWidth',5),1:NClust)
            arrayfun(@(x) text([-500 -500],[ClustIndsMid(x) ClustIndsMid(x)],num2str(x),...
                'FontSize',10),1:NClust)
            axis equal
            axis off
            
            varargout{1}=gcf;
            obj.FigParams.FormatAxes(gca);
            
        end
        function  varargout=PlotMotifClusts(obj,Motifs,TemplateMotif,CoreMotifs,clust_ind,TsneX,PcaX,varargin)   %plots clusters for motifs
            fprintf('\nPlotMotifClusts: Plotting Motif Clusters and Example Motifs ...')
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %% setup parameters
            Nclust=length(unique(clust_ind));
            MotifCol=distinguishable_colors(Nclust);
            Ts=1/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            NexmMotif=14; %Number of example motifs
            Ncol=NexmMotif+1;
            NtimCoreMtf=size(CoreMotifs{1},2);
            NrowInd=4; %indidual motifs with example
            NcolInd=4;  % individual motifs with example
            
            varargout(1:2+Nclust)=obj.FigParams.RenderFigure(Nclust+2,[0 0 0.8 0.8]); % figures
            
            Legends=[];
            %%% plt scatter in tsne space
            figure(varargout{1});
            if ~isempty(TsneX)
                for i=1:Nclust
                    ind=find(clust_ind==i);
                    subplot(131); hold on % plot Tsne here
                    plt1=arrayfun(@(x) plot(TsneX.Val2D(x,1),TsneX.Val2D(x,2),...
                        '.','color',MotifCol(i,:),'MarkerSize',7),ind,'UniformOutput',0);
                    
                    subplot(132); hold on % plot Tsne in 3 dim here
                    arrayfun(@(x) plot3(TsneX.Val3D(x,1),TsneX.Val3D(x,2),TsneX.Val3D(x,3),...
                        '.','color',MotifCol(i,:),'MarkerSize',7),ind,'UniformOutput',0);
                    
                    subplot(133);hold on % plot PCA here
                    plt2=arrayfun(@(x) plot(PcaX(x,1),PcaX(x,2),'.',...
                        'color',MotifCol(i,:),'MarkerSize',7),ind,'UniformOutput',0);
                    Legends=[Legends plt1{1}];
                    
                end
                subplot(131)
                xlabel('TSNE1');ylabel('TSNE2');axis square
                subplot(132)
                xlabel('TSNE1');ylabel('TSNE2');ylabel('TSNE3');axis square
                subplot(133)
                xlabel('PCA1');ylabel('PCA2');axis square
                legend(Legends,arrayfun(@(x) ['Clust' num2str(x)] ,1:Nclust,'UniformOutput',0),...
                    'Location','bestoutside')
                title('clustering of motifs')
            end
            
            %% plot example motifs in one plot for all and seperate plots for each motif
            k=1;
            for i=1:Nclust%SortOccDistInd
                ind=find(clust_ind==i);
                
                % plot core motif
                %              figure(varargout{3});
                %              subplot(Nclust,Ncol,(k-1)*Ncol+1);hold on
                %              helperCWTTimeFreqPlot(CoreMotifs{i},(1:NtimCoreMtf)*Ts,AnalysisData.cwt_f,'justplot1',['Core mtf clust' num2str(i)],'Time','Frequency(Hz)',obj.LogScale)
                %              axis off
                
                figure(varargout{1+i});
                subplot(NrowInd,NcolInd,1);hold on
                helperCWTTimeFreqPlot(CoreMotifs{i},(1:NtimCoreMtf)*Ts,AnalysisData.cwt_f,'justplot1',['Core mtf clust' num2str(i)],'Time','Frequency(Hz)',obj.LogScale)
                axis square
                
                % plot template motif
                figure(varargout{1+i});
                subplot(NrowInd,NcolInd,2);hold on
                helperCWTTimeFreqPlot(TemplateMotif{i},(1:NtimCoreMtf)*Ts,AnalysisData.cwt_f,'justplot1',['Template mtf clust' num2str(i)],'Time','Frequency(Hz)',obj.LogScale)
                axis square
                % plot some example motifs
                ExmpMotfInd=randsample(ind,NexmMotif);
                for ne=1:NexmMotif
                    %                  figure(varargout{3});subplot(Nclust,Ncol,(k-1)*Ncol+ne+1);hold on
                    ThisExampleMotif=reshape(Motifs(ExmpMotfInd(ne),:),AnalysisData.SizeW(1),AnalysisData.SizeW(2)*3);
                    %                  helperCWTTimeFreqPlot(ThisExampleMotif,(1:3*AnalysisData.SizeW(2))*Ts,AnalysisData.cwt_f,'justplot1',['exm mtf' num2str(ne)],'Time','Frequency(Hz)',obj.LogScale)
                    %                  axis off
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    figure(varargout{1+i});subplot(NrowInd,NcolInd,2+ne);hold on
                    helperCWTTimeFreqPlot(ThisExampleMotif,(1:3*AnalysisData.SizeW(2))*Ts,AnalysisData.cwt_f,'justplot1',['exm mtf' num2str(ne)],'Time','Frequency(Hz)',obj.LogScale)
                    axis square
                end
                k=k+1;
            end
            %% plot distribution of power and freq of occurance of each cluster
            figure(varargout{2+Nclust})
            subplot(121);hold on  % plot dist of power
            arrayfun(@(x) plot(AnalysisData.cwt_f,mean(CoreMotifs{x},2),'linewidth',5,...
                'color',MotifCol(x,:)),1:length(CoreMotifs))
            xlabel('freq');ylabel('Time Avg Power')
            %% plot freq of occurance
            subplot(122);hold on  % plot dist of power
            [OccDist]=hist(clust_ind,1:Nclust); %dist of occurance
            OccDist=OccDist/sum(OccDist);
            arrayfun(@(x) bar(x,OccDist(x),'facecolor',MotifCol(x,:)),1:Nclust)
            axis tight
            xlabel('Cluster #')
            ylabel('prob of occurence')
            
        end
        
        function [clust_ind_all]=GeneralizeMotifClustering(obj,Motifs,RandMotifInd,clust_ind,CoreMotifs,SizeW,varargin)
            % Generalizes clustered motifs to withhelp data, helps
            % increasing speed of clustring
            %inputs:
            %Motifs: all of the motifs
            %NRandMotif: set of motifs that are used in clustring
            %clust_ind:   cluster ind of random motifs
            %CoreMotifs: cell of core  motifs
            % SizeW: SizeW
            % different methods can be used eg.g prototype or maximum
            % likelihood, we will implemet that later
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            % prepare some vars
            NMotifs=size(Motifs,1);
            clust_ind_all=nan*ones(NMotifs,1);
            clust_ind_all(RandMotifInd)=clust_ind;
            NRemainMotifs=setdiff(1:NMotifs,RandMotifInd);
            RemainMotifs=Motifs(NRemainMotifs,:);
            % calculate distance of remaining motifs to each core motif remaining motifs
            [RemainDist] = cellfun(@(x) TemporalCrossCorrDistance1D(x,RemainMotifs,SizeW,0),CoreMotifs,'UniformOutput',0);
            RemainDist=cell2mat(RemainDist);
            
            [~,clust_ind_remain]=max(RemainDist,[],2);
            clust_ind_all(NRemainMotifs)=clust_ind_remain;
        end
        
        
        function  TemporalXCorrTensor_onClust(obj,Motifs,PairNum,XCorrFileName,SizeW,varargin) % calcultes X by X atrix of temporal xcorr of motifs runs on cluster
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % check if we already have this file so just skip it
            [PATH,TargFile]=fileparts(XCorrFileName);
            SaveFileName=[PATH AnalysisOpts.FS TargFile '_'  num2str(PairNum) AnalysisOpts.MotifAnalysis.XCorrTensorMethod '.mat'];
            if exist(SaveFileName,'file') & ~obj.OverWrite
                fprintf('\n File for pair %i exists...',PairNum)
                return
            end
            
            fprintf('\nComputing Distances...')
            tic
            %%  first check if we have the random set of motifs we want to look at
            if ~obj.ManData.IsVarExistinFile(XCorrFileName,'MotifPairChunk')
                Nmotifs=AnalysisOpts.MotifAnalysis.NSampledMotifs; % we randomly take X motifs from this pool and do all of the analaysis on it
                N_orgMotifs=size(Motifs,1); % number of original motifs
                SampledMotifsInd=randsample(1:N_orgMotifs,Nmotifs);
                MotifPairs=nchoosek(1:Nmotifs,2);  % take pairs of motifs
                NMotifPairs=size(MotifPairs,1);
                NSteps=AnalysisOpts.MotifAnalysis.NchunksXCorrTensor+1;
                %% chop-off Motif pairs in NSteps pieces
                Steps=NMotifPairs/NSteps;
                MotifPairs_Step=floor([1:Steps:NMotifPairs]); MotifPairs_Step=[MotifPairs_Step(1:end-1) NMotifPairs];
                MotifPairChunk=[MotifPairs_Step(1:end-1)' MotifPairs_Step(2:end)'];
                MotifPairChunk(2:end,1)=MotifPairChunk(2:end,1)+1;
                save(XCorrFileName,'MotifPairChunk','SampledMotifsInd','MotifPairs','-append')
            else
                load(XCorrFileName,'MotifPairChunk','SampledMotifsInd','MotifPairs');
            end
            if ~isempty(PairNum)
                SampledMotifs=Motifs(SampledMotifsInd,:);
                ThisPair=MotifPairChunk(PairNum,:);
                NThisPair=ThisPair(2)-ThisPair(1)+1;
                k=1;
                for x=ThisPair(1):ThisPair(2)
                    i=MotifPairs(x,1);j=MotifPairs(x,2);
                    if strcmpi(AnalysisOpts.MotifAnalysis.XCorrTensorMethod,'DTW') % if we are using dtw
                        eval(['[DistMat_' num2str(PairNum) '(' num2str(k) '),ShiftMat_' num2str(PairNum) '(' num2str(k) '),EuDistMat_' num2str(PairNum) '(' num2str(k) ')]= TemporalCrossCorrDistance1DParLoopDTW(i,j,SampledMotifs,SizeW);']);
                    elseif strcmpi(AnalysisOpts.MotifAnalysis.XCorrTensorMethod,'')
                        eval(['[DistMat_' num2str(PairNum) '(' num2str(k) '),ShiftMat_' num2str(PairNum) '(' num2str(k) '),EuDistMat_' num2str(PairNum) '(' num2str(k) ')]= TemporalCrossCorrDistance1DParLoop(i,j,SampledMotifs,SizeW);']);
                    end
                    k=k+1;
                end
                
                if AnalysisOpts.SaveData
                    
                    save(SaveFileName,['DistMat_' num2str(PairNum)],['ShiftMat_' num2str(PairNum)],['EuDistMat_' num2str(PairNum)],'-v7.3')
                    
                    % matfile route
                    % MatFile=matfile(XCorrFileName,'Writable',true); % connect to file as matfile
                    % eval(['MatFile.DistMat_'  num2str(PairNum) '=DistMat_' num2str(PairNum) ';']);
                    % eval(['MatFile.ShiftMat_' num2str(PairNum) '=ShiftMat_' num2str(PairNum) ';']);
                    % clear MatFile
                    
                    fprintf('\n finished calculating distances ...')
                end
                toc
            end
        end
        function  ReshapeXCorrArea_Old(obj,XCorrFileName,varargin) % OLD version reshapes all of the chunks of pairs calculated in TemporalXCorrTensor_onClust
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            load(XCorrFileName)
            Nmotifs=size(Motifs,1);clearvars Motifs
            MotifPairs=nchoosek(1:Nmotifs,2);  % take pairs of motifs
            NMotifPairs=size(MotifPairs,1);
            NSteps=AnalysisOpts.NchunksXCorrTensor+1;
            %% chop-off Motif pairs in NSteps pieces
            Steps=NMotifPairs/NSteps;
            MotifPairs_Step=floor([1:Steps:NMotifPairs]); MotifPairs_Step=[MotifPairs_Step(1:end-1) NMotifPairs];
            MotifPairChunk=[MotifPairs_Step(1:end-1)' MotifPairs_Step(2:end)'];
            %  MotifPairChunk(2:end,1)=MotifPairChunk(2:end,1)+1;
            %% concatinate all of the distmats that are in the workspace
            DistMat=[];ShiftMat=[];
            for m=1:AnalysisOpts.NchunksXCorrTensor
                ThisPair=MotifPairChunk(m,:);
                x=ThisPair(1):ThisPair(2);
                NThisPair=length(x);
                if exist(['DistMat_' num2str(m)],'var')
                    if m~=1
                        eval(['DistMat=cat(2,DistMat,DistMat_' num2str(m) '(2:end));']);
                        eval(['ShiftMat=cat(2,ShiftMat,ShiftMat_' num2str(m) '(2:end));']);
                        eval(['IQ(m)=double(length(DistMat_' num2str(m) ')== NThisPair);']);
                    else
                        eval(['DistMat=cat(2,DistMat,DistMat_' num2str(m) ');']);
                        eval(['ShiftMat=cat(2,ShiftMat,ShiftMat_' num2str(m) ');']);
                        eval(['IQ(m)=double(length(DistMat_' num2str(m) ')== NThisPair);']);
                    end
                    eval(['clearvars DistMat_' num2str(m)]);
                else % if the file doesn't yet ext replace it with NAN for now
                    IQ(m)=NaN;
                    NanMat   =nan*ones(1,NThisPair-1);
                    DistMat  =cat(2,DistMat,NanMat);
                    ShiftMat =cat(2,ShiftMat,NanMat);
                end
            end
            %           m=1;
            %             while exist(['DistMat_' num2str(m)],'var')
            %                 eval(['DistMat=cat(2,DistMat,DistMat_' num2str(m) ');']);
            %                 eval(['ShiftMat=cat(2,ShiftMat,ShiftMat_' num2str(m) ');']);
            %                 m=m+1;
            %             end
            %%% put them in similarity matrix now
            DistMat=obj.ManData.ReshapeSquareMatrix(MotifPairs,DistMat);
            ShiftMat=obj.ManData.ReshapeSquareMatrix(MotifPairs,ShiftMat);
            % append them to the file we ahave
            save(XCorrFileName,'DistMat','ShiftMat','-append')
            
        end
        function  ReshapeXCorrArea(obj,XCorrFileName,varargin) % reshapes all of the chunks of pairs calculated in TemporalXCorrTensor_onClust
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % define the name of the DistMat based on the method used to
            % calculate it
            DistMatVarName=['DistMat' AnalysisOpts.MotifAnalysis.XCorrTensorMethod];% xcorr distance
            EuDistMatVarName=['EuDistMat' AnalysisOpts.MotifAnalysis.XCorrTensorMethod]; % correlation distance
            
            if obj.ManData.IsVarExistinFile(XCorrFileName,DistMatVarName) & obj.ManData.IsVarExistinFile(XCorrFileName,EuDistMatVarName) & ~obj.OverWrite
                return
            end
            [PATH,TargFile]=fileparts(XCorrFileName);
            %% check if all of these files exist otherwise throw an error
            FileNameFunc=@(m) [PATH AnalysisOpts.FS TargFile '_'  num2str(m) AnalysisOpts.MotifAnalysis.XCorrTensorMethod '.mat'];
            ExistXcorrs=arrayfun(@(x) exist(FileNameFunc(x),'file'),1:AnalysisOpts.MotifAnalysis.NchunksXCorrTensor);
            if find(~ExistXcorrs);error('There are unprocesses pair chunks ...');end
            %% concatinate all of the distmats that are in the workspace
            DistMat=[];EuDistMat=[];ShiftMat=[];
            for m=1:AnalysisOpts.MotifAnalysis.NchunksXCorrTensor
                FileName=[PATH AnalysisOpts.FS TargFile '_'  num2str(m) AnalysisOpts.MotifAnalysis.XCorrTensorMethod '.mat'];
                
                if exist(FileName,'file')
                    load(FileName);
                    
                    eval(['DistMat=cat(2,DistMat,DistMat_' num2str(m) ');']);
                    eval(['EuDistMat=cat(2,EuDistMat,EuDistMat_' num2str(m) ');']);
                    eval(['ShiftMat=cat(2,ShiftMat,ShiftMat_' num2str(m) ');']);
                    eval(['clearvars DistMat_' num2str(m) ' EuDistMat_' num2str(m) ' ShiftMat_' num2str(m)]);
                    IQ(m)=1;
                else % if the file doesn't yet ext replace it with NAN for now
                    warning(['File for pairchunk ' num2str(m) ' does not exist ... replacing with NAN']);
                    IQ(m)=NaN;
                    NanMat   =nan*ones(1,NThisPair-1);
                    DistMat  =cat(2,DistMat,NanMat);
                    EuDistMat=cat(2,EuDistMat,NanMat);
                    ShiftMat =cat(2,ShiftMat,NanMat);
                end
            end
            %%% put them in similarity matrix now
            load(XCorrFileName,'MotifPairs')
            DistMat=obj.ManData.ReshapeSquareMatrix(MotifPairs,DistMat);
            EuDistMat=obj.ManData.ReshapeSquareMatrix(MotifPairs,EuDistMat);
            ShiftMat=obj.ManData.ReshapeSquareMatrix(MotifPairs,ShiftMat);
            % copy this to a new varaible according to the method we use to
            % calculate it
            eval(['DistMat' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '=DistMat;']);
            eval(['EuDistMat' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '=EuDistMat;']);
            eval(['ShiftMat' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '=ShiftMat;']);
            eval(['IQ' AnalysisOpts.MotifAnalysis.XCorrTensorMethod '=IQ;']);
            
            % append them to the file we ahave
            % save(XCorrFileName,'DistMat','ShiftMat','IQ','-append')
            txt=sprintf('save(XCorrFileName,''DistMat%s'',''EuDistMat%s'',''ShiftMat%s'',''IQ%s'',''-append'')',...
                AnalysisOpts.MotifAnalysis.XCorrTensorMethod,AnalysisOpts.MotifAnalysis.XCorrTensorMethod,...
                AnalysisOpts.MotifAnalysis.XCorrTensorMethod,AnalysisOpts.MotifAnalysis.XCorrTensorMethod);
            eval(txt);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function  varargout=ClassifyMotifs(obj,MotifData,varargin)   % classifies motifs with classification algorithm
            % data is the W metrix for each neuron stacked in 3 Dim
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % normalize each otif to 0 and 1
            % take the Ws and put them into a huge 3 dim matrix
            Ws=arrayfun(@(x) MotifData{AnalysisData.Identifier(x,1)}.AnalysisData.MotifSpecs_Disc{AnalysisData.Identifier(x,3)}.W,...
                1:size(AnalysisData.Identifier,1),'UniformOutput',0);
            W=obj.ConcatinateWs(Ws);
            % take the empty Ws
            SumW=sum(W,1);GoodW=SumW~=0;
            W=W(:,SumW~=0);
            AnalysisData.IdentifierH=NaN*ones(size(AnalysisData.IdentifierW,1),1);
            AnalysisData.IdentifierW=AnalysisData.IdentifierW(SumW~=0,:);
            AnalysisData.NormW=obj.ManData.Normalize_0_1(W');  % Normalize Ws
            %% get some params
            AnalysisData.f=MotifData{1}.AnalysisData.cwt_f;
            AnalysisData.SizeW=[size(Ws{1},1) size(Ws{1},3)];
            AnalysisData.COr=@(x,y) corr(x(:),y(:)); % setup correlation function
            %% cluster the data
            clustFig=obj.FigParams.RenderFigure(1,[]);
            %    if ~isfield(AnalysisData,'TsneX') % if we aready have this don't do it again
            if strcmpi(obj.MotifClassficationMethod,'phenograph')
                %% if we don't have distances saved compute them
                SNW=size(AnalysisData.NormW);
                if SNW(1)<SNW(2)
                    NormW=transpose(downsample(AnalysisData.NormW',ceil(SNW(2)/SNW(1))));
                    AnalysisData.TsneX=tsne(NormW);  % TSNE
                    [~,AnalysisData.PCA]=pca(NormW);
                else
                    AnalysisData.TsneX=tsne(AnalysisData.NormW);
                    [~,AnalysisData.PCA]=pca(AnalysisData.NormW);
                end
                [AnalysisData.clust_ind,ovr_Q] = Phenograph(AnalysisData.NormW, zeros(size(AnalysisData.NormW,1),1),...
                    'DistanceMetric',@TemporalCrossCorrDistance1D); % phenograph clustering
                
            elseif strcmpi(obj.MotifClassficationMethod,'PhenoCluster')
                % first compute distance between members
                XCorrFileName=['MotifsClust_Dist_' AnalysisOpts.Area_2look{1} '_Rec' AnalysisOpts.RecDate '.mat'];
                if ~exist([AnalysisOpts.DataSavePath XCorrFileName],'file')
                    obj.TemporalXCorrTensor(AnalysisData.NormW)
                end
                Distance=load([AnalysisOpts.DataSavePath XCorrFileName]);
                X=triu(Distance.DistMat)+triu(Distance.DistMat)';
                AnalysisData.MotifSimilarityMat=X;
                AnalysisData.clust_ind=PhenoCluster(X,'k',AnalysisOpts.Phenograph_K); % cluster with phenograph
                AnalysisData.clust_ind=AnalysisData.clust_ind(:,1); % take the level 1 for now
                
                
            elseif strcmpi(obj.MotifClassficationMethod,'SOM')
                % first compute distance between members
                [~,~,XCorrFileName]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisOpts.AnalysisPathName,AnalysisOpts.Animal,AnalysisOpts.RecDate,obj.ChNum,['_MotifDist' AnalysisOpts.ExtraParamTxtWrite],'ext','.mat');
                
                if ~exist(XCorrFileName,'file')
                    obj.TemporalXCorrTensor(AnalysisData.NormW);
                end
                ClustInfo=load(XCorrFileName);
                X=zeros(size(ClustInfo.DistMat));X=X+triu(ClustInfo.DistMat,1)+triu(ClustInfo.DistMat,1)';
                AnalysisData.MotifSimilarityMat=triu(ClustInfo.DistMat)+triu(ClustInfo.DistMat)';
                % perform clustring
                % [SOMClsut,AnalysisData.clust_ind,AnalysisData.TsneX,AnalysisData.PCA]=SOMCluster(X,'method','matlab','verbose',1,'clustringMethod',AnalysisOpts.SOMClustringMethod);
                if ~isfield(ClustInfo,'clust_ind_pheno')   %
                    [clust_ind,clust_ind_pheno]=SOMCluster(X,'method','matlab','verbose',0,'clustringMethod',AnalysisOpts.SOMClustringMethod);
                    save(XCorrFileName,'clust_ind','clust_ind_pheno','-append')
                    AnalysisData.clust_ind=clust_ind_pheno;
                else
                    AnalysisData.clust_ind=ClustInfo.clust_ind_pheno;
                end
                
            elseif strcmpi(obj.MotifClassficationMethod,'EncoderPhenograph')
                % build autoencoder, find features and feed them to
                % clusterign algorithm
                AutoEncFeature=BuildAutoEncoderMotifs(AnalysisData.NormW,obj.ClassificaitonAutoEncLayerSiz);
                Features=AutoEncFeature{obj.ClassificaitonAutoEncTargLayer+1}';
                [AnalysisData.clust_ind,ovr_Q] = Phenograph(Features, ...
                    zeros(size(Features,1),1),'DistanceMetric','euclidean'); % phenograph clustering
                % get the TSNE
                AnalysisData.TsneX=tsne(Features);  % TSNE
                [~,AnalysisData.PCA]=pca(Features);
                
            end
            %     end
            %% work on getting out the Hs now
            NCh=length(AnalysisData.Ch);
            NH=size(MotifData{1}.AnalysisData.MotifSpecs_Disc{1}.H,1); % how many Hz
            NtimH=size(MotifData{1}.AnalysisData.MotifSpecs_Disc{1}.H,2);
            NTimChunks=length(MotifData{1}.AnalysisData.MotifSpecs_Disc);
            Nclust=length(unique(AnalysisData.clust_ind));
            col=distinguishable_colors(Nclust);
            
            AnalysisData.IdentifierH(GoodW)=AnalysisData.clust_ind;
            AnalysisData.IdentifierH=arrayfun(@(x) reshape(AnalysisData.IdentifierH(NTimChunks*NH*(x-1)+1:NTimChunks*NH*x),NH,NTimChunks),1:NCh,'UniformOutput',0);
            
            %% find motif templates and align everything to that
            obj.FindTemplateMotif(AnalysisData.clust_ind,AnalysisData.NormW)
            AnalysisData.NormW=AnalysisData.NormW_Aligned;
            AnalysisData.clust_ind=AnalysisData.clust_ind_Aligned;
            
            %% take TSNE and PCA of Motifs
            if ~isfield(ClustInfo,'TsneX')   %
                TsneX=tsne(AnalysisData.NormW);
                [~,PCA]=pca(AnalysisData.NormW);
                save(XCorrFileName,'TsneX','PCA','-append')
                AnalysisData.TsneX=TsneX;
                AnalysisData.PCA=PCA;
            else
                AnalysisData.TsneX=ClustInfo.TsneX;
                AnalysisData.PCA=ClustInfo.PCA;
            end
            %             AnalysisData.TsneX=tsne(AnalysisData.NormW);
            %             [~,AnalysisData.PCA]=pca(AnalysisData.NormW);
            
            % plot Hs of clusters
            for ch=1:NCh
                Hs{ch}=cell2mat(arrayfun(@(x) MotifData{ch}.AnalysisData.MotifSpecs_Disc{x}.H,...
                    1:NTimChunks,'UniformOutput',false));
                Clust_Hs{ch}=cell2mat(arrayfun(@(x) AnalysisData.IdentifierH{ch}(:,x).*ones(NH,NtimH),...
                    1:NTimChunks,'UniformOutput',0));
                ThisClustH{ch}=zeros(Nclust,NtimH*NTimChunks);
                for k=1:NH
                    for i=1:Nclust
                        ind=Clust_Hs{ch}(k,:)==i;
                        temp=zeros(1,NtimH*NTimChunks);
                        temp(ind)=Hs{ch}(k,ind);
                        ThisClustH{ch}(i,:)=ThisClustH{ch}(i,:)+temp;
                    end
                end
                varargout{ch}=figure;
                %                      for x=1:Nclust
                %                          subplot(Nclust,1,x);
                %                          plot(ThisClustH{ch}(x,:),'color',col(x,:));
                %                          title(['H for cluster]' num2str(x) ' Ch' num2str(ch)])
                %                          axis tight
                %                      end
                % plot(find(ind),Hs{ch}(k,ind),'color',col(i,:))
            end
            varargout{ch+1}=clustFig{1};
            AnalysisData.ClustH=ThisClustH;
        end
        function  [DistMat,ShiftMat]=TemporalXCorrTensor(obj,Motifs,varargin) % calcultes X by X atrix of temporal xcorr of motifs
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %%
            Nmotifs=size(Motifs,1);
            MotifPairs=nchoosek(1:Nmotifs,2);  % take pairs of motifs
            NMotifPairs=size(MotifPairs,1);
            DistMat=nan(1,NMotifPairs);
            ShiftMat=nan(1,NMotifPairs);
            i=MotifPairs(:,1);j=MotifPairs(:,2);
            SizeW=AnalysisData.SizeW;
            if ~AnalysisOpts.RunonCluster
                parfor x=1:NMotifPairs
                    [DistMat(x),ShiftMat(x)]= TemporalCrossCorrDistance1DParLoop(i(x),j(x),Motifs,SizeW);
                    if mod(x,floor(NMotifPairs/100))==0
                        fprintf('\n Pair #%i %i Percent ELT',x,x*100/NMotifPairs)
                    end
                end
            else
                for x=1:NMotifPairs
                    
                    [DistMat(x),ShiftMat(x)]= TemporalCrossCorrDistance1DParLoop(i(x),j(x),Motifs,SizeW);
                    if mod(x,floor(NMotifPairs/100))==0
                        fprintf('\n Pair #%i %i Percent ELT',x,x*100/NMotifPairs)
                    end
                end
            end
            %%% put them in similarity matrix now
            DistMat=obj.ManData.ReshapeSquareMatrix(MotifPairs,DistMat);
            ShiftMat=obj.ManData.ReshapeSquareMatrix(MotifPairs,ShiftMat);
            
            if AnalysisOpts.SaveData
                [~,~,XCorrFileName]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,AnalysisOpts.AnalysisPathName,AnalysisOpts.Animal,AnalysisOpts.RecDate,obj.ChNum,['_MotifDist'  AnalysisOpts.ExtraParamTxtWrite],'ext','.mat');
                save(XCorrFileName,'DistMat','ShiftMat','-v7.3')
            end
            
        end
        
        function  [CoreMotifs,TemplateMotif,NormW_Aligned,clust_ind_Aligned,clust_ind_notAligned,TestCoreMotifs]=FindTemplateMotif(obj,Clusters,Motifs,SimMat,varargin)  % find the template cluster in a motif
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %%
            fprintf('\n Finding Template Motif...')
            Clusts=[unique(Clusters)]';
            AnalysisData.NormW_Aligned=[];
            AnalysisData.clust_ind_Aligned=[];
            AnalysisData.clust_ind_notAligned=[]; % how do we match this to orginal inds
              for thiscl=Clusts
                ClstInd=find(Clusters==thiscl);
                ThsClstMtfs=Motifs(ClstInd,:);
                NclustElem=length(ClstInd);      % number of elements of cluster
                
                % find similarity of each motif of this clust to the rest of members and find the remplate motif
                MeanMembDist=arrayfun(@(x) mean(SimMat(x,setdiff(ClstInd',x))),ClstInd)'; % mean distance of each member to the rest
                [sortMotifsSim,sortMotifsSimInd]=sort(MeanMembDist,'descend'); % sort the mean similarities
                [tempMotifInd]=sortMotifsSimInd(1);
                % are we puting the best motif as template or choosing
                % random
                if strcmpi(AnalysisOpts.MotifAnalysis.TemplateMotif,'best')
                    TempMotif=tempMotifInd;   %randsample([1:NclustElem],1);
                elseif strcmpi(AnalysisOpts.MotifAnalysis.TemplateMotif,'random')
                    TempMotif=randsample([1:NclustElem],1);
                end
                
                if strcmpi(AnalysisOpts.MotifAnalysis.AlignMotifs,'xcorr')
                    % align motifs to template motif using max temp corr
                    NewMotif_corctRS=obj.AlignMotifs(ThsClstMtfs(TempMotif,:),ThsClstMtfs);
                    NewMotif_corct=arrayfun(@(x) reshape(NewMotif_corctRS(x,:),AnalysisData.SizeW(1),AnalysisData.SizeW(2)*3),...
                        1:size(NewMotif_corctRS,1),'UniformOutput' ,0);
                    NewMotif_corct=obj.ManData.ReshapeCell2Mat(NewMotif_corct,3);
                    CoreCommunity{thiscl}=NewMotif_corct(:,:,sortMotifsSimInd(1:floor(NclustElem*obj.CoreMotifsPerc)));
                    CoreMotifs{thiscl}= mean(CoreCommunity{thiscl},3);
                    TemplateMotif{thiscl}=NewMotif_corct(:,:,TempMotif);
                elseif strcmpi(AnalysisOpts.MotifAnalysis.AlignMotifs,'dtw')
                    %   Align motifs using dtw
                    [~,NewMotif_corct,CorrAfterDTW]=obj.AlignMotifsDTW(ThsClstMtfs(TempMotif,:),ThsClstMtfs);
                    CoreCommunity{thiscl}=NewMotif_corct(:,:,sortMotifsSimInd(1:floor(NclustElem*obj.CoreMotifsPerc)));
                    CoreMotifs{thiscl}=nanmean(CoreCommunity{thiscl},3);
                    TemplateMotif{thiscl}=NewMotif_corct(:,:,TempMotif);
                end
                %                 figure
                %                 helperCWTTimeFreqPlot(CoreMotifs{thiscl},1:size(CoreMotifs{thiscl},2),1:size(CoreMotifs{thiscl},1),'justplot1',[''],'Time','Frequency(Hz)',0)
                %
                %                 [sortMotifsDTW,sortMotifsDTWSimInd]=sort(CorrAfterDTW,'descend'); % sort the mean similarities
                %                 for i=0.1:0.05:1
                %                     figure;obj.CoreMotifsPerc=i;
                %                     MeanClusterDTW{thiscl}=nanmean(NewMotifDTW_padd(:,:,sortMotifsDTWSimInd(1:floor(NclustElem*obj.CoreMotifsPerc))),3);
                %                     %  MeanClusterDTW{thiscl}=nanmean(NewMotifDTW_padd(:,:,sortMotifsDTWSimInd(end-floor(NclustElem*obj.CoreMotifsPerc:end))),3);
                %                     helperCWTTimeFreqPlot(MeanClusterDTW{thiscl},1:size(MeanClusterDTW{thiscl},2),1:size(MeanClusterDTW{thiscl},1),'justplot1',[num2str(i) ],'Time','Frequency(Hz)',0);
                %                 end
                
                % generate core motifs with different percent of data so we
                % can seeit
                for PercCM=[0.1 0.5 1]
                    temp=NewMotif_corct(:,:,sortMotifsSimInd(1:floor(NclustElem*PercCM)));
                    TestCoreMotifs.(['P' strrep(num2str(PercCM),'.','')] ){thiscl}=mean(temp,3);
                end
                %
                AnalysisData.NormW_Aligned=cat(1,AnalysisData.NormW_Aligned, NewMotif_corctRS);
                AnalysisData.clust_ind_Aligned=cat(1,AnalysisData.clust_ind_Aligned,thiscl*ones(NclustElem,1));
                AnalysisData.clust_ind_notAligned=cat(1,AnalysisData.clust_ind_notAligned,ClstInd);
            end
            AnalysisData.CoreMotifs=CoreMotifs;
            NormW_Aligned=AnalysisData.NormW_Aligned;
            clust_ind_Aligned=AnalysisData.clust_ind_Aligned;
            clust_ind_notAligned=AnalysisData.clust_ind_notAligned;
        end
        function  NewMotif_corct=AlignMotifs(obj,TempMotif,AllMotifs,varargin)  % Aligns motifs to template Motif
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            % align motifs to template with max xcorr
            
            [~,~,NewMotif_corct] = TemporalCrossCorr(TempMotif,AllMotifs);
            
            %             PlotDetails=0;
            %             SizeW=AnalysisData.SizeW;
            %             LMotf=SizeW(2); % number of samples for motif
            %             TempMotf=reshape(TempMotif,SizeW);
            %             TempMotf_pad=padarray(TempMotf,[0,LMotf]);
            %             COr=@(x,y) corr(x(:),y(:));
            %             ShiftRng=-LMotf:LMotf;
            %
            %             for j=1:size(AllMotifs)
            %
            %                 ThisMotf=reshape(AllMotifs(j,:),SizeW);
            %                 ThisMotf_pad=padarray(ThisMotf,[0,LMotf]);
            %                 crr=arrayfun(@(x) COr(TempMotf_pad,circshift(ThisMotf_pad,x,2)),ShiftRng,'UniformOutput',1);
            %                 [ssr,snd] = max(crr(:));
            %                 Shift=ShiftRng(snd); % how many to shift
            %
            %                 % shift this motif to match template motif
            %                 NewMotif=circshift(ThisMotf_pad,Shift,2);
            %                 NewMotif_corct(j,:)=reshape(NewMotif,[1 size(NewMotif,1)*size(NewMotif,2)]);
            %
            %                 if PlotDetails
            %                     subplot(141)
            %                     helperCWTTimeFreqPlot(TempMotf_pad,1:size(TempMotf_pad,2),AnalysisData.cwt_f,'justplot1',['PSD'],'Time','Frequency(Hz)',0)
            %                     subplot(142)
            %                     helperCWTTimeFreqPlot(ThisMotf_pad,1:size(ThisMotf_pad,2),AnalysisData.cwt_f,'justplot1',['PSD'],'Time','Frequency(Hz)',0)
            %                     subplot(143)
            %                     plot(crr(:))
            %                     title('Cross-Correlation')
            %                     hold on
            %                     plot(snd,ssr,'or')
            %                     hold off
            %                     text(snd*1.05,ssr,'Maximum')
            %                     subplot(144)
            %                     helperCWTTimeFreqPlot(NewMotif,1:size(NewMotif,2),AnalysisData.cwt_f,'justplot1',['PSD'],'Time','Frequency(Hz)',0)
            %                 end
            %             end
            %
        end
        function  [NewMotif_corctDTW,NewMotifDTW_padd,CorrAfter]=AlignMotifsDTW(obj,TempMotif,AllMotifs,varargin)  % Aligns motifs to template Motif using dynamic time warping
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %%  global AnalysisData
            PlotDetails=0;
            SizeW=AnalysisData.SizeW;
            LMotf=SizeW(2); % number of samples for motif
            TempMotf=reshape(TempMotif,SizeW);
            TempMotf_pad=padarray(TempMotf,[0,LMotf]);
            
            for j=1:size(AllMotifs,1)
                
                ThisMotf=reshape(AllMotifs(j,:),SizeW);
                ThisMotf_pad=padarray(ThisMotf,[0,LMotf]);
                
                [dist,ix,iy]=dtw(ThisMotf,TempMotf);
                ThisMotf_warp=ThisMotf(:,ix);
                TempMotf_warp=TempMotf(:,iy);
                NewMotif_corctDTW{j}= ThisMotf_warp;
                CorrBefore=corr2(TempMotf,ThisMotf);
                CorrAfter(j)=corr2(TempMotf_warp,ThisMotf_warp); % correlation after dtw
                if PlotDetails
                    
                    subplot(221)
                    helperCWTTimeFreqPlot(TempMotf,1:size(TempMotf,2),1:size(TempMotf,1),'justplot1',['TemoMotif ' num2str(CorrBefore,2) ],'Time','Frequency(Hz)',0)
                    subplot(222)
                    helperCWTTimeFreqPlot(TempMotf_warp,1:size(TempMotf_warp,2),1:size(TempMotf_warp,1),'justplot1',['TempMotif warped ' num2str(CorrAfter,2)],'Time','Frequency(Hz)',0)
                    subplot(223)
                    helperCWTTimeFreqPlot(ThisMotf,1:size(ThisMotf,2),1:size(ThisMotf,1),'justplot1',['Motif 2'],'Time','Frequency(Hz)',0)
                    subplot(224)
                    helperCWTTimeFreqPlot(ThisMotf_warp,1:size(ThisMotf_warp,2),1:size(ThisMotf_warp,1),'justplot1',['Motif 2 warped'],'Time','Frequency(Hz)',0)
                    pause
                end
            end
            SizeMtfs=cellfun(@(x) size(x,2),NewMotif_corctDTW);
            MaxSize=max(SizeMtfs);
            NewMotif_padd=arrayfun(@(x) padarray(NewMotif_corctDTW{x},[0 MaxSize-SizeMtfs(x)],NaN,'post'),1:size(AllMotifs,1),'UniformOutput',0);
            NewMotifDTW_padd=obj.ManData.ReshapeCell2Mat(NewMotif_padd,3);
            Mean_NewMotif_padd=nanmean(NewMotifDTW_padd,3);
            %  helperCWTTimeFreqPlot(Mean_NewMotif_padd,1:size(Mean_NewMotif_padd,2),1:size(Mean_NewMotif_padd,1),'justplot1',[''],'Time','Frequency(Hz)',0)
            
        end
        
        function  Motifs=RemoveNoloadingMotifs(~,Motifs)% removes motifs that have no loadings
            % loadings=Motifs.loadings;
            SumMotifs=sum(sum(Motifs.W,3),1);
            
            OkMotifs=SumMotifs>0;
            Motifs.W=Motifs.W(:,OkMotifs,:);
            Motifs.H=Motifs.H(OkMotifs,:);
            % Motifs.loadings=Motifs.loadings(OkMotifs);
            if isfield(Motifs,'PEVind')
                Motifs.PEVind= Motifs.PEVind(OkMotifs);
            end
        end
        function  varargout=PlotMotifClustsH(obj,TrialTimes,BlockSpec,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            varargout=[];
            
            NBlks=length(BlockSpec.Rule);
            
            for ch=1:length(AnalysisData.Ch)
                for b=1:NBlks
                    ThisTrialTimes=BlockSpec.ThisBlkTrialTimes{b};
                    H_blk{b}=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(...
                        AnalysisData.ClustH{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,...
                        'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                    H_blk_ro(:,:,b)=obj.ManData.Cell2Mat(H_blk{b}); % reorganized H blck
                end
                Time=-1:1/obj.Fs:(0.5-0.01); %% what was atually in the experiment
                for i=1:length(AnalysisData.MeanCluster);W(:,i,:)=AnalysisData.MeanCluster{i};end
                % take sum example blocks and mean block
                ExpBlkNum=1;
                H_exmp=arrayfun(@(x) H_blk_ro(:,:,x),ExpBlkNum,'UniformOutput',0);
                BlkTrialNumber=arrayfun(@(x) BlockSpec.ThisBlkTrialNum(x),ExpBlkNum,'UniformOutput',1);
                H_exmp=[H_exmp {mean(H_blk_ro,3)}]; % add mean block
                BlkTrialNumber=[BlkTrialNumber BlkTrialNumber(1)];
                ExpBlkNum=[ExpBlkNum NaN];
                
                Wfigs=cell(1,length(AnalysisData.MeanCluster));
                [Wfigs{:}]=obj.PlotW_H_Trial(W,H_exmp,AnalysisData.f,Time,BlkTrialNumber,BlockSpec.TrialOrder,ExpBlkNum,'ChNum',ch) ; % plot Ws
                
                varargout=[varargout Wfigs];
            end
            
        end
        function  [MotifSpecs_Refit]=RefitMotifsGeneral(obj,ChNum,ThisDate,CoreMotifs,MaxTime,SizeW,varargin)
            % general function refits motifs to each channels's wavelet to get the H
            % weightings
            % Maxtime is duration of random chuck taken from data to refit.
            % Maxtime=0 if all of the signal needs to be taken
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % if we already have it don't recalculate it
            if isempty(CoreMotifs);CoreMotifs=AnalysisData.CoreMotifs;end
            CoreMotifs=obj.ManData.ReshapeCell2Mat(CoreMotifs,32);
            % trim motifs so it matches the size we want
            CoreMotifs=obj.ManData.TrimMatrix(CoreMotifs,SizeW(2));
            NMotifs=size(CoreMotifs,2);LMotifs=size(CoreMotifs,3);
            ChCount=1;
            for Ch=ChNum % load each channels's Motif fit data and Original
                % start processing
                tic
                fprintf('\nRefitting to held-out data CH %i',Ch)
                % load motifs data of this channel
                AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;% save current rec date
                AnalysisOpts.RecDate=obj.ManData.DetermineDate(ThisDate);
                AnalysisOpts.CurrentCh=Ch;
                [~,~,MotifsFullPath]=obj.ManData.GetFileName('Motif_Discover',AnalysisOpts.ExtraParamTxtRead);
                MotifData=load(MotifsFullPath);
                [AnalysisData.cwt,AnalysisData.cwt_f]=obj.TlBxFuncs.getWaveletData(Ch,0);
                AnalysisData.cwt=AnalysisData.cwt{1};
                
                % now adjust the TF data to what had when we were fitting
                % the data
                FsWaveTarg=MotifData.AnalysisOpts.MotifAnalysis.FsWaveTarg;
                FsWave=MotifData.AnalysisOpts.MotifAnalysis.FsWave;
                AnalysisData.cwt=transpose(downsample(AnalysisData.cwt',FsWave/FsWaveTarg));
                AnalysisData.Time_cwt=0:1/FsWaveTarg:(size(AnalysisData.cwt,2)-1)/FsWaveTarg;
                % get a random chunk of data
                if MaxTime~=0 % put max time to zero if you want the whole signal
                    TimeStr=randsample(AnalysisData.Time_cwt(AnalysisData.Time_cwt>=0 & AnalysisData.Time_cwt<(AnalysisData.Time_cwt(end)-2*MaxTime-1/FsWaveTarg)),1);
                    % timing for discovering the motifs
                    TimeInd=AnalysisData.Time_cwt>=TimeStr & AnalysisData.Time_cwt<(MaxTime+TimeStr);
                else
                    TimeInd=1:length(AnalysisData.Time_cwt);
                end
                % now match the frequencies
                flimitind=AnalysisData.cwt_f>=MotifData.AnalysisData.cwt_f(1) & AnalysisData.cwt_f<=MotifData.AnalysisData.cwt_f(end);
                AnalysisData.cwt_f=AnalysisData.cwt_f(flimitind);
                AnalysisData.cwt=AnalysisData.cwt(flimitind,TimeInd);
                AnalysisData.Time_cwt=AnalysisData.Time_cwt(TimeInd);
                % now refit the core motifs(or any other set of motifs)
                [MotifSpecs_Refit{ChCount}]=obj.PlotMotifs(AnalysisData.cwt,AnalysisData.cwt_f,[],...
                    'SavePlot',0,'SaveData',0,'ShowPlot',0,...
                    'DownSampleFactor',1,'ExtraString','','LogScale',0,'PowerMethod' ,MotifData.AnalysisOpts.MotifAnalysis.PowerMethod,...
                    'maxiter',AnalysisOpts.MotifAnalysis.maxiterRefit,'K',NMotifs,...
                    'L',LMotifs,'lambda',0,'lambdaOrthoH',MotifData.AnalysisOpts.MotifAnalysis.lambdaOrthoH,...
                    'lambdaOrthoW',MotifData.AnalysisOpts.MotifAnalysis.lambdaOrthoW,'lambdaL1W',MotifData.AnalysisOpts.MotifAnalysis.lambdaL1W,...
                    'lambdaL1H',MotifData.AnalysisOpts.MotifAnalysis.lambdaL1H,...
                    'W_init',CoreMotifs,'W_fixed',1,'useWupdate',0);
                % if we need this chunk's raw data as well
                % get raw data
                [AnalysisData.RawData,FsRaw]=obj.TlBxFuncs.getRawEphysData(Ch,0);AnalysisData.RawData=AnalysisData.RawData{1};
                AnalysisData.Time_raw=0:1/FsRaw:(size(AnalysisData.RawData,2)-1)/FsRaw;
                if MaxTime~=0 % put max time to zero if you want the whole signal
                    % timing for discovering the motifs
                    RawTimeInd=AnalysisData.Time_raw>=TimeStr & AnalysisData.Time_raw<(MaxTime+TimeStr);
                else
                    RawTimeInd=1:length(AnalysisData.Time_raw);
                end
                AnalysisData.RawData=AnalysisData.RawData(RawTimeInd);
                % save this Chs data
                MotifSpecs_Refit=MotifSpecs_Refit{ChCount};
                ChCount=ChCount+1;
                toc
            end
            AnalysisOpts.RecDate=AnalysisOpts.CurrentRecDate; % revert back to previous rec date
            
        end
        function  [MotifSpecs_Refit]=RefitMotifs(obj,ChNum,CoreMotifs,SizeW,XCorrFileName,varargin)
            % refits motifs to each channels's wavelet to get the H
            % weightings
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % if we already have it don't recalculate it
            if obj.ManData.IsVarExistinFile(XCorrFileName,'MotifSpecs_Refit') & ~obj.OverWrite
                load(XCorrFileName,'MotifSpecs_Refit');
                return
            end
            if isempty(ChNum);ChNum=AnalysisData.Ch;end
            %% prepare Core motifs
            %    MatFile=matfile(XCorrFileName,'Writable',true);
            [PATH,TargFile]=fileparts(XCorrFileName);
            
            if isempty(CoreMotifs);CoreMotifs=AnalysisData.CoreMotifs;end
            CoreMotifs=obj.ManData.ReshapeCell2Mat(CoreMotifs,32);
            % trim motifs so it matches the size we want
            CoreMotifs=obj.ManData.TrimMatrix(CoreMotifs,SizeW(2));
            NMotifs=size(CoreMotifs,2);LMotifs=size(CoreMotifs,3);
            ChCount=1;
            for Ch=ChNum % load each channels's Motif fit data and Original
                % if we have this Chs file don't recalculate it
                ChMotifFile=[PATH AnalysisOpts.FS TargFile '_rf' num2str(Ch) '.mat'];
                if exist(ChMotifFile,'file') & ~obj.OverWrite % if there is then concatinate them all
                    for ThisCh= AnalysisData.Ch
                        fprintf('\nLoading refit motifs from Channel %i',ThisCh)
                        ThisChMotifFile=[PATH AnalysisOpts.FS TargFile '_rf' num2str(ThisCh) '.mat'];
                        ThisMotifSpecs_Refit=load(ThisChMotifFile);
                        IndCh=find(ThisCh==AnalysisData.Ch);
                        MotifSpecs_Refit{IndCh}=ThisMotifSpecs_Refit.MotifSpecs_Refit;
                        % now get the cwt from the first channel
                        if IndCh==1
                            AnalysisOpts.CurrentCh=ThisCh;
                            [~,~,MotifsFullPath]=obj.ManData.GetFileName('Motif_Discover',AnalysisOpts.ExtraParamTxtRead);
                            % [~,~,MotifsFullPath]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.DataSavePath,'Motif_Discover',AnalysisOpts.Animal,AnalysisOpts.RecDate,Ch,['_L' num2str(AnalysisOpts.MotifAnalysis.L_ms) 'T' num2str(AnalysisOpts.MotifAnalysis.MaxTime)]);
                            MotifData=load(MotifsFullPath);
                            AnalysisData.cwt_f=MotifData.AnalysisData.cwt_f; % load frequency
                        end
                    end
                    return
                end
                % start processing
                tic
                fprintf('\nRefitting to held-out data CH %i',Ch)
                % load motifs data of this channel
                [~,ThisDate]=fileparts(PATH);
                AnalysisOpts.CurrentRecDate=AnalysisOpts.RecDate;% save current rec date
                AnalysisOpts.RecDate=ThisDate; % revert to this date
                AnalysisOpts.CurrentCh=Ch;
                [~,~,MotifsFullPath]=obj.ManData.GetFileName('Motif_Discover',AnalysisOpts.ExtraParamTxtRead);
                MotifData=load(MotifsFullPath);
                %obj.TlBxFuncs.getEphysData(Ch); % load wavelet data
                [AnalysisData.cwt,AnalysisData.cwt_f]=obj.TlBxFuncs.getWaveletData(Ch,0);
                AnalysisData.cwt=AnalysisData.cwt{1};
                % now adjust the TF data to what had when we were fitting
                % the data
                FsWaveTarg=MotifData.AnalysisOpts.MotifAnalysis.FsWaveTarg;
                FsWave=MotifData.AnalysisOpts.MotifAnalysis.FsWave;
                AnalysisData.cwt=transpose(downsample(AnalysisData.cwt',FsWave/FsWaveTarg));
                AnalysisData.Time_cwt=0:1/FsWaveTarg:(size(AnalysisData.cwt,2)-1)/FsWaveTarg;
                % now match the frequencies
                flimitind=AnalysisData.cwt_f>=MotifData.AnalysisData.cwt_f(1) & AnalysisData.cwt_f<=MotifData.AnalysisData.cwt_f(end);
                AnalysisData.cwt_f=AnalysisData.cwt_f(flimitind);
                AnalysisData.cwt=AnalysisData.cwt(flimitind,:);
                
                % now refit the core motifs(or any other set of motifs)
                [MotifSpecs_Refit{ChCount}]=obj.PlotMotifs(AnalysisData.cwt,AnalysisData.cwt_f,[],...
                    'SavePlot',0,'SaveData',0,'ShowPlot',0,...
                    'DownSampleFactor',1,'ExtraString','','LogScale',0,'PowerMethod' ,MotifData.AnalysisOpts.MotifAnalysis.PowerMethod,...
                    'maxiter',AnalysisOpts.MotifAnalysis.maxiterRefit,'K',NMotifs,...
                    'L',LMotifs,'lambda',0,'lambdaOrthoH',MotifData.AnalysisOpts.MotifAnalysis.lambdaOrthoH,...
                    'lambdaOrthoW',MotifData.AnalysisOpts.MotifAnalysis.lambdaOrthoW,'lambdaL1W',MotifData.AnalysisOpts.MotifAnalysis.lambdaL1W,'lambdaL1H',1,...
                    'W_init',CoreMotifs,'W_fixed',1,'useWupdate',0);
                
                % save this Chs data
                MotifSpecs_Refit=MotifSpecs_Refit{ChCount};
                save(ChMotifFile,'MotifSpecs_Refit')
                %  eval(['MatFile.MotifSpecs_Refit_' num2str(Ch) '=MotifSpecs_Refit{ChCount};']);
                
                ChCount=ChCount+1;
                toc
            end
            %   save(XCorrFileName,'MotifSpecs_Refit','-append')
            AnalysisOpts.RecDate=AnalysisOpts.CurrentRecDate; % revert back to previous rec date
            
        end
        function  [W,H,Lags,varargout]=ShiftWH2Zero(~,W,H,TH)
            
            NWs=size(W,2);
            % TH=0.25; % percent of max value that we consider as start time
            Ncol=ceil(NWs/5);
            Nraw=5;
            Nch=length(H);
            varargout{1}=figure;
            for i=1:NWs
                MeanThisMotif=mean(squeeze(W(:,i,:)),1);
                Maxval=max(MeanThisMotif);
                Lags(i)=find(MeanThisMotif>=Maxval*TH,1,'first');
                %
                subplot(Nraw,Ncol,i)
                imagesc(squeeze(W(:,i,:)));
                set(gca,'YDir','normal');
                v=axis;
                hold on
                plot([Lags(i) Lags(i)],[v(3) v(4)],'r')
                % plot(MeanThisMotif/Maxval,'-r')
                title(num2str(i))
            end
            % now shift H and W values
            for ch=1:Nch
                TempH{ch}=transpose(cell2mat(arrayfun(@(x) transpose(circshift(H{ch}(x,:),Lags(x))),1:NWs,'Uniformoutput',0)));
            end
            H=TempH;
            
        end
        function  [varargout]=PlotCoreMotifs(obj,CoreMotifs,Nmotif2plot,TrimMotifs,varargin) % plots the core motifs
            fprintf('\nPlotCoreMotifs: Plotting Core Motifs')
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
                       
            
            if nargin<3;TrimMotifs=1;end
            opts.TrimMotifs=TrimMotifs;
            opts.NormalizeMotif=1;
            
            if isempty(CoreMotifs)
                CoreMotifs=AnalysisData.CoreMotifs;
            end
            if ~iscell(CoreMotifs) % if it is a 3 dim matrix then turn it into a cell
                CoreMotifs=arrayfun(@(x) squeeze(CoreMotifs(:,x,:)),1:size(CoreMotifs,2),'UniformOutput',0);
            end
            % if we need to trim the motifs
            if opts.TrimMotifs
                CoreMotifs=obj.ManData.TrimMatrix(CoreMotifs,AnalysisData.SizeW(2));
            end
            
            if ~isempty(Nmotif2plot)
                CoreMotifs=CoreMotifs(Nmotif2plot);
            end
            if opts.NormalizeMotif
                CoreMotifs=obj.ManData.NormalizeMotifs(CoreMotifs);
            end
            NWs=length(CoreMotifs);
            NTim=size(CoreMotifs{1},2);
            
            FsWaveTarg=AnalysisOpts.MotifAnalysis.FsWaveTarg;
            if NWs>9;Ncol=ceil(NWs/3);Nraw=3;else;Ncol=3;Nraw=3;end
            % if NWs>4;Nraw=4;else;Nraw=NWs;end
            if   obj.NewFig
                obj.FigParams.RenderFigure(1,[]);
            end

            for i=1:NWs
                if length(Nmotif2plot)~=1
                    subplot(Nraw,Ncol,i)
                end
                hold on
                if strcmp(obj.FigTitle,'#')
                    Title=[obj.FigTitle  num2str(i)];
                else
                    Title=obj.FigTitle;
                end
                helperCWTTimeFreqPlot(CoreMotifs{i},(1:NTim)/FsWaveTarg,AnalysisData.cwt_f,'Justplot1',[Title],'Time(s)','Freq(Hz)',obj.LogScale);
                % plot lag if we need it
                if ~isempty(obj.FigTitle)
                    MeanThisMotif=mean(CoreMotifs{i},1);
                    Maxval=max(MeanThisMotif);
                    Lags(i)=find(MeanThisMotif>=Maxval*AnalysisOpts.MotifAnalysis.ShiftWH2ZeroTH,1,'first');
                    v=axis;
                    plot(Lags(i)/FsWaveTarg,v(3),'*r')
                end
                axis square
            end
            varargout{1}=gcf;
        end
        function  [varargout]=PlotMotifStas(obj,Motifs)
        end
        function Perf=CalBhvPerf(~,TrialTimes,navg)
            global AnalysisOpts
            CorrectInd=strcmp(AnalysisOpts.TrialTimesFields,'ResponseError');%% we don't use this here but still good to have it
            Bhv=[TrialTimes(:,CorrectInd)]';
            sm_kern = ones(1, navg);
            sm_kern = sm_kern./sum(sm_kern);
            Perf=convn(Bhv, sm_kern, 'same');
        end
        
        function  varargout=PlotMotifClustsH_Refit_Area(obj,W,H,TrialTimes,BlockSpec,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            varargout=[];
            
            AnalysisData.Time=-1:1/obj.Fs:(0.5-0.01); %% what was atually in the experiment
            %% grab and plot Reward data for all of trails and channels
            % [H_tot_Reward]=obj.GrabHTimeTrialEpisode(H,TrialTimes,'GIVE_REWARD');
            %   [H_tot_FixOn] =obj.GrabHTimeTrialEpisode(H,TrialTimes,'FIX_ON');
            %             load H_data H_tot_Reward
            %             % plot average H for rewards for each channel and each area
            %             obj.PlotW_H_EachW_Condition(W,H_tot_Reward,AnalysisData.cwt_f,AnalysisData.Time,TrialTimes,AnalysisData.ChArea,'NumRewards',0:6);
            
            %% grab Color data for all of trails and channels
            [H_tot_TargOn]=obj.GrabHTimeTrialEpisode(H,TrialTimes,'TARGET_ON');ColList=unique(TrialTimes(:,7));
            %             % plot color data now
            ColorFigs=obj.PlotW_H_EachW_Condition(W,H_tot_TargOn,AnalysisData.cwt_f,AnalysisData.Time,TrialTimes,AnalysisData.ChArea,'BestColor',ColList);
            
            
            % grab H values based on 35 before and after switch
            [H_exmp_tot,BlkTrialNumber,ExpBlkNum,Bhv_blk_Perf]=obj.GrabHTimeTrial(H,BlockSpec);
            
            obj.PlotW_H_Trial_Condition(W,H_exmp_tot.Reward,AnalysisData.cwt_f,AnalysisData.Time,BlockSpec,AnalysisData.ChArea,[],[]);
            
            
            
            
            Wfigs=cell(1,length(AnalysisData.MeanCluster));
            %% area based Analysis
            % plot W-H Area
            [Wfigs{:}]=obj.PlotW_H_Trial_Area(W,H_exmp_tot.FixOn(:,WBlks),AnalysisData.cwt_f,Time,BlkTrialNumber(WBlks),BlockSpec.TrialOrder,ExpBlkNum(WBlks),AnalysisData.ChArea,'ChNum',NaN) ; % plot Ws
            
            % get the H data for all of the trials
            
            
            %           HRaw=obj.ManData.ReshapeCell2Mat(H_exmp_tot,42); % reshape all channels into a big matrix
            %                     Nblks=length(ExpBlkNum);
            %                     AreaInds=[1 4 5];
            %                     for i=AreaInds % go over areas and take average for each block
            %                         for b=1:Nblks
            %                             Ind=find(i==AreaInds);
            %                             ThisAreaNeu=AnalysisData.ChArea==i;
            %                             H_Area{Ind,b}=mean(HRaw{b}(:,:,:,ThisAreaNeu),4);
            %                         end
            %                     end
            %                    Fig_Xcorr=cell(1,length(AnalysisData.MeanCluster)); % look at PFC and LIP for now
            %                    [XcorrPairBlck,T_XcorrPairBlck]=arrayfun(@(x) obj.CalXcorrHChPair(H_Area{1,x},H_Area{2,x}),1:size(H_Area,2),'UniformOutput',0);
            %                    T_XcorrPairBlck=T_XcorrPairBlck{1}{1};
            %                    [Fig_Xcorr{:}]=obj.PlotXcorrHChPair(XcorrPairBlck,T_XcorrPairBlck,BlockSpec.TrialOrder,ExpBlkNum,Bhv_blk_Perf) ;
            %
            varargout=[varargout Wfigs ];
            
        end
        function  varargout=PlotMotifClustsH_Refit_Channel(obj,W,H,TrialTimes,BlockSpec,varargin)
            % plots refit of Cluster H for each channel
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            varargout=[];
            
            NBlks=length(BlockSpec.Rule);
            H_exmp_tot=[];
            for ch=1:length(AnalysisData.Ch)
                H_blk=[];H_blk_ro=[];
                for b=1:NBlks
                    ThisTrialTimes=BlockSpec.ThisBlkTrialTimes{b};
                    H_blk{b}=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(...
                        H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,...
                        'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                    H_blk_ro(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_blk{b},3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    % calculate performance
                    if ch==1
                        Bhv_blk_Perf(b,:)=obj.CalBhvPerf(ThisTrialTimes,10);
                    end
                end
                Time=-1:1/obj.Fs:(0.5-0.01); %% what was atually in the experiment
                % take sum example blocks and mean block
                ExpBlkNum=1:3;
                H_exmp=arrayfun(@(x) H_blk_ro(:,:,:,x),ExpBlkNum,'UniformOutput',0);
                BlkTrialNumber=arrayfun(@(x) BlockSpec.ThisBlkTrialNum(x),ExpBlkNum,'UniformOutput',1);
                H_exmp=[H_exmp {mean(H_blk_ro,4)}]; % add mean block
                BlkTrialNumber=[BlkTrialNumber BlkTrialNumber(1)];
                ExpBlkNum=[ExpBlkNum NaN];
                
                H_exmp_tot=[H_exmp_tot;H_exmp];
            end
            Wfigs=cell(1,length(AnalysisData.MeanCluster));
            
            %% Neuron based analysis
            %% do cross correlation analysis on this pair of neurons for each block
            [Wfigs{:}]=obj.PlotW_H_Trial(W,H_exmp_tot,AnalysisData.cwt_f,Time,BlkTrialNumber,BlockSpec.TrialOrder,ExpBlkNum,'ChNum',ch) ; % plot Ws
            Fig_Xcorr=cell(1,length(AnalysisData.MeanCluster));
            [XcorrPairBlck,T_XcorrPairBlck]=arrayfun(@(x) obj.CalXcorrHChPair(H_exmp_tot{1,x},H_exmp_tot{2,x}),1:size(H_exmp_tot,2),'UniformOutput',0);
            T_XcorrPairBlck=T_XcorrPairBlck{1}{1};
            [Fig_Xcorr{:}]=obj.PlotXcorrHChPair(XcorrPairBlck,T_XcorrPairBlck,BlockSpec.TrialOrder,ExpBlkNum,Bhv_blk_Perf) ;
            
            varargout=[varargout Wfigs Fig_Xcorr];
            
        end
        function  varargout=PlotMotifClustsH_Refit_Flow(obj,W,H,TrialTimes,BlockSpec,varargin)
            % plots refit of Cluster H for all of the channels looking at
            % flow of information
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            obj.UseDataPointer=0; % use data pointer to load data
            varargout=[];
            Nch=length(AnalysisData.Ch);
            NBlks=length(BlockSpec.Rule);
            Blks2Look=1:NBlks; % which blocks are we looking at
            ChPairs=nchoosek(1:Nch,2);
            % get H data first
            [H_exmp_tot,BlkTrialNumber,ExpBlkNum,Bhv_blk_Perf]=obj.GrabHTimeTrial(H,BlockSpec,'Fs',AnalysisOpts.NeuralAnalysis.FsWaveTarg);
            obj.ManData.SaveVar([],H_exmp_tot,'H_exmp_tot','TrialData')
            % get raw data now
            [RawData,FsLFP,DataPointer]=obj.TlBxFuncs.getRawEphysData(AnalysisData.Ch,obj.UseDataPointer);
            % grab Raw data
            [RawData_exmp_tot]=obj.GrabHTimeTrial(RawData,BlockSpec,'Fs',FsLFP,'UseDataPointer',obj.UseDataPointer,'DataPointerVar','RawData');
            obj.ManData.SaveVar([],RawData_exmp_tot,'RawData_exmp_tot','TrialData')
            clear RawData
            % grab wavelet data
            %    [WavData_exmp_tot]=obj.GrabHTimeTrial(RawData,BlockSpec,'Fs',FsLFP,'UseDataPointer',obj.UseDataPointer,'DataPointerVar','Wavelet_Linear');
            %     obj.ManData.SaveVar([],WavData_exmp_tot,'WavData_exmp_tot','TrialData')
            
            for Chp=1:size(ChPairs,1)
                Chp1=ChPairs(Chp,1);Chp2=ChPairs(Chp,2);
                % calculte cross correlation for motifs for all of the
                % channels
                for blk=Blks2Look
                    blkInd=find(blk==Blks2Look);
                    [XcorrH{Chp,blkInd}]= obj.CalXcorrHChPair(H_exmp_tot.TargetOn{Chp1,blk},H_exmp_tot.TargetOn{Chp2,blk},'Navg',5,'Nshuffle',10);
                end
                
                % calculate coherence for LFP signal now
                RawData1=RawData_exmp_tot.TargetOn(Chp1,Blks2Look);RawData2=RawData_exmp_tot.TargetOn(Chp2,Blks2Look);
                % [mscohPair,cpsdPair,wavcohPair]=obj.CalCoherenceChPair(RawData1,RawData2,FsLFP,10);
                [mscohPair{Chp}]=obj.CalCoherenceBlockChPair(RawData1,RawData2,FsLFP,10);
            end
            obj.ManData.SaveVar([],mscohPair,'mscohPair','CohData');
            obj.ManData.SaveVar([],XcorrH,'XcorrH','CohData');
            
            %            Wfigs=cell(1,length(AnalysisData.MeanCluster));
            %
            %            %% Neuron based analysis
            %            %% do cross correlation analysis on this pair of neurons for each block
            %            ChPairs=nchoosek(1:Nch,2);
            %            for i=1:NBlks2look
            %                tic
            %               fprintf('\nCalculating Xcorr for block %i ...',i)
            %               [XcorrPairBlck{i},XcorrPairBlckAvg{i},T_XcorrPairBlck]=arrayfun(@(x) obj.CalXcorrHChPair(H_exmp_tot{ChPairs(x,1),i},H_exmp_tot{ChPairs(x,2),i},'Navg',10,'Nshuffle',10),1:size(ChPairs,1),'UniformOutput',0);
            %            %    MaxXcorr{i}=obj.RearrangXcorrHchPair(XcorrPairBlck{i},T_XcorrPairBlck{1}{1},ChPairs);
            %               toc
            %            end
            %            T_XcorrPairBlck=T_XcorrPairBlck{1}{1};
            %            % aggregate different blocks across each % we dont'
            %            % include the last block
            %            XcorrPairBlck_All=obj.ManData.AggregateCellVals(XcorrPairBlck(1:end-1));
            %            XcorrPairBlckAvg_All=obj.ManData.AggregateCellVals(XcorrPairBlckAvg(1:end-1));
            %            %% now rearrange Xcorr
            %            Col=distinguishable_colors(length(ExpBlkNum));
            %            for b=[1:NBlks2look]
            %                MaxXcorr{b}=obj.RearrangXcorrHchPair(XcorrPairBlckAvg{b},T_XcorrPairBlck{1}{1},ChPairs);
            %                for w=1:nW
            %                    for i=1:BlkTrialNumber(1)
            %                        Q(i)=CommunityStructureLouvian(MaxXcorr{b}{w}(:,:,i),0);
            %                       % title(['Trl' num2str(i) ',Q:' num2str(Q(i))]);
            %                      %  drawnow;pause(0);
            %                    end
            %                    subplot(4,4,w);hold on
            %                    plot(BlockSpec.TrialOrder,Q,'color',Col(b,:))
            %                end
            %            end
            %            %%
            %            Fig_Xcorr=cell(1,length(AnalysisData.MeanCluster));
            %           [Fig_Xcorr{:}]=obj.PlotXcorrHChPair(XcorrPairBlck,T_XcorrPairBlck,BlockSpec.TrialOrder,ExpBlkNum,Bhv_blk_Perf) ;
            %            varargout=[varargout Wfigs Fig_Xcorr];
            
        end
        function  varargout=Clust_PlotMotifClustsH_Refit_Flow(obj,H,BlockSpec,ChPairNum,varargin)
            %Runs on cluster
            % plots refit of Cluster H for all of the channels looking at
            % flow of information
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            %             % grab all of the trial data
            %            [H_exmp_tot,RawData_exmp_tot,Wavelet_exmp_tot,f_Linear]=obj.GrabAllTrialData(H,BlockSpec);
            %             obj.CalCoherenceAllChPairs(BlockSpec,H_exmp_tot,RawData_exmp_tot,Wavelet_exmp_tot)
            %
            %            Wfigs=cell(1,length(AnalysisData.MeanCluster));
            %
            %            %% Neuron based analysis
            %            %% do cross correlation analysis on this pair of neurons for each block
            %            ChPairs=nchoosek(1:Nch,2);
            %            for i=1:NBlks2look
            %                tic
            %               fprintf('\nCalculating Xcorr for block %i ...',i)
            %               [XcorrPairBlck{i},XcorrPairBlckAvg{i},T_XcorrPairBlck]=arrayfun(@(x) obj.CalXcorrHChPair(H_exmp_tot{ChPairs(x,1),i},H_exmp_tot{ChPairs(x,2),i},'Navg',10,'Nshuffle',10),1:size(ChPairs,1),'UniformOutput',0);
            %            %    MaxXcorr{i}=obj.RearrangXcorrHchPair(XcorrPairBlck{i},T_XcorrPairBlck{1}{1},ChPairs);
            %               toc
            %            end
            %            T_XcorrPairBlck=T_XcorrPairBlck{1}{1};
            %            % aggregate different blocks across each % we dont'
            %            % include the last block
            %            XcorrPairBlck_All=obj.ManData.AggregateCellVals(XcorrPairBlck(1:end-1));
            %            XcorrPairBlckAvg_All=obj.ManData.AggregateCellVals(XcorrPairBlckAvg(1:end-1));
            %            %% now rearrange Xcorr
            %            Col=distinguishable_colors(length(ExpBlkNum));
            %            for b=[1:NBlks2look]
            %                MaxXcorr{b}=obj.RearrangXcorrHchPair(XcorrPairBlckAvg{b},T_XcorrPairBlck{1}{1},ChPairs);
            %                for w=1:nW
            %                    for i=1:BlkTrialNumber(1)
            %                        Q(i)=CommunityStructureLouvian(MaxXcorr{b}{w}(:,:,i),0);
            %                       % title(['Trl' num2str(i) ',Q:' num2str(Q(i))]);
            %                      %  drawnow;pause(0);
            %                    end
            %                    subplot(4,4,w);hold on
            %                    plot(BlockSpec.TrialOrder,Q,'color',Col(b,:))
            %                end
            %            end
            %            %%
            %            Fig_Xcorr=cell(1,length(AnalysisData.MeanCluster));
            %           [Fig_Xcorr{:}]=obj.PlotXcorrHChPair(XcorrPairBlck,T_XcorrPairBlck,BlockSpec.TrialOrder,ExpBlkNum,Bhv_blk_Perf) ;
            %            varargout=[varargout Wfigs Fig_Xcorr];
        end
        function  [RawData_exmp_tot,Wavelet_exmp_tot,f_Linear]=GrabAllTrialData(obj,BlockSpec,varargin) % grabs raw data and wavelet data 
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            obj.UseDataPointer=0; % use data pointer to load data
            AnalysisOpts.CurrentCh='';%AnalysisOpts.Ch(end);
            
            % this is to take the raw data for now we don't want it 
%             if ~obj.ManData.IsVarExistinFileAnalysis([],'RawData_exmp_tot','TrialData') | obj.OverWrite
%                 % get raw data now
%                 [RawData,FsLFP,DataPointer]=obj.TlBxFuncs.getRawEphysData(AnalysisData.Ch,obj.UseDataPointer);
%                 % grab Raw data
%                 [RawData_exmp_tot]=obj.GrabHTimeTrial(RawData,BlockSpec,'Fs',FsLFP,'UseDataPointer',obj.UseDataPointer,'DataPointerVar','RawData');
%                 obj.ManData.SaveVar([],RawData_exmp_tot,'RawData_exmp_tot','TrialData')
%                 clear RawData
%             else
%                 RawData_exmp_tot=obj.ManData.LoadVar([],'RawData_exmp_tot','TrialData',obj.UseDataPointer);
%             end
            
            % grab wavelet data
            if ~obj.ManData.IsVarExistinFileAnalysis(obj.AnalysisPathName,'Wavelet_exmp_tot','TrialData_Wave') | obj.OverWrite
                [WaveletData,f_Linear,FsWave,DataPointer]=obj.TlBxFuncs.getWaveletData(AnalysisData.Ch,obj.UseDataPointer);
                [Wavelet_exmp_tot]=obj.GrabHTimeTrial(WaveletData,BlockSpec,'Fs',FsWave,'UseDataPointer',obj.UseDataPointer,'DataPointerVar','Wavelet_Linear');
                AnalysisOpts.CurrentCh='';
                obj.ManData.DeleteFile(obj.AnalysisPathName,'TrialData_Wave',1);% delete previous file
                obj.ManData.SaveVar(obj.AnalysisPathName,Wavelet_exmp_tot,'Wavelet_exmp_tot','TrialData_Wave')
                obj.ManData.SaveVar(obj.AnalysisPathName,FsWave,'FsWave','TrialData_Wave')
                obj.ManData.SaveVar(obj.AnalysisPathName,f_Linear,'f_Linear','TrialData_Wave')
                obj.ManData.SaveVar(obj.AnalysisPathName,BlockSpec,'BlockSpec','TrialData_Wave')
                clear WaveletData
            else
                Wavelet_exmp_tot=obj.ManData.LoadVar(obj.AnalysisPathName,'Wavelet_exmp_tot','TrialData_Wave',obj.UseDataPointer);
                FsWave=obj.ManData.LoadVar(obj.AnalysisPathName,'FsWave','TrialData_Wave',obj.UseDataPointer);
                f_Linear=obj.ManData.LoadVar(obj.AnalysisPathName,'f_Linear','TrialData_Wave',obj.UseDataPointer);
            end
        end
        function  [H_exmp_tot]=GrabAllHTrialData(obj,H,BlockSpec,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            obj.UseDataPointer=0; % use data pointer to load data
            AnalysisOpts.CurrentCh='';%AnalysisOpts.Ch(end);
            
            if ~obj.ManData.IsVarExistinFileAnalysis([],'H_exmp_tot','HTrialData') | obj.OverWrite
                if isempty(H) % H doesn't esit so we aren't doing anlaysis on this
                    H_exmp_tot=[];
                else
                    % get H data first
                    [H_exmp_tot,BlkTrialNumber,ExpBlkNum,Bhv_blk_Perf]=obj.GrabHTimeTrial(H,BlockSpec,'Fs',AnalysisOpts.MotifAnalysis.FsWaveTarg);
                    obj.ManData.DeleteFile([],'HTrialData',1);% delete previous file
                    obj.ManData.SaveVar([],H_exmp_tot,'H_exmp_tot','HTrialData')
                    obj.ManData.SaveVar([],BlockSpec,'BlockSpec','HTrialData')
                end
            else
                H_exmp_tot=obj.ManData.LoadVar([],'H_exmp_tot','HTrialData',obj.UseDataPointer);
            end        
        end
        
        function varargout=PlotCoherenceStatsPerTrl(obj,H_exmp_tot,RawData,varargin) % calculates coherence stats for single trials
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            NWs=size(H_exmp_tot.TargetOn{1,1},1);
            NTrl=size(H_exmp_tot.TargetOn{1,1},3);
            MaxFreq=80; %Hz what is the max freq we want to look at
            FsLFP=1000; % sampling frequency of LFP
            FsWave=50;
            NCol=3;NRow=2;
            NCol2=11;NRow2=4;
            AnalysisData.RawTime=-1:1/FsLFP:(0.5-1/FsLFP);
            AnalysisData.Time   =-1:1/FsWave:(0.5-1/FsWave);
            MotfCol=distinguishable_colors(NWs);
            varargout=[];
            % first let look at xcorr of Hs for this trial
            trl=1;
            [XcorrH,XcorrHavg,T_XcorrH]= obj.CalXcorrHChPair(H_exmp_tot.TargetOn{1,6},H_exmp_tot.TargetOn{2,6},'Navg',5,'Nshuffle',10);
            % calculate coherence for LFP signal now
            RawData1=RawData.TargetOn{1,6};RawData2=RawData.TargetOn{2,6};
            %  [mscohPair,cpsdPair,wavcohPair]=obj.CalCoherenceChPair(RawData1,RawData2,FsLFP,10);
            [mscohPair]=obj.CalCoherenceChPair(RawData1,RawData2,FsLFP,5);
            for trl=1:NTrl
                ThisTrlFigs=obj.FigParams.RenderFigure(2,[]);
                %plot the raw data
                figure(ThisTrlFigs{1})
                subplot(NRow,NCol,1);hold on
                plot(AnalysisData.RawTime,RawData1(1,:,trl),'b','linewidth',obj.FigParams.line_width);
                plot(AnalysisData.RawTime,RawData2(1,:,trl),'r','linewidth',obj.FigParams.line_width);
                xlabel('Time Stim Onset');ylabel('uV');title('LFP data');
                % plot ms coherence with phase
                coh_f=mscohPair.f;
                f_ind=coh_f<MaxFreq;
                subplot(NRow,NCol,2)
                plot(coh_f(f_ind),mscohPair.cohr{trl}(f_ind),'linewidth',obj.FigParams.line_width)
                xlabel('Freq');ylabel('Coherence');title('Mag Sqr Coh');
                % plot phase of corss-spectrum
                %                 cpsd_f=cpsdPair.f{1};
                %                 cpsdPhase=angle(cpsdPair.cpsd{1}{trl})/pi;
                %                 cpsd_f_ind=cpsd_f<MaxFreq;
                coh_phase=mscohPair.phase{trl}(f_ind);
                subplot(NRow,NCol,3)
                plot(coh_f(f_ind),coh_phase,'linewidth',obj.FigParams.line_width)
                xlabel('Freq');ylabel('X spect ang(rad)');title('Coh Phase');
                % plot wavelet coherence during trial
                %                 subplot(NRow,NCol,4)
                %                 wavcoh_f=wavcohPair.f{1};
                %                 wavcoh_f_ind=wavcoh_f<MaxFreq;
                %                 thiswavcohPair=wavcohPair.wavcoh{1}{trl};
                %                 helperPlotCoherence(thiswavcohPair(wavcoh_f_ind,:),AnalysisData.RawTime,wavcoh_f(wavcoh_f_ind),...
                %                 wavcohPair.coi{1}{trl},'Time from Stim Onset','Wave Coh')
                % now calculate the value of motif xcorr at lag zero and report
                % it here
                ZeroLagInd=T_XcorrH{1}==0;
                XcorrZeroLag=arrayfun(@(x) XcorrH{x}(ZeroLagInd,trl),1:NWs);
                subplot(NRow,NCol,5);hold on
                arrayfun(@(x) bar(x,XcorrZeroLag(x),'FaceColor',MotfCol(x,:)),1:NWs)
                xlabel('Motif #');ylabel('corr');title('Zerolag Xcorr')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                % start the second figure with details on each motif
                figure(ThisTrlFigs{2})
                % plot xcorr for each of the motifs now so we make sure
                % everything is consistent
                for w=1:NWs
                    subplot(NRow2,NCol2,w)
                    plot(T_XcorrH{1},XcorrH{w}(:,trl),'color',MotfCol(w,:),'linewidth',obj.FigParams.line_width)
                    ylim([-0.2 1]);ylim([T_XcorrH{1}(1),T_XcorrH{1}(end)])
                    xlabel('Time');ylabel('XCorr');title(['Xcorr mtf ' num2str(w)])
                    
                    subplot(NRow2,NCol2,w+1)
                    helperCWTTimeFreqPlot(AnalysisData.CoreMotifs{w} ,(1:size(AnalysisData.CoreMotifs{w},2))/FsWave,AnalysisData.cwt_f,'justplot1',['Mtf ' num2str(w)],'Time(s)','Freq',obj.LogScale);
                end
                
                varargout=[varargout ThisTrlFigs];
            end
            
        end
        function  varargout=CalCoherenceStatsRec(obj,H,BlockSpec,ChPairNum,varargin) % calcultes coherence for recording
            
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            % grab all of the trial data
            [RawData_exmp_tot,Wavelet_exmp_tot,f_Linear]=obj.GrabAllTrialData(H,BlockSpec);
            obj.CalCoherenceAllChPairs(BlockSpec,H_exmp_tot,RawData_exmp_tot,Wavelet_exmp_tot)
            
        end
        function  varargout=CalHCoherenceStatsRec(obj,H,BlockSpec,ChPairNum,varargin) % calcultes coherence for recording
            
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            if AnalysisOpts.MotifAnalysis.ProcessStep==1
                % grab all of the trial data
                [H_exmp_tot]=obj.GrabAllHTrialData(H,BlockSpec,'OverWrite',1);
            elseif AnalysisOpts.MotifAnalysis.ProcessStep==2
                % load data and calculate coherence
                [H_exmp_tot]=obj.GrabAllHTrialData(H,BlockSpec,'OverWrite',0);
                obj.CalHCoherenceAllChPairs(BlockSpec,H_exmp_tot)
            end
        end
        
        function varargout=PlotCoherenceStatsRec(obj,DateNums,varargin) % Plots coherence stats for different recording dates
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % grab data from all of the recordings
            [RecData,RecSpecs]=obj.ManData.LoadDataFromRec(DateNums,'','wavcohPair','','CohData','All',0);
            %% rearrange data based on recording and pair area
            PairIndsAll.FrontParietPairs=[];PairIndsAll.FEFPairs=[]; % initialize pairs
            PairIndsAll.PFCPairs=[];PairIndsAll.LIPPairs=[];         % initialize pairs
            PairIndsAll.AllPairs=[];PairIndsAll.PairArea=[];wavcohPairAll=[];
            for Rec=1:length(RecData)
                ChPairs=nchoosek(1:length(RecSpecs(Rec).ChsSet),2);
                ChAr=RecSpecs(Rec).ChannelArea(RecSpecs(Rec).ChsSet);
                NPairs=size(ChPairs,1);
                PairArea=[];
                PairArea(:,1)=ChAr(ChPairs(:,1));
                PairArea(:,2)=ChAr(ChPairs(:,2));
                
                ThisRecAllPairs=1:NPairs;
                ThisRecFrontParietPairs=find(((PairArea(:,1)==1 | PairArea(:,1)==4) & PairArea(:,2)==5) |((PairArea(:,2)==1 | PairArea(:,2)==4) & PairArea(:,1)==5))';
                ThisRecFEFPairs=find((PairArea(:,1)==4) & PairArea(:,2)==4)';
                ThisRecPFCPairs=find((PairArea(:,1)==1) & PairArea(:,2)==1)';
                ThisRecLIPPairs=find((PairArea(:,1)==5) & PairArea(:,2)==5)';
                
                PairIndsAll.FrontParietPairs=cat(2,PairIndsAll.FrontParietPairs,ThisRecFrontParietPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.FEFPairs=cat(2,PairIndsAll.FEFPairs,ThisRecFEFPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.PFCPairs=cat(2,PairIndsAll.PFCPairs,ThisRecPFCPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.LIPPairs=cat(2,PairIndsAll.LIPPairs,ThisRecLIPPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.AllPairs=cat(2,PairIndsAll.AllPairs,ThisRecAllPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.PairArea=cat(1,PairIndsAll.PairArea,PairArea);
                % load all of the data into a big matrix
                wavcohPairAll=cat(2,wavcohPairAll,RecData(Rec).out);
            end
            %% now you can just run the analysis similar to what you would have done with area
            obj.PlotCoherenceStatsArea(PairIndsAll,wavcohPairAll)
            
        end
        %% traditional coehrence and selectivity analysis 
        function varargout=PlotSelectivityStatsRec(obj,DateNums,varargin) % Plots selectivity stats for different recording dates
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % grab data from all of the recordings
            AnalysisOpts.CurrentCh='';
            [RecData,RecSpecs]=obj.ManData.LoadDataFromRec(DateNums,obj.AnalysisPathName,'Wavelet_exmp_tot','','TrialData_Wave',[],0);
            AnalysisOpts.CurrentCh='';
            % load BlockSpec from first recording
         %   temp=obj.ManData.LoadDataFromRec(DateNums(1),obj.AnalysisPathName,'BlockSpec','','TrialData_Wave','',0);
          %  AnalysisData.BlockSpec=temp.out;
            AnalysisData.RecSpecs=RecSpecs;
            %% rearrange data based on recording and  area
            ChsIndsAll.PFCChs=[]; % initialize Chs
            ChsIndsAll.FEFChs=[];ChsIndsAll.LIPChs=[];
            ChsIndsAll.AllChs=[];HTrialAll=[];ChsIndsAll.RecInd=[];
            for Rec=1:length(RecData)
                ChAr=RecSpecs(Rec).ChannelArea(RecSpecs(Rec).ChsSet);
                
                ThisRecPFCChs=find(ChAr==1)';
                ThisRecFEFChs=find(ChAr==4)';
                ThisRecLIPChs=find(ChAr==5)';
                ThisRecAllChs=1:length(ChAr);
                
                ChsIndsAll.PFCChs=cat(2,ChsIndsAll.PFCChs,ThisRecPFCChs+length(ChsIndsAll.AllChs));
                ChsIndsAll.FEFChs=cat(2,ChsIndsAll.FEFChs,ThisRecFEFChs+length(ChsIndsAll.AllChs));
                ChsIndsAll.LIPChs=cat(2,ChsIndsAll.LIPChs,ThisRecLIPChs+length(ChsIndsAll.AllChs));
                ChsIndsAll.AllChs=cat(2,ChsIndsAll.AllChs,ThisRecAllChs+length(ChsIndsAll.AllChs));
                % save the indexes of recording date
                ChsIndsAll.RecInd=cat(2,ChsIndsAll.RecInd,Rec*ones(1,length(ChAr)));
                % load all of the data into a big matrix
                HTrialAll=cat(2,HTrialAll,RecData(Rec).out);
            end
            % clear extra data to save space
            clear RecData
            %% now you can just run the analysis similar to what you would have done with area
            % obj.PlotHSelectivityStatsAreaBhvModel(ChsIndsAll,HTrialAll)
            obj.PlotCorrStatsAreaBhvModel(ChsIndsAll,HTrialAll)
        end
        function varargout=PlotCorrStatsAreaBhvModel(obj,ChsIndsAll,HTrialAll,varargin) % has BhvModel-correlation stats for traditional analysis
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % set-up options
            opts.TimePeriodName={'FixOn'};
            opts.AreaToPlot= {'PFCChs','FEFChs','LIPChs'};
            opts.AnalysisFields={'precision_ind','glmfit_Bias','glmfit_RPE','glmfit_precision','glmfit_value_for_choice_ro_mean'};%'value_for_choice_Info','RPECorr',
            opts.AnalysisFieldsTitle={'Precision','glm Bias','glm RPE','glm precision','glm value for choice'};%'Value','RPE'
            opts.FreqBox={[0 6]   ,[7 21]  ,[22 34] ,[35 55] ,[60 100]}; % frequencies to look for power
            opts.TimeBox={[-0.5 0],[-0.5 0],[-0.5 0],[-0.5 0],[-0.5 0]}; % time to look for power
          % opts.TimeBox={[0 0.2],[0 0.2],[0 0.2],[0 0.2],[0 0.2]};
            opts.WantedTrialTim=[-0.7 0.2]; 
            opts.Time=obj.GenerateTimeAxis(-1,0.5,AnalysisOpts.MotifAnalysis.FsWave);
            opts.WantedTimeInd=opts.Time>=opts.WantedTrialTim(1) & opts.Time<opts.WantedTrialTim(2);
            opts.WantedTime=opts.Time(opts.WantedTimeInd);
            opts.TrialOrder=AnalysisData.BlockSpec.BhvModelInfo.BlockSpec.TrialOrder;
            opts.NavgTimeEpoch=10;
            opts.NRow=length(opts.AreaToPlot);opts.NCol=length(opts.FreqBox)+1;
            
            FEFChs=ChsIndsAll.FEFChs;
            PFCChs=ChsIndsAll.PFCChs;
            LIPChs=ChsIndsAll.LIPChs;
            AllChs=ChsIndsAll.AllChs;
             %% prepare all of the data ( put the data for each channel into one cell)
            for f=opts.TimePeriodName % loop on time period
                HChDataAll.(f{1})  =obj.ManData.ReshapeCell2Mat(HTrialAll,52,f{1});
                for ChInd=1:length(AllChs) % loop on channel
                    fprintf('\nProcessing info for ch:%i time period:%s',ChInd,f{1});
                    ThisChRec=ChsIndsAll.RecInd(ChInd);
                    ThisRecSwchBlks=AnalysisData.RecSpecs(ThisChRec).BlockSpec.BhvModelInfo.model_outputs.Total_Switch;
                    %    SwchBlks=find(ThisRecSwchBlks);NoSwchBlks=find(ThisRecSwchBlks==0);
                    Blks=1:length(ThisRecSwchBlks);
                    % get the info for blocks where there is a switch
                    [CorrVals.(f{1}){ChInd}]=obj.CalCorrBlockChBhvModel_trad(HChDataAll.(f{1}){ChInd},...
                        AnalysisOpts.CohAnalysis.NAvgCoh,opts.Time,AnalysisData.RecSpecs(ThisChRec).TrialTimes,...
                        AnalysisData.RecSpecs(ThisChRec).BlockSpec.BhvModelInfo.BlockSpec,Blks);
                end                        
            end
             
          %  load('SelectivityAnalysis_Beaker_181115_TradCorrModelExampleFig.mat');
            % save this var for now
       %     obj.ManData.SaveVar('SelectivityAnalysis',CorrVals,'CorrVals','TradCorrModel');
       %     obj.ManData.SaveVar('SelectivityAnalysis',ChsIndsAll,'ChsIndsAll','TradCorrModel');
            
            % reorganize into full matrix
            CorrVals_ro.(f{1})=obj.ManData.ReshapeCellStruct2Mat(CorrVals.(f{1}),'glmfit_value_for_choice','precision_ind','glmfit_Bias','glmfit_RPE','glmfit_precision'); 
            TempMat=squeeze(CorrVals_ro.(f{1}).glmfit_value_for_choice);
            CorrVals_ro.(f{1}).glmfit_value_for_choice_ro=arrayfun(@(x) obj.ManData.ReshapeCell2Mat(TempMat(x,:),3),1:size(TempMat,1),'UniformOutput',0);
            CorrVals_ro.(f{1}).glmfit_value_for_choice_ro_mean=mean(obj.ManData.ReshapeCell2Mat(CorrVals_ro.(f{1}).glmfit_value_for_choice_ro,4),4);
            % plot the results
            FieldCOL=distinguishable_colors(length(opts.AnalysisFields));
            NRow=2;NCol=length(opts.AreaToPlot);
            FigNum=1;% plot Xcorr values first
            for field=opts.AnalysisFields
                f=find(strcmp(field{1},opts.AnalysisFields));
                for TimPeriod=1:length(opts.TimePeriodName)
                    ThisTimPrd=opts.TimePeriodName{TimPeriod};
                    CorrFig(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                 %  PowFig(FigNum) =obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]); % power fig                    
                    for i=1:length(opts.AreaToPlot) % loop on areas to plot
                           eval(['ThisChs=' opts.AreaToPlot{i} ';']);                                                        
                           %% plot correlation to model variables
                           figure(CorrFig{FigNum})
                           subplot(NRow,NCol,i)
                           FieldTitle=opts.AnalysisFieldsTitle{f};
                           % take average for these channels 
                           if ~iscell(CorrVals_ro.(ThisTimPrd).(field{1}))
                               ThisCorrMat=CorrVals_ro.(ThisTimPrd).(field{1})(:,opts.WantedTimeInd,ThisChs);
                           else
                               
                           end
                           % perform statistical test on each value
                           [SigTest,p_SigTest]=obj.ManData.StatTeston3Dmat(ThisCorrMat,0,'right');
                           MeanCorr=nanmean(ThisCorrMat,3);
                           % MeanCorr=imgaussfilt(MeanCorr,5);
                           % plot the results in an image
                           helperCWTTimeFreqPlot(MeanCorr.*SigTest,opts.WantedTime,AnalysisData.cwt_f,'justplot1',...
                               [FieldTitle ' ' opts.AreaToPlot{i} ',n=' num2str(length(ThisChs))] ,...
                               ['Time from ' ThisTimPrd] ,'Frequency(Hz)',0)
                           colorbar
                           axis square
                           % plot data without any sig test
                           subplot(NRow,NCol,NCol+i)
                           helperCWTTimeFreqPlot(MeanCorr,opts.WantedTime,AnalysisData.cwt_f,'justplot1',...
                               [FieldTitle ' ' opts.AreaToPlot{i} ',n=' num2str(length(ThisChs))] ,...
                               ['Time from ' ThisTimPrd] ,'Frequency(Hz)',0)
                           colorbar
                           axis square                                                      
                           %% plot average power for different frequency
                           %% bands as well
%                            if ~isempty(ThisChs)
%                                 figure(PowFig{FigNum})
%                                 obj.PlotPowerAroundSwitch(CorrVals.(ThisTimPrd)(ThisChs),ChsIndsAll,ThisChs,opts,ThisTimPrd,...
%                                     ['Norm Power ' opts.AreaToPlot{i}],PowFig{FigNum},i)
%                             end          
                   end
                  %  legend(opts.AreaToPlot,'Location','best','FontSize',5);    
                  FigNum=FigNum+1;
                end
            end
            % save the figures
            AnalysisOpts.CurrentCh='';
            [~,~,AnalysisData.SelectivityFig]=obj.ManData.GetFileName('SelectivityAnalysis','TradCorrModelGLMFixOn','SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.SaveFigSeries([],AnalysisData.SelectivityFig,[CorrFig],'SaveEachFrame',1);
            varargout=[CorrFig ];           
        end
        function [varargout]=PlotPowerAroundSwitch(obj,PowerData,ChsIndsAll,ThisChs,opts,AxisTitle,Title,Figh,StrRow)
            %plots change of power during leaning for differnet channels
            %and areas
            global AnalysisOpts AnalysisData
            
            opts.Col={'r','b','k','m','g','m'};
            
            %% preppare parameters
            TimeAxis=opts.Time;
            FreqAxis=AnalysisData.cwt_f;
            %% prepare power data 
            % normalize power data (freq,time,trial,blk,ch)
            % PowerData=PowerData.*FreqAxis';
            % normalize to max power for each electrode now
            %NormPowerData=arrayfun(@(x) PowerData(:,:,:,:,x)/max(PowerData(:,:,:,:,x),[],'all'),1:size(PowerData,5),'UniformOutput',0);
            %NormPowerData=obj.ManData.ReshapeCell2Mat(NormPowerData,5);
            %% normalize power by 1/f and max
%             NormPowerData=arrayfun(@(x) PowerData{x}.CwtPower.*FreqAxis',1:length(PowerData),...
%                 'UniformOutput',0);
%             NormPowerData=arrayfun(@(x) NormPowerData{x}/max(NormPowerData{x},[],'all'),1:length(PowerData),...
%                 'UniformOutput',0);
%             AvgNormPowerCh=obj.ManData.ReshapeCell2Mat(cellfun(@(x) nanmean(nanmean(x,4),3),NormPowerData,'UniformOutput',0),3);            
%             AvgNormPowerData=nanmean(AvgNormPowerCh,3);
            %% normalize by percent modulation from the mean power .*FreqAxis' 
            fnormPowerData=arrayfun(@(x) PowerData{x}.CwtPower.*FreqAxis,1:length(PowerData),...
                 'UniformOutput',0);  % normalize by f
            AvgNormPowerCh=obj.ManData.ReshapeCell2Mat(arrayfun(@(x) nanmean(nanmean(nanmean(fnormPowerData{x},4),3),2),...
                1:length(PowerData),'UniformOutput',0),3);            
             NormPowerData=arrayfun(@(x) (fnormPowerData{x}-repmat(AvgNormPowerCh(:,:,x),[1 size(fnormPowerData{x},2) size(fnormPowerData{x},3) size(fnormPowerData{x},4)]))./(repmat(AvgNormPowerCh(:,:,x),[1 size(fnormPowerData{x},2) size(fnormPowerData{x},3) size(fnormPowerData{x},4)])),...
                1:length(PowerData),'UniformOutput',0);
            AvgNormPowerData=mean(obj.ManData.ReshapeCell2Mat(arrayfun(@(x) nanmean(nanmean(NormPowerData{x},4),3),...
                1:length(PowerData),'UniformOutput',0),3),3);
            %% get the switch data for each channel and cut data into switch and non switch 
            ThisRecSwchBlks=arrayfun(@(x) AnalysisData.RecSpecs(ChsIndsAll.RecInd(x)).BlockSpec.BhvModelInfo.model_outputs.Total_Switch,ThisChs,'UniformOutput',0);
            NormPowerDataSwch=arrayfun(@(x) NormPowerData{x}(:,:,:,find(ThisRecSwchBlks{x})),1:length(ThisChs),'UniformOutput',0);
            NormPowerDataNoSwch=arrayfun(@(x) NormPowerData{x}(:,:,:,find(ThisRecSwchBlks{x}==0)),1:length(ThisChs),'UniformOutput',0);
            %% compute average in the time boxed period
            NBoxs=length(opts.TimeBox);
            for box=1:NBoxs
                ThisTimeInd{box}=TimeAxis>=opts.TimeBox{box}(1) & TimeAxis<=opts.TimeBox{box}(2);
                ThisFreqInd{box}=FreqAxis>=opts.FreqBox{box}(1) & FreqAxis<=opts.FreqBox{box}(2);
                %TimBoxTrls{box}=NormPowerData(ThisFreqInd{box},ThisTimeInd{box},:,:,:);
                %AvgTimBoxTrls{box}=squeeze(mean(mean(TimBoxTrls{box},2),1));
                TimBoxTrls.Swch{box}=cell2mat(arrayfun(@(x) squeeze(nanmean(nanmean(nanmean(NormPowerDataSwch{x}(ThisFreqInd{box},ThisTimeInd{box},:,:),4),2),1)),...
                    1:length(ThisChs),'UniformOutput',0))';
                TimBoxTrls.NoSwch{box}=cell2mat(arrayfun(@(x) squeeze(nanmean(nanmean(nanmean(NormPowerDataNoSwch{x}(ThisFreqInd{box},ThisTimeInd{box},:,:),4),2),1)),...
                    1:length(ThisChs),'UniformOutput',0))';               
               % AvgTimBoxTrls{box}=squeeze(nanmean(TimBoxTrls{box},[1 2 4]))';
            end           
            %% figure plot for each pair all of the information across blcoks
            subplot(opts.NRow,opts.NCol,1+opts.NCol*(StrRow-1))
            % plot average wave coherence for all of the trials and pairs
            % to find time boxing
            helperCWTTimeFreqPlot(AvgNormPowerData,TimeAxis,FreqAxis,'image',[Title ],['Time(s) from ' AxisTitle],'Freq',0);
            colorbar
            % plot the rectangles now
            hold on
            arrayfun(@(x) rectangle('Position',[opts.TimeBox{x}(1),opts.FreqBox{x}(1),opts.TimeBox{x}(2)-opts.TimeBox{x}(1),opts.FreqBox{x}(2)-opts.FreqBox{x}(1)],...
                'EdgeColor',opts.Col{x},'LineWidth',3),1:NBoxs);
            % chi=get(gca, 'Children');
            % Reverse the stacking order so that the patch overlays the line
            % set(gca, 'Children',flipud(chi));
            %% plot avg of each time box with respect to trials
            opts.TrialOrderNew=opts.TrialOrder(1:end-opts.NavgTimeEpoch)+opts.NavgTimeEpoch-1;
            for box=1:NBoxs
                subplot(opts.NRow,opts.NCol,1+box+opts.NCol*(StrRow-1))
                Title=sprintf('%i-%iHz',opts.FreqBox{box}(1),opts.FreqBox{box}(2));
                MeanBoxSwch=movmean(TimBoxTrls.Swch{box},opts.NavgTimeEpoch,2);
                MeanBoxNoSwch=movmean(TimBoxTrls.NoSwch{box},opts.NavgTimeEpoch,2);
              %  obj.FigParams.PlotMeanStd(opts.TrialOrder,movmean(AvgTimBoxTrls{box},opts.NavgTimeEpoch,2),[],'Trl to switch','mean norm power',opts.Col{box},0,Title)
                obj.FigParams.PlotMeanStd(opts.TrialOrderNew,MeanBoxSwch(:,1:length(opts.TrialOrderNew)),[],'Trl to switch','mean norm power',opts.Col{box},1,Title)              
                hold on
                obj.FigParams.PlotMeanStd(opts.TrialOrderNew,MeanBoxNoSwch(:,1:length(opts.TrialOrderNew)),[],'Trl to switch','mean norm power',1,1,Title,'--')              
             %   v=axis;
             %   plot([-opts.NavgTimeEpoch -opts.NavgTimeEpoch],[v(3) v(4)],'r')
            end
            legend({'HardSwtch','SoftSwtch'},'Location','best');
            varargout{1}=Figh;
        end

        function [Vals]=CalCorrBlockChBhvModel(obj,HChData,NTrlAvg,Time,TrialTimes,BlockSpec,WantedBlks,varargin)% includes Bhv Model
            % calculates correlation of values of power or coherence to
            % different values of the model 
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            opts.NormalizeH=0; % devides H for each motif to it's max
            opts.TimPridLength=0.2; % lenght of time period we are looking after or before an event
            GrabTrialData=@(x,y)  [nan*ones(length(find(x<1)),size(y,2)) ;y(x(x>=1),:)]; % function to get data from trial times
            % NTrlAvg number of trials to be averaged for coherence stats
            HChData=HChData(:,:,:,WantedBlks);
            NMotif=size(HChData,1);
            Ntim=size(HChData,2); % time points
            NTrl=size(HChData,3);
            NBlks=size(HChData,4); % how many blocks we want to average% last blcok is the average of all (should be fixed)
            %% attention removing last block
            HChData=HChData(:,:,:,1:NBlks);
            Vals=[];
            for m=1:NMotif
                ThisMotif_HChData=squeeze(HChData(m,:,:,:)); % take all of the H channel data only for this motif
                if opts.NormalizeH % if we are normalizing H values
                    ThisMotif_HChData=ThisMotif_HChData/max(ThisMotif_HChData(:));
                end
                ThisMotif_HChDataAll=reshape(ThisMotif_HChData,[size(ThisMotif_HChData,1) size(ThisMotif_HChData,2)*size(ThisMotif_HChData,3)]);
                % loop in time points and trials calculate correlation to
                % RPE, confidence and value_for Choice
                RPE=cell2mat(BlockSpec.RPE(WantedBlks)');
                precision_ind=cell2mat(BlockSpec.precision_ind(WantedBlks)');
                value_for_choice=cell2mat(BlockSpec.down_sampled_Value_for_choice4(WantedBlks))';
                % cal info for RPE
                [Vals.RPECorr(m,:),Vals.RPECorr_p(m,:)]=corr(ThisMotif_HChDataAll',RPE,'rows','pairwise');
                Vals.RPEInfo(m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',RPE)',1:Ntim,'UniformOutput',0));
                % cal info for precision_ind
                [Vals.precision_ind(m,:),Vals.precision_ind_p(m,:)]=corr(ThisMotif_HChDataAll',precision_ind,'rows','pairwise');
                Vals.precision_ind_Info(m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',precision_ind)',1:Ntim,'UniformOutput',0));
                % cal info for belief
                for k=1:4
                    Vals.value_for_choice(k,m,:)=corr(ThisMotif_HChDataAll',value_for_choice(:,k),'rows','pairwise');
                    Vals.value_for_choice_Info(k,m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',value_for_choice(:,k))',1:Ntim,'UniformOutput',0));
                end
            end
            
        end
        function [Vals]=CalCorrBlockChBhvModel_trad(obj,ChData,NTrlAvg,Time,TrialTimes,BlockSpec,WantedBlks,varargin)% includes Bhv Model
            % calculates correlation of model values for using wavelet and
            % traditional methods
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            opts.NormalizeH=0; % devides H for each motif to it's max
            opts.TimPridLength=0.2; % lenght of time period we are looking after or before an event
            GrabTrialData=@(x,y)  [nan*ones(length(find(x<1)),size(y,2)) ;y(x(x>=1),:)]; % function to get data from trial times
            % NTrlAvg number of trials to be averaged for coherence stats
            ChData=ChData(:,:,:,WantedBlks);
            NFreq=size(ChData,1);
            Ntim=size(ChData,2); % time points
            NTrl=size(ChData,3);
            NBlks=size(ChData,4); % how many blocks we want to average% last blcok is the average of all (should be fixed)
            %% attention removing last block
            ChData=ChData(:,:,:,1:NBlks);
            % calculate power of data
            ChData_Pow=ChData.*conj(ChData);
            % take the norm mean power of data 
            ChData_Pow=obj.ManData.NormMeanPower(ChData_Pow,[]);%AnalysisData.cwt_f);
            %% loop in time points and trials calculate correlation to
            %% RPE, confidence and value_for Choice
            RPE=cell2mat(BlockSpec.RPE(WantedBlks)');
            precision_ind=cell2mat(BlockSpec.precision_ind(WantedBlks)');
            value_for_choice=cell2mat(BlockSpec.down_sampled_Value_for_choice4(WantedBlks))';
            value_for_choice_full=cell2mat(BlockSpec.Value_for_choice(WantedBlks))';
            % take the shifted values as well 
            RPE_sh=cell2mat(cellfun(@(x) x(1:end-1)',BlockSpec.RPE(WantedBlks),'UniformOutput',0))';
            precision_ind_sh=cell2mat(cellfun(@(x) x(1:end-1)',BlockSpec.precision_ind(WantedBlks),'UniformOutput',0))';
            value_for_choice_sh=cell2mat(cellfun(@(x) x(:,1:end-1),BlockSpec.down_sampled_Value_for_choice4(WantedBlks),'UniformOutput',0))';
           value_for_choice_full_sh=cell2mat(cellfun(@(x) x(:,1:end-1),BlockSpec.Value_for_choice(WantedBlks),'UniformOutput',0))';
            Vals=[];            
            for f=1:NFreq % loop on frequency 
                for t=1:Ntim % loop on time
                    
                    ThisTimeFreq=squeeze(ChData_Pow(f,t,2:end,:)); % take all of the H channel data only for this motif
                    ThisTimeFreq=reshape(ThisTimeFreq,[size(ThisTimeFreq,1)*size(ThisTimeFreq,2) 1]);
                    
                    %  cal info for RPE
                    %  [Vals.RPECorr(f,t),Vals.RPECorr_p(f,t)]=corr(ThisTimeFreq,RPE,'rows','pairwise');
                    %  Vals.RPEInfo(m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',RPE)',1:Ntim,'UniformOutput',0));
                    %  cal info for precision_ind
                    [Vals.precision_ind(f,t),Vals.precision_ind_p(f,t)]=corr(ThisTimeFreq,precision_ind_sh,'rows','pairwise');
                    %  Vals.precision_ind_Info(m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',precision_ind)',1:Ntim,'UniformOutput',0));
                    
                    %  cal info for belief
                    %  for k=1:4
                    %     Vals.value_for_choice(k,f,t)=corr(ThisTimeFreq,value_for_choice(:,k),'rows','pairwise');
                    %     Vals.value_for_choice_Info(k,m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',value_for_choice(:,k))',1:Ntim,'UniformOutput',0));
                    %  end
                    %% build a glm model and correlate it to time/freq
                    GlmFit(f,t,:)=glmfit([RPE_sh,precision_ind_sh,value_for_choice_sh],ThisTimeFreq);                    
                end
            end
            Vals.glmfit_Bias=GlmFit(:,:,1);
            Vals.glmfit_RPE=GlmFit(:,:,2);
            Vals.glmfit_precision=GlmFit(:,:,3);
            Vals.glmfit_value_for_choice=arrayfun(@(x) GlmFit(:,:,x),4:7,'UniformOutput',0);
            % output the power as well
            Vals.CwtPower=ChData_Pow;
        end
        %% Motifs coherence and selectivity analysis
        function varargout=PlotHCoherenceStatsRec(obj,DateNums,Animal,varargin) % Plots coherence stats for different recording dates
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % grab data from all of the recordings
            AnalysisOpts.CurrentCh='';
            [RecData,RecSpecs]=obj.ManData.LoadDataFromRec(DateNums,Animal,'','XcorrH','','HCohData','All',0);
            AnalysisOpts.CurrentCh='';
            % load BlockSpec and TXcorr from first recording
            temp=obj.ManData.LoadDataFromRec(DateNums(1),Animal,'','BlockSpec','','HTrialData','',0);
            AnalysisData.BlockSpec=temp.out;
            AnalysisOpts.CurrentCh='';
            temp=obj.ManData.LoadDataFromRec(DateNums(1),Animal,'','TXcorr','','HCohData_1','',0);
            AnalysisData.TXcorr=temp.out;
            %% rearrange data based on recording and pair area
            PairIndsAll.FrontParietPairs=[];PairIndsAll.FEFPairs=[]; % initialize pairs
            PairIndsAll.PFCPairs=[];PairIndsAll.LIPPairs=[];         % initialize pairs
            PairIndsAll.AllPairs=[];PairIndsAll.PairArea=[];XcorrPairAll=[];
            for Rec=1:length(RecData)
                ChPairs=nchoosek(1:length(RecSpecs(Rec).ChsSet),2);
                ChAr=RecSpecs(Rec).ChannelArea(RecSpecs(Rec).ChsSet);
                NPairs=size(ChPairs,1);
                PairArea=[];
                PairArea(:,1)=ChAr(ChPairs(:,1));
                PairArea(:,2)=ChAr(ChPairs(:,2));
                
                ThisRecAllPairs=1:NPairs;
                ThisRecFrontParietPairs=find(((PairArea(:,1)==1 | PairArea(:,1)==4) & PairArea(:,2)==5) |((PairArea(:,2)==1 | PairArea(:,2)==4) & PairArea(:,1)==5))';
                ThisRecFEFPairs=find((PairArea(:,1)==4) & PairArea(:,2)==4)';
                ThisRecPFCPairs=find((PairArea(:,1)==1) & PairArea(:,2)==1)';
                ThisRecLIPPairs=find((PairArea(:,1)==5) & PairArea(:,2)==5)';
                
                PairIndsAll.FrontParietPairs=cat(2,PairIndsAll.FrontParietPairs,ThisRecFrontParietPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.FEFPairs=cat(2,PairIndsAll.FEFPairs,ThisRecFEFPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.PFCPairs=cat(2,PairIndsAll.PFCPairs,ThisRecPFCPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.LIPPairs=cat(2,PairIndsAll.LIPPairs,ThisRecLIPPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.AllPairs=cat(2,PairIndsAll.AllPairs,ThisRecAllPairs+length(PairIndsAll.AllPairs));
                PairIndsAll.PairArea=cat(1,PairIndsAll.PairArea,PairArea);
                % load all of the data into a big matrix
                XcorrPairAll=cat(2,XcorrPairAll,RecData(Rec).out);
            end
            %% now you can just run the analysis similar to what you would have done with area
            obj.PlotHCoherenceStatsArea(PairIndsAll,XcorrPairAll)            
        end
        function varargout=PlotHSelectivityStatsRec(obj,DateNums,Animal,varargin) % Plots selectivity stats for different recording dates
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % grab data from all of the recordings
            AnalysisOpts.CurrentCh='';
            [RecData,RecSpecs]=obj.ManData.LoadDataFromRec(DateNums,'','H_exmp_tot','','HTrialData',[],0);
            AnalysisOpts.CurrentCh='';
            % load BlockSpec from first recording
            temp=obj.ManData.LoadDataFromRec(DateNums(1),'','BlockSpec','','HTrialData','',0);
            AnalysisData.BlockSpec=temp.out;
            AnalysisData.RecSpecs=RecSpecs;
            %% rearrange data based on recording and  area
            ChsIndsAll.PFCChs=[]; % initialize Chs
            ChsIndsAll.FEFChs=[];ChsIndsAll.LIPChs=[];
            ChsIndsAll.AllChs=[];HTrialAll=[];ChsIndsAll.RecInd=[];
            for Rec=1:length(RecData)
                ChAr=RecSpecs(Rec).ChannelArea(RecSpecs(Rec).ChsSet);
                
                ThisRecPFCChs=find(ChAr==1)';
                ThisRecFEFChs=find(ChAr==4)';
                ThisRecLIPChs=find(ChAr==5)';
                ThisRecAllChs=1:length(ChAr);
                
                ChsIndsAll.PFCChs=cat(2,ChsIndsAll.PFCChs,ThisRecPFCChs+length(ChsIndsAll.AllChs));
                ChsIndsAll.FEFChs=cat(2,ChsIndsAll.FEFChs,ThisRecFEFChs+length(ChsIndsAll.AllChs));
                ChsIndsAll.LIPChs=cat(2,ChsIndsAll.LIPChs,ThisRecLIPChs+length(ChsIndsAll.AllChs));
                ChsIndsAll.AllChs=cat(2,ChsIndsAll.AllChs,ThisRecAllChs+length(ChsIndsAll.AllChs));
                % save the indexes of recording date
                ChsIndsAll.RecInd=cat(2,ChsIndsAll.RecInd,Rec*ones(1,length(ChAr)));
                % load all of the data into a big matrix
                HTrialAll=cat(2,HTrialAll,RecData(Rec).out);
            end
            % clear extra data to save space
            clear RecData
            %% now you can just run the analysis similar to what you would have done with area
            % obj.PlotHSelectivityStatsAreaBhvModel(ChsIndsAll,HTrialAll)
            obj.PlotHSCorrStatsAreaBhvModel(ChsIndsAll,HTrialAll)
        end
        function varargout=PlotHSelectivityStatsAreaBhvModel(obj,ChsIndsAll,HTrialAll,varargin) % has BhvModel-Selectivity stats for Motifs
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            FEFChs=ChsIndsAll.FEFChs;
            PFCChs=ChsIndsAll.PFCChs;
            LIPChs=ChsIndsAll.LIPChs;
            AllChs=ChsIndsAll.AllChs;
            NCoreMotifs=length(AnalysisData.CoreMotifs);
            %% prepare all of the data ( put the data for each channel into one cell)
            HChDataAll.Reward=obj.ManData.ReshapeCell2Mat(HTrialAll,52,'TargetOn');
            %% precalculate all of the selectivities we care about
            Ts=1/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            Time=-1:Ts:(0.5-Ts); %% what was atually in the experiment
            Vals=[];
            for ChInd=1:length(AllChs)
                ThisChRec=ChsIndsAll.RecInd(ChInd);
                ThisRecSwchBlks=AnalysisData.RecSpecs(ThisChRec).BlockSpec.BhvModelInfo.model_outputs.Total_Switch;
                SwchBlks=find(ThisRecSwchBlks);NoSwchBlks=find(ThisRecSwchBlks==0);
                % get the info for blocks where there is a switch
                
                %                 [CorrVals{ChInd}]=obj.CalHCorrBlockChBhvModel(HChDataAll.Reward{ChInd},...
                %                     AnalysisOpts.CohAnalysis.NAvgCoh,Time,AnalysisData.RecSpecs(ThisChRec).TrialTimes,...
                %                     AnalysisData.RecSpecs(ThisChRec).BlockSpec.BhvModelInfo.BlockSpec,SwchBlks);
                %
                [Vals.SwchBlks{ChInd},Params.SwchBlks]=obj.CalHSelectivityBlockChBhvModel(HChDataAll.Reward{ChInd},...
                    AnalysisOpts.CohAnalysis.NAvgCoh,Time,AnalysisData.RecSpecs(ThisChRec).TrialTimes,...
                    AnalysisData.RecSpecs(ThisChRec).BlockSpec.BhvModelInfo.BlockSpec,SwchBlks);
                %get the info for blocks where there is no switch
                [Vals.NoSwchBlks{ChInd},Params.NoSwchBlks]=obj.CalHSelectivityBlockChBhvModel(HChDataAll.Reward{ChInd},...
                    AnalysisOpts.CohAnalysis.NAvgCoh,Time,AnalysisData.RecSpecs(ThisChRec).TrialTimes,...
                    AnalysisData.RecSpecs(ThisChRec).BlockSpec.BhvModelInfo.BlockSpec,NoSwchBlks);
            end
            %              %% grab Color data for all of trails and channels
            %             % plot color data now
            %             ColorFigs=obj.PlotW_H_EachW_Condition(W,H_tot_TargOn,AnalysisData.cwt_f,AnalysisData.Time,TrialTimes,AnalysisData.ChArea,'BestColor',ColList);
            %             % grab H values based on 35 before and after switch
            %             [H_exmp_tot,BlkTrialNumber,ExpBlkNum,Bhv_blk_Perf]=obj.GrabHTimeTrial(H,BlockSpec);
            %             obj.PlotW_H_Trial_Condition(W,H_exmp_tot.Reward,AnalysisData.cwt_f,AnalysisData.Time,BlockSpec,AnalysisData.ChArea,[],[]);
            %             Wfigs=cell(1,length(AnalysisData.MeanCluster));
            %             %% area based Analysis
            %             % plot W-H Area
            %             [Wfigs{:}]=obj.PlotW_H_Trial_Area(W,H_exmp_tot.FixOn(:,WBlks),AnalysisData.cwt_f,Time,BlkTrialNumber(WBlks),BlockSpec.TrialOrder,ExpBlkNum(WBlks),AnalysisData.ChArea,'ChNum',NaN) ; % plot Ws
            %% calculate some of the vars we need
            opts.Ntrls2Switch2Look=10; % How many trials around the switch we want to look at
            trlSeq=AnalysisData.BlockSpec.TrialOrder;
            opts.SwitchTrl=find(trlSeq==0);
            opts.SwitchTrlAvg=find(Params.SwchBlks.TrlAvgInd(:,2)==opts.SwitchTrl);
            opts.NAvg=opts.SwitchTrl-opts.SwitchTrlAvg+1;
            opts.trlSeqAvg=trlSeq(1:end-opts.NAvg+1);
            % get the new time seq average
            opts.TrlSeqAvg=trlSeq(Params.SwchBlks.TrlAvgInd(:,2));
            opts.TimSeqAvg=Time(Params.SwchBlks.TimAvgInd(:,2));
            NCol=6;
            WantedTime=Time>=0 & Time<=0.2;
            %  opts.trl2lookInd=(SwitchTrl-opts.Ntrls2Switch2Look):(SwitchTrl+opts.Ntrls2Switch2Look-1-4); % trials to look at
            % opts.trlSeq2look=trlSeq;%-opts.Ntrls2Switch2Look:opts.Ntrls2Switch2Look;
            TimePeriodName={'Reward'};
            AreaToPlot={'PFCChs','FEFChs','LIPChs','AllChs'};
            AnalysisFields={'TimeAvg','TrlAvg'};
            FieldCOL=distinguishable_colors(length(AnalysisFields));
            NRow=length(AreaToPlot);
            FigNum=1;% plot Xcorr values first
            for TimPeriod=1:length(TimePeriodName)
                ThisTimPrd=TimePeriodName{TimPeriod};
                for m=1:NCoreMotifs
                    Figs(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    obj.PlotCoreMotifs([],m,1); % plot this motif
                    Figs(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    for i=1:length(AreaToPlot)
                        eval(['ThisChs=' AreaToPlot{i} ';']);
                        % plot fix on
                        for f=1%:length(AnalysisFields)
                            
                            temp=arrayfun(@(x) Vals.SwchBlks{x}(m).AvgTimTrl,ThisChs,'UniformOutput',0);
                            AvgHChs.TimeTrlAvg=mean(obj.ManData.ReshapeCell2Mat(temp,3),3);
                            % time and trial for no switch
                            temp=arrayfun(@(x) Vals.NoSwchBlks{x}(m).AvgTimTrl,ThisChs,'UniformOutput',0);
                            AvgHChs.TimeTrlAvgNoSwch=mean(obj.ManData.ReshapeCell2Mat(temp,3),3);
                            % average for time
                            AvgHChs.AvgTim=transpose(cell2mat(arrayfun(@(x) Vals.SwchBlks{x}(m).AvgTim',ThisChs,'UniformOutput',0)));
                            
                            %% with switch blocks
                            % average pre trial
                            AvgHChs.AvgTrlPre=transpose(cell2mat(arrayfun(@(x) Vals.SwchBlks{x}(m).AvgTrlPre',ThisChs,'UniformOutput',0)));
                            % average post trial
                            AvgHChs.AvgTrlPst=transpose(cell2mat(arrayfun(@(x) Vals.SwchBlks{x}(m).AvgTrlPst',ThisChs,'UniformOutput',0)));
                            
                            %% noswitch blocks
                            AvgHChs.AvgTrlPreNoSwch=transpose(cell2mat(arrayfun(@(x) Vals.NoSwchBlks{x}(m).AvgTrlPre',ThisChs,'UniformOutput',0)));
                            % average post trial
                            AvgHChs.AvgTrlPstNoSwch=transpose(cell2mat(arrayfun(@(x) Vals.NoSwchBlks{x}(m).AvgTrlPst',ThisChs,'UniformOutput',0)));
                            
                            %  AvgHChs.TimeAvg=cell2mat(arrayfun(@(x) mean(mean(HChDataAll.(ThisTimPrd){x}(m,:,:,:),3),4)',ThisChs,'UniformOutput',0))';
                            %  AvgHChs.TrlAvg=cell2mat(arrayfun(@(x) squeeze(mean(mean(HChDataAll.(ThisTimPrd){x}(m,WantedTime,:,:),2),4)),ThisChs,'UniformOutput',0))';
                            % avegrage rule selectivity
                            temp=arrayfun(@(x) Vals.SwchBlks{x}(m).RuleSelctivity,ThisChs,'UniformOutput',0);
                            AvgHChs.RuleSelectivity=mean(obj.ManData.ReshapeCell2Mat(temp,3),3);
                            % avegrage reward selectivity
                            temp=arrayfun(@(x) Vals.SwchBlks{x}(m).RewardSelctivity,ThisChs,'UniformOutput',0);
                            AvgHChs.RewardSelectivity=mean(obj.ManData.ReshapeCell2Mat(temp,3),3);
                            
                            % Plot Time average
                            subplot(NRow,NCol,1+(i-1)*NCol)
                            hold on
                            obj.FigParams.PlotMeanStd(opts.TimSeqAvg,AvgHChs.AvgTim,[],'Time to Stim',AnalysisFields{f},FieldCOL(f,:),1,'')
                            title(['Avg H over time' AreaToPlot{i}]);xlabel('Time');ylabel('Norm H')
                            % plot Trial Average pre
                            subplot(NRow,NCol,2+(i-1)*NCol)
                            hold on
                            obj.FigParams.PlotMeanStd(opts.TrlSeqAvg,AvgHChs.AvgTrlPre,[],'trl2swtch',AnalysisFields{f},FieldCOL(f,:),1,'')
                            title('Avg H Pre Stim On');xlabel('Trl 2 swtch');ylabel('Norm H')
                            % plot no switch
                            obj.FigParams.PlotMeanStd(opts.TrlSeqAvg,AvgHChs.AvgTrlPreNoSwch,[],'trl2swtch',AnalysisFields{f},1-FieldCOL(f,:),1,'')
                            title('Avg H Pre Stim On');xlabel('Trl 2 swtch');ylabel('Norm H')
                            % plot Trial Average post
                            subplot(NRow,NCol,3+(i-1)*NCol)
                            hold on
                            obj.FigParams.PlotMeanStd(opts.TrlSeqAvg,AvgHChs.AvgTrlPst,[],'trl2swtch',AnalysisFields{f},FieldCOL(f,:),1,'')
                            title('Avg H Post Stim On');xlabel('Trl 2 swtch');ylabel('Norm H')
                            % plot No switch
                            obj.FigParams.PlotMeanStd(opts.TrlSeqAvg,AvgHChs.AvgTrlPstNoSwch,[],'trl2swtch',AnalysisFields{f},1-FieldCOL(f,:),1,'')
                            title('Avg H Post Stim On');xlabel('Trl 2 swtch');ylabel('Norm H')
                            
                            % plot Time and Trial Average
                            subplot(NRow,NCol,4+(i-1)*NCol)
                            hold on
                            imagesc(opts.TimSeqAvg,opts.TrlSeqAvg, AvgHChs.TimeTrlAvg');
                            colorbar;title('Average H');xlabel('Time');ylabel('Trial to Switch')
                            
                            % plot Time and Trial Average no switch
                            subplot(NRow,NCol,5+(i-1)*NCol)
                            hold on
                            imagesc(opts.TimSeqAvg,opts.TrlSeqAvg, AvgHChs.TimeTrlAvgNoSwch');
                            colorbar;title('Average H');xlabel('Time');ylabel('Trial to Switch(no Switch Blks)')
                            
                            %                             % plot rule selectivity
                            %                             subplot(NRow,NCol,5+(i-1)*NCol)
                            %                             hold on
                            %                             imagesc(opts.TimSeqAvg,opts.TrlSeqAvg,AvgHChs.RuleSelectivity');
                            %                             colorbar;title('Rule Selectivity');xlabel('Time');ylabel('Trial to Switch')
                            %                             % plot reward selectivity
                            subplot(NRow,NCol,6+(i-1)*NCol)
                            hold on
                            imagesc(opts.TimSeqAvg,opts.TrlSeqAvg,AvgHChs.RewardSelectivity');
                            colorbar;title('Reward Selectivity');xlabel('Time');ylabel('Trial to Switch')
                            %  obj.FigParams.PlotMeanStd(1:75,AvgHChs.(AnalysisFields{f}),[],'trl2swtch',AnalysisFields{f},FieldCOL(f,:),1,'')
                            %  v=axis;
                            %  plot([-opts.NAvg -opts.NAvg]+1,[v(3) v(4)],'r')
                            %  title([PairsToPlot{i} TimePeriodName{TimPeriod} num2str(m)])
                        end
                    end
                end
            end
            AnalysisOpts.CurrentCh='';
            [~,~,AnalysisData.SelectivityFig]=obj.ManData.GetFileName([],'SelectivityTargOn3','SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.SaveFigSeries([],AnalysisData.SelectivityFig,[Figs])
            
            varargout=[Figs ];
            
            
            %                     figure(2)
            %                     subplot(5,6,sp)
            %                     hold on
            %                     [~,PCAMeanZeroLag]=pca(AvgPairsFixOn.MeanZeroLag');
            %                     hold on;arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',COL(x,:),'MarkerSize',10),1:46)
            %                     hold on;arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
            %
            %                      figure(3)
            %                     subplot(5,6,sp)
            %                     hold on
            %                     [~,PCAMeanZeroLag]=pca(AvgPairsFixOn.MaxXcorrLag');
            %                     hold on;arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',COL(x,:),'MarkerSize',10),1:46)
            %                     hold on;arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
            
            
            
            %                 % plot target on now
            %                 Fig_TargOn{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                 AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 2}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' TargOn ' num2str(sp)])
            %                 end
            %
            %                 % plot resp on
            %                 Fig_RespOn{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                  AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 3}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' RespOn ' num2str(sp)])
            %                 end
            %                 % plot reward
            %                 Fig_Reward{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                  AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 4}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' Reward ' num2str(sp)])
            
            
            
            
            %  obj.ManData.SaveVar([],AvgTimBoxTrls,'AvgTimBoxTrls','AvgTimBoxTrlsData');
            
            
        end
        function varargout=PlotHSCorrStatsAreaBhvModel(obj,ChsIndsAll,HTrialAll,varargin) % has BhvModel-correlation stats for Motifs
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % set-up options
            opts.TimePeriodName={'Reward','TargetOn'};
            opts.AreaToPlot= {'PFCChs','FEFChs','LIPChs'};
            opts.AnalysisFields={'value_for_choice_Info','RPEInfo','precision_ind_Info'};
            opts.AnalysisFieldsTitle={'Value','RPE','Precision'};
            Ts=1/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            Time=-1:Ts:(0.5-Ts); %% what was atually in the experiment
            
            FEFChs=ChsIndsAll.FEFChs;
            PFCChs=ChsIndsAll.PFCChs;
            LIPChs=ChsIndsAll.LIPChs;
            AllChs=ChsIndsAll.AllChs;
            NCoreMotifs=length(AnalysisData.CoreMotifs);
            %% prepare all of the data ( put the data for each channel into one cell)
            for f=opts.TimePeriodName % loop on time period
                HChDataAll.(f{1})  =obj.ManData.ReshapeCell2Mat(HTrialAll,52,f{1});
                for ChInd=1:length(AllChs) % loop on channel
                    fprintf('\nProcessing info for ch:%i time period:%s',ChInd,f{1});
                    ThisChRec=ChsIndsAll.RecInd(ChInd);
                    ThisRecSwchBlks=AnalysisData.RecSpecs(ThisChRec).BlockSpec.BhvModelInfo.model_outputs.Total_Switch;
                    %    SwchBlks=find(ThisRecSwchBlks);NoSwchBlks=find(ThisRecSwchBlks==0);
                    Blks=1:length(ThisRecSwchBlks);
                    % get the info for blocks where there is a switch
                    [CorrVals.(f{1}){ChInd}]=obj.CalHCorrBlockChBhvModel(HChDataAll.(f{1}){ChInd},...
                        AnalysisOpts.CohAnalysis.NAvgCoh,Time,AnalysisData.RecSpecs(ThisChRec).TrialTimes,...
                        AnalysisData.RecSpecs(ThisChRec).BlockSpec.BhvModelInfo.BlockSpec,Blks);
                end
                %                 % reorganize into full matrix
                CorrVals_ro.(f{1})=obj.ManData.ReshapeCellStruct2Mat(CorrVals.(f{1}));
            end
            
            FieldCOL=distinguishable_colors(length(opts.AnalysisFields));
            NRow=5;NCol=NCoreMotifs/5;
            FigNum=1;% plot Xcorr values first
            for field=opts.AnalysisFields
                f=find(strcmp(field{1},opts.AnalysisFields));
                for TimPeriod=1:length(opts.TimePeriodName)
                    ThisTimPrd=opts.TimePeriodName{TimPeriod};
                    Figs(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    
                    for i=1:length(opts.AreaToPlot) % loop on areas to plot
                        eval(['ThisChs=' opts.AreaToPlot{i} ';']);
                        for m=1:NCoreMotifs       % loop on core motifs
                            % plot information with null distribution
                            subplot(NRow,NCol,m)
                            FieldTitle=opts.AnalysisFieldsTitle{f};
                            if contains(FieldTitle,'Value')
                                Info=squeeze(mean(CorrVals_ro.(ThisTimPrd).(field{1})(4,m,1,:,ThisChs),1))';
                                InfoNull=squeeze(mean(CorrVals_ro.(ThisTimPrd).(field{1})(4,m,2:end,:,ThisChs),3))';
                                
                            else
                                Info=squeeze(CorrVals_ro.(ThisTimPrd).(field{1})(m,1,:,ThisChs))';
                                InfoNull=squeeze(mean(CorrVals_ro.(ThisTimPrd).(field{1})(m,2:end,:,ThisChs),2))';
                            end
                            hold on
                            obj.FigParams.PlotMeanStd(Time,movmean(Info,5,2),[],['Time to ' ThisTimPrd],'I(bits)',i,3,[FieldTitle ',' num2str(m)])
                            % plot null information now
                            obj.FigParams.PlotMeanStd(Time,InfoNull,[],['Time to ' ThisTimPrd],'I(bits)',i,3,[FieldTitle ',' num2str(m)],'--')
                            
                            %    correlation
                            %                              subplot(NRow,NCol,m)
                            %                             Info=squeeze(CorrVals_ro.(ThisTimPrd).(field{1})(m,:,ThisChs))';
                            %                              hold on
                            %                             obj.FigParams.PlotMeanStd(Time,Info,[],['Time to ' ThisTimPrd],'corr',i,1,[field{1} ',' num2str(m)])
                        end
                    end
                    legend(opts.AreaToPlot,'Location','best','FontSize',5);
                    
                end
            end
            
            AnalysisOpts.CurrentCh='';
            [~,~,AnalysisData.SelectivityFig]=obj.ManData.GetFileName([],'InformationModelWNull','SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.SaveFigSeries([],AnalysisData.SelectivityFig,[Figs])
            varargout=[Figs ];
            
            
            %                     figure(2)
            %                     subplot(5,6,sp)
            %                     hold on
            %                     [~,PCAMeanZeroLag]=pca(AvgPairsFixOn.MeanZeroLag');
            %                     hold on;arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',COL(x,:),'MarkerSize',10),1:46)
            %                     hold on;arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
            %
            %                      figure(3)
            %                     subplot(5,6,sp)
            %                     hold on
            %                     [~,PCAMeanZeroLag]=pca(AvgPairsFixOn.MaxXcorrLag');
            %                     hold on;arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',COL(x,:),'MarkerSize',10),1:46)
            %                     hold on;arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
            
            
            
            %                 % plot target on now
            %                 Fig_TargOn{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                 AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 2}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' TargOn ' num2str(sp)])
            %                 end
            %
            %                 % plot resp on
            %                 Fig_RespOn{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                  AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 3}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' RespOn ' num2str(sp)])
            %                 end
            %                 % plot reward
            %                 Fig_Reward{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                  AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 4}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' Reward ' num2str(sp)])
            
            
            
            
            %  obj.ManData.SaveVar([],AvgTimBoxTrls,'AvgTimBoxTrls','AvgTimBoxTrlsData');
            
            
        end        
        function varargout=PlotHSelectivityStatsArea(obj,ChsIndsAll,HTrialAll,varargin) % plots coherence stats between areas
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            FEFChs=ChsIndsAll.FEFChs;
            PFCChs=ChsIndsAll.PFCChs;
            LIPChs=ChsIndsAll.LIPChs;
            AllChs=ChsIndsAll.AllChs;
            NCoreMotifs=length(AnalysisData.CoreMotifs);
            %% prepare all of the data ( put the data for each channel into one cell)
            HChDataAll.Reward=obj.ManData.ReshapeCell2Mat(HTrialAll,52,'TargetOn');
            %% precalculate all of the selectivities we care about
            Ts=1/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            Time=-1:Ts:(0.5-Ts); %% what was atually in the experiment
            Vals=[];
            for ChInd=1:length(AllChs)
                ThisChRec=ChsIndsAll.RecInd(ChInd);
                [Vals{ChInd},Params]=obj.CalHSelectivityBlockChBhvModel(HChDataAll.Reward{ChInd},...
                    AnalysisOpts.CohAnalysis.NAvgCoh,Time,AnalysisData.RecSpecs(ThisChRec).TrialTimes,...
                    AnalysisData.RecSpecs(ThisChRec).BlockSpec);
            end
            %              %% grab Color data for all of trails and channels
            %             % plot color data now
            %             ColorFigs=obj.PlotW_H_EachW_Condition(W,H_tot_TargOn,AnalysisData.cwt_f,AnalysisData.Time,TrialTimes,AnalysisData.ChArea,'BestColor',ColList);
            %             % grab H values based on 35 before and after switch
            %             [H_exmp_tot,BlkTrialNumber,ExpBlkNum,Bhv_blk_Perf]=obj.GrabHTimeTrial(H,BlockSpec);
            %             obj.PlotW_H_Trial_Condition(W,H_exmp_tot.Reward,AnalysisData.cwt_f,AnalysisData.Time,BlockSpec,AnalysisData.ChArea,[],[]);
            %             Wfigs=cell(1,length(AnalysisData.MeanCluster));
            %             %% area based Analysis
            %             % plot W-H Area
            %             [Wfigs{:}]=obj.PlotW_H_Trial_Area(W,H_exmp_tot.FixOn(:,WBlks),AnalysisData.cwt_f,Time,BlkTrialNumber(WBlks),BlockSpec.TrialOrder,ExpBlkNum(WBlks),AnalysisData.ChArea,'ChNum',NaN) ; % plot Ws
            %% calculate some of the vars we need
            opts.Ntrls2Switch2Look=10; % How many trials around the switch we want to look at
            trlSeq=AnalysisData.BlockSpec.TrialOrder;
            opts.SwitchTrl=find(trlSeq==0);
            opts.SwitchTrlAvg=find(Params.TrlAvgInd(:,2)==opts.SwitchTrl);
            opts.NAvg=opts.SwitchTrl-opts.SwitchTrlAvg+1;
            opts.trlSeqAvg=trlSeq(1:end-opts.NAvg+1);
            % get the new time seq average
            opts.TrlSeqAvg=trlSeq(Params.TrlAvgInd(:,2));
            opts.TimSeqAvg=Time(Params.TimAvgInd(:,2));
            NCol=6;
            WantedTime=Time>=0 & Time<=0.2;
            %  opts.trl2lookInd=(SwitchTrl-opts.Ntrls2Switch2Look):(SwitchTrl+opts.Ntrls2Switch2Look-1-4); % trials to look at
            % opts.trlSeq2look=trlSeq;%-opts.Ntrls2Switch2Look:opts.Ntrls2Switch2Look;
            TimePeriodName={'Reward'};
            AreaToPlot={'PFCChs','FEFChs','LIPChs','AllChs'};
            AnalysisFields={'TimeAvg','TrlAvg'};
            FieldCOL=distinguishable_colors(length(AnalysisFields));
            NRow=length(AreaToPlot);
            FigNum=1;% plot Xcorr values first
            for TimPeriod=1:length(TimePeriodName)
                ThisTimPrd=TimePeriodName{TimPeriod};
                for m=1:NCoreMotifs
                    Figs(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    obj.PlotCoreMotifs([],m); % plot this motif
                    Figs(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    for i=1:length(AreaToPlot)
                        eval(['ThisChs=' AreaToPlot{i} ';']);
                        % plot fix on
                        for f=1%:length(AnalysisFields)
                            
                            temp=arrayfun(@(x) Vals{x}(m).AvgTimTrl,ThisChs,'UniformOutput',0);
                            AvgHChs.TimeTrlAvg=mean(obj.ManData.ReshapeCell2Mat(temp,3),3);
                            % average for time
                            AvgHChs.AvgTim=transpose(cell2mat(arrayfun(@(x) Vals{x}(m).AvgTim',ThisChs,'UniformOutput',0)));
                            
                            % average pre trial
                            AvgHChs.AvgTrlPre=transpose(cell2mat(arrayfun(@(x) Vals{x}(m).AvgTrlPre',ThisChs,'UniformOutput',0)));
                            % average post trial
                            AvgHChs.AvgTrlPst=transpose(cell2mat(arrayfun(@(x) Vals{x}(m).AvgTrlPst',ThisChs,'UniformOutput',0)));
                            
                            %  AvgHChs.TimeAvg=cell2mat(arrayfun(@(x) mean(mean(HChDataAll.(ThisTimPrd){x}(m,:,:,:),3),4)',ThisChs,'UniformOutput',0))';
                            %  AvgHChs.TrlAvg=cell2mat(arrayfun(@(x) squeeze(mean(mean(HChDataAll.(ThisTimPrd){x}(m,WantedTime,:,:),2),4)),ThisChs,'UniformOutput',0))';
                            % avegrage rule selectivity
                            temp=arrayfun(@(x) Vals{x}(m).RuleSelctivity,ThisChs,'UniformOutput',0);
                            AvgHChs.RuleSelectivity=mean(obj.ManData.ReshapeCell2Mat(temp,3),3);
                            % avegrage reward selectivity
                            temp=arrayfun(@(x) Vals{x}(m).RewardSelctivity,ThisChs,'UniformOutput',0);
                            AvgHChs.RewardSelectivity=mean(obj.ManData.ReshapeCell2Mat(temp,3),3);
                            
                            % Plot Time average
                            subplot(NRow,NCol,1+(i-1)*NCol)
                            hold on
                            obj.FigParams.PlotMeanStd(opts.TimSeqAvg,AvgHChs.AvgTim,[],'Time to Stim',AnalysisFields{f},FieldCOL(f,:),1,'')
                            title(['Avg H over time' AreaToPlot{i}]);xlabel('Time');ylabel('Norm H')
                            % plot Trial Average pre
                            subplot(NRow,NCol,2+(i-1)*NCol)
                            hold on
                            obj.FigParams.PlotMeanStd(opts.TrlSeqAvg,AvgHChs.AvgTrlPre,[],'trl2swtch',AnalysisFields{f},FieldCOL(f,:),1,'')
                            title('Avg H Pre Stim On');xlabel('Trl 2 swtch');ylabel('Norm H')
                            % plot Trial Average post
                            subplot(NRow,NCol,3+(i-1)*NCol)
                            hold on
                            obj.FigParams.PlotMeanStd(opts.TrlSeqAvg,AvgHChs.AvgTrlPst,[],'trl2swtch',AnalysisFields{f},FieldCOL(f,:),1,'')
                            title('Avg H Post Stim On');xlabel('Trl 2 swtch');ylabel('Norm H')
                            
                            
                            % plot Time and Trial Average
                            subplot(NRow,NCol,4+(i-1)*NCol)
                            hold on
                            imagesc(opts.TimSeqAvg,opts.TrlSeqAvg, AvgHChs.TimeTrlAvg');
                            colorbar;title('Average H');xlabel('Time');ylabel('Trial to Switch')
                            % plot rule selectivity
                            subplot(NRow,NCol,5+(i-1)*NCol)
                            hold on
                            imagesc(opts.TimSeqAvg,opts.TrlSeqAvg,AvgHChs.RuleSelectivity');
                            colorbar;title('Rule Selectivity');xlabel('Time');ylabel('Trial to Switch')
                            % plot reward selectivity
                            subplot(NRow,NCol,6+(i-1)*NCol)
                            hold on
                            imagesc(opts.TimSeqAvg,opts.TrlSeqAvg,AvgHChs.RewardSelectivity');
                            colorbar;title('Reward Selectivity');xlabel('Time');ylabel('Trial to Switch')
                            %  obj.FigParams.PlotMeanStd(1:75,AvgHChs.(AnalysisFields{f}),[],'trl2swtch',AnalysisFields{f},FieldCOL(f,:),1,'')
                            %  v=axis;
                            %  plot([-opts.NAvg -opts.NAvg]+1,[v(3) v(4)],'r')
                            %  title([PairsToPlot{i} TimePeriodName{TimPeriod} num2str(m)])
                        end
                    end
                end
            end
            AnalysisOpts.CurrentCh='';
            [~,~,AnalysisData.SelectivityFig]=obj.ManData.GetFileName([],'SelectivityTargOn3','SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.SaveFigSeries([],AnalysisData.SelectivityFig,[Figs])
            
            varargout=[Figs ];
            
            
            %                     figure(2)
            %                     subplot(5,6,sp)
            %                     hold on
            %                     [~,PCAMeanZeroLag]=pca(AvgPairsFixOn.MeanZeroLag');
            %                     hold on;arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',COL(x,:),'MarkerSize',10),1:46)
            %                     hold on;arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
            %
            %                      figure(3)
            %                     subplot(5,6,sp)
            %                     hold on
            %                     [~,PCAMeanZeroLag]=pca(AvgPairsFixOn.MaxXcorrLag');
            %                     hold on;arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',COL(x,:),'MarkerSize',10),1:46)
            %                     hold on;arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
            
            
            
            %                 % plot target on now
            %                 Fig_TargOn{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                 AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 2}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' TargOn ' num2str(sp)])
            %                 end
            %
            %                 % plot resp on
            %                 Fig_RespOn{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                  AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 3}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' RespOn ' num2str(sp)])
            %                 end
            %                 % plot reward
            %                 Fig_Reward{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                  AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 4}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' Reward ' num2str(sp)])
            
            
            
            
            %  obj.ManData.SaveVar([],AvgTimBoxTrls,'AvgTimBoxTrls','AvgTimBoxTrlsData');
            
            
        end
        function [Vals,Params]=CalHSelectivityBlockChBhvModel(obj,HChData,NTrlAvg,Time,TrialTimes,BlockSpec,WantedBlks,varargin)% includes Bhv Model
            % includes the trials across all of the blocks to calculate
            % selectvity
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            opts.NormalizeH=0; % devides H for each motif to it's max
            opts.TimPridLength=0.2; % lenght of time period we are looking after or before an event
            GrabTrialData=@(x,y)  [nan*ones(length(find(x<1)),size(y,2)) ;y(x(x>=1),:)]; % function to get data from trial times
            % NTrlAvg number of trials to be averaged for coherence stats
            HChData=HChData(:,:,:,WantedBlks);
            NMotif=size(HChData,1);
            Ntim=size(HChData,2); % time points
            NTrl=size(HChData,3);
            NBlks=size(HChData,4); % how many blocks we want to average% last blcok is the average of all (should be fixed)
            %% attention removing last block
            HChData=HChData(:,:,:,1:NBlks);
            %%
            NTrlBlk=NTrlAvg*NBlks; % number of  Trials * Num blocks
            NTimAvg=10;% is 10*1/50=200ms time averaging
            % generate the inds for averaging acorss trials
            TrlAvgInd=obj.ManData.GenMovAvgInds(NTrl,NTrlAvg,'valid');
            Params.TrlAxis=TrlAvgInd(:,2)';
            Params.TrlAvgInd=TrlAvgInd;
            % generate the inds for averaging acorss time
            TimAvgInd=obj.ManData.GenMovAvgInds(Ntim,NTimAvg,'valid');
            Params.TimeAxis=Time(TimAvgInd(:,2));
            Params.TimAvgInd=TimAvgInd;
            % generate Inds for averaging across specific period
            FloorTime=int8(Time*100); %********attention double check this ***
            FloorPreEventInds=int8(opts.TimPridLength*100);
            Params.PreEventInds=FloorTime>=-FloorPreEventInds & FloorTime<0;
            Params.PstEventInds=FloorTime>=0 & FloorTime<FloorPreEventInds;
            %
            for m=1:NMotif
                ThisMotif_HChData=squeeze(HChData(m,:,:,:)); % take all of the H channel data only for this motif
                if opts.NormalizeH % if we are normalizing H values
                    ThisMotif_HChData=ThisMotif_HChData/max(ThisMotif_HChData(:));
                    % z-score values across trials
                    %                     ThisMotif_HChDataReshape=reshape(ThisMotif_HChData,...
                    %                         [size(ThisMotif_HChData,1) size(ThisMotif_HChData,2)*size(ThisMotif_HChData,3)]);
                    %                     z_ThisMotif_HChDataReshape=zscore(ThisMotif_HChDataReshape,0,1);
                    %                     ThisMotif_HChData=reshape(z_ThisMotif_HChDataReshape,[size(ThisMotif_HChData,1) size(ThisMotif_HChData,2) size(ThisMotif_HChData,3)]);
                end
                ThisMotif_HChDataAll=reshape(ThisMotif_HChData,[size(ThisMotif_HChData,1) size(ThisMotif_HChData,2)*size(ThisMotif_HChData,3)]);
                
                for Tims=1:size(TimAvgInd)
                    ThisTims=TimAvgInd(Tims,:);
                    % calculte average overtime
                    temp=ThisMotif_HChDataAll(ThisTims(1):ThisTims(2),:);
                    Vals(m).AvgTim(Tims)=nanmean(nanmean(temp,1),2);
                    
                    for Trls=1:size(TrlAvgInd,1)
                        ThisTrls=TrlAvgInd(Trls,:);
                        ThisMotifH=[];BlkTrls=[]; ThisMotifH_Pre=[];ThisMotifH_Pst=[];
                        for blk=1:NBlks
                            ThisMotifH=[ThisMotifH,ThisMotif_HChData(ThisTims(1):ThisTims(2),ThisTrls(1):ThisTrls(2),blk)];
                            % get the information about stimulus and reward
                            % for each trial
                            BlkTrls=[BlkTrls,BlockSpec.ThisBlkTrials{WantedBlks(blk)}(ThisTrls(1):ThisTrls(2))];
                            % get time specific averages
                            if Tims==1
                                ThisMotifH_Pre=[ThisMotifH_Pre,ThisMotif_HChData(Params.PreEventInds,ThisTrls(1):ThisTrls(2),blk)];
                                ThisMotifH_Pst=[ThisMotifH_Pst,ThisMotif_HChData(Params.PstEventInds,ThisTrls(1):ThisTrls(2),blk)];
                            end
                        end
                        % get the number of Rewards
                        RwdInd=strcmp(AnalysisOpts.TrialTimesFields,'NumRewards');
                        NumRewards=GrabTrialData(BlkTrls,TrialTimes(:,RwdInd))';
                        % get the current rule of the trial
                        BstClrInd=strcmp(AnalysisOpts.TrialTimesFields,'BestColor');
                        CurrClrRule=GrabTrialData(BlkTrls,TrialTimes(:,BstClrInd))';
                        %% calculate information about Current color Rule
                        %  CurrClrRule_bin=obj.ManData.BinColorSpace(4,CurrClrRule(Trls,:));
                        MeanHThisTrl=nanmean(ThisMotifH,1);
                        %  I=obj.ManData.CalculateInformation(CurrClrRule_bin,MeanHThisTrl);
                        % calculate Selectivity
                        Vals(m).RuleSelctivity(Tims,Trls)=obj.ManData.CalSelectivity(CurrClrRule,MeanHThisTrl);
                        Vals(m).RewardSelctivity(Tims,Trls)=obj.ManData.CalSelectivity(NumRewards,MeanHThisTrl);
                        
                        %% take average across Time and Trial
                        Vals(m).AvgTimTrl(Tims,Trls)=nanmean(nanmean(ThisMotifH,1),2);
                        if Tims==1
                            Vals(m).AvgTrlPre(Trls)= nanmean(nanmean(ThisMotifH_Pre,1),2);
                            Vals(m).AvgTrlPst(Trls)= nanmean(nanmean(ThisMotifH_Pst,1),2);
                        end
                    end
                end
            end
            
        end
        function [Vals,Params]=CalHSelectivityBlockCh(obj,HChData,NTrlAvg,Time,TrialTimes,BlockSpec,varargin)
            % includes the trials across all of the blocks to calculate
            % selectvity
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            opts.NormalizeH=1; % devides H for each motif to it's max
            opts.TimPridLength=0.2; % lenght of time period we are looking after or before an event
            % NTrlAvg number of trials to be averaged for coherence stats
            NMotif=size(HChData,1);
            Ntim=size(HChData,2); % time points
            NTrl=size(HChData,3);
            NBlks=size(HChData,4)-1; % how many blocks we want to average% last blcok is the average of all (should be fixed)
            %% attention removing last block
            HChData=HChData(:,:,:,1:NBlks);
            %%
            NTrlBlk=NTrlAvg*NBlks; % number of  Trials * Num blocks
            NTimAvg=10;% is 10*1/50=200ms time averaging
            % generate the inds for averaging acorss trials
            TrlAvgInd=obj.ManData.GenMovAvgInds(NTrl,NTrlAvg,'valid');
            Params.TrlAxis=TrlAvgInd(:,2)';
            Params.TrlAvgInd=TrlAvgInd;
            % generate the inds for averaging acorss time
            TimAvgInd=obj.ManData.GenMovAvgInds(Ntim,NTimAvg,'valid');
            Params.TimeAxis=Time(TimAvgInd(:,2));
            Params.TimAvgInd=TimAvgInd;
            % generate Inds for averaging across specific period
            FloorTime=int8(Time*100); %********attention double check this ***
            FloorPreEventInds=int8(opts.TimPridLength*100);
            Params.PreEventInds=FloorTime>=-FloorPreEventInds & FloorTime<0;
            Params.PstEventInds=FloorTime>=0 & FloorTime<FloorPreEventInds;
            %
            for m=1:NMotif
                ThisMotif_HChData=squeeze(HChData(m,:,:,:)); % take all of the H channel data only for this motif
                if opts.NormalizeH % if we are normalizing H values
                    ThisMotif_HChData=ThisMotif_HChData/max(ThisMotif_HChData(:));
                    % z-score values across trials
                    %                     ThisMotif_HChDataReshape=reshape(ThisMotif_HChData,...
                    %                         [size(ThisMotif_HChData,1) size(ThisMotif_HChData,2)*size(ThisMotif_HChData,3)]);
                    %                     z_ThisMotif_HChDataReshape=zscore(ThisMotif_HChDataReshape,0,1);
                    %                     ThisMotif_HChData=reshape(z_ThisMotif_HChDataReshape,[size(ThisMotif_HChData,1) size(ThisMotif_HChData,2) size(ThisMotif_HChData,3)]);
                end
                ThisMotif_HChDataAll=reshape(ThisMotif_HChData,[size(ThisMotif_HChData,1) size(ThisMotif_HChData,2)*size(ThisMotif_HChData,3)]);
                
                for Tims=1:size(TimAvgInd)
                    ThisTims=TimAvgInd(Tims,:);
                    % calculte average overtime
                    temp=ThisMotif_HChDataAll(ThisTims(1):ThisTims(2),:);
                    Vals(m).AvgTim(Tims)=nanmean(nanmean(temp,1),2);
                    
                    for Trls=1:size(TrlAvgInd,1)
                        ThisTrls=TrlAvgInd(Trls,:);
                        ThisMotifH=[];BlkTrls=[]; ThisMotifH_Pre=[];ThisMotifH_Pst=[];
                        for blk=1:NBlks
                            ThisMotifH=[ThisMotifH,ThisMotif_HChData(ThisTims(1):ThisTims(2),ThisTrls(1):ThisTrls(2),blk)];
                            % get the information about stimulus and reward
                            % for each trial
                            BlkTrls=[BlkTrls,BlockSpec.ThisBlkTrials{blk}(ThisTrls(1):ThisTrls(2))];
                            % get time specific averages
                            if Tims==1
                                ThisMotifH_Pre=[ThisMotifH_Pre,ThisMotif_HChData(Params.PreEventInds,ThisTrls(1):ThisTrls(2),blk)];
                                ThisMotifH_Pst=[ThisMotifH_Pst,ThisMotif_HChData(Params.PstEventInds,ThisTrls(1):ThisTrls(2),blk)];
                            end
                        end
                        % get the number of Rewards
                        RwdInd=strcmp(AnalysisOpts.TrialTimesFields,'NumRewards');
                        NumRewards=TrialTimes(BlkTrls,RwdInd)';
                        % get the current rule of the trial
                        BstClrInd=strcmp(AnalysisOpts.TrialTimesFields,'BestColor');
                        CurrClrRule=TrialTimes(BlkTrls,BstClrInd)';
                        %% calculate information about Current color Rule
                        %  CurrClrRule_bin=obj.ManData.BinColorSpace(4,CurrClrRule(Trls,:));
                        MeanHThisTrl=nanmean(ThisMotifH,1);
                        %  I=obj.ManData.CalculateInformation(CurrClrRule_bin,MeanHThisTrl);
                        % calculate Selectivity
                        Vals(m).RuleSelctivity(Tims,Trls)=obj.ManData.CalSelectivity(CurrClrRule,MeanHThisTrl);
                        Vals(m).RewardSelctivity(Tims,Trls)=obj.ManData.CalSelectivity(NumRewards,MeanHThisTrl);
                        
                        %% take average across Time and Trial
                        Vals(m).AvgTimTrl(Tims,Trls)=nanmean(nanmean(ThisMotifH,1),2);
                        if Tims==1
                            Vals(m).AvgTrlPre(Trls)= nanmean(nanmean(ThisMotifH_Pre,1),2);
                            Vals(m).AvgTrlPst(Trls)= nanmean(nanmean(ThisMotifH_Pst,1),2);
                        end
                    end
                end
            end
            
        end
        function [Vals]=CalHCorrBlockChBhvModel(obj,HChData,NTrlAvg,Time,TrialTimes,BlockSpec,WantedBlks,varargin)% includes Bhv Model
            % calculates correlation of H values for each motif to
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            opts.NormalizeH=0; % devides H for each motif to it's max
            opts.TimPridLength=0.2; % lenght of time period we are looking after or before an event
            GrabTrialData=@(x,y)  [nan*ones(length(find(x<1)),size(y,2)) ;y(x(x>=1),:)]; % function to get data from trial times
            % NTrlAvg number of trials to be averaged for coherence stats
            HChData=HChData(:,:,:,WantedBlks);
            NMotif=size(HChData,1);
            Ntim=size(HChData,2); % time points
            NTrl=size(HChData,3);
            NBlks=size(HChData,4); % how many blocks we want to average% last blcok is the average of all (should be fixed)
            %% attention removing last block
            HChData=HChData(:,:,:,1:NBlks);
            Vals=[];
            for m=1:NMotif
                ThisMotif_HChData=squeeze(HChData(m,:,:,:)); % take all of the H channel data only for this motif
                if opts.NormalizeH % if we are normalizing H values
                    ThisMotif_HChData=ThisMotif_HChData/max(ThisMotif_HChData(:));
                end
                ThisMotif_HChDataAll=reshape(ThisMotif_HChData,[size(ThisMotif_HChData,1) size(ThisMotif_HChData,2)*size(ThisMotif_HChData,3)]);
                % loop in time points and trials calculate correlation to
                % RPE, confidence and value_for Choice
                RPE=cell2mat(BlockSpec.RPE(WantedBlks)');
                precision_ind=cell2mat(BlockSpec.precision_ind(WantedBlks)');
                value_for_choice=cell2mat(BlockSpec.down_sampled_Value_for_choice4(WantedBlks))';
                % cal info for RPE
                [Vals.RPECorr(m,:),Vals.RPECorr_p(m,:)]=corr(ThisMotif_HChDataAll',RPE,'rows','pairwise');
                Vals.RPEInfo(m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',RPE)',1:Ntim,'UniformOutput',0));
                % cal info for precision_ind
                [Vals.precision_ind(m,:),Vals.precision_ind_p(m,:)]=corr(ThisMotif_HChDataAll',precision_ind,'rows','pairwise');
                Vals.precision_ind_Info(m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',precision_ind)',1:Ntim,'UniformOutput',0));
                % cal info for belief
                for k=1:4
                    Vals.value_for_choice(k,m,:)=corr(ThisMotif_HChDataAll',value_for_choice(:,k),'rows','pairwise');
                    Vals.value_for_choice_Info(k,m,:,:)=cell2mat(arrayfun(@(x) obj.ManData.CalculateInformation(ThisMotif_HChDataAll(x,:)',value_for_choice(:,k))',1:Ntim,'UniformOutput',0));
                end
            end
            
        end

        function [Vals,Params]=CalHCoherenceBlockChPair(obj,XcorrPair,NTrlAvg,TimPeriod,varargin)
            % includes the trials across all of the blocks to calculate
            % coherence
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % Ch1 and Ch2 are raw data values  organized in 1*Time*Trial
            % NTrlAvg number of trials to be averaged for coherence stats
            RefXcorr=XcorrPair{1, 1}; %take one ref data point
            NMotif=length(RefXcorr);
            Ntim=size(RefXcorr{1},1); % time points of Xcorr
            NTrl=size(RefXcorr{1},2);
            NBlks=size(XcorrPair,2); % how many blocks we want to average
            NTrlBlk=NTrlAvg*NBlks; % number of  Trials * Num blocks
            % generate the inds for averaging acorss trials
            AvgInd=obj.ManData.GenMovAvgInds(NTrl,NTrlAvg,'valid');
            ZeroLag=AnalysisData.TXcorr==0;
            for m=1:NMotif
                for Trls=1:size(AvgInd,1)
                    ThisTrls=AvgInd(Trls,:);ThisMotifXcorr=[];
                    for blk=1:NBlks
                        ThisMotifXcorr=[ThisMotifXcorr,XcorrPair{TimPeriod,blk}{m}(:,ThisTrls(1):ThisTrls(2))];
                    end
                    temp=mean(ThisMotifXcorr,2);
                    Vals(m).MeanZeroLag(Trls)=temp(ZeroLag);
                    [Vals(m).MaxXcorr(Trls,:),Vals(m).MaxXcorrLag(Trls,:)]=max(ThisMotifXcorr,[],1);
                    % calculate Zscore of Xcorr for each trial
                    MeanXcorrTrl=mean(ThisMotifXcorr,1);StdXcorrTrl=std(ThisMotifXcorr,1,1);
                    ZeroLagXcorrTrl=ThisMotifXcorr(ZeroLag,:);MaxXcorrTrl=max(ThisMotifXcorr,[],1);
                    Vals(m).ZscoreZeroLag(Trls)=abs((ZeroLagXcorrTrl-MeanXcorrTrl)/StdXcorrTrl);
                    Vals(m).ZscoreMaxXcorr(Trls)=abs((MaxXcorrTrl-MeanXcorrTrl)/StdXcorrTrl);
                end
            end
            Params.AvgInd=AvgInd;
            Params.ZeroLag=ZeroLag;
        end
        function varargout=PlotHCoherenceStatsArea(obj,PairIndsAll,XcorrPairAll,varargin) % plots coherence stats between areas
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            FrontParietPairs=PairIndsAll.FrontParietPairs;
            FEFPairs=PairIndsAll.FEFPairs;
            PFCPairs=PairIndsAll.PFCPairs;
            LIPPairs=PairIndsAll.LIPPairs;
            AllPairs=PairIndsAll.AllPairs;
            NCoreMotifs=length(AnalysisData.CoreMotifs);
            % load all of the pairs first (this may take a while)
            opts.UseZscore=1;
            Vals=[];
            for TimPeriod=1:2 % which time period we are looking at(Fix on etc etc)
                for PairInd=1:length(XcorrPairAll)
                    %  tempBlk(:,:,Blk)= XcorrPairAll{1,PairInd}{TimPid,Blk}{1,Motif};
                    [Vals{TimPeriod,PairInd},Params]=obj.CalHCoherenceBlockChPair(XcorrPairAll{1,PairInd},AnalysisOpts.CohAnalysis.NAvgCoh,TimPeriod);
                end
            end
            obj.ManData.SaveVar([],Vals,'Vals','XcorrValues','WantedDate','ALL');
            obj.ManData.SaveVar([],Params,'Params','XcorrValues','WantedDate','ALL');
            obj.ManData.SaveVar([],PairIndsAll,'PairIndsAll','XcorrValues','WantedDate','ALL');
            %% prepare ms coherence
            %% calculate some of the vars we need
            opts.Ntrls2Switch2Look=10; % How many trials around the switch we want to look at
            trlSeq=AnalysisData.BlockSpec.TrialOrder;
            opts.SwitchTrl=find(trlSeq==0);
            opts.SwitchTrlAvg=find(Params.AvgInd(:,2)==opts.SwitchTrl);
            opts.NAvg=opts.SwitchTrl-opts.SwitchTrlAvg+1;
            opts.trlSeqAvg=trlSeq(1:end-opts.NAvg+1);
            %  opts.trl2lookInd=(SwitchTrl-opts.Ntrls2Switch2Look):(SwitchTrl+opts.Ntrls2Switch2Look-1-4); % trials to look at
            % opts.trlSeq2look=trlSeq;%-opts.Ntrls2Switch2Look:opts.Ntrls2Switch2Look;
            TimePeriodName={'TargOn','Reward'};
            PairsToPlot={'FrontParietPairs','FEFPairs','PFCPairs','LIPPairs','AllPairs'};
            XcorrFields={'MeanZeroLag','MaxXcorrLag','MaxXcorr','ZscoreZeroLag','ZscoreMaxXcorr'};
            FieldCOL=distinguishable_colors(length(XcorrFields));
            NCol=length(XcorrFields);NRow=length(PairsToPlot);
            FigNum=1;% plot Xcorr values first
            for TimPeriod=1:length(TimePeriodName)
                for m=1:NCoreMotifs
                    Figs(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    obj.PlotCoreMotifs([],m); % plot this motif
                    Figs(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    for i=1:length(PairsToPlot)
                        eval(['ThisPairs=' PairsToPlot{i} ';']);
                        % plot fix on
                        for f=1:length(XcorrFields)
                            AvgXcorrPairs.MeanZeroLag=cell2mat(arrayfun(@(x) Vals{TimPeriod,x }(m).MeanZeroLag',ThisPairs,'UniformOutput',0))';
                            AvgXcorrPairs.MaxXcorrLag=cell2mat(arrayfun(@(x) mean(AnalysisData.TXcorr(Vals{TimPeriod,x }(m).MaxXcorrLag),2),ThisPairs,'UniformOutput',0))'*1000/AnalysisOpts.MotifAnalysis.FsWaveTarg;
                            AvgXcorrPairs.MaxXcorr=cell2mat(arrayfun(@(x) mean(Vals{TimPeriod,x }(m).MaxXcorr,2),ThisPairs,'UniformOutput',0))';
                            AvgXcorrPairs.ZscoreZeroLag=cell2mat(arrayfun(@(x) Vals{TimPeriod,x }(m).ZscoreZeroLag',ThisPairs,'UniformOutput',0))';
                            AvgXcorrPairs.ZscoreMaxXcorr=cell2mat(arrayfun(@(x) Vals{TimPeriod,x }(m).ZscoreMaxXcorr',ThisPairs,'UniformOutput',0))';
                            
                            subplot(NRow,NCol,f+(i-1)*NCol)
                            hold on
                            obj.FigParams.PlotMeanStd(opts.trlSeqAvg,AvgXcorrPairs.(XcorrFields{f}),[],'trl2swtch',XcorrFields{f},FieldCOL(f,:),0,'')
                            v=axis;
                            plot([-opts.NAvg -opts.NAvg]+1,[v(3) v(4)],'r')
                            title([PairsToPlot{i} TimePeriodName{TimPeriod} num2str(m)])
                        end
                    end
                end
            end
            AnalysisOpts.CurrentCh='';
            [~,~,AnalysisData.XcorrValsFig]=obj.ManData.GetFileName([],'XcorrVals','SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.SaveFigSeries([],AnalysisData.XcorrValsFig,[Figs])
            close all
            %% plot Projected values into PC space
            FigNum=1;% plot Xcorr values first
            TrlCOL=colormap(jet(length(opts.trlSeqAvg)));
            for TimPeriod=1:length(TimePeriodName)
                for m=1:NCoreMotifs
                    FigsPC(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    obj.PlotCoreMotifs([],m); % plot this motif
                    FigsPC(FigNum)=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                    FigNum=FigNum+1;
                    for i=1:length(PairsToPlot)
                        eval(['ThisPairs=' PairsToPlot{i} ';']);
                        % plot fix on
                        for f=1:length(XcorrFields)
                            AvgXcorrPairs.MeanZeroLag=cell2mat(arrayfun(@(x) Vals{TimPeriod,x }(m).MeanZeroLag',ThisPairs,'UniformOutput',0))';
                            AvgXcorrPairs.MaxXcorrLag=cell2mat(arrayfun(@(x) mean(AnalysisData.TXcorr(Vals{TimPeriod,x }(m).MaxXcorrLag),2),ThisPairs,'UniformOutput',0))'*1000/AnalysisOpts.MotifAnalysis.FsWaveTarg;
                            AvgXcorrPairs.MaxXcorr=cell2mat(arrayfun(@(x) mean(Vals{TimPeriod,x }(m).MaxXcorr,2),ThisPairs,'UniformOutput',0))';
                            AvgXcorrPairs.ZscoreZeroLag=cell2mat(arrayfun(@(x) Vals{TimPeriod,x }(m).ZscoreZeroLag',ThisPairs,'UniformOutput',0))';
                            AvgXcorrPairs.ZscoreMaxXcorr=cell2mat(arrayfun(@(x) Vals{TimPeriod,x }(m).ZscoreMaxXcorr',ThisPairs,'UniformOutput',0))';
                            
                            subplot(NRow,NCol,f+(i-1)*NCol)
                            hold on
                            [~,PCAMeanZeroLag]=pca(AvgXcorrPairs.(XcorrFields{f})');
                            arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',TrlCOL(x,:),'MarkerSize',10),1:46)
                            arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
                            xlabel(['PC1 ' XcorrFields{f}]);ylabel(['PC2 ' XcorrFields{f}])
                            title([PairsToPlot{i} TimePeriodName{TimPeriod} num2str(m)])
                        end
                    end
                end
            end
            AnalysisOpts.CurrentCh='';
            [~,~,AnalysisData.XcorrValsFig]=obj.ManData.GetFileName([],'XcorrPCReward','SaveInResults',1,'WantedDate','ALL');
            obj.FigParams.SaveFigSeries([],AnalysisData.XcorrValsFig,[FigsPC])
            
            varargout=[Figs FigsPC];
            
            
            %                     figure(2)
            %                     subplot(5,6,sp)
            %                     hold on
            %                     [~,PCAMeanZeroLag]=pca(AvgPairsFixOn.MeanZeroLag');
            %                     hold on;arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',COL(x,:),'MarkerSize',10),1:46)
            %                     hold on;arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
            %
            %                      figure(3)
            %                     subplot(5,6,sp)
            %                     hold on
            %                     [~,PCAMeanZeroLag]=pca(AvgPairsFixOn.MaxXcorrLag');
            %                     hold on;arrayfun(@(x) plot(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),'.','color',COL(x,:),'MarkerSize',10),1:46)
            %                     hold on;arrayfun(@(x) text(PCAMeanZeroLag(x,1),PCAMeanZeroLag(x,2),num2str(opts.trlSeqAvg(x))),1:46)
            
            
            
            %                 % plot target on now
            %                 Fig_TargOn{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                 AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 2}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' TargOn ' num2str(sp)])
            %                 end
            %
            %                 % plot resp on
            %                 Fig_RespOn{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                  AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 3}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' RespOn ' num2str(sp)])
            %                 end
            %                 % plot reward
            %                 Fig_Reward{i}=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                  AvgPairsFixOn=arrayfun(@(x) mean(obj.ManData.ReshapeCell2Mat(XCorrPair{1, 4}(ThisPairs,x),3),3),1:NCoreMotifs,'UniformOutput',0);
            %                 for sp=1:NCoreMotifs
            %                     subplot(5,6,sp)
            %                     plot(-10:9,movmean(max(AvgPairsFixOn{sp},[],1),5));
            %                     xlabel('trl to switch') ;
            %                     ylabel('max xcorr');
            %                     title([PairsToPlot{i} ' Reward ' num2str(sp)])
            
            
            
            
            %  obj.ManData.SaveVar([],AvgTimBoxTrls,'AvgTimBoxTrls','AvgTimBoxTrlsData');
            
            
        end
        function varargout=PlotCoherenceStatsArea(obj,PairIndsAll,wavcohPairAll,varargin) % plots coherence stats between areas
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            FrontParietPairs=PairIndsAll.FrontParietPairs;
            FEFPairs=PairIndsAll.FEFPairs;
            PFCPairs=PairIndsAll.PFCPairs;
            LIPPairs=PairIndsAll.LIPPairs;
            AllPairs=PairIndsAll.AllPairs;
            
            % load all of the pairs first (this may take a while)
            NavgTimeEpochSet=[5 5]; % how many averages per time epoch
            TimeEpochAvgInd=1; % we have for now [1=avg 1 trl target on, 2=5 target on, 3=1 Respon 4=5 Respon]
            opts.NavgTimeEpoch=NavgTimeEpochSet(TimeEpochAvgInd);
            opts.Ntrls2Switch2Look=10; % How many trials around the switch we want to look at
            opts.UseZscore=1;
            %             % load ms coherence data
            %             mscohPairAll=obj.ManData.LoadVarSeries('','mscohPair','','CohData',AllPairs);
            %             mscohPair.coh=cellfun(@(x) x{TimeEpochAvgInd,end},mscohPairAll,'UniformOutput',0);  % esing 'end' to incluse all of the blocks
            %             mscohPair.f_mscoh=mscohPair.coh{1}.f;
            
            % load wavelet coherence data
            % Load data first
            %   wavcohPairAll=obj.ManData.LoadVarSeries('','wavcohPair','','CohData',AllPairs);
            if opts.UseZscore
                wavcohPair_FixOn =cellfun(@(x) obj.ManData.ReshapeCell2Mat(x{1,end}.cohrZscore,3),wavcohPairAll,'UniformOutput',0);
                wavcohPair_TargOn=cellfun(@(x) obj.ManData.ReshapeCell2Mat(x{2,end}.cohrZscore,3),wavcohPairAll,'UniformOutput',0);
                wavcohPair_RespOn=cellfun(@(x) obj.ManData.ReshapeCell2Mat(x{3,end}.cohrZscore,3),wavcohPairAll,'UniformOutput',0);
                wavcohPair_Reward=cellfun(@(x) obj.ManData.ReshapeCell2Mat(x{4,end}.cohrZscore,3),wavcohPairAll,'UniformOutput',0);
            else
                wavcohPair_FixOn=cellfun(@(x) obj.ManData.ReshapeCell2Mat(x{1,end}.cohr,3),wavcohPairAll,'UniformOutput',0);
                wavcohPair_TargOn=cellfun(@(x) obj.ManData.ReshapeCell2Mat(x{2,end}.cohr,3),wavcohPairAll,'UniformOutput',0);
                wavcohPair_RespOn=cellfun(@(x) obj.ManData.ReshapeCell2Mat(x{3,end}.cohr,3),wavcohPairAll,'UniformOutput',0);
                wavcohPair_Reward=cellfun(@(x) obj.ManData.ReshapeCell2Mat(x{4,end}.cohr,3),wavcohPairAll,'UniformOutput',0);
            end
            % get the fix on data
            wavcohPairStacked_FixOn.coh=obj.ManData.ReshapeCell2Mat(wavcohPair_FixOn,4); % stack pairs on 4th dim
            wavcohPairStacked_FixOn.f_wavcoh=obj.ManData.LoadVar('','f_Linear','TrialData_Wave',0);
            % get the target on data
            wavcohPairStacked_TargOn.coh=obj.ManData.ReshapeCell2Mat(wavcohPair_TargOn,4); % stack pairs on 4th dim
            wavcohPairStacked_TargOn.f_wavcoh=obj.ManData.LoadVar('','f_Linear','TrialData_Wave',0);
            % get the response on data
            wavcohPairStacked_RespOn.coh=obj.ManData.ReshapeCell2Mat(wavcohPair_RespOn,4); % stack pairs on 4th dim
            wavcohPairStacked_RespOn.f_wavcoh=obj.ManData.LoadVar('','f_Linear','TrialData_Wave',0);
            % get the reward data
            wavcohPairStacked_Reward.coh=obj.ManData.ReshapeCell2Mat(wavcohPair_Reward,4); % stack pairs on 4th dim
            wavcohPairStacked_Reward.f_wavcoh=obj.ManData.LoadVar('','f_Linear','TrialData_Wave',0);
            %% prepare ms coherence
            %% calculate some of the vars we need
            NTrl=size(wavcohPairStacked_TargOn.coh,3); % has been changed beacuse of averaging
            trlSeq=-10:5;%-34:(35-(70-NTrl));
            SwitchTrl=find(trlSeq==0);
            opts.trl2lookInd=(SwitchTrl-opts.Ntrls2Switch2Look):(SwitchTrl+opts.Ntrls2Switch2Look-1-4); % trials to look at
            opts.trlSeq2look=trlSeq;%-opts.Ntrls2Switch2Look:opts.Ntrls2Switch2Look;
            
            % now run this for a bunch of different conditions and save off
            % the plots
            TimeBoxFixOn ={[-0.6 0],[-0.6 0],[-0.6 0],[0 0.5],[0 0.5]};FreqBoxFixOn ={[17 35],[0 6],[7 15],[17 35],[35 55]};
            TimeBoxTargOn={[-0.6 0],[-0.6 0],[-0.6 0],[0 0.5],[0 0.5]};FreqBoxTargOn={[17 35],[0 6],[7 15],[17 35],[35 55]};
            TimeBoxRespOn={[-0.6 0],[-0.6 0],[-0.6 0],[0 0.5],[0 0.5]};FreqBoxRespOn={[17 35],[0 6],[7 15],[17 35],[35 55]};
            TimeBoxReward={[-0.6 0],[-0.6 0],[-0.6 0],[0 0.5],[0 0.5]};FreqBoxReward={[17 35],[0 6],[7 15],[17 35],[35 55]};
            
            PairsToPlot={'FrontParietPairs','FEFPairs','PFCPairs','LIPPairs','AllPairs'};
            for i=1:length(PairsToPlot)
                eval(['ThisPairs=' PairsToPlot{i} ';']);
                % plot fix on
                Fig_FixOn=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                opts.TimeBox=TimeBoxFixOn;opts.FreqBox=FreqBoxFixOn;
                [wavcohPairStacked(i).FixOn,AvgTimBoxTrls(i).FixOn,Figs_FixOn{i}]=obj.PlotCoherencePairs(wavcohPairStacked_FixOn,ThisPairs,opts,PairsToPlot{i},'FixOn',Fig_FixOn,1);
                
                % plot target on now
                Fig_TargOn=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                opts.TimeBox=TimeBoxTargOn;opts.FreqBox=FreqBoxTargOn;
                [wavcohPairStacked(i).TargOn,AvgTimBoxTrls(i).TargOn,Figs_TargOn{i}]=obj.PlotCoherencePairs(wavcohPairStacked_TargOn,ThisPairs,opts,PairsToPlot{i},'TargOn',Fig_TargOn,1);
                
                % plot resp on
                Fig_RespOn=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                opts.TimeBox=TimeBoxRespOn;opts.FreqBox=FreqBoxRespOn;
                [wavcohPairStacked(i).TargOn,AvgTimBoxTrls(i).FixOn,Figs_RespOn{i}]=obj.PlotCoherencePairs(wavcohPairStacked_RespOn,ThisPairs,opts,PairsToPlot{i},'RespOn',Fig_RespOn,1);
                
                % plot reward
                Fig_Reward=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
                opts.TimeBox=TimeBoxReward;opts.FreqBox=FreqBoxReward;
                [wavcohPairStacked(i).Reward,AvgTimBoxTrls(i).Reward,Figs_Reward{i}]=obj.PlotCoherencePairs(wavcohPairStacked_Reward,ThisPairs,opts,PairsToPlot{i},'Reward',Fig_Reward,1);
            end
            obj.ManData.SaveVar([],AvgTimBoxTrls,'AvgTimBoxTrls','AvgTimBoxTrlsData');
            
            %             % plot all of the individual pairs as well now
            %             for i=1:NPairs
            %                 Figh=obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %                 opts.TimeBox=TimeBoxRespOn;opts.FreqBox=FreqBoxRespOn;
            %                 obj.PlotCoherencePairs(wavcohPairStacked_RespOn,i,opts,['Pair' num2str(i)],'RespOn',Figh,1);
            %                 % plot target on now
            %                 opts.TimeBox=TimeBoxTargOn;opts.FreqBox=FreqBoxTargOn;
            %                 IndPairFigs{i}=obj.PlotCoherencePairs(wavcohPairStacked_TargOn,i,opts,['Pair' num2str(i)],'TargOn',Figh,2);
            %             end
            AnalysisOpts.AnalysisPathName='data';
            [FigsSaveFileName]=GenerateFileName(AnalysisOpts.FS,AnalysisOpts.ResultsSavePath,AnalysisOpts.AnalysisPathName,AnalysisOpts.Animal,AnalysisOpts.RecDate,'',['AllPairsCoh']);
            obj.FigParams.SaveFigSeries(FigsSaveFileName,AnalysisOpts.ResultsSavePath,[Figs_FixOn Figs_TargOn Figs_RespOn Figs_Reward ])
            
            varargout=[Figs_FixOn Figs_TargOn Figs_RespOn Figs_Reward ];
            if 0
                for b=1:4
                    figure
                    % cluster these values now
                    X=indAllCohMat(f_ind_set{b},:,:);
                    X1=squeeze(mean(X,1))';
                    TsneX1=tsne(X1);
                    Clust=Phenograph(X1);
                    Col=distinguishable_colors(length(unique(Clust)));
                    % subplot(4,2,1+2*(b-1))
                    subplot(121)
                    hold on
                    arrayfun(@(x) plot(TsneX1(Clust==x,1),TsneX1(Clust==x,2),'.','color',Col(x,:)),1:length(unique(Clust)));
                    xlabel('Tsne1');ylabel('Tsne2');title(['Clusting of pairs ' F_vals{b}]);
                    subplot(122)
                    hold on;arrayfun(@(x) plot(trlSeq,mean(X1(Clust==x,:),1),'color',Col(x,:)),1:length(unique(Clust)));
                    xlabel('Trl to switch');ylabel('Z-coh');title(['mean z-coh of pairs ' F_vals{b}]);
                end
                %%% look at the Xcorr now
                T_Xcorr=T_XcorrH{1, 1}{1, 1};
                ZerolagInd=find(T_Xcorr==0);
                for mot=1:22
                    for b=1:13
                        for p=FrotoParaietalPairs'
                            ind_p=find(p==FrotoParaietalPairs);
                            XcorrMatAll{mot}(:,:,ind_p,b)=XcorrHavg{p,b}{1,mot};
                        end
                    end
                end
                %% Xcorr for H vals analysis
                for mot=1:22
                    Figs{14+mot}=figure;
                    set(gcf,'Units','normalized','Position',[0 0 0.9 0.9])
                    MeanBlkXcorr{mot}=mean(XcorrMatAll{mot},4);
                    MeanPairXcorr{mot}=mean(XcorrMatAll{mot},3);
                    MeanBlkPairXcorr{mot}=mean(MeanBlkXcorr{mot},3);
                    subplot(3,5,1)
                    % imagesc(MeanBlkPairXcorr{mot})
                    helperCWTTimeFreqPlot(AnalysisData.CoreMotifs{mot} ,(1:size(AnalysisData.CoreMotifs{mot},2))/50,AnalysisData.cwt_f,'justplot1',['Mtf ' num2str(mot)],'Time(s)','Freq',0);
                    subplot(3,5,2);hold on
                    plot(-34:35,MeanBlkPairXcorr{mot}(ZerolagInd,:))
                    plot(-34:35,movmean(MeanBlkPairXcorr{mot}(ZerolagInd,:),5),'r')
                    xlabel('Trl 2 switch');ylabel('Xcorr')
                    title('Avg Xcorr for all PFC-LIP pairs all Blks')
                    for i=1:13
                        subplot(3,5,2+i)
                        plot(-34:35,squeeze(MeanPairXcorr{mot}(ZerolagInd,:,i)))
                        xlabel('Trl 2 switch');ylabel('Xcorr')
                        title(['Avg Xcorr Blk' num2str(i)])
                    end
                end
            end
            
        end
        function [wavcohPairStacked,AvgTimBoxTrls,varargout]=PlotCoherencePairs(obj,wavcohPair,ThisPairs,opts,Title,AxisTitle,Figh,StrRow)
            global AnalysisOpts AnalysisData
            
            % define params
            opts.FsWave=AnalysisOpts.MotifAnalysis.FsWave; % tomake it easier
            opts.NPairs=length(ThisPairs);
            if opts.NPairs>1
                opts.NRow=3;
            else
                opts.NRow=2;
            end
            opts.NCol=1+length(opts.TimeBox);
            %% define time and frequency boxing
            %   opts.TimeBox={[-0.6 0],[-0.6 0],[0.2 0.5]};
            %   opts.FreqBox={[12 30],[0 6],[12 30]};
            opts.Col={'r','b','k','m','g','m'};
            %% %%%%%%%%%%%%%%%%% wavelet coherence analysis
            %% prepare wavelet coherence
            wavcohPairStacked=wavcohPair.coh(:,:,opts.trl2lookInd,ThisPairs); % limit analysis to the trials around switch
            AvgWavcohPairs=mean(wavcohPairStacked,4); % avg acorss pairs keep trials
            AvgWavcohTrials=squeeze(mean(wavcohPairStacked,3)); % avg acorss trls keep pairs
            AvgWavcohTrialsPairs=mean(AvgWavcohPairs,3); % avg across trls and pairs
            %% preppare parameters
            TimeAxis=obj.GenerateTimeAxis(-1,0.5,opts.FsWave);
            FreqAxis=wavcohPair.f_wavcoh;
            
            %% compute average in the time boxed period
            NBoxs=length(opts.TimeBox);
            for box=1:NBoxs
                ThisTimeInd{box}=TimeAxis>=opts.TimeBox{box}(1) & TimeAxis<=opts.TimeBox{box}(2);
                ThisFreqInd{box}=FreqAxis>=opts.FreqBox{box}(1) & FreqAxis<=opts.FreqBox{box}(2);
                TimBoxTrls{box}=wavcohPairStacked(ThisFreqInd{box},ThisTimeInd{box},:,:);
                AvgTimBoxTrls{box}=squeeze(mean(mean(TimBoxTrls{box},2),1))';
            end
            %% save off some vars for future analysis
            
            %% plot everythign now
            varargout=Figh;%obj.FigParams.RenderFigure(1,[0 0 0.9 0.9]);
            %% figure plot for each pair all of the information across blcoks
            subplot(opts.NRow,opts.NCol,1+opts.NCol*(StrRow-1))
            % plot average wave coherence for all of the trials and pairs
            % to find time boxing
            helperCWTTimeFreqPlot(AvgWavcohTrialsPairs,TimeAxis,FreqAxis,'image',[Title ' Coh'],['Time(s) from ' AxisTitle],'Freq',0);
            colorbar
            % plot the rectangles now
            hold on
            arrayfun(@(x) rectangle('Position',[opts.TimeBox{x}(1),opts.FreqBox{x}(1),opts.TimeBox{x}(2)-opts.TimeBox{x}(1),opts.FreqBox{x}(2)-opts.FreqBox{x}(1)],...
                'EdgeColor',opts.Col{x},'LineWidth',3),1:NBoxs);
            % chi=get(gca, 'Children');
            % Reverse the stacking order so that the patch overlays the line
            % set(gca, 'Children',flipud(chi));
            
            %% plot avg of each time box with respect to trials
            for box=1:NBoxs
                subplot(opts.NRow,opts.NCol,1+box+opts.NCol*(StrRow-1))
                Title=sprintf('%i-%iHz',opts.FreqBox{box}(1),opts.FreqBox{box}(2));
                obj.FigParams.PlotMeanStd(opts.trlSeq2look,AvgTimBoxTrls{box},[],'Trl to switch','mean Z-scored Coh',opts.Col{box},0,Title)
                hold on
                v=axis;
                plot([-opts.NavgTimeEpoch -opts.NavgTimeEpoch],[v(3) v(4)],'r')
            end
            %% now plot clustring results for each of these conditions if the number of pairs is more than 1
            if opts.NPairs>1
                arrayfun(@(x) obj.ClusterCoherencePairs(AvgTimBoxTrls{x},opts,x,StrRow),1:NBoxs);
            end
            
            
            if 0
                %% plot ms coh pairs
                % plot average of all of the pairs
                f_ind_all=f>0 & f<80;
                indAllCohMatAvgTrl=squeeze(mean(indAllCohMat,2));
                PlotMeanStd(f(f_ind_all),indAllCohMatAvgTrl(f_ind_all,:)',[],['Freq'],['Mean Z-Coherence'],'b',1,'')
                % plot wavelet coherence
                indAllCohMat=[];
                indAllCohMatShuff=[];
                
                
                for j=1:length(Pairs)
                    for i=1:NTrl
                        indAllCohMat(:,i,j)=ThismscohPair{j}.cohrZscore{i};
                        % TempAll=zscore([(mscohPair{j}.Shuffcohr{i}); (mscohPair{j}.cohr{i}) ]);
                        % indAllCohMatShuff(:,i,j)=real(mean(TempAll(1:end-1,:),1));
                    end
                end
                
                f_ind_set{1}=f>2 & f<=6;
                f_ind_set{2}=f>=10 & f<=30;
                f_ind_set{3}=f>=17 & f<=25;
                f_ind_set{4}=f>=0 & f<=150;
                F_vals={'0-10Hz','10-25Hz','35-80Hz','100-150Hz'};
                Figs{i}=figure;
                
                for b=1:4 % b represents bad
                    f_ind=f_ind_set{b};
                    indAllCohMatWant=indAllCohMat(f_ind,:,:); % take the frequency band we are intrested
                    
                    indAllCohMatWantMeanFreq=squeeze(mean(indAllCohMatWant,1));
                    
                    RealindAllCohMatWant=mean(indAllCohMatWant,3);
                    % RealindAllCohMatWantShuff=mean(indAllCohMatShuff(f_ind,:,Pair,Blks),3);
                    
                    subplot(4,2,1+2*(b-1))
                    imagesc(trlSeq,f(f_ind),movmean(RealindAllCohMatWant,1,2));
                    xlabel('Trl to switch');ylabel('Freq');
                    title(['FrontoParietal Z-scored Coh ' F_vals{b}] )
                    colorbar
                    subplot(4,2,2+2*(b-1))
                    PlotMeanStd(trlSeq,indAllCohMatWantMeanFreq',[],[],[],'b',1,'')
                    hold on
                    %  plot(trlSeq,mean(RealindAllCohMatWantShuff,1),'r')
                    xlabel('Trl to switch');ylabel('mean Z-scored Coh');
                    title('Avg FrontoParietal Coh')
                    v=axis;
                    plot([-10 -10],[v(3) v(4)],'r')
                end
                
            end
        end
        function  ClusterCoherencePairs(obj,X,opts,Box,StrRow,varargin) % clusters coherence pairs and plots them
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            if isempty(opts);opts.NRow=2;opts.NCol=2;StrRow=0;Box=0;
                opts.trlSeq2look=-10:5;opts.NavgTimeEpoch=5;end
            
            opts.DimRed='PC'; % can be 'TSNE' or 'umap'
            %% Calculate TSNE and then cluster with phenograph
            %  Xmap=X-mean(X,2);%mapminmax(X,0,1); % map everything between 0 and 1 to pick up only on the shape
            % Xmap=mapminmax(X,0,1);
            Xmap=X;
            if strcmpi(opts.DimRed,'PC')
                [~,XDimRc]=pca(Xmap);
            elseif strcmpi(opts.DimRed,'TSNE')
                [XDimRc]=tsne(Xmap);
            elseif strcmpi(opts.DimRed,'umap')
                XDimRc=run_umap(Xmap,'metric','correlation','method','Java','verbose','none');
            end
            %Clust=Phenograph(XDimRc,'DistanceMetric','correlation');
            Clust=Phenograph(XDimRc);
            
            NClust=length(unique(Clust));
            Col=distinguishable_colors(NClust);
            % plot the results where we want them to be
            
            subplot(opts.NRow,opts.NCol,1+Box+opts.NCol*(StrRow))
            hold on;
            arrayfun(@(x) obj.FigParams.PlotMeanStd(opts.trlSeq2look,X(Clust==x,:),[],'Trl to switch',...
                'Z-coh',Col(x,:),3,'mean z-coh of pairs'),1:NClust);
            v=axis;
            plot([-opts.NavgTimeEpoch -opts.NavgTimeEpoch],[v(3) v(4)],'r')
            
            subplot(opts.NRow,opts.NCol,1+Box+opts.NCol*(StrRow+1))
            hold on
            arrayfun(@(x) obj.FigParams.Plot(XDimRc(Clust==x,1),XDimRc(Clust==x,2),Col(x,:),...
                [opts.DimRed ' 1'],[opts.DimRed ' 2'],'Clusting ','p_line_style','.'),1:NClust);
        end
        function  TimeAxis=GenerateTimeAxis(~,StrTim,StpTim,Fs) % generates a time axis
            Ts=1/Fs;
            TimeAxis=StrTim:Ts:(StpTim-Ts);
            
        end
        function  [H_exmp_tot,BlkTrialNumber,ExpBlkNum,Bhv_blk_Perf]=GrabHTimeTrialOld(obj,H,BlockSpec,varargin)  
            % grabs the H data centered on different episodes of the task
            % has all blocks mean at the end
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            NBlks=length(BlockSpec.Rule);
            % take  example blocks and mean block
            ExpBlkNum=[1:NBlks NaN];
            BlkTrialNumber=arrayfun(@(x) BlockSpec.ThisBlkTrialNum(x),1:NBlks,'UniformOutput',1);
            BlkTrialNumber=[BlkTrialNumber BlkTrialNumber(1)];
            % grab data now
            H_exmp_tot.FixOn=[]; H_exmp_tot.TargetOn=[];H_exmp_tot.Reward=[];H_exmp_tot.RespOn=[];
            for ch=1:length(H)
                H_blk.FixOn=[];H_blk.TargetOn=[];H_blk.Reward=[];H_blk.RespOn=[]; % initialize H blocks
                for b=1:NBlks
                    ThisTrialTimes=BlockSpec.ThisBlkTrialTimes{b};
                    %% grab data centered on Fixation
                    AnalysisOpts.TrialTiming.StartFieldName='FIX_ON';
                    H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,'UseDataPointer',obj.UseDataPointer,'DataPointerVar',obj.DataPointerVar,...
                        'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                    H_blk.FixOn(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    
                    %% grab data centered on stimulus presentaion
                    AnalysisOpts.TrialTiming.StartFieldName='TARGET_ON';
                    H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,'UseDataPointer',obj.UseDataPointer,'DataPointerVar',obj.DataPointerVar,...
                        'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                    H_blk.TargetOn(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    
                    %% grab data centered on Reward Time
                    AnalysisOpts.TrialTiming.StartFieldName='GIVE_REWARD';
                    H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,'UseDataPointer',obj.UseDataPointer,'DataPointerVar',obj.DataPointerVar,...
                        'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                    H_blk.Reward(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    
                    %% grab data centered on Response Time
                    AnalysisOpts.TrialTiming.StartFieldName='RESPONSE_ON';
                    H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,'UseDataPointer',obj.UseDataPointer,'DataPointerVar',obj.DataPointerVar,...
                        'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                    H_blk.RespOn(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    
                    % calculate performance
                    if ch==1
                        Bhv_blk_Perf(b,:)=obj.CalBhvPerf(ThisTrialTimes,10);
                    end
                end
                
                %% re-organize H values now
                % fix on
                H_exmp=[arrayfun(@(x) H_blk.FixOn(:,:,:,x),1:NBlks,'UniformOutput',0) {mean(H_blk.FixOn,4)}];
                H_exmp_tot.FixOn=[H_exmp_tot.FixOn;H_exmp];
                % target on
                H_exmp=[arrayfun(@(x) H_blk.TargetOn(:,:,:,x),1:NBlks,'UniformOutput',0) {mean(H_blk.TargetOn,4)}];
                H_exmp_tot.TargetOn=[H_exmp_tot.TargetOn;H_exmp];
                % reward on
                H_exmp=[arrayfun(@(x) H_blk.Reward(:,:,:,x),1:NBlks,'UniformOutput',0) {mean(H_blk.Reward,4)}];
                H_exmp_tot.Reward=[H_exmp_tot.Reward;H_exmp];
                % Response on
                H_exmp=[arrayfun(@(x) H_blk.RespOn(:,:,:,x),1:NBlks,'UniformOutput',0) {mean(H_blk.RespOn,4)}];
                H_exmp_tot.RespOn=[H_exmp_tot.RespOn;H_exmp];
            end
            
        end
        function  [H_exmp_tot,BlkTrialNumber,ExpBlkNum,Bhv_blk_Perf]=GrabHTimeTrial(obj,H,BlockSpec,varargin)  
            % grabs the H data centered on different episodes of the task
            % doesn't have all blocks mean
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            WantedStartFieldName=cellfun(@(x) sum(strcmpi(x,obj.StartFieldName)),{'FIX_ON','TARGET_ON','GIVE_REWARD','RESPONSE_ON'});
            
            NBlks=length(BlockSpec.Rule);
            % take  example blocks and mean block
            ExpBlkNum=1:NBlks;
            BlkTrialNumber=arrayfun(@(x) BlockSpec.ThisBlkTrialNum(x),1:NBlks,'UniformOutput',1);
            % grab data now
            H_exmp_tot.FixOn=[]; H_exmp_tot.TargetOn=[];H_exmp_tot.Reward=[];H_exmp_tot.RespOn=[];
            for ch=1:length(H)
                H_blk.FixOn=[];H_blk.TargetOn=[];H_blk.Reward=[];H_blk.RespOn=[]; % initialize H blocks
                for b=1:NBlks
                    ThisTrialTimes=BlockSpec.ThisBlkTrialTimes{b};
                    %% grab data centered on Fixation
                    if WantedStartFieldName(1)                   
                        AnalysisOpts.TrialTiming.StartFieldName='FIX_ON';
                        H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,'UseDataPointer',obj.UseDataPointer,'DataPointerVar',obj.DataPointerVar,...
                            'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                        H_blk.FixOn(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    end
                    %% grab data centered on stimulus presentaion
                    if WantedStartFieldName(2)
                        AnalysisOpts.TrialTiming.StartFieldName='TARGET_ON';
                        H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,'UseDataPointer',obj.UseDataPointer,'DataPointerVar',obj.DataPointerVar,...
                            'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                        H_blk.TargetOn(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    end
                    %% grab data centered on Reward Time
                    if WantedStartFieldName(3)
                        AnalysisOpts.TrialTiming.StartFieldName='GIVE_REWARD';
                        H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,'UseDataPointer',obj.UseDataPointer,'DataPointerVar',obj.DataPointerVar,...
                            'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                        H_blk.Reward(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    end
                    %% grab data centered on Response Time
                    if WantedStartFieldName(4)
                        AnalysisOpts.TrialTiming.StartFieldName='RESPONSE_ON';
                        H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,'UseDataPointer',obj.UseDataPointer,'DataPointerVar',obj.DataPointerVar,...
                            'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum(b),'UniformOutput',0);
                        H_blk.RespOn(:,:,:,b)=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
                    end
                    % calculate performance
                    if ch==1
                        Bhv_blk_Perf(b,:)=obj.CalBhvPerf(ThisTrialTimes,10);
                    end
                end
                
                %% re-organize H values now
                % fix on
                if WantedStartFieldName(1)
                    H_exmp=arrayfun(@(x) H_blk.FixOn(:,:,:,x),1:NBlks,'UniformOutput',0);
                    H_exmp_tot.FixOn=[H_exmp_tot.FixOn;H_exmp];
                end                
                % target on
                if WantedStartFieldName(2)
                    H_exmp=arrayfun(@(x) H_blk.TargetOn(:,:,:,x),1:NBlks,'UniformOutput',0);
                    H_exmp_tot.TargetOn=[H_exmp_tot.TargetOn;H_exmp];
                end
                % reward on
                if WantedStartFieldName(3)
                    H_exmp=arrayfun(@(x) H_blk.Reward(:,:,:,x),1:NBlks,'UniformOutput',0);
                    H_exmp_tot.Reward=[H_exmp_tot.Reward;H_exmp];
                end
                % Response on
                if WantedStartFieldName(4)
                    H_exmp=arrayfun(@(x) H_blk.RespOn(:,:,:,x),1:NBlks,'UniformOutput',0);
                    H_exmp_tot.RespOn=[H_exmp_tot.RespOn;H_exmp];
                end
            end           
        end
        function  [H_tot]=GrabHTimeTrialEpisode(obj,H,TrialTimes,FieldName,varargin)  % grabs the H data centered on specific episodes of the task
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            % grab data now
            ThisTrialTimes=TrialTimes;
            Ntrls=size(ThisTrialTimes,1);
            AnalysisOpts.TrialTiming.StartFieldName=FieldName;
            Nch=length(AnalysisData.Ch);
            H_tot=cell(Nch,1);
            for ch=1:Nch
                %% grab data centered on Fixation
                fprintf('Fetching H for channel %',AnalysisData.Ch(ch));
                H_temp=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(H{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,...
                    'Fs',obj.Fs,'AnalysisType','Wavelet_Trial_H'),1:Ntrls,'UniformOutput',0);
                H_tot{ch,1}=obj.ManData.ReshapeCell2Mat(H_temp,3); % 1st dim: Cluster, 2nd: Time,3rd: Trial, 4th Block
            end
        end
        function MaxXcorr=RearrangXcorrHchPair(obj,XcorrPairBlck,T_XcorrPairBlck,ChPairs)
            % rearranges Xcorr matrix so we have similarity matrix for each
            % motif across trials
            Nmotifs=length(XcorrPairBlck{1});
            Ntrl=size(XcorrPairBlck{1,1}{1, 1},2);
            NChPair=size(ChPairs,1);
            NCh=length(unique(ChPairs(:)));
            ZerolagInd=(T_XcorrPairBlck==0);
            X=zeros(size(NCh));
            for m=1:Nmotifs
                for trl=1:Ntrl
                    % fetch the value of xcorr for this motif and trial at
                    % zero-lag
                    temp=arrayfun(@(x) mean(XcorrPairBlck{x}{m}(ZerolagInd,trl),3),1:NChPair);
                    temp=obj.ManData.ReshapeSquareMatrix(ChPairs,temp);
                    temp=X+triu(temp,1)+triu(temp,1)';
                    MaxXcorr{m}(:,:,trl)=temp;
                end
            end
        end
        function [XcorrPair,XcorrPairAvg,T_XcorrPair,MeanShuffAll]=CalXcorrHChPair(obj,H1,H2,varargin)  % calculates xcorr bewtween pair of channels across trails and blocks
            % H1 and H2 are H values  organized in Motif*Time*Trial
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            NMotif=size(H1,1);
            Ntim=size(H1,2);
            NTrl=size(H1,3);
            Nshuff=obj.Nshuffle; % repeat shuffling 50 times
            
            [~,TXcorr]=xcorr(ones(Ntim,1),ones(Ntim,1));
            TXcorr=TXcorr/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            IndTXcorr=(TXcorr>=-AnalysisOpts.MotifAnalysis.XcorrTimInt & TXcorr <=AnalysisOpts.MotifAnalysis.XcorrTimInt);
            for m=1:NMotif
                ThisMotifH1=squeeze(H1(m,:,:));ThisMotifH2=squeeze(H2(m,:,:));
                % ThisMotifH1=ThisMotifH1-mean(ThisMotifH1,1);
                % ThisMotifH2=ThisMotifH2-mean(ThisMotifH2,1);
                
                %  compute shuffle xcorr
                for nsh=1:Nshuff
                    RndOrd1=randperm(NTrl,NTrl);RndOrd2=randperm(NTrl,NTrl);
                    XcorrPair_shuff=cell2mat(arrayfun(@(x) xcorr(ThisMotifH1(:,RndOrd1(x)),ThisMotifH2(:,RndOrd2(x)),'coef'),1:NTrl,'UniformOutput',0));
                    MeanShuff{m}(:,nsh)=mean(XcorrPair_shuff,2);
                end
                MeanShuffAll{m}=mean(MeanShuff{m},2);
                %
                [XcorrPair{m},T_XcorrPair]=(arrayfun(@(x) xcorr(ThisMotifH1(:,x),ThisMotifH2(:,x),'coef'),1:NTrl,'UniformOutput',0));
                XcorrPair{m}=cell2mat(XcorrPair{m})-repmat(MeanShuffAll{m},[1 NTrl]);
                XcorrPair{m}=XcorrPair{m}(IndTXcorr,:); % take only the time period we care about
                XcorrPairAvg{m}=movmean(XcorrPair{m},obj.Navg,2); % take average across trials as well
            end
            T_XcorrPair=T_XcorrPair{1}(IndTXcorr);
            %%
            %             % chop based on Trials
            %             TrlStr=[1 floor(NTrl/3) 2*floor(NTrl/3) NTrl];
            %             TrlIntrvl=[TrlStr(1:end-1)' TrlStr(2:end)'];
            %             for  m=1:NMotif
            %                 MeanXcorr=arrayfun(@(x) mean(XcorrPair{m}(:,TrlIntrvl(x,1):TrlIntrvl(x,2)),2),1:3,'UniformOutput',0);
            %             end
            %
        end
        function [mscohPair,cpsdPair,wavcohPair]=CalCoherenceChPairOld(obj,Raw1,Raw2,Fs,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % Ch1 and Ch2 are raw data values  organized in 1*Time*Trial
            opts.NHamming=200;
            opts.Noverlap=199;
            
            NMotif=size(Raw1,1);
            Ntim=size(Raw1,2);
            NTrl=size(Raw1,3);
            for m=1:NMotif
                ThisMotifRaw1=squeeze(Raw1(m,:,:));ThisMotifRaw2=squeeze(Raw2(m,:,:));
                % calculate magnitude squared  coherence
                [mscohPair.coh{m},mscohPair.f]=(arrayfun(@(x) mscohere(ThisMotifRaw1(:,x),ThisMotifRaw2(:,x),hamming(opts.NHamming),opts.Noverlap,[],Fs),1:NTrl,'UniformOutput',0));
                % calculate cross spectrum
                [cpsdPair.cpsd{m},cpsdPair.f]=(arrayfun(@(x) cpsd(ThisMotifRaw1(:,x),ThisMotifRaw2(:,x),hamming(opts.NHamming),opts.Noverlap,[],Fs),1:NTrl,'UniformOutput',0));
                % calculate wavlet coherence
                [wavcohPair.wavcoh{m},wavcohPair.cs{m},wavcohPair.f,wavcohPair.coi{m}]=(arrayfun(@(x) wcoherence(ThisMotifRaw1(:,x),ThisMotifRaw2(:,x),Fs),1:NTrl,'UniformOutput',0));
            end
            mscohPair.f=mscohPair.f{1};cpsdPair.f=cpsdPair.f{1};wavcohPair.f=wavcohPair.f{1};
        end
        function CalCoherenceAllChPairs(obj,BlockSpec,RawData_exmp_tot,Wavelet_exmp_tot,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            Nch=length(AnalysisData.Ch);
            NBlks=length(BlockSpec.Rule);
            % Blks2LookHCohr=1:NBlks; % which blocks are we looking at for H coherence
            % Blks2LookCohr=mat2cell([1:NBlks 1:NBlks],1,[ones(1,NBlks) NBlks]); % blcoks for rawdata coherence; get the blocks one by one and then all of them together
            Blks2LookCohr=mat2cell(1:NBlks,1,NBlks); % Get all of the blocks blcoks for rawdata coherence; get the blocks one by one and then all of them together
            
            ChPairs=nchoosek(1:Nch,2);
            AnalysisData.TimeFreq.NavgCoh=5; % how many trials to average for coherence analysis
            
            if strcmpi(AnalysisOpts.PairNum,'ALL') % if we want it all the sweep across
                AnalysisOpts.PairNum=1:size(ChPairs,1);
            end
            
            for Chp=AnalysisOpts.PairNum
                Chp1=ChPairs(Chp,1);Chp2=ChPairs(Chp,2);
                fprintf('\nCalculating Coherence for Pair %i',Chp)
                
                %% calculate raw coherence for all of the blocks
                for blk=1:length(Blks2LookCohr)
                    Blks2Look=Blks2LookCohr{blk}; % list of blocks to look
                    if 0 % no ms coehernce for now
                        %% calculate coherence for LFP signal now with 1 and 5 averaging
                        % target on as reference
                        RawData1=RawData_exmp_tot.TargetOn(Chp1,Blks2Look);RawData2=RawData_exmp_tot.TargetOn(Chp2,Blks2Look);
                        % [mscohPair,cpsdPair,wavcohPair]=obj.CalCoherenceChPair(RawData1,RawData2,FsLFP,10);
                        %  eval(['mscohPair_' num2str(Chp) '{1,' num2str(blk) '}=obj.CalCoherenceBlockChPair(RawData1,RawData2,AnalysisOpts.MotifAnalysis.FsLFP,1);']);
                        eval(['mscohPair_' num2str(Chp) '{1,' num2str(blk) '}=obj.CalCoherenceBlockChPair(RawData1,RawData2,AnalysisOpts.MotifAnalysis.FsLFP,AnalysisData.TimeFreq.NavgCoh);']);
                        % eval(['mscohPair_' num2str(Chp) '{3,' num2str(blk) '}=obj.CalCoherenceBlockChPair(RawData1,RawData2,AnalysisOpts.MotifAnalysis.FsLFP,10);']);
                        % response on as reference
                        RawData1=RawData_exmp_tot.RespOn(Chp1,Blks2Look);RawData2=RawData_exmp_tot.RespOn(Chp2,Blks2Look);
                        % eval(['mscohPair_' num2str(Chp) '{3,' num2str(blk) '}=obj.CalCoherenceBlockChPair(RawData1,RawData2,AnalysisOpts.MotifAnalysis.FsLFP,1);']);
                        eval(['mscohPair_' num2str(Chp) '{2,' num2str(blk) '}=obj.CalCoherenceBlockChPair(RawData1,RawData2,AnalysisOpts.MotifAnalysis.FsLFP,AnalysisData.TimeFreq.NavgCoh);']);
                        
                    end
                    %% calculate wavelet coherence for LFP signal with 1 and 5 averaging
                    % Fix on
                    WaveData1=Wavelet_exmp_tot.FixOn(Chp1,Blks2Look);WaveData2=Wavelet_exmp_tot.FixOn(Chp2,Blks2Look);
                    eval(['wavcohPair_' num2str(Chp) '{1,' num2str(blk) '}=obj.CalWaveCoherenceBlockChPair(WaveData1,WaveData2,AnalysisOpts.MotifAnalysis.FsWave,AnalysisData.TimeFreq.NavgCoh);']);
                    
                    % target on
                    WaveData1=Wavelet_exmp_tot.TargetOn(Chp1,Blks2Look);WaveData2=Wavelet_exmp_tot.TargetOn(Chp2,Blks2Look);
                    eval(['wavcohPair_' num2str(Chp) '{2,' num2str(blk) '}=obj.CalWaveCoherenceBlockChPair(WaveData1,WaveData2,AnalysisOpts.MotifAnalysis.FsWave,AnalysisData.TimeFreq.NavgCoh);']);
                    
                    % response on
                    WaveData1=Wavelet_exmp_tot.RespOn(Chp1,Blks2Look);WaveData2=Wavelet_exmp_tot.RespOn(Chp2,Blks2Look);
                    eval(['wavcohPair_' num2str(Chp) '{3,' num2str(blk) '}=obj.CalWaveCoherenceBlockChPair(WaveData1,WaveData2,AnalysisOpts.MotifAnalysis.FsWave,AnalysisData.TimeFreq.NavgCoh);']);
                    
                    % Reward on
                    WaveData1=Wavelet_exmp_tot.Reward(Chp1,Blks2Look);WaveData2=Wavelet_exmp_tot.Reward(Chp2,Blks2Look);
                    eval(['wavcohPair_' num2str(Chp) '{4,' num2str(blk) '}=obj.CalWaveCoherenceBlockChPair(WaveData1,WaveData2,AnalysisOpts.MotifAnalysis.FsWave,AnalysisData.TimeFreq.NavgCoh);']);
                    
                end
                %  eval(['obj.ManData.SaveVar([],mscohPair_' num2str(Chp) ',''mscohPair_' num2str(Chp) ''',''CohData_' num2str(Chp) ''');']);
                
                eval(['obj.ManData.SaveVar([],wavcohPair_' num2str(Chp) ',''wavcohPair_' num2str(Chp) ''',''CohData_' num2str(Chp) ''');']);
            end
            %% calculate PSD for each channel now if we are on the first pair
            %             if strcmpi(AnalysisOpts.PairNum,'ALL')
            %                 for Ch=AnalysisData.Ch
            %                     ChInd=find(Ch==AnalysisData.Ch);
            %                     Blks2Look=Blks2LookCohr{end};
            %                     RawData=RawData_exmp_tot.TargetOn(Ch,Blks2Look);
            %                     PSD_CH{ChInd}=obj.CalPSDBlockCh(RawData,AnalysisOpts.MotifAnalysis.FsLFP,AnalysisData.TimeFreq.NavgCoh); % average for 5 trials for this
            %                 end
            %                 obj.ManData.SaveVar([],PSD_CH,'PSD_CH','TrialData');
            %             end
        end
        function CalHCoherenceAllChPairs(obj,BlockSpec,H_exmp_tot,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            Nch=length(AnalysisData.Ch);
            NBlks=length(BlockSpec.Rule);
            Blks2LookHCohr=1:NBlks; % which blocks are we looking at for H coherence
            % Blks2LookCohr=mat2cell([1:NBlks 1:NBlks],1,[ones(1,NBlks) NBlks]); % blcoks for rawdata coherence; get the blocks one by one and then all of them together
            Blks2LookCohr=mat2cell(1:NBlks,1,NBlks); % Get all of the blocks blcoks for rawdata coherence; get the blocks one by one and then all of them together
            ChPairs=nchoosek(1:Nch,2);
            AnalysisData.TimeFreq.NavgCoh=5; % how many trials to average for coherence analysis
            
            if strcmpi(AnalysisOpts.PairNum,'ALL') % if we want it all the sweep across
                AnalysisOpts.PairNum=1:size(ChPairs,1);
            end
            
            for Chp=AnalysisOpts.PairNum
                Chp1=ChPairs(Chp,1);Chp2=ChPairs(Chp,2);
                fprintf('\nCalculating Coherence for Pair %i',Chp)
                %% calculate cross correlation for motifs for all of the
                %% channels
                if ~isempty(H_exmp_tot)
                    for blk=Blks2LookHCohr
                        blkInd=find(blk==Blks2LookHCohr);
                        %                         % Cal for FixOn
                        %                         H1=H_exmp_tot.FixOn{Chp1,blk};H2=H_exmp_tot.FixOn{Chp2,blk};
                        %                         eval(['[XcorrH_' num2str(Chp) '{1,blkInd},~,TXcorr]= obj.CalXcorrHChPair(H1,H2,''Navg'',obj.Navg,''Nshuffle'',obj.Nshuffle);']);
                        
                        % Cal for TargetOn
                        H1=H_exmp_tot.TargetOn{Chp1,blk};H2=H_exmp_tot.TargetOn{Chp2,blk};
                        eval(['[XcorrH_' num2str(Chp) '{1,blkInd},~,TXcorr]= obj.CalXcorrHChPair(H1,H2,''Navg'',obj.Navg,''Nshuffle'',obj.Nshuffle);']);
                        
                        %                          % Cal for RespOn
                        %                         H1=H_exmp_tot.RespOn{Chp1,blk};H2=H_exmp_tot.RespOn{Chp2,blk};
                        %                         eval(['XcorrH_' num2str(Chp) '{3,blkInd}= obj.CalXcorrHChPair(H1,H2,''Navg'',obj.Navg,''Nshuffle'',obj.Nshuffle);']);
                        
                        % Cal for Reward
                        H1=H_exmp_tot.Reward{Chp1,blk};H2=H_exmp_tot.Reward{Chp2,blk};
                        eval(['XcorrH_' num2str(Chp) '{2,blkInd}= obj.CalXcorrHChPair(H1,H2,''Navg'',obj.Navg,''Nshuffle'',obj.Nshuffle);']);
                    end
                end
                obj.ManData.DeleteFile([],['HCohData_' num2str(Chp)],1);% delete previous file
                eval(['obj.ManData.SaveVar([],XcorrH_' num2str(Chp) ',''XcorrH_' num2str(Chp) ''',''HCohData_' num2str(Chp) ''');']);
                eval(['obj.ManData.SaveVar([],TXcorr,''TXcorr'',''HCohData_' num2str(Chp) ''');']);
                
            end
            
        end
        
        function [mscohPair,wavcohPair]=CalCoherenceChPair(obj,Raw1,Raw2,Fs,NTrlAvg,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % Ch1 and Ch2 are raw data values  organized in 1*Time*Trial
            % NTrlAvg number of trials to be averaged for coherence stats
            
            opts.UseHann=1; % if use hanning window for FFT
            NMotif=size(Raw1,1);
            Ntim=size(Raw1,2);
            NTrl=size(Raw1,3);
            % generate the inds for averaging acorss trials
            AvgInd=obj.ManData.GenMovAvgInds(NTrl,NTrlAvg,'same');
            
            for m=1:NMotif
                for Trls=1:size(AvgInd,1)
                    ThisTrls=AvgInd(Trls,:);
                    if (ThisTrls(2)-ThisTrls(1))>0
                        ThisMotifRaw1=[squeeze(Raw1(m,:,ThisTrls(1):ThisTrls(2)))]';ThisMotifRaw2=[squeeze(Raw2(m,:,ThisTrls(1):ThisTrls(2)))]';
                    elseif (ThisTrls(2)-ThisTrls(1))==0
                        ThisMotifRaw1=[squeeze(Raw1(m,:,ThisTrls(1):ThisTrls(2)))];ThisMotifRaw2=[squeeze(Raw2(m,:,ThisTrls(1):ThisTrls(2)))];
                    end
                    % calculate magnitude squared coherence
                    [mscohPair.cohr{Trls},mscohPair.phase{Trls},mscohPair.f]=obj.TimFreq.CalCoherence(ThisMotifRaw1,ThisMotifRaw2,Fs,opts.UseHann,0);
                end
            end
            wavcohPair=[];
        end
        function [mscohPairFinal]=CalCoherenceBlockChPair(obj,Raw1,Raw2,Fs,NTrlAvg,varargin)
            % includes the trials across all of the blocks to calculate
            % coherence
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % Ch1 and Ch2 are raw data values  organized in 1*Time*Trial
            % NTrlAvg number of trials to be averaged for coherence stats
            
            opts.UseHann=obj.UseHann; % if use hanning window for FFT
            NMotif=size(Raw1{1},1);
            Ntim=size(Raw1{1},2);
            NTrl=size(Raw1{1},3);
            NBlks=size(Raw1,2); % how many blocks we want to average
            NTrlBlk=NTrlAvg*NBlks; % number of  Trials * Num blocks
            % generate the inds for averaging acorss trials
            AvgInd=obj.ManData.GenMovAvgInds(NTrl,NTrlAvg,'valid');
            
            for m=1:NMotif
                for Trls=1:size(AvgInd,1)
                    ThisTrls=AvgInd(Trls,:);ThisMotifRaw1=[];ThisMotifRaw2=[];
                    for blk=1:NBlks
                        if (ThisTrls(2)-ThisTrls(1))>0
                            ThisMotifRaw1=[ThisMotifRaw1;transpose(squeeze(Raw1{blk}(m,:,ThisTrls(1):ThisTrls(2))))];ThisMotifRaw2=[ThisMotifRaw2;transpose(squeeze(Raw2{blk}(m,:,ThisTrls(1):ThisTrls(2))))];
                        elseif (ThisTrls(2)-ThisTrls(1))==0
                            ThisMotifRaw1=[ThisMotifRaw1;squeeze(Raw1{blk}(m,:,ThisTrls(1):ThisTrls(2)))];ThisMotifRaw2=[ThisMotifRaw2;squeeze(Raw2{blk}(m,:,ThisTrls(1):ThisTrls(2)))];
                        end
                    end
                    % calculate magnitude squared coherence
                    [mscohPair.cohr{Trls},mscohPair.phase{Trls},mscohPair.f]=obj.TimFreq.CalCoherence(ThisMotifRaw1,ThisMotifRaw2,Fs,opts.UseHann,0);
                    % calculate the shuffle of this coherence
                    [Shuffcohr,Shuffphase]=arrayfun (@(x) obj.TimFreq.CalCoherence(ThisMotifRaw1(randperm(NTrlBlk),:),ThisMotifRaw2(randperm(NTrlBlk),:),Fs,opts.UseHann,0),1:obj.Nshuffle,'UniformOutput',0);
                    mscohPair.Shuffcohr{Trls}= obj.ManData.ReshapeCell2Mat(Shuffcohr,2);
                    mscohPair.Shuffphase{Trls}=obj.ManData.ReshapeCell2Mat(Shuffphase,2);
                    % now calculate z-score of the true value
                    allCohZScore=zscore([mscohPair.Shuffcohr{Trls} ;mscohPair.cohr{Trls}]);
                    mscohPair.cohrZscore{Trls}=allCohZScore(end,:);
                    % take the mean of the shuffle and store that only
                    mscohPair.Shuffcohr{Trls}=mean(mscohPair.Shuffcohr{Trls},1);
                    mscohPair.Shuffphase{Trls}=mean(mscohPair.Shuffphase{Trls},1);
                end
            end
            % only keep the variables we care about
            mscohPairFinal.cohr=mscohPair.cohr;
            mscohPairFinal.cohrZscore=mscohPair.cohrZscore;
            mscohPairFinal.f=mscohPair.f;
            
        end
        function [wavcohPairFinal]=CalWaveCoherenceBlockChPair(obj,WaveData1,WaveData2,FsWave,NTrlAvg,varargin)
            % calculates wavelet coherence
            % includes the trials across all of the blocks to calculate
            % coherence
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % Ch1 and Ch2 are raw data values  organized in 1*Time*Trial
            % NTrlAvg number of trials to be averaged for coherence stats
            Nfreq=size(WaveData1{1},1);
            Ntim=size(WaveData1{1},2);
            NTrl=size(WaveData1{1},3);
            NBlks=size(WaveData1,2); % how many blocks we want to average
            NTrlBlk=NTrlAvg*NBlks; % number of  Trials * Num blocks
            % generate the inds for averaging acorss trials
            AvgInd=obj.ManData.GenMovAvgInds(NTrl,NTrlAvg,'valid');
            
            for Trls=1:size(AvgInd,1)
                ThisTrls=AvgInd(Trls,:);ThisWavData1=[];ThisWavData2=[];
                for blk=1:NBlks  % add trials from all of the blocks here now
                    ThisWavData1=cat(3,ThisWavData1,WaveData1{blk}(:,:,ThisTrls(1):ThisTrls(2)));
                    ThisWavData2=cat(3,ThisWavData2,WaveData2{blk}(:,:,ThisTrls(1):ThisTrls(2)));
                end
                
                % calculate wavelet coherence
                [wavcohPair.cohr{Trls},wavcohPair.phase{Trls}]=obj.TimFreq.CalWaveCoherence(ThisWavData1,ThisWavData2,FsWave,0);
                
                % calculate the shuffle of this coherence
                [Shuffcohr,Shuffphase]=arrayfun (@(x) obj.TimFreq.CalWaveCoherence(ThisWavData1(:,:,randperm(NTrlBlk)),ThisWavData2(:,:,randperm(NTrlBlk)),FsWave,0),1:obj.Nshuffle,'UniformOutput',0);
                wavcohPair.Shuffcohr{Trls}= obj.ManData.ReshapeCell2Mat(Shuffcohr,3);
                wavcohPair.Shuffphase{Trls}=obj.ManData.ReshapeCell2Mat(Shuffphase,3);
                
                % now calculate z-score of the true value
                allCohZScore=zscore(cat(3,wavcohPair.Shuffcohr{Trls},wavcohPair.cohr{Trls}),0,3); % zscore along 3rd dim
                wavcohPair.cohrZscore{Trls}=allCohZScore(:,:,end);
                
                % only store the mean of the shuffle
                wavcohPair.Shuffcohr{Trls}=mean(wavcohPair.Shuffcohr{Trls},3);
                wavcohPair.Shuffphase{Trls}=mean(wavcohPair.Shuffphase{Trls},3);
            end
            % keep the variables we care about
            wavcohPairFinal.cohr=wavcohPair.cohr;
            wavcohPairFinal.cohrZscore=wavcohPair.cohrZscore;
            
        end
        function  PSD=CalPSDBlockCh(obj,RawData,Fs,NTrlAvg,varargin) % calculates PSD across trials and blocks
            % includes the trials across all of the blocks to calculate
            % PSD
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % RawData are raw data values  organized in 1*Time*Trial
            % NTrlAvg number of trials to be averaged for coherence stats
            
            opts.UseHann=obj.UseHann; % if use hanning window for FFT
            NMotif=size(RawData{1},1);
            Ntim=size(RawData{1},2);
            NTrl=size(RawData{1},3);
            NBlks=size(RawData,2); % how many blocks we want to average
            NTrlBlk=NTrlAvg*NBlks; % number of  Trials * Num blocks
            % generate the inds for averaging acorss trials
            AvgInd=obj.ManData.GenMovAvgInds(NTrl,NTrlAvg,'valid');
            
            
            for m=1:NMotif
                for Trls=1:size(AvgInd,1)
                    ThisTrls=AvgInd(Trls,:);ThisRaw1=[];
                    for blk=1:NBlks
                        if (ThisTrls(2)-ThisTrls(1))>0
                            ThisRaw1=[ThisRaw1;transpose(squeeze(RawData{blk}(m,:,ThisTrls(1):ThisTrls(2))))];
                        elseif (ThisTrls(2)-ThisTrls(1))==0
                            ThisRaw1=[ThisRaw1;squeeze(RawData{blk}(m,:,ThisTrls(1):ThisTrls(2)))];
                        end
                    end
                    % calculate PSD for this trial sets
                    [PSD.mean_PSD{Trls},PSD.f_PSD,PSD.df,PSD.fNQ]=obj.TimFreq.CalPSD(ThisRaw1,Fs,opts.UseHann);
                end
            end
            
        end
        function varargout=PlotXcorrHChPair(obj,XcorrPairBlck,T_XcorrPair,TrialOrder,BlkNum,Bhv_blk_Perf,varargin) % plots H cross correlation between pairs of channels
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            Navg=obj.Navg;
            NMotifs=length(XcorrPairBlck{1});
            varargout=obj.FigParams.RenderFigure(NMotifs,[]);
            NExmp=size(XcorrPairBlck,2);
            for m=1:NMotifs
                figure(varargout{m})
                % plot imagesc of trials and xcorr
                for i=1:NExmp
                    ThisXcorr=movmean(XcorrPairBlck{i}{m},obj.Navg,2);
                    [MaxCorr,MaxCorrInd]=max(ThisXcorr,[],1);
                    % plot image of xcorr
                    subplot(4,NExmp,i)
                    helperCWTTimeFreqPlot(ThisXcorr,TrialOrder,T_XcorrPair,...
                        'image',['Block' num2str(BlkNum(i)) ],...
                        'Trial','Time(s)',obj.LogScale)
                    
                    % plot Max X corr
                    subplot(4,NExmp,NExmp+i);hold on
                    plot(TrialOrder,MaxCorr)
                    xlabel('Trial');ylabel('Max xcorr')
                    % plot ind of max corr
                    subplot(4,NExmp,2*NExmp+i)
                    plot(TrialOrder,T_XcorrPair(MaxCorrInd))
                    xlabel('Trial');ylabel('Max time')
                    
                    % plot behavior and report correlation to Max value
                    subplot(4,NExmp,3*NExmp+i);hold on
                    plot(TrialOrder,Bhv_blk_Perf(i,:),'r')
                    [a,p]=corr(MaxCorr',Bhv_blk_Perf(i,:)');
                    xlabel('Trial');ylabel('bhv Perf')
                    title(['Corr' num2str(a,3),'p' num2str(p,3)])
                end
            end
            
        end
        function  varargout=PlotInfo_H(obj,TrialTimes,BlockSpec,varargin) % discovers what is the information in Hs
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            varargout=[];
            figure;hold on
            Col=distinguishable_colors(length(AnalysisData.Ch));
            % set up trial times
            ThisTrialTimes=BlockSpec.ThisBlkTrialTimes{1};
            %   ThisTrialTimes(:,1:5)=ThisTrialTimes(:,1:5)-BlockSpec.ThisBlkTrialTimes{1}(1);
            for ch=1:length(AnalysisData.Ch)
                H=arrayfun(@(x)  obj.TrialFunc.GrabDataTimeTrial(...
                    AnalysisData.ClustH{ch} ,ThisTrialTimes(x,:),[],'WaveletDownSampleFactor',1,...
                    'Fs',obj.Fs,...
                    'AnalysisType','Wavelet_Trial_H'),1:BlockSpec.ThisBlkTrialNum,'UniformOutput',0);
                H=obj.ManData.Cell2Mat(H);
                H=obj.ManData.ReshapeTrials(H,BlockSpec.ThisBlkTrialNum(1)); % now reshape these trials
                BlockTrialTimes=cellfun(@(x) TrialTimes(x,:),BlockSpec.ThisBlkTrials(1),'UniformOutput',0);
                %% train classifer on each Time Point
                Time=-1:1/obj.Fs:(0.5-0.01); %% what was atually in the experiment
                NTim=length(Time);NHs=size(H,1);
                Class_TrlTims=cell2mat(BlockTrialTimes');
                Data_Hs=arrayfun(@(x) cell2mat(H(x,:)'),1:NHs,'UniformOutput',0); % organize Hs
                
                Perf_Avgh=obj.H_Info_SVM(Class_TrlTims,Data_Hs,Time,'avgh',BlockSpec);
                Perf_Trial=obj.H_Info_SVM(Class_TrlTims,Data_Hs,Time,'trial',BlockSpec);
                Perf_Rule   =obj.H_Info_SVM(Class_TrlTims,Data_Hs,Time,'rule');
                Perf_Axis   =obj.H_Info_SVM(Class_TrlTims,Data_Hs,Time,'axis');
                Perf_Feature=obj.H_Info_SVM(Class_TrlTims,Data_Hs,Time,'feature');
                [varargout]=obj.Plot_H_Info(Perf_Avgh,Perf_Trial,Perf_Rule,Perf_Axis,Perf_Feature);
                
            end
        end
        function  [PerfTest,score,varargout]=SVM_Binary(obj,Class,Data,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            NdataPnt=size(Data,1);
            NTrainSamp=floor(NdataPnt*obj.ClassTrainFrac);
            TrainingInds=randsample(1:NdataPnt,NTrainSamp);
            TestInds=setdiff(1:NdataPnt,TrainingInds);
            DataTrain=Data(TrainingInds,:);ClassTrain=Class(TrainingInds);  % train data
            DataTest =Data(TestInds,:);    ClassTest=Class(TestInds);     % test data
            
            %Train the SVM Classifier
            cl = fitcsvm(DataTrain,ClassTrain,'KernelFunction','linear');
            % predict test
            [label,score] = predict(cl,DataTest);
            
            CorrPredict=label==ClassTest;
            PerfTest=sum(CorrPredict)/length(CorrPredict);
            
        end
        function  [PerfTest,ScoreTest,varargout]=SVM_Binary_Time(obj,Class,Data_Hs,Time,NHs,WantedTrls,varargin)
            global AnalysisOpts  AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            NTim=length(Time); if isempty(Time);NTim=1;end
            if isempty(WantedTrls);WantedTrls=1:size(Class,1);end
            for t=1:NTim
                fprintf('\nTime point:%i',t)
                Data=cell2mat(arrayfun(@(x) Data_Hs{x}(WantedTrls,t),NHs,'UniformOutput',0));
                [PerfTest_1{t},ScoreTest{t}]=arrayfun(@(x) obj.SVM_Binary(Class,Data),1:obj.Nrep,'UniformOutput',0);
            end
            PerfTest.Mean=arrayfun(@(x) mean(cell2mat(PerfTest_1{x})),1:NTim);
            PerfTest.STD =arrayfun(@(x) std(cell2mat(PerfTest_1{x})),1:NTim);
            
            
        end
        function  [PerfTest,ScoreTest,varargout]=H_Info_SVM(obj,Class,Data_Hs,Time,Feature,varargin) %cals H info about features
            global AnalysisOpts  AnalysisData
            if ~isempty(varargin);BlockSpec=varargin{1};end
            NHs=1:7;
            RuleInd=strcmp(AnalysisOpts.TrialTimesFields,'CONDITION_NUM_OFFSET'); %% we don't use this here but still good to have it
            switch  Feature
                case 'axis'
                    ClassAxis=Class(:,RuleInd);
                    ClassAxis(ClassAxis==3)=1;  % look at axis
                    
                    [PerfTest,ScoreTest]=SVM_Binary_Time(obj,ClassAxis,Data_Hs,Time,NHs,[])  ;
                    
                case 'feature'
                    ClassFeature=Class(:,RuleInd);
                    ClassFeature(ClassFeature==3)=2; % Look at feature
                    
                    [PerfTest,ScoreTest]=SVM_Binary_Time(obj,ClassFeature,Data_Hs,Time,NHs,[])   ;
                    
                case 'rule'   % look at rule
                    RuleInd=strcmp(AnalysisOpts.TrialTimesFields,'CONDITION_NUM_OFFSET'); %% we don't use this here but still good to have it
                    Class=Class(:,RuleInd);
                    %%  compare rule 1 and 2
                    ThisRuleTrls=(Class==1 | Class==2);
                    ClassRule=Class(ThisRuleTrls);
                    [PerfTest.R12,ScoreTest.R12]=obj.SVM_Binary_Time(ClassRule,Data_Hs,Time,NHs,ThisRuleTrls);
                    %%  Compare Rule 1 and 3
                    ThisRuleTrls=(Class==1 | Class==3);
                    ClassRule=Class(ThisRuleTrls);
                    [PerfTest.R13,ScoreTest.R13]=obj.SVM_Binary_Time(ClassRule,Data_Hs,Time,NHs,ThisRuleTrls);
                    %%  Compare Rule 2 and 3
                    ThisRuleTrls=(Class==2 | Class==3);
                    ClassRule=Class(ThisRuleTrls);
                    [PerfTest.R23,ScoreTest.R23]=obj.SVM_Binary_Time(ClassRule,Data_Hs,Time,NHs,ThisRuleTrls);
                    
                case 'trial'   % go forward and back in trials
                    NTrl=10;NMaxTrl=50;Nsteps=NMaxTrl-NTrl+1;
                    RuleTrls=arrayfun(@(x) obj.FindRuleTrials(x,BlockSpec,NTrl,NMaxTrl),1:3,'UniformOutput',0);
                    RuleComb=[1 2;1 3;2 3];
                    for i=1:3
                        Rules=RuleComb(i,:);
                        for s=1:Nsteps  % loop on time points
                            ThisRuleTrls=[RuleTrls{Rules(1)}.Forward{s} RuleTrls{Rules(2)}.Forward{s}];
                            LR(1)=length(RuleTrls{Rules(1)}.Forward{s});
                            LR(2)=length(RuleTrls{Rules(2)}.Forward{s});
                            ClassRule=[Rules(1)*ones(LR(1),1);Rules(2)*ones(LR(2),1)];
                            temp=obj.SVM_Binary_Time(ClassRule,Data_Hs,Time,NHs,ThisRuleTrls);
                            PerfTest{i}(s,:)=temp.Mean;
                        end
                    end
                    ScoreTest=[];
                case 'avgh'  % plots average H across time and trials
                    NTrl=1;NMaxTrl=30;
                    Rule=Class(:,RuleInd);
                    for i=1:3
                        RuleInd=find(Rule==i);
                        %% Time Avg
                        PerfTest.TimeAvg{i}=arrayfun(@(x) mean(Data_Hs{x}(RuleInd,:),1),NHs,'UniformOutput',0);
                        %% Trial Avg
                        RuleTrls=obj.FindRuleTrials(i,BlockSpec,NTrl,NMaxTrl);
                        for j=1:NMaxTrl
                            ForwardTrls= cell2mat(RuleTrls.ForwardVert{j});
                            BakwardTrls= cell2mat(RuleTrls.BakwardVert{j});
                            PerfTest.TrialAvgFrw{i}(j,:)=cell2mat(arrayfun(@(x) mean(mean(Data_Hs{x}(ForwardTrls,:),1)),NHs,'UniformOutput',0));
                            PerfTest.TrialAvgBak{i}(j,:)=cell2mat(arrayfun(@(x) mean(mean(Data_Hs{x}(BakwardTrls,:),1)),NHs,'UniformOutput',0));
                        end
                    end
                    
            end
        end
        function  [varargout]=Plot_H_Info(obj,Perf_Avgh,Perf_Trial,Perf_Rule,Perf_Axis,Perf_Feature) % plots results from H_Info
            Col=distinguishable_colors(3); % Rule colors
            Time=-1:.01:(0.5-0.01); %% what was atually in the experiment
            %% Plot Average H
            if ~isempty(Perf_Avgh)
                varargout{1}=figure;set(gcf,'Units','Normalized','Position',[0 0 0.8 0.8]);
                NH=length(Perf_Avgh.TimeAvg{1});
                NTim=length(Perf_Avgh.TimeAvg{1}{1});
                %% avg Time
                for i=1:NH
                    subplot(1,NH,i)
                    arrayfun(@(x) PlotMeanStd(Time,Perf_Avgh.TimeAvg{x}{i},...
                        [],'Time from Stim','H',Col(x,:)),1:3,'UniformOutput',0)
                    title(['Avg H(' num2str(i) ') in Time'])
                    axis square
                end
                %% Avg Trial
                varargout{2}=figure;set(gcf,'Units','Normalized','Position',[0 0 0.8 0.8]);
                NTrls=size(Perf_Avgh.TrialAvgBak{1},1);
                for i=1:NH
                    subplot(2,NH,i) % to switch
                    arrayfun(@(x) PlotMeanStd(-1*NTrls+1:0,Perf_Avgh.TrialAvgBak{x}(:,i),...
                        [],'Trial to switch','H',Col(x,:)),1:3,'UniformOutput',0)
                    title(['Avg H(' num2str(i) ') to switch'])
                    
                    subplot(2,NH,i+NH) % from swith
                    arrayfun(@(x) PlotMeanStd(0:NTrls-1,Perf_Avgh.TrialAvgFrw{x}(:,i),...
                        [],'Trial from switch','H',Col(x,:)),1:3,'UniformOutput',0)
                    title(['Avg H(' num2str(i) ') from switch'])
                end
            end
            %% plot SVM performance for rule
            RuleComb=[1 2;1 3;2 3];
            varargout{3}=figure;set(gcf,'Units','Normalized','Position',[0 0 0.8 0.8]);
            
            subplot(1,3,1)
            PlotMeanStd(Time,Perf_Rule.R12.Mean,Perf_Rule.R12.STD,'Time from Stim','SVM Perf Rule1&2',Col(1,:))
            subplot(1,3,2)
            PlotMeanStd(Time,Perf_Rule.R13.Mean,Perf_Rule.R13.STD,'Time from Stim','SVM Perf Rule1&3',Col(2,:))
            subplot(1,3,3)
            PlotMeanStd(Time,Perf_Rule.R23.Mean,Perf_Rule.R23.STD,'Time from Stim','SVM Perf Rule2&3',Col(3,:))
            %% Plot SVM Performance for feature
            varargout{4}=figure;set(gcf,'Units','Normalized','Position',[0 0 0.8 0.8]);
            subplot(121)
            PlotMeanStd(Time,Perf_Feature.Mean,Perf_Feature.STD,'Time from Stim','SVM Perf Feature',Col(1,:))
            subplot(122)
            PlotMeanStd(Time,Perf_Axis.Mean,Perf_Axis.STD,'Time from Stim','SVM Perf Axis',Col(2,:))
            %% Plot SMV Performance over Trial and Time
            varargout{5}=figure;set(gcf,'Units','Normalized','Position',[0 0 0.8 0.8]);
            for i=1:3
                subplot(1,3,i)
                helperCWTTimeFreqPlot(Perf_Trial{i},Time,1:size(Perf_Trial{i},1),...
                    'justplot2',['SVM Perf Rule:' num2str(i)],...
                    'Time(s)','Trials from Switch',obj.LogScale)
            end
            
        end
        function  [varargout]=PlotExampleHs(obj,X_Wave,X_Raw,W,H,f,Time,varargin) % plots random Hs and looks at reconstruction
            global AnalysisOpts   AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            varargout{1}=obj.FigParams.RenderFigure(1,[]);
            % take a random time
            X_Wave=obj.PrepapreCWT(X_Wave,'PowerMethod',AnalysisOpts.PowerMethod,'DownSampleFactor',1,'freq',f); % prepare CWT data for processing
            NRaw=4;NCol=3;
            NTim=size(H,2);NW=size(W,2);
            Lt=AnalysisOpts.MotifAnalysis.L;
            RandTimPnt=randsample([Lt:NTim-Lt],1);
            T=RandTimPnt-Lt:RandTimPnt+Lt-1;
            %% show raw data
            H_fig{1}=subplot(NRaw,NCol,1);
            % plot(Time(T),X_Raw(T),'r','linewidth',1.5);
            plot(X_Raw,'r','linewidth',1.5);
            xlabel('Time(s)')
            title('RawData');
            obj.FigParams.FormatAxes(H_fig{1})
            %% show PSD
            H_fig{2}=subplot(NRaw,NCol,2);
            helperCWTTimeFreqPlot(X_Wave(:,T),Time(T),f,'justplot1',['CWT CH' num2str(obj.ChNum)],...
                'Time(s)','Frequency',obj.LogScale);
            obj.FigParams.FormatAxes(H_fig{2})
            %% show H value
            H_fig{3}=subplot(NRaw,NCol,3);
            plot(Time(T),H(:,T),'linewidth',1.5);
            legend(arrayfun(@(x) ['W' num2str(x)],1:NW,'UniformOutput',0),'Location','northeastoutside')
            xlabel('Time(s)')
            title('H value');
            obj.FigParams.FormatAxes(H_fig{3})
            %% Reconstruct the signal now
            Rconst_X=helper.reconstruct(W,H(:,T));
            H_fig{4}=subplot(NRaw,NCol,4);
            helperCWTTimeFreqPlot(Rconst_X,Time(T),f,'justplot1','Reconstrcution',...
                'Time(s)','Frequency',obj.LogScale);
            obj.FigParams.FormatAxes(H_fig{4})
            %% Plot PEV for this example
            H_fig{5}=subplot(NRaw,NCol,5);
            MotifData{1}.W=W;MotifData{1}.H=H(:,T);X{1}=X_Wave(:,T);
            obj.PlotMotifPEV(X,MotifData,'PEVmethod','var')
            %% Plot Ws now
            k=1;
            for i=1:NW
                ThisW=squeeze(W(:,i,:));
                if sum(ThisW(:))~=0
                    H_fig{5}=subplot(NRaw,NCol,5+k);
                    helperCWTTimeFreqPlot(ThisW,(0:Lt-1)*10,f,'justplot1',['W ' num2str(i)],...
                        'Time(ms)','Frequency',obj.LogScale);
                    k=k+1;
                end
            end
            
        end
        function  varargout=CrosFreqAnalysis(obj,X_Wave,X_Raw,W,H,loadings,f,Time,MotifNum,PerctileTH,varargin) % Plots H triggered PSD and spiking
            global AnalysisOpts AnalysisData
            %PerctileTH the percentile where we find Th threshold for detecting H triggered events
            
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            opts.ThBinSize=0.5; % in ms size of bin to find maximum of data
            opts.PACfreqRange=[30 50]; % frequency range to plot PAC
            opts.Lwind_ms=200; % how many ms to look around max power of each event to find zero phase
            opts.LwindPow_ms=500;% % how many ms to look around mean max power  find zero phase
            opts.FreqThSet= [0.2]; %0.2 this val is multiplied to the max of meanW to find point of crossing
            opts.Tepoch_ms=300;% takes T ms tine around the max power (T/2 each side)
            opts.WantedMotifs=MotifNum; % which motifs to look at
            opts.FindCenterFreq=1; % find center frequency for phase proving frequency
            opts.FreqW=1; % width of frequency filtering for lower freq
            % Prepare varables
            Lt=AnalysisOpts.MotifAnalysis.L;
            FsLFP=AnalysisOpts.MotifAnalysis.FsLFP;
            FsWaveTarg=AnalysisOpts.MotifAnalysis.FsWaveTarg;
            [W,H]=obj.DelRedunWH(W,H,loadings); % remove zero loadings
            NTim=size(H,2);NW=size(W,2);
            
            X_Wave=obj.PrepapreCWT(X_Wave,'PowerMethod',AnalysisOpts.MotifAnalysis.PowerMethod,'DownSampleFactor',1,'freq',f); % prepare CWT data for processing
            %% calculte distribution of peak values of H and find the threshold for each motif
            [Th]=arrayfun(@(x) obj.CalDistPeakVals(H(x,:),opts.ThBinSize*FsWaveTarg,PerctileTH),1:NW);
            
            %%  Now extract H related signal from Rawdata, PSD and spiking data
            [HTrigData.PSD,Tim_PSD,HTrigData.PSDind,HTrigData.NEvent]=obj.ExtractHrelatedSig(Lt,Th,H,X_Wave,3);
            [HTrigData.Raw,Tim_Raw,HTrigData.Rawind]=obj.ExtractHrelatedSig(Lt,Th,H,X_Raw,2);
            HTrigData.Time=obj.ExtractHrelatedSig(Lt,Th,H,Time,2);  % chunk the times we have
            
            ValidData=find(arrayfun(@(x) ~isempty(HTrigData.PSD{x}),1:length(HTrigData.PSD)));
            T_PSD=Tim_PSD/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            T_Raw=Tim_Raw/AnalysisOpts.MotifAnalysis.FsLFP;
            T_W=(1:Lt)/AnalysisOpts.MotifAnalysis.FsWaveTarg;
            
            %%  Step 1: find the frequency of intrest for each Motif select one of them and then filter data in seperate bands
            varargout(1)=obj.FigParams.RenderFigure(1,[]); % plot core motifs
            obj.PlotCoreMotifs(W,[],'FigTitle','#')
            
            varargout(2)=obj.FigParams.RenderFigure(1,[]); % show where we are filtering the motifs
            SumW=sum(W,3); % take mean across time
            SumW=SumW./max(SumW,[],1);
            WantedMotifs= [ opts.WantedMotifs ];%11
            NWantedMotifs=length(WantedMotifs);
            FiltData=[];
            for i=WantedMotifs
                %% find Freq where there is crossing from FreqTh of Max value
                k=find(i==WantedMotifs); % index
                FreqTh=opts.FreqThSet(k);
                ThisW=SumW(:,i);
                ThThisW=FreqTh*max(ThisW(:));
                PassTh=find(diff(ThisW>=ThThisW))+1;
                if opts.FindCenterFreq
                    [~,CentFreqInd]=max(ThisW(PassTh));
                    CentFreq=f(PassTh(CentFreqInd));
                    % add 1 Hz for each band
                    FiltData(i).BandFreq(1)=CentFreq-opts.FreqW;
                    FiltData(i).BandFreq(2)=CentFreq+opts.FreqW;
                else
                    FiltData(i).BandFreq=f(PassTh);  % These are band frequencies
                    if mod(length(PassTh),2)~=0;FiltData(i).BandFreq=[0.1 FiltData(i).BandFreq];end
                    Nf=length(FiltData(i).BandFreq);
                    % add 1 Hz for each band
                    FiltData(i).BandFreq(1:2:Nf)=FiltData(i).BandFreq(1:2:Nf)-opts.FreqW;
                    FiltData(i).BandFreq(2:2:Nf)=FiltData(i).BandFreq(2:2:Nf)+opts.FreqW;
                    FiltData(i).BandFreq(FiltData(i).BandFreq<=0)=0.5;
                end
                %% now plot the data
                subplot(2,NWantedMotifs,k);hold on
                obj.FigParams.Plot(f,ThisW,'b','Freq','Power','Sum of power in Freq')
                arrayfun(@(x) plot(x,ThThisW,'r*'),FiltData(i).BandFreq);
                axis square
                subplot(2,NWantedMotifs,k+NWantedMotifs)
                obj.PlotCoreMotifs(W,i,'FigTitle','#')
                axis square
            end
            %% filter all of the data around these frequencies now
            for i=WantedMotifs
                [FiltData(i).FiltRaw]=arrayfun(@(x) obj.FilterFunc.BandPassFilter(X_Raw,...
                    AnalysisOpts.MotifAnalysis.FsLFP,'PassBand',FiltData(i).BandFreq(x:x+1)),...
                    find(mod(1:length(FiltData(i).BandFreq),2)~=0),...
                    'uniformoutput',0);
                
                [FiltData(i).HTrigFilt]=cellfun(@(x) obj.ExtractHrelatedSig(Lt,Th,H,x,2,i),FiltData(i).FiltRaw,...
                    'UniformOutput',0);
            end
            %% plot results now
            varargout(3)=obj.FigParams.RenderFigure(1,[]);
            for i=WantedMotifs
                NRow=2;NCol=4;
                bin=0.2;
                T_Raw_Step=bin*AnalysisOpts.MotifAnalysis.FsLFP;
                T_Raw_inds=1:T_Raw_Step:length(T_Raw);
                T_Raw_Tim=T_Raw(T_Raw_inds(1:end-1));
                %% align everything to the trough of lower frequency
                i=WantedMotifs;
                N=HTrigData.NEvent(i); % number of events
                Col=distinguishable_colors(N);
                % plot all of the trials and then plot the grand average
                subplot(321);hold on
                arrayfun(@(x) obj.FigParams.Plot(T_Raw,FiltData(i).HTrigFilt{1}{i}(x,:),Col(x,:),[],[],[]),1:N)
                % plot mean signal
                obj.FigParams.Plot(T_Raw,mean(FiltData(i).HTrigFilt{1}{i},1),'k','Time(s)','Volts','Non Aligned Events','p_line_width',5)
                
                % plot angle of hilbert transform
                subplot(322);hold on%obj.ManData.CalAngle
                AngleHilbert=cell2mat(arrayfun(@(x) angle(hilbert(FiltData(i).HTrigFilt{1}{i}(x,:)))',1:N,'UniformOutput',0))';
                arrayfun(@(x) obj.FigParams.Plot(T_Raw,AngleHilbert(x,:),Col(x,:),'Time(s)','Phase(Rad)','Instantanous Phase'),1:N)
                
                %             hold on; subplot(323);
                %             Tbin=cell2mat(arrayfun(@(x) [x;x+100],1:900,'UniformOutput',0));
                %             for k=1:N
                %             RMS(k,:)=arrayfun(@(x) rms(FiltData(i).HTrigFilt{1}{i}(k,Tbin(1,x):Tbin(2,x))),1:900);end
                %             hold on;arrayfun(@(x)  plot(T_Raw(1:900),RMS(x,:),'color',Col(x,:)),1:6)
                % plot power
                subplot(323);hold on
                obj.FigParams.Plot(Time,H(i,:),'b',[],[],[])
                v=axis;
                obj.FigParams.Plot([v(1) v(2)],[Th(i) Th(i)],'r','Time(s)','H',['H value motif ' num2str(i)] )
                
                
                subplot(324);%obj.ManData.CalAngle
                PWRsig=cell2mat(arrayfun(@(x)  transpose(abs(hilbert(FiltData(i).HTrigFilt{1}{i}(x,:))).^2),1:N,'UniformOutput',0))';
                MeanPWR=mean(PWRsig,1);
                hold on;
                arrayfun(@(x) obj.FigParams.Plot(T_Raw,PWRsig(x,:),Col(x,:),[],[],[]),1:N)
                % plot mean signal
                obj.FigParams.Plot(T_Raw,MeanPWR,'k','Time(s)','Power','Power of signal','p_line_width',5)
                
                % find max of the PWRsig and then take period around it and align
                % the min phase to it
                % look around 100 ms
                Lwind=opts.Lwind_ms*(1000/FsLFP);
                LwindPow=opts.LwindPow_ms*(1000/FsLFP);
                LSamp=length(MeanPWR);
                [~,PwrInd]=sort(MeanPWR,'descend');
                MaxPwrIndMean=PwrInd(find(PwrInd>Lwind & PwrInd<(LSamp-Lwind),1,'first'));
                % find max power of each event around this mean max power
                [~,MaxPwrInd]=arrayfun(@(x) max(PWRsig(x,MaxPwrIndMean-LwindPow:MaxPwrIndMean+LwindPow)),1:N);
                MaxPwrInd=MaxPwrInd+MaxPwrIndMean-LwindPow;
                % [~,MaxPwrInd]=max(MeanPWR);
                [~,indMinPhas]=arrayfun(@(x) min(AngleHilbert(x,MaxPwrInd(x)-Lwind:MaxPwrInd(x)+Lwind)),1:N);
                % plot the min phases to watch
                subplot(322)
                v=axis;
                arrayfun(@(x) plot(T_Raw(indMinPhas(x)+MaxPwrInd(x)-Lwind),v(3),'*','color',Col(x,:)),1:N);
                % shift everything towards the max power
                
                % take Tepoch around the max power (Tepoch/2 each side)
                Tepoch=opts.Tepoch_ms*(1000/FsLFP);
                TepochTim=(-Tepoch:Tepoch)/FsLFP;
                TroughInd= arrayfun(@(x) MaxPwrInd(x)-opts.Lwind_ms+indMinPhas(x),1:N); % what the alighned index
                ShifftSig=cell2mat(arrayfun(@(x) FiltData(i).HTrigFilt{1}{i}(x,(TroughInd(x)-Tepoch):(TroughInd(x)+Tepoch))',...
                    1:N,'UniformOutput',0))';
                ShifftSigRaw=cell2mat(arrayfun(@(x) HTrigData.Raw{i}(x,(TroughInd(x)-Tepoch):(TroughInd(x)+Tepoch))',...
                    1:N,'UniformOutput',0))';
                
                %  ShifftSig=cell2mat(arrayfun(@(x) circshift(FiltData(i).HTrigFilt{1}{i}(x,:),100-indMinPhas(x))',1:6,'UniformOutput',0));
                %  ShifftSigRaw=cell2mat(arrayfun(@(x) circshift(HTrigData.Raw{i}(x,:),100-indMinPhas(x))',1:6,'UniformOutput',0));
                %WantedInd=MaxPwrInd-2*Lwind:MaxPwrInd+2*Lwind;
                %WantedTim=-2*Lwind:2*Lwind;
                
                subplot(325);hold on
                arrayfun(@(x) obj.FigParams.Plot(TepochTim,ShifftSig(x,:),Col(x,:),[],[],[]),1:N)
                obj.FigParams.Plot(TepochTim,mean(ShifftSig,1),'k',['Time(s)'],['volts'],['Average Aligned Signal'],'p_line_width',5);
                
                % calculte time frequency representation for each filtered
                % epoch and average across them
                %  [TFR,f_epoch]=arrayfun(@(x) obj.TimFreq.CalTFRhilbert(ShifftSigRaw(x,:),FsLFP,[30 50],2,4),1:N,'UniformOutput',0);
                
                [TFR,f_epoch]=arrayfun(@(x) cwt(ShifftSigRaw(x,:),'amor',FsLFP,'VoicesPerOctave',20,...
                    'FrequencyLimits',opts.PACfreqRange ),1:N,'UniformOutput',0);
                f_epoch=f_epoch{1};
                TFR=cellfun(@(x) obj.PrepapreCWT(x,'PowerMethod',AnalysisOpts.MotifAnalysis.PowerMethod,'DownSampleFactor',1,'freq',f_epoch),TFR,'UniformOutput',0); % prepare CWT data for processing
                
                TFR=obj.ManData.ReshapeCell2Mat(TFR,3);
                %   TFR=TFR(:,WantedInd,:);
                subplot(326);hold on
                helperCWTTimeFreqPlot(mean(TFR,3),TepochTim,f_epoch,'image',[],'','',obj.LogScale);  axis square
                v=axis;
                MeanSigSupIm=mapminmax(mean(ShifftSig,1),v(3),v(4)); % superimposed mean signal
                obj.FigParams.Plot(TepochTim,MeanSigSupIm,'k',['Time(s)'],['Freq' ],['Avg PSD n=' num2str(N)],'p_line_width',1);
            end
            %% show raw data
            %             varargout(4)=obj.FigParams.RenderFigure(1,[]);
            %             CF=CrossFreqCopling;
            %             CF.CalComodulationCWT(X_Raw(1:30000) ,FsLFP)
            % %
            %            for i=WantedMotifs%length(ValidData)
            %                 varargout(find(i==WantedMotifs)+1)=obj.FigParams.RenderFigure(1,[]);
            %                 Ntrl=size(HTrigData.PSD{i},3);
            %                 exmp=1:Ntrl;
            %                 H_fig{1}=subplot(NRow,NCol,1);
            %                 plot(T_Raw,mean(FiltData(i).HTrigFilt{1}{i}(exmp,:) ,1))
            %                 hold on
            %                 plot(T_Raw,mean(FiltData(i).HTrigFilt{3}{i}(exmp,:) ,1),'r')
            %                 title('Mean LFP')
            %                 xlabel('time(s)');ylabel('Voltage')
            %                 H_fig{2}=subplot(NRow,NCol,2);
            %                 helperCWTTimeFreqPlot(mean(HTrigData.PSD{i}(:,:,exmp),3),T_PSD,f,'image',...
            %                     ['Avg PSD n=' num2str(size(HTrigData.PSD{i},3))],...
            %                     'Time(s)','Frequency',obj.LogScale);  axis square
            %                 now calculate it over time
            %                 calculate selectivity index for each mean trial and plot
            %                 distribution of that
            %                 Data1=mean(FiltData(i).HTrigFilt{3}{i},1);
            %                 Data2=mean(FiltData(i).HTrigFilt{1}{i},1);
            %                 [SIm_mean,SIpRad_mean ]=CalSynchronizationIndex(Data1,Data2,AnalysisOpts.MotifAnalysis.FsLFP);
            %                 calculate across time
            %                 for Ttrl=1:length(T_Raw_inds)-1
            %                     [SIm_tim_mean(Ttrl),SIpRad_tim_mean(Ttrl)]=CalSynchronizationIndex(Data1(T_Raw_inds(Ttrl):T_Raw_inds(Ttrl+1)),...
            %                         Data2(T_Raw_inds(Ttrl):T_Raw_inds(Ttrl+1)),AnalysisOpts.MotifAnalysis.FsLFP);
            %                 end
            %
            %                 H_fig{3}=subplot(NRow,NCol,3);hold on
            %                 plot(T_Raw_Tim,SIm_tim_mean);
            %                  xlabel('Time')
            %                 ylabel('SI index')
            %                 title('SI index for mean signal')
            %
            %                 H_fig{4}=subplot(NRow,NCol,4);hold on
            %                 plot(T_Raw_Tim,SIpRad_tim_mean,'r')
            %                 xlabel('Time')
            %                 ylabel('Radians')
            %                 title('Phase SI index for mean signal')
            %
            %                 %%%%%%%%%%%%%%%%%%%
            %                 exmp=1;
            %                 H_fig{5}=subplot(NRow,NCol,5);
            %                 plot(T_Raw,mean(FiltData(i).HTrigFilt{1}{i}(exmp,:) ,1))
            %                 hold on
            %                 plot(T_Raw,mean(FiltData(i).HTrigFilt{2}{i}(exmp,:) ,1),'r')
            %
            %                 xlabel('time(s)');ylabel('Voltage')
            %                 H_fig{6}=subplot(NRow,NCol,6);
            %                 helperCWTTimeFreqPlot(mean(HTrigData.PSD{i}(:,:,exmp),3),T_PSD,f,'image',...
            %                     ['Example PSD'],...
            %                     'Time(s)','Frequency',obj.LogScale);  axis square
            %
            %                hold on
            %                plot(mean(FiltData(i).HTrigFilt{1, 1}{1, i}(exmp,:) ,1))
            %
            %                 calculate selectivity index for each trial and plot
            %                 distribution of that
            %                  for trl=1:Ntrl
            %                      [SIm(trl),SIpRad(trl)]=CalSynchronizationIndex(FiltData(i).HTrigFilt{2}{i}(trl,:),FiltData(i).HTrigFilt{1}{i}(trl,:),AnalysisOpts.MotifAnalysis.FsLFP);
            %                      calculate across time
            %                      for Ttrl=1:length(T_Raw_inds)-1
            %                          [SIm_tim(trl,Ttrl),SIpRad_tim(trl,Ttrl)]=CalSynchronizationIndex(FiltData(i).HTrigFilt{2}{i}(trl,T_Raw_inds(Ttrl):T_Raw_inds(Ttrl+1)),...
            %                              FiltData(i).HTrigFilt{1}{i}(trl,T_Raw_inds(Ttrl):T_Raw_inds(Ttrl+1)),AnalysisOpts.MotifAnalysis.FsLFP);
            %                      end
            %                  end
            %                  H_fig{7}=subplot(NRow,NCol,7);hold on
            %                  obj.FigParams.PlotMeanStd(T_Raw_Tim,SIm_tim,[],'','','b',1,'')
            %                   xlabel('Time')
            %                  ylabel('SI index')
            %                  title('SI index across trials')
            %
            %                  H_fig{8}=subplot(NRow,NCol,8);hold on
            %                   obj.FigParams.PlotMeanStd(T_Raw_Tim,SIpRad_tim,[],'','','r',1,'')
            %                  xlabel('Time')
            %                  ylabel('Radians')
            %                  title('Phase of SI index across trials')
            %
            %                 if 0
            %                 PlotMeanStd(T_Raw,HTrigData.Raw{ValidData(i)},[],'Time','V','b',1)
            %                 title('Avg LFP');axis square
            %                 obj.FigParams.FormatAxes(H_fig{1});
            %                 % show PSD
            %                 H_fig{2}=subplot(NRow,NCol,2);
            %                 helperCWTTimeFreqPlot(mean(HTrigData.PSD{ValidData(i)},3),T_PSD,f,'justplot1',...
            %                     ['Avg PSD n=' num2str(size(HTrigData.PSD{ValidData(i)},3))],...
            %                     'Time(s)','Frequency',obj.LogScale); axis square
            %                 obj.FigParams.FormatAxes(H_fig{2});
            %                 % Show corresponding W
            %                 H_fig{3}=subplot(NRow,NCol,3);
            %                 helperCWTTimeFreqPlot(squeeze(W(:,ValidData(i),:)),T_W,f,'justplot1',['W'],...
            %                     'Time(s)','Frequency',obj.LogScale); axis square
            %                 obj.FigParams.FormatAxes(H_fig{3});
            %
            %                     % plot PSTH
            %                     calculate PSTH for all the neurons of this recording
            %                     [PSTH,Raster,SpkCount] = cellfun(@(x) obj.NeuAnaFunc.CalRasterPSTH(x,HTrigData.Time{ValidData(i)}),SpkData.ts,'UniformOutput',0);
            %                     find this channel and plot it only for this
            %                     NNeu=length(SpkData.Ch);
            %                     Col=distinguishable_colors(NNeu);
            %                     ThisChInd=find(SpkData.Ch==obj.ChNum);
            %                     [CatPSTH,CatRaster,Color]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
            %                     H_fig{4}=subplot(NRow,NCol,4);
            %                     arrayfun(@(x) obj.NeuAnaFunc.PlotPSTH(PSTH{x},Col(x,:),'PSTH for this CH',1,H_fig{4}),ThisChInd,'uniformoutput',0) ;
            %                     axis square
            %                     % plot raster
            %                     H_fig{5}=subplot(NRow,NCol,5);
            %                     obj.NeuAnaFunc.PlotRaster(CatRaster,Color,'Raster This Ch',H_fig{5});axis square
            %                     % plot PSTH for All neurons
            %                     H_fig{6}=subplot(NRow,NCol,6);
            %                     ThisChInd=1:NNeu;
            %                     [CatPSTH,~,~]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
            %                     obj.NeuAnaFunc.PlotPSTH(CatPSTH,'b','PSTH All Chs',0,H_fig{6});
            %                     axis square
            %                     % Plot raster for 5 example Chs
            %                                      H_fig{7}=subplot(NRow,NCol,7);
            %                                      ThisChInd=randsample(1:NNeu,5);
            %                                      [~,CatRaster,Color]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),Col(ThisChInd,:));
            %                                      obj.NeuAnaFunc.PlotRaster(CatRaster,Color,'Raster Exmp Ch',H_fig{7});axis square
            %                                      axis square
            %                     % Plot PSTH for PFC Neurons
            %                     H_fig{7}=subplot(NRow,NCol,7);
            %                     ThisChInd=find(SpkData.ChArea==1);
            %                     [CatPSTH,~,~]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
            %                     obj.NeuAnaFunc.PlotPSTH(CatPSTH,'b','PSTH PFC Chs',0,H_fig{7});axis square
            %                     % Plot PSTH for FEF Neurons
            %                     H_fig{8}=subplot(NRow,NCol,8);
            %                     ThisChInd=find(SpkData.ChArea==4);
            %                     [CatPSTH,~,~]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
            %                     obj.NeuAnaFunc.PlotPSTH(CatPSTH,'b','PSTH FEF Chs',0,H_fig{8});axis square
            %                     % Plot PSTH for PFC Neurons
            %                     H_fig{9}=subplot(NRow,NCol,9);
            %                     ThisChInd=find(SpkData.ChArea==5);
            %                     [CatPSTH,~,~]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
            %                     obj.NeuAnaFunc.PlotPSTH(CatPSTH,'b','PSTH LIP Chs',0,H_fig{9});axis square
            %                     if 1
            %                         % Do Network Analysis
            %                         H_fig{10}=subplot(NRow,NCol,10);
            %                         A=cell2mat(SpkCount);Corrmat=corr(A,A);
            %                         Corrmat=triu(Corrmat,1)+[triu(Corrmat,1)]';
            %                         NaNCol=~isnan(Corrmat(:,1));
            %                         Corrmat=Corrmat(NaNCol,NaNCol);
            %                         [M,Q]=community_louvain(Corrmat,1,[],'negative_asym');
            %                         [~,X] = sort(M);
            %                         imagesc(Corrmat(X,X));
            %                         colorbar
            %                         title(['# clust' num2str(length(unique(M)))]);
            %                         % show Network
            %                         H_fig{11}=subplot(NRow,NCol,[11 12]);
            %                         NPFC=length(find(SpkData.ChArea==1));NFEF=length(find(SpkData.ChArea==4));
            %                         NLIP=length(find(SpkData.ChArea==5));
            %                         AreaCol=[repmat([1 0 0],NPFC,1); repmat([0 1 0],NFEF,1); repmat([0 0 1],NLIP,1)];
            %                         AreaCol=AreaCol(NaNCol,:);
            %                         GraphCorrMat_m(Corrmat,AreaCol,0,5)
            %                     end
            %                 end
            %             end
            %            for i=WantedMotifs
            %    CF.CalComodulation(reshape(HTrigData.Raw{i},1,size(HTrigData.Raw{i},1)*size(HTrigData.Raw{i},2)),1000)
            %  CF.CalComodulationCWT(ShifftSigRaw(1,:),FsLFP)
            %end
            
        end
        function  [Th,varargout]=CalDistPeakVals(obj,Data,BinSiz,PerctileTH,varargin)
            % calculates distribution of peak values for data and returns
            % the threshold of PerctileTH percentile of data
            % BinSiz size of th ebins in samples
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            opts.ShowPlot=0; % Shows an example plot
            
            if PerctileTH<1;PerctileTH=PerctileTH*100;end
            Samp=1:size(Data,2);
            Bins=obj.ManData.BinData(BinSiz,[],1,size(Data,2));
            MaxVals=arrayfun(@(x) unique(max(Data(Samp>=Bins(x) & Samp<Bins(x+1)))),1:length(Bins)-1);
            Th = prctile(MaxVals,PerctileTH);
            
            if opts.ShowPlot
                varargout(1)=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]); % create figures
                [N,edges] = hist(MaxVals);
                
                subplot(121)
                obj.FigParams.BarPlot(edges,N,'b',['Val'],['Probability'],['Histogram of Max values']);
                hold on;v=axis;plot([Th Th],[v(3) v(4)],'r');
                subplot(122)
                obj.FigParams.Plot(1:length(Data),Data,'b',[],[],[]);
                hold on
                v=axis;
                obj.FigParams.Plot([v(1) v(2)],[Th Th],'r','Sample','Val','Data with Th');
            end
        end
        function  varargout=HTriggeredAnalysis(obj,X_Wave,freq,X_Raw,MotifOnset,SpkData,Lt,varargin) % Plots H triggered PSD and spiking
            global AnalysisOpts AnalysisData
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % Lt time we want to look at before and after motif onset 
            % Prepare varables
            FsImaging=AnalysisOpts.MotifAnalysis.FsImaging;
            FsWave=AnalysisOpts.MotifAnalysis.FsWave;
            FsRaw=AnalysisOpts.MotifAnalysis.FsLFP;
             
            X_Wave=obj.PrepapreCWT(X_Wave,'PowerMethod',AnalysisOpts.MotifAnalysis.PowerMethod,'DownSampleFactor',1,'freq',freq); % prepare CWT data for processing
            %%  Now extract H related signal from Rawdata, PSD and spiking data
            [HTrigData.PSD,Tim_PSD]=obj.ExtractMotifOnsetSig(Lt,MotifOnset,FsImaging,X_Wave,FsWave,3);
            [HTrigData.Raw,Tim_Raw]=obj.ExtractMotifOnsetSig(Lt,MotifOnset,FsImaging,X_Raw,FsRaw,2);

           % HTrigData.Time=obj.ExtractHrelatedSig(Lt,Th,H,Time,2);  % chunk the times we have
            
            ValidData=find(arrayfun(@(x) ~isempty(HTrigData.PSD{x}),1:length(HTrigData.PSD)));
            T_PSD=Tim_PSD/FsWave;
            T_Raw=Tim_Raw/FsRaw;
         %   T_W=(1:Lt)/AnalysisOpts.MotifAnalysis.FsWave;
            
         NRow=3;NCol=5;

         %% show raw data
         %             for i=1:length(ValidData)
         %                 varargout(i)=obj.FigParams.RenderFigure(1,[]);
         %                 H_fig{1}=subplot(NRow,NCol,1);
         %                 PlotMeanStd(T_Raw,HTrigData.Raw{ValidData(i)},[],'Time','V','b',1)
         %                 title('Avg LFP');axis square
         %                 obj.FigParams.FormatAxes(H_fig{1});
         %                 %% show PSD
         %                 H_fig{2}=subplot(NRow,NCol,2);
         %                 helperCWTTimeFreqPlot(mean(HTrigData.PSD{ValidData(i)},3),T_PSD,freq,'justplot1',...
         %                     ['Avg PSD n=' num2str(size(HTrigData.PSD{ValidData(i)},3))],...
         %                     'Time(s)','Frequency',obj.LogScale); axis square
         %                 obj.FigParams.FormatAxes(H_fig{2});
         %                 %% Show corresponding W
         %                 H_fig{3}=subplot(NRow,NCol,3);
         %                 helperCWTTimeFreqPlot(squeeze(W(:,ValidData(i),:)),T_W,f,'justplot1',['W'],...
         %                     'Time(s)','Frequency',obj.LogScale); axis square
         %                 obj.FigParams.FormatAxes(H_fig{3});
         %                 if 1
         %                     %% plot PSTH
         %                     % calculate PSTH for all the neurons of this recording
         %                     [PSTH,Raster,SpkCount] = cellfun(@(x) obj.NeuAnaFunc.CalRasterPSTH(x,HTrigData.Time{ValidData(i)}),SpkData.ts,'UniformOutput',0);
         %                     % find this channel and plot it only for this
         %                     NNeu=length(SpkData.Ch);
         %                     Col=distinguishable_colors(NNeu);
         %                     ThisChInd=find(SpkData.Ch==obj.ChNum);
         %                     [CatPSTH,CatRaster,Color]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
         %                     H_fig{4}=subplot(NRow,NCol,4);
         %                     arrayfun(@(x) obj.NeuAnaFunc.PlotPSTH(PSTH{x},Col(x,:),'PSTH for this CH',1,H_fig{4}),ThisChInd,'uniformoutput',0) ;
         %                     axis square
         %                     %% plot raster
         %                     H_fig{5}=subplot(NRow,NCol,5);
         %                     obj.NeuAnaFunc.PlotRaster(CatRaster,Color,'Raster This Ch',H_fig{5});axis square
         %                     %% plot PSTH for All neurons
         %                     H_fig{6}=subplot(NRow,NCol,6);
         %                     ThisChInd=1:NNeu;
         %                     [CatPSTH,~,~]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
         %                     obj.NeuAnaFunc.PlotPSTH(CatPSTH,'b','PSTH All Chs',0,H_fig{6});
         %                     axis square
         %                     %% Plot raster for 5 example Chs
         %                     %                  H_fig{7}=subplot(NRow,NCol,7);
         %                     %                  ThisChInd=randsample(1:NNeu,5);
         %                     %                  [~,CatRaster,Color]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),Col(ThisChInd,:));
         %                     %                  obj.NeuAnaFunc.PlotRaster(CatRaster,Color,'Raster Exmp Ch',H_fig{7});axis square
         %                     %                  axis square
         %                     %% Plot PSTH for PFC Neurons
         %                     H_fig{7}=subplot(NRow,NCol,7);
         %                     ThisChInd=find(SpkData.ChArea==1);
         %                     [CatPSTH,~,~]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
         %                     obj.NeuAnaFunc.PlotPSTH(CatPSTH,'b','PSTH PFC Chs',0,H_fig{7});axis square
         %                     %% Plot PSTH for FEF Neurons
         %                     H_fig{8}=subplot(NRow,NCol,8);
         %                     ThisChInd=find(SpkData.ChArea==4);
         %                     [CatPSTH,~,~]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
         %                     obj.NeuAnaFunc.PlotPSTH(CatPSTH,'b','PSTH FEF Chs',0,H_fig{8});axis square
         %                     %% Plot PSTH for PFC Neurons
         %                     H_fig{9}=subplot(NRow,NCol,9);
         %                     ThisChInd=find(SpkData.ChArea==5);
         %                     [CatPSTH,~,~]=obj.NeuAnaFunc.CatRasterPSTH(PSTH(ThisChInd),Raster(ThisChInd),T_Raw,Col(ThisChInd,:));
         %                     obj.NeuAnaFunc.PlotPSTH(CatPSTH,'b','PSTH LIP Chs',0,H_fig{9});axis square
         %                     if 1
         %                         %% Do Network Analysis
         %                         H_fig{10}=subplot(NRow,NCol,10);
         %                         A=cell2mat(SpkCount);Corrmat=corr(A,A);
         %                         Corrmat=triu(Corrmat,1)+[triu(Corrmat,1)]';
         %                         NaNCol=~isnan(Corrmat(:,1));
         %                         Corrmat=Corrmat(NaNCol,NaNCol);
         %                         [M,Q]=community_louvain(Corrmat,1,[],'negative_asym');
         %                         [~,X] = sort(M);
         %                         imagesc(Corrmat(X,X));
         %                         colorbar
         %                         title(['# clust' num2str(length(unique(M)))]);
         %                         %% show Network
         %                         H_fig{11}=subplot(NRow,NCol,[11 12]);
         %                         NPFC=length(find(SpkData.ChArea==1));NFEF=length(find(SpkData.ChArea==4));
         %                         NLIP=length(find(SpkData.ChArea==5));
         %                         AreaCol=[repmat([1 0 0],NPFC,1); repmat([0 1 0],NFEF,1); repmat([0 0 1],NLIP,1)];
         %                         AreaCol=AreaCol(NaNCol,:);
         %                         GraphCorrMat_m(Corrmat,AreaCol,0,5)
         %                     end
         %                 end
         %             end
         %
         varargout=obj.FigParams.RenderFigure(1,[]);
         FreqLimit=[20 80];
         IndF=(freq>=FreqLimit(1) & freq<=FreqLimit(2));
         for i=1:length(ValidData)
             H_fig{1}=subplot(NRow,NCol,i);
             PSD=mean(HTrigData.PSD{ValidData(i)},3);
             if i==1
                 helperCWTTimeFreqPlot(PSD(IndF,:),T_PSD,freq(IndF),'justplot1',...
                     ['Ch' num2str(AnalysisOpts.CurrentCh) 'Pr' num2str(AnalysisOpts.ProbeNum) 'Dep' num2str(AnalysisOpts.CurrentChDepth) 'Mot' num2str(i) 'n=' num2str(size(HTrigData.PSD{ValidData(i)},3))],...
                     'Time(s)','Frequency',obj.LogScale); axis square
             else
                 helperCWTTimeFreqPlot(PSD(IndF,:),T_PSD,freq(IndF),'justplot1',...
                     ['Mot' num2str(i) 'n=' num2str(size(HTrigData.PSD{ValidData(i)},3))],...
                     'Time(s)','Frequency',obj.LogScale); axis square
             end
             obj.FigParams.FormatAxes(H_fig{1});
         end
        end
        function [out,Tim,HTrigInd,NEvent]=ExtractMotifOnsetSig(obj,Lt,MotifOnset,FsImaging,Data,FsData,Dim,varargin) 
            % extract the signal in Data that is triggered after the onset
            % of the motif for the duration of Lt before and after the
            % onset 
            FsFact=(FsData/FsImaging);% what is the sampling difference between Imaging and Data
            LtImaging=Lt*FsImaging;LtData=Lt*FsData; % convert Lt into samples 
            opts.IncludePreEvent=1; % includes the data pre event we want to look at
            
            
            if ~isempty(varargin) % ind defines which motif we want to look at
                ind=varargin{1};
            else
                ind=1:length(MotifOnset);
            end
            
            for i=ind  % LOOP on different Hs
               
                H_trig_tim=MotifOnset{i};
                Ntim=H_trig_tim(end); 
                H_trig_tim=H_trig_tim(H_trig_tim<(Ntim-2*LtImaging) & H_trig_tim>(2*LtImaging));
                if ~isempty(H_trig_tim)
                    if  opts.IncludePreEvent==1
                        HTrigData=arrayfun(@(x) Data(:,floor((x-LtImaging)*FsFact):floor((x+LtImaging)*FsFact-1)),H_trig_tim,'UniformOutput',0); % take the data
%                     elseif opts.IncludePreEvent==0
%                         HTrigData=arrayfun(@(x) Data(:,x*FsFact:(x+Lt)*FsFact-1),H_trig_tim,'UniformOutput',0);
%                     elseif opts.IncludePreEvent==2
%                         HTrigData=arrayfun(@(x) Data(:,(x-Lt)*FsFact:(x+2*Lt)*FsFact-1),H_trig_tim,'UniformOutput',0); % take the data
                    end
                    HTrigInd{i}=arrayfun(@(x) x*FsFact,H_trig_tim,'UniformOutput',1); % index of H trigger in time
                    out{i}=obj.ManData.ReshapeCell2Mat(HTrigData,Dim);
                end
                NEvent(i)=length(HTrigData);
            end
            if  opts.IncludePreEvent==1
                Tim=-floor(LtImaging*FsFact):floor((LtImaging*FsFact)-1);
%             elseif opts.IncludePreEvent==0
%                 Tim=0:(Lt)*FsFact-1;
%             elseif opts.IncludePreEvent==2
%                 Tim=-(Lt*FsFact):(2*Lt*FsFact)-1;
            end
            
        end

        function [out,Tim,HTrigInd,NEvent]=ExtractHrelatedSig(obj,Lt,Th,H,Data,Dim,varargin)  % extract the H triggered signals
            opts.IncludePreEvent=2; % includes the data pre event we want to look at
            
            Ntim=size(H,2); % Nsampels
            
            if ~isempty(varargin) % Hind defines which H we want to look at
                Hind=varargin{1};
            else
                Hind=1:size(H,1);
            end
            % Th  threshould for H to be triggred , data can be PSD, LFP
            % or spiking.
            FsFact=ceil(size(Data,2)/size(H,2)); % what is the sampling difference between H and Data
            for i=Hind  % LOOP on different Hs
                %Hs=mapminmax(H(i,:),0,1);
                Hs_Th=H(i,:)>=Th(i); % restrict H to values that are bigger than Th;
                
                H_trig_tim=find(diff(Hs_Th)>0)+1;  % find times that H is triggered
                H_trig_tim=H_trig_tim(H_trig_tim<(Ntim-2*Lt-5) & H_trig_tim>(2*Lt+5));
                if ~isempty(H_trig_tim)
                    if  opts.IncludePreEvent==1
                        HTrigData=arrayfun(@(x) Data(:,(x-Lt)*FsFact:(x+Lt)*FsFact-1),H_trig_tim,'UniformOutput',0); % take the data
                    elseif opts.IncludePreEvent==0
                        HTrigData=arrayfun(@(x) Data(:,x*FsFact:(x+Lt)*FsFact-1),H_trig_tim,'UniformOutput',0);
                    elseif opts.IncludePreEvent==2
                        HTrigData=arrayfun(@(x) Data(:,(x-Lt)*FsFact:(x+2*Lt)*FsFact-1),H_trig_tim,'UniformOutput',0); % take the data
                    end
                    HTrigInd{i}=arrayfun(@(x) x*FsFact,H_trig_tim,'UniformOutput',1); % index of H trigger in time
                    out{i}=obj.ManData.ReshapeCell2Mat(HTrigData,Dim);
                end
                NEvent(i)=length(HTrigData);
            end
            if  opts.IncludePreEvent==1
                Tim=-(Lt*FsFact):(Lt*FsFact)-1;
            elseif opts.IncludePreEvent==0
                Tim=0:(Lt)*FsFact-1;
            elseif opts.IncludePreEvent==2
                Tim=-(Lt*FsFact):(2*Lt*FsFact)-1;
            end
            
        end
        function [W,H]=DelRedunWH(~,W,H,loadings) % removes W and H with 0 loadings
            %             NonZeroLoadings=loadings>0;
            %             H=H(NonZeroLoadings,:);
            %             W=W(:,NonZeroLoadings,:);
            SumW=arrayfun(@(x) squeeze(sum(sum(W(:,x,:)))),1:size(W,2));
            NonZeroLoadings=SumW>0;
            H=H(NonZeroLoadings,:);
            W=W(:,NonZeroLoadings,:);
        end
        
        function  RuleTrls=FindRuleTrials(obj,Rule,BlockSpec,NTrl,NMaxTrl)  % finds trials from beg and end of rule blocks
            RuleInd=find(BlockSpec.Rule==Rule & BlockSpec.ThisBlkTrialNum>=NMaxTrl);
            Trls2Look=[1:NMaxTrl-NTrl+1;NTrl:NMaxTrl];
            for l=1:length(Trls2Look)
                RuleTrls.Forward{l}=cell2mat(arrayfun(@(x) BlockSpec.ThisBlkTrials{x}(Trls2Look(1,l):Trls2Look(2,l)),RuleInd,'UniformOutput',0));
                RuleTrls.Bakward{l}=cell2mat(arrayfun(@(x) BlockSpec.ThisBlkTrials{x}(end-Trls2Look(2,l):end-Trls2Look(1,l)),RuleInd,'UniformOutput',0));
                RuleTrls.ForwardVert{l}=(arrayfun(@(x) BlockSpec.ThisBlkTrials{x}(Trls2Look(1,l):Trls2Look(2,l)),RuleInd,'UniformOutput',0));
                RuleTrls.BakwardVert{l}=(arrayfun(@(x) BlockSpec.ThisBlkTrials{x}(end-Trls2Look(2,l):end-Trls2Look(1,l)),RuleInd,'UniformOutput',0));
            end
        end
        function  [MotifPEV_ind,MotifPEV_all]=CalMotifPEV(obj,X,W,H,varargin) %  Cal explained variance for ind motifs
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            % X is the original signal
            %loop on Ws
            if ~isreal(X)
                X=obj.PrepapreCWT(X); % prepare CWT data for processing
            end
            NW=size(W,2);
            for i=1:NW
                Rec_X=helper.reconstruct(W(:,i,:),H(i,:));
                MotifPEV_ind(i,1)=obj.PEV(X,Rec_X);
            end
            % calculte PEV for all
            Rec_X=helper.reconstruct(W,H);
            MotifPEV_all=obj.PEV(X,Rec_X);
        end
        
        function  out=PEV(obj,X,Y)  % calculate PEVe
            if strcmpi(obj.PEVmethod,'var')
                % X  oiginal and Y model estimated
                residuals=X-Y;
                out = 1-nanvar(residuals(:))./nanvar(X(:));
            elseif strcmpi(obj.PEVmethod,'pls') % uses regression
                [XL,yl,~,~,beta,out] = plsregress(Y',X',1);
                out=100*out(2);
            end
        end
        function  varargout=PlotMotifPEV(obj,X,MotifData,varargin)  % plot collective explained variance
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            NStep=length(X);
            PEV=arrayfun(@(x) obj.CalMotifPEV(X{x},MotifData{x}.W,MotifData{x}.H),1:NStep,'UniformOutput',0);
            SortPEV=cell2mat(cellfun(@(x) sort(x,'ascend'),PEV,'UniformOutput',0));
            CumSum_SortPEV=cumsum(SortPEV,1);
            varargout{1}=gcf;
            plot(CumSum_SortPEV);hold on
            plot(nanmean(CumSum_SortPEV,2),'k','linewidth',5)
            xlabel('Motif ordered by PEV')
            ylabel('Cumulative sum of PEV')
            title(['Motfis PEV for CH' num2str(obj.ChNum)])
            
        end
        function [idx, crit,varargout] = InternallyValidateWs(obj,X,MotifSpecs,Ws,Hs,criterion,hp,varargin)
             global AnalysisOpts
             obj=obj.ParseParams(varargin) ; %%Process optional inputs
           
            %Camden MacDowell -timeless
            %Ws and Hs are cell arrays from multiple fits
            if ~isreal(X)
                X=obj.PrepapreCWT(X,'PowerMethod',obj.PowerMethod,'DownSampleFactor',obj.DownSampleFactor); % prepare CWT data for processing
            end
            
            N = numel(MotifSpecs);
            crit = zeros(1,N);
            for i = 1:N   %loop through fits
                if ~isempty(MotifSpecs)
                   ThisMotifs=obj.RemoveNoloadingMotifs(MotifSpecs{i});
                   W=ThisMotifs.W;
                   H=ThisMotifs.H;
                else
                    W=Ws{i};
                    H=Hs{i};
                end
                
                k = size(W,2);
                
                switch criterion
                    case 'AIC' %Residuals are normally distributed so we can approximate AIC=nlog(?^2Z)+2k,
                        Xhat = obj.reconstruct(W,H);
                        r_var = nanvar(X(:)-Xhat(:));
                        crit(i) = numel(X(:))*log(r_var)+(2*k);
                        targetfun = @min;
                        %         figure; histogram(X(:)-Xhat(:));
                    case 'BIC'
                        Xhat = obj.reconstruct(W,H);
                        r_var = nanvar(X(:)-Xhat(:));
                        crit(i) = k*log(numel(X)) - (2*log(r_var));
                        targetfun = @min;
                    case 'PEV'
                        Xhat = obj.reconstruct(W,H);
                        crit(i) = obj.PEV(X,Xhat);
                        targetfun = @max;
                    case 'DISS' % calculate dissimilarity between fits 
                        OtherMotfs=setdiff(1:N,i);
                        crit(i)=mean(arrayfun(@(x) helper.DISSX(MotifSpecs{i}.H,MotifSpecs{i}.W,MotifSpecs{x}.H,MotifSpecs{x}.W),OtherMotfs));
                        targetfun = @max;
                    case 'COST'
                        [crit(i),~,~] = helper.get_seqNMF_cost(X,MotifSpecs{i}.W,MotifSpecs{i}.H);
                        targetfun = @min;
                    case 'NMOTIF' % number of motifs
                        ThisMotifSpec=obj.RemoveNoloadingMotifs(MotifSpecs{i});
                        crit(i)=size(ThisMotifSpec.H,1);
                        targetfun = @min;
                    otherwise
                        error('unknown crossvalidation criterion');
                end %switch                
            end %fit loop
            
            [~, idx] = targetfun(crit); %get the index with the best value
            
            %get the best fit (lowest
            if obj.verbose
                if isempty(hp)
                    varargout(1)=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]); % create figures
                else
                    subplot(hp)
                end
                hold on; plot(crit,'linewidth',2,'color',[0.5 0.5 0.5],'marker','o');
                plot(idx,crit(idx),'marker','x','linestyle','none','color','r','markersize',20,'linewidth',2);
                ylabel(sprintf('%s (AU)',criterion)); xlabel('fit iteration');
            end
            
        end %function
        function varargout=CharctrizeValidateWs(obj,X,MotifSpecs,varargin) % plots information about repetitions of motifs and compares them
            % shows some examples of repetitions 
             global AnalysisOpts
             obj=obj.ParseParams(varargin) ; %%Process optional inputs
             % define options
             opts.CriteriaList={'AIC','BIC','PEV','COST','NMOTIF'};
             % calculate vars
             N=size(MotifSpecs,2); % number of repetitions
             Nepoc=size(MotifSpecs,1); % number of epocs
             Nexm=3; % number of example epocs 
             Ncol_cri=numel(opts.CriteriaList);Nrow_cri=Nexm+3; % show criteria resutls
             Ncol_exm=7;Nrow_exm=N; % show one example of these repetitions 
             exmInd=randsample(Nepoc,Nexm)';
             
             % calculate InternallyValidateWs first 
             for cri=opts.CriteriaList
                 [idx.(cri{1}), critrion.(cri{1})] = arrayfun(@(x) obj.InternallyValidateWs(X{x},MotifSpecs(x,:),[],[],cri{1},[],...
                     'verbose',0,'PowerMethod',AnalysisOpts.MotifAnalysis.PowerMethod,'DownSampleFactor',1),1:Nepoc,'UniformOutput',0);
             end
             %% Plot results 
             % plot charactristics first
             varargout(1)=obj.FigParams.RenderFigure(1,[0 0 0.8 0.8]); % create figures             
             for i=1:Nrow_cri
                 for ncri=1:Ncol_cri
                     hp=subplot(Nrow_cri,Ncol_cri,(i-1)*Ncol_cri+ncri);
                     thiscrit=opts.CriteriaList{ncri};
                     if i<=Nexm % plot examples
                         obj.InternallyValidateWs(X{i},MotifSpecs(i,:),[],[],thiscrit,hp,'verbose',1,...
                             'PowerMethod',AnalysisOpts.MotifAnalysis.PowerMethod,'DownSampleFactor',1);
                     elseif i==Nexm+1  % plot all of them
                         [critrionsort.(thiscrit),critrionsortind.(thiscrit)]=sort(obj.ManData.ReshapeCell2Mat(critrion.(thiscrit),2),2);
                         obj.FigParams.PlotMeanStd(1:N,critrionsort.(thiscrit),[],'fit iteration',sprintf('%s (AU)',thiscrit),...
                             1,1,sprintf('Avg %s',thiscrit)) 
                     elseif i==Nexm+2 % plot differences
                         crit_rs=obj.ManData.ReshapeCell2Mat(critrion.(thiscrit),2);
                         MinMaxdiff=max(crit_rs,[],2)-min(crit_rs,[],2);
                         obj.FigParams.PlotMeanStd(1,MinMaxdiff,[],'',sprintf( '%s (AU)',thiscrit),...
                             1,2,sprintf('Max-Min %s',thiscrit)) 
                     elseif i==Nexm+3 % plot std
                         crit_rs=obj.ManData.ReshapeCell2Mat(critrion.(thiscrit),2);
                         Var=var(crit_rs,1,2);
                         obj.FigParams.PlotMeanStd(1,Var,[],'',sprintf( '%s (AU)',thiscrit),...
                             1,2,sprintf('Var %s',thiscrit))                          
                     end
                 end
             end
             
             % plot examples of motifs for runs
             NReps=5;MaxMotifs=7;
             ExmMtfFigs=arrayfun(@(x) obj.PlotValidateMotifs(MotifSpecs(x,critrionsortind.PEV(x,end-NReps+1:end)),MaxMotifs,NReps,'NewFig',1),exmInd,'UniformOutput',0);
        end
        function varargout=PlotValidateMotifs(obj,MotifSpecs,MaxMotifs,NReps,varargin) % plots motifs from number of iterations 
             global AnalysisOpts AnalysisData
             obj=obj.ParseParams(varargin) ; %%Process optional inputs
            
            opts.MaxMotifs=MaxMotifs; % maximum number of motifs shown for each repetition
            opts.NormalizeMotif=1; % should we normalize the motifs
            %MotifSpecs is cell array where each repetition is stored in a
            %cell 
            N=NReps;%numel(MotifSpecs);
            FsWaveTarg=AnalysisOpts.MotifAnalysis.FsWaveTarg;
            Ncol=opts.MaxMotifs;Nraw=N; 
            if  obj.NewFig
                obj.FigParams.RenderFigure(1,[]);
            end
            for n=1:N % loop on each
                ThisMotifSpec=obj.RemoveNoloadingMotifs(MotifSpecs{n});
                ThisMotifs=ThisMotifSpec.W;
                
                if ~iscell(ThisMotifs) % if it is a 3 dim matrix then turn it into a cell
                    ThisMotifs=arrayfun(@(x) squeeze(ThisMotifs(:,x,:)),1:size(ThisMotifs,2),'UniformOutput',0);
                end
                if opts.NormalizeMotif
                    ThisMotifs=obj.ManData.NormalizeMotifs(ThisMotifs);
                end                
                NWs=length(ThisMotifs);
                NTim=size(ThisMotifs{1},2);
                               
                for i=1:min(opts.MaxMotifs,NWs)                  
                    subplot(Nraw,Ncol,i+(n-1)*Ncol)                     
                    hold on                     
                    Title=sprintf('Rep:%i Mtf:%i',n,i);                                        
                   % helperCWTTimeFreqPlot(ThisMotifs{i},(1:NTim)/FsWaveTarg,AnalysisData.cwt_f,'image',Title,'Time(s)','Freq(Hz)',obj.LogScale);                   
                    axis square
                    helperCWTTimeFreqPlot(ThisMotifs{i},(1:NTim)/FsWaveTarg,AnalysisData.cwt_f,'image',Title,'','',obj.LogScale);                                    
               %    axis off
                end
            end
            varargout{1}=gcf;
        end
        
        function X_hat = reconstruct(~,W,H)
            % ------------------------------------------------------------------------
            % USAGE: X_hat = helper.reconstruct(W,H)
            % ------------------------------------------------------------------------
            % INPUTS
            % W:      W is a NxKxL tensor which gives the neuron basis
            %         functions which are used for the reconstructions. The L'th NxK slice
            %         of W is the neural basis set for a lag of L.
            %
            % H:      H is a KxT matrix which gives timecourses for each factor
            % ------------------------------------------------------------------------
            % OUTPUTS
            % X_hat:  The reconstruction X_hat = W (*) H;
            % ------------------------------------------------------------------------
            % Emily Mackevicius and Andrew Bahle
            
            [N,K,L] = size(W);
            [~,T] = size(H);
            
            % zeropad by L
            H = [zeros(K,L),H,zeros(K,L)];
            T = T+2*L;
            X_hat = zeros(N,T);
            
            for tau = 1:L % go through every offset from 1:tau
                %X_hat = X_hat + W(:, :, tau) * circshift(H,tau-1,2);
                X_hat = X_hat + W(:, :, tau) * circshift(H,[0,tau-1]);
            end
            
            % undo zer0padding
            X_hat = X_hat(:,(L+1):(end-L));
            
        end
        function [varargout]= ShowRawDataTF(obj,RawData,CWTdata,f_cwt,TrialTimes,FsRaw,FsWave,FsWaveTarg,TimStp,varargin)  % brings examples of raw data to see if we can find motifs by eye
            global AnalysisOpts
            obj=obj.ParseParams(varargin) ; %%Process optional inputs
            %% Prepare data
            flimit=[0 80]; % limits of Freq we want to look at
            Navg=AnalysisOpts.MotifAnalysis.Navg;
            %%
            flimitind=f_cwt>=flimit(1) & f_cwt<=flimit(2);
            f_cwt=f_cwt(flimitind);
            CWTdata=CWTdata(flimitind,:);
            CWTdata=movmean(CWTdata,obj.Navg,2);
            if size(f_cwt,1)==1;f_cwt=f_cwt';end
            if size(RawData,1)==1;RawData=RawData';end
            CWTdata=transpose(downsample(CWTdata',FsWave/FsWaveTarg));
            FsWave=FsWaveTarg;
            %%
            % filter raw data in different bands
            ffilt=[0.5 50 100 150 200]; % list of frequency we want to look at the signal in
            [FiltRaw,BandSet]=obj.FilterFunc.BandPassFilter(RawData,FsRaw,ffilt); % filter raw data
            % FiltRawDS=cellfun(@(x) downsample(x,Fs/FsLFP),FiltRaw,'UniformOutput',0);
            %calculate power of data
            PwrMethodSet={'Power','NormPower','NormMedianPower','Angle'};
            CwtPower=cellfun(@(x) obj.PrepapreCWT(CWTdata,'PowerMethod',x,'DownSampleFactor',1,'freq',f_cwt),PwrMethodSet,'uniformoutput',0);
            % get some initial information about what we have
            NSampCWT=size(CWTdata,2);NSampRaw=size(RawData,1);
            TimeRaw=0:1/FsRaw:(NSampRaw-1)/FsRaw; % build time vector
            TimeCWT=0:1/FsWave:(NSampCWT-1)/FsWave;
            %% initilize plots
            NRaw=4;Ncol=4;
            varargout{1}=obj.FigParams.RenderFigure(1,[]);
            NTimeSet=length(0:TimStp:TimeRaw(end));
            for i=1:NTimeSet
                %% now plot everything in a plot to look at
                nTS=0;
                for TS=[TimStp TimStp/2]
                    %  TimStp=5; % 1 second time steps to look at
                    TimeSet=0:TimStp:TimeRaw(end);
                    % TimeSet=[TimeSet(1:end-1)' TimeSet(2:end)'];
                    TimeSet=[TimeSet(1:end-1)' (TimeSet(1:end-1)+TS)'];
                    
                    ThisTimeCWT=TimeCWT(TimeCWT>=TimeSet(i,1) & TimeCWT<=TimeSet(i,2));
                    ThisTimeCWTInd=(TimeCWT>=TimeSet(i,1) & TimeCWT<=TimeSet(i,2));
                    
                    ThisTimeRaw=TimeRaw(TimeRaw>=TimeSet(i,1) & TimeRaw<=TimeSet(i,2));
                    ThisTimeRawInd=(TimeRaw>=TimeSet(i,1) & TimeRaw<=TimeSet(i,2));
                    %% start ploting stuff
                    % plot TF plot first
                    for k=1:length(PwrMethodSet)
                        h{nTS+1,k}=subplot(NRaw,Ncol,k+NRaw*(Ncol/2)*nTS);cla
                        Pow=CwtPower{k}(:,ThisTimeCWTInd);
                        helperCWTTimeFreqPlot(Pow ,ThisTimeCWT,f_cwt,'justplot1',['CWT CH ' num2str(obj.ChNum) PwrMethodSet{k}],...
                            'Time(s)','Frequency',obj.LogScale);
                    end
                    %% now plot the raw data
                    for j=1:size(BandSet,1)
                        % plot the raw data now
                        h{nTS+1,j+k}=subplot(NRaw,Ncol,j+k+NRaw*(Ncol/2)*nTS);
                        X=FiltRaw(j,ThisTimeRawInd);
                        plot(ThisTimeRaw,X,'k','LineWidth',1)
                        xlabel('Time(s)');ylabel('volts')
                        title([num2str(BandSet(j,1)) '-' num2str(BandSet(j,2))])
                        axis tight
                    end
                    % take the trail times and stamp fix,stim and reward time
                    cellfun(@(x) obj.PlotTrialTimes(x,TrialTimes,1,ThisTimeRaw),h(nTS+1,:),'UniformOutput',0);
                    nTS=nTS+1;
                end
                pause
            end
        end
    end
end

