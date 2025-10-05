% process all of the data of the classifier with different methods
function [Hsum,TargetSps_sum,NStgTxt,Ext,RawDataFig]=RunAllStattestsonData(ClassifierOpts,Metric,AvgDim2,Dim,WantedNStg,CorrMethodSet,Hsum,TargetSps_sum)
global AnalysisOpts

warning off
%close all
NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
fp=fig_params;
Prmt=@(x) permute(x,[1 3 2]); % permute
Time=-0.6:0.01:0.55;
opts.PlotSummaryStatSpecs=0; % summary plot for STD pval ....
opts.PlotSummaryStatCorrMethod=1;% summery plot for each correlation method e.g. Kendall Spearman

%TimeIndSet=fp.LimitTimeAxis4Classifier(Time,ClassifierOpts);
%TimeInd=TimeIndSet{Dim};

 %ClassifierOpts=NeuAna.DefineClassifierTestOptions(TaskName,'CalShuffTrlOrder',1);
%% load data
% the data should be already loaded and passed along to this
%% preprocess the data
DimNum=[0 1 3];
[~,TargetFactorDim]=arrayfun(@(dim) NeuAna.getClassifierDimInfo(ClassifierOpts,dim),[0 1 3],'UniformOutput',0);
FileNameSpecs=@(x) [ClassifierOpts.Name '_' TargetFactorDim{x}];
ThisFileNameSpecs=FileNameSpecs(Dim);

PATH=AnalysisOpts.CodeTestDataPath;

LoadFileName=[ThisFileNameSpecs '_D' num2str(Dim(1)) 'ShData'];
% prepare data
Data=load([PATH LoadFileName],[Metric '_Obsv'],[Metric '_Shuff']);
if ~(isfield(Data,[Metric '_Obsv']) & isfield(Data,[Metric '_Shuff']))
    fprintf('\nWe dont have the data to process this condition')
    return
end
Observed=ManData.SmoothData(Data.([Metric '_Obsv']){1},NeuAna.WidthSmoothing,'SmoothingMethod',NeuAna.SmoothingMethod,'DimSmoothing',2);
Shuffle=ManData.SmoothData(Data.([Metric '_Shuff']){1},NeuAna.WidthSmoothing,'SmoothingMethod',NeuAna.SmoothingMethod,'DimSmoothing',2);
NTim=size(Observed,2);

if strcmp(Metric,'Accuracy')
    MetricName=TargetFactorDim{Dim};
else  
    ThisFileNameSpecs=  [ClassifierOpts.Name '_' Metric];
    MetricName=Metric;
end
fprintf('\nLoading data of classifier results from %s',LoadFileName)

if WantedNStg==0 | size(Observed,1)<WantedNStg
    TrlRng=[1:size(Observed,1)]';NStg=size(Observed,1);NStgTxt='';
else
    TrlRng=[1:WantedNStg]';NStg=WantedNStg;NStgTxt=['Nstg' num2str(NStg)];
end

% define averaging time period 
if strcmp(MetricName,'ShapeCat') | strcmp(MetricName,'ColorCat') | ...
        strcmp(MetricName,'ColorAvgDist') | strcmp(MetricName,'ShapeAvgDist') | strcmp(MetricName,'CPI')
    AvgTimePeriod=[0.1 0.3];
elseif strcmp(MetricName,'ResponseLoc')
    AvgTimePeriod=[0.2 0.4];
elseif strcmp(MetricName,'Rule')
    AvgTimePeriod=[-0.4 0];
end
TimeIndAvg=Time>=AvgTimePeriod(1) & Time<=AvgTimePeriod(2);
Observed_TimAvg=mean(Data.([Metric '_Obsv']){1}(:,TimeIndAvg),2);
Shuffle_TimAvg=mean(Data.([Metric '_Shuff']){1}(:,TimeIndAvg,:),2);


% if we don't want the default then define it 
if AvgDim2>0 & AvgDim2~=1
    NeuAna.WidthSmoothingDim2=AvgDim2;
end
% avergae across dim2?
if AvgDim2
    Observed=ManData.SmoothData(Observed,NeuAna.WidthSmoothingDim2,'SmoothingMethod','movmean','DimSmoothing',1);
    Shuffle=ManData.SmoothData(Shuffle,NeuAna.WidthSmoothingDim2,'SmoothingMethod','movmean','DimSmoothing',1);
    % smooth also average data in Dim 2 
    Observed_TimAvg=ManData.SmoothData(Observed_TimAvg,NeuAna.WidthSmoothingDim2,'SmoothingMethod','movmean','DimSmoothing',1);
    Shuffle_TimAvg=ManData.SmoothData(Shuffle_TimAvg,NeuAna.WidthSmoothingDim2,'SmoothingMethod','movmean','DimSmoothing',1);   
end
Shuffle=permute(Shuffle,[3 1 2]);Shuffle_TimAvg=permute(Shuffle_TimAvg,[3 1 2]);
NrepShuff=size(Shuffle,1);

%% start plotting the data
% plot observed data and shuffle
RepSpace=@(x) strrep(x,'_',' ');
col=copper(NStg);

%CorrMethodSet={'Kendall','Spearman','3PW_p','3PW_slop','TheilSen','Regression','Sigmoid','ARIMAlag2','ARIMAlag4','ARIMAlag6'};
if isempty(CorrMethodSet)
    CorrMethodSet={'Kendall','Spearman','TheilSen','Regression','Sigmoid','3PW_p','3PW_slop','TFPW_Y_p','TFPW_Y_slop',...
        'TFPW_WS_p','TFPW_WS_slop','PW_p','PW_slop'};
end

nMethods=length(CorrMethodSet);

for nThisCorrMethod=1:nMethods
    CorrMethod=CorrMethodSet{nThisCorrMethod};
    PlotData=1;
    fprintf('\n calculating trend correlation using %s ...',CorrMethod)
    %% loop on different methods and calculate statistics
    if strcmp(CorrMethod,'Kendall')
        CorrObsv=arrayfun(@(x) corr(TrlRng,Observed(1:NStg,x),'type','Kendall'),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(Shuffle(rep,1:NStg,x))','type','Kendall'),1:NTim)',1:NrepShuff,'UniformOutput',0))';
         % do the same for avergae time 
        CorrObsv_TimAvg=arrayfun(@(x) corr(TrlRng,Observed_TimAvg(1:NStg,x),'type','Kendall'),1);
        CorrShuff_TimAvg=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(Shuffle_TimAvg(rep,1:NStg,x))','type','Kendall'),1)',1:NrepShuff,'UniformOutput',0))';
    elseif strcmp(CorrMethod,'GPTKendall')
        CorrObsv=arrayfun(@(x) mannKendallTestAC(Observed(1:NStg,x),0.05),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) mannKendallTestAC(squeeze(Shuffle(rep,1:NStg,x)),0.05),1:NTim)',1:NrepShuff,'UniformOutput',0))';    
    elseif strcmp(CorrMethod,'TheilSen')
        CorrObsv=arrayfun(@(x) TheilSen([TrlRng,Observed(1:NStg,x)]),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) TheilSen([TrlRng,squeeze(Shuffle(rep,1:NStg,x))']),1:NTim)',1:NrepShuff,'UniformOutput',0))';
        % do the same for avergae time   
        CorrObsv_TimAvg=arrayfun(@(x) TheilSen([TrlRng,Observed_TimAvg(1:NStg,x)]),1);
        CorrShuff_TimAvg=cell2mat(arrayfun(@(rep) arrayfun(@(x) TheilSen([TrlRng,squeeze(Shuffle_TimAvg(rep,1:NStg,x))']),1)',1:NrepShuff,'UniformOutput',0))';
  
    elseif strcmp(CorrMethod,'Spearman')
        CorrObsv=arrayfun(@(x) corr(TrlRng,Observed(1:NStg,x),'type','Spearman'),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(Shuffle(rep,1:NStg,x))','type','Spearman'),1:NTim)',1:NrepShuff,'UniformOutput',0))';
        % do the same for avergae time
        CorrObsv_TimAvg=arrayfun(@(x) corr(TrlRng,Observed_TimAvg(1:NStg,x),'type','Spearman'),1);
        CorrShuff_TimAvg=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(Shuffle_TimAvg(rep,1:NStg,x))','type','Spearman'),1)',1:NrepShuff,'UniformOutput',0))';
 
    elseif strcmp(CorrMethod,'Regression')
        CorrObsv=arrayfun(@(x) LinearRegression(TrlRng,Observed(1:NStg,x)),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) LinearRegression(TrlRng,squeeze(Shuffle(rep,1:NStg,x))),1:NTim)',1:NrepShuff,'UniformOutput',0))';
        % do the same for avergae time
        CorrObsv_TimAvg=arrayfun(@(x) LinearRegression(TrlRng,Observed_TimAvg(1:NStg,x)),1);
        CorrShuff_TimAvg=cell2mat(arrayfun(@(rep) arrayfun(@(x) LinearRegression(TrlRng,squeeze(Shuffle_TimAvg(rep,1:NStg,x))),1)',1:NrepShuff,'UniformOutput',0))';

    elseif strcmp(CorrMethod,'MMK')
        significance_value_tau = AnalysisOpts.Modified_MannKendall_significance_value_tau;
        significance_value_ac = AnalysisOpts.Modified_MannKendall_significance_value_ac;
        CorrObsv=arrayfun(@(x) Modified_MannKendall_test(TrlRng, Observed(1:NStg,x)', significance_value_tau, significance_value_ac),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) Modified_MannKendall_test(TrlRng, squeeze(Shuffle(rep,1:NStg,x))', significance_value_tau, significance_value_ac),1:NTim)',1:NrepShuff,'UniformOutput',0))';
    elseif strcmp(CorrMethod,'Sigmoid')
        [SlopeObsv,~,~,DifMinMaxObsv,NormDifMinMaxObsv,NormDifMinObsv]=arrayfun(@(x) fitAndAdjustSigmoidLsq(TrlRng,Observed(1:NStg,x),0,num2str(Time(x))),1:NTim,'UniformOutput',0);
        for rep=1:NrepShuff
            [SlopeShuff(rep,:),~,~,DifMinMaxShuff(rep,:),NormDifMinMaxShuff(rep,:),NormDifMinShuff(rep,:)]=arrayfun(@(x) fitAndAdjustSigmoidLsq(TrlRng,squeeze(Shuffle(rep,1:NStg,x))',0,num2str(Time(x))),1:NTim,'UniformOutput',0);
        end
        CorrObsv=cell2mat(DifMinMaxObsv);
        CorrShuff=cell2mat(DifMinMaxShuff);
        % do the same for avergae time
        [~,~,~,DifMinMaxObsv_TimAvg]=arrayfun(@(x) fitAndAdjustSigmoidLsq(TrlRng,Observed_TimAvg(1:NStg,x),0,num2str(Time(x))),1,'UniformOutput',0);
        for rep=1:NrepShuff
            [~,~,~,DifMinMaxShuff_TimAvg(rep,:)]=arrayfun(@(x) fitAndAdjustSigmoidLsq(TrlRng,squeeze(Shuffle_TimAvg(rep,1:NStg,x))',0,num2str(Time(x))),1,'UniformOutput',0);
        end
        CorrObsv_TimAvg=cell2mat(DifMinMaxObsv_TimAvg);
        CorrShuff_TimAvg=cell2mat(DifMinMaxShuff_TimAvg);
    elseif contains(CorrMethod,'ARIMA')
        [OutVar,~,SaveVarName]=FitStatMetric2Data(str2double(CorrMethod(end)),'ARIMA',AvgDim2,NStg,NStgTxt,LoadFileName,Observed,Shuffle,MetricName);
        CorrObsv=OutVar.(SaveVarName).Observed_a;
        CorrShuff=OutVar.(SaveVarName).Shuffle_a;
        % do the same for avergae time
        [OutVar_TimAvg,~,SaveVarName_TimAvg]=FitStatMetric2Data(str2double(CorrMethod(end)),'ARIMA',AvgDim2,NStg,NStgTxt,LoadFileName,Observed_TimAvg,Shuffle_TimAvg,[MetricName '_TimAvg']);
        CorrObsv_TimAvg=OutVar_TimAvg.(SaveVarName_TimAvg).Observed_a;
        CorrShuff_TimAvg=OutVar_TimAvg.(SaveVarName_TimAvg).Shuffle_a;
        
    elseif strcmp(CorrMethod,'3PW_p') | strcmp(CorrMethod,'PW_p') |  strcmp(CorrMethod,'TFPW_Y_p') | strcmp(CorrMethod,'TFPW_WS_p')
        [OutVar,~,SaveVarName]=FitStatMetric2Data(2,CorrMethod,AvgDim2,NStg,NStgTxt,LoadFileName,Observed,Shuffle,MetricName);        
        CorrObsv=OutVar.(SaveVarName).Observed_p;
        CorrShuff=OutVar.(SaveVarName).Shuffle_p;
        % do the same for avergae time       
        [OutVar_TimAvg,~,SaveVarName_TimAvg]=FitStatMetric2Data(2,CorrMethod,AvgDim2,NStg,NStgTxt,LoadFileName,Observed_TimAvg,Shuffle_TimAvg,[MetricName '_TimAvg']);        
        CorrObsv_TimAvg=OutVar_TimAvg.(SaveVarName_TimAvg).Observed_p;
        CorrShuff_TimAvg=OutVar_TimAvg.(SaveVarName_TimAvg).Shuffle_p;
      
    elseif strcmp(CorrMethod,'3PW_slop')  | strcmp(CorrMethod,'PW_slop') |  strcmp(CorrMethod,'TFPW_Y_slop') | strcmp(CorrMethod,'TFPW_WS_slop')
        CorrObsv=OutVar.(SaveVarName).Observed_a;
        CorrShuff=OutVar.(SaveVarName).Shuffle_a;
        CorrObsv_TimAvg=OutVar_TimAvg.(SaveVarName_TimAvg).Observed_a;
        CorrShuff_TimAvg=OutVar_TimAvg.(SaveVarName_TimAvg).Shuffle_a;
    end

    if PlotData      
        % plot with cluster correction
%         if strcmp(TargetFactorDim{Dim},'Rule')
%             TimeInd=Time>=-0.4 & Time<=0;
%         else
%             TimeInd=Time>=-0.2 & Time<=0.6;
%         end
        TimeInd=logical(ones(1,NTim));
        %% plot mean of observed and shuffle
        if nThisCorrMethod==1
            RawDataFig=figure;
            subplot(131);hold on;arrayfun(@(x) plot(Time,Observed(x,1:NTim),'color',col(x,:)),TrlRng');
            xlabel('Time');title('Observed')
            
            MeanShuff=squeeze(mean(Shuffle,1));
            subplot(132);hold on;arrayfun(@(x) plot(Time,MeanShuff(x,1:NTim),'color',col(x,:)),TrlRng');
            xlabel('Time');title('Mean Shuffle')
            % plot average data as well
            subplot(133);hold on;plot(1:NStg,Observed_TimAvg(1:NStg));
            xlabel('Time');title('Mean Shuffle')
            sgtitle(['Shuffle and observed for ' RepSpace(ThisFileNameSpecs)])
        end
        %% plot stat data
        AnalysisOpts.Time=Time(TimeInd);
        ManData.ClusterMassCorrection_permutationTwoTail(Prmt(CorrShuff(:,TimeInd)),Prmt(CorrObsv(1,TimeInd)),0.1,1,'ShowClustCorrectionPlot',1);
        ClusterCorrectionTJB(CorrObsv(1,TimeInd)',CorrShuff(:,TimeInd)',{[],[2,4,7]},'Time',Time(TimeInd),'ThresholdPercentage',15,'DoubleSided',1);
        ClusterCorrectionTJB(CorrObsv(1,TimeInd)',CorrShuff(:,TimeInd)',{[],[2,4,8]},'Time',Time(TimeInd),'ThresholdPercentage',24,'DoubleSided',1);
        subplot(2,4,2);cla
        %ClusterCorrectionTJB(CorrObsv(1,TimeInd)',CorrShuff(:,TimeInd)',{[],[2,4,2]},'Time',Time(TimeInd),'ThresholdPercentage',35,'DoubleSided',1);
        % above line was changed on 09/30/2025
        ClusterCorrectionTJB(CorrObsv(1,TimeInd)',CorrShuff(:,TimeInd)',{[],[2,4,2]},'Time',Time(TimeInd),'ThresholdPercentage',30,'DoubleSided',1);
        

        % replace the subplot for distribution with pvalue for average time 
        subplot(2,4,3);cla
        ManData.PlotPvals(CorrObsv_TimAvg,CorrShuff_TimAvg,0);

        sgtitle({RepSpace(ThisFileNameSpecs);CorrMethod})
        SourceFig=gcf;
        %  fp.FormatAllAxesFig(gcf)

        %% create a summery plot for each of the stat specs for each metric so we can look at all of them together
        if opts.PlotSummaryStatSpecs
            SourceSpSet={[2,4,6],[2,4,4],[2,4,1],[2,4,7],[2,4,8],[2,4,3]}; % referint to subplots inside the figure
            SourceSpSetName={'STD','uncorrpval','Cluster','ClustTJBTHR15','ClusterTJBTHR24','pvalAvgTime'};
            if nThisCorrMethod==1
                h=fp.RenderFigure(length(SourceSpSet),[]);
                [Htot,TargetSps]=cellfun(@(x) fp.RenderSubplots([],[],x,nMethods),h,'UniformOutput',0);
            end

            % loop on figures and methods
            for tsp=1:length(SourceSpSet)
                figure(SourceFig);
                SourceSp=subplot(SourceSpSet{tsp}(1),SourceSpSet{tsp}(2),SourceSpSet{tsp}(3)); % choose the relevant Sp
                fp.CopySubplot(SourceSp,TargetSps{tsp}(nThisCorrMethod)); % copy into target
                fp.AppendFigTitles(TargetSps{tsp}(nThisCorrMethod),{CorrMethod},'AppendTitles',1)
                figure(Htot{tsp})
                sgtitle(RepSpace(ThisFileNameSpecs))
            end
        end
        %% create a summary figure for each of the factors for each correlation method
        if opts.PlotSummaryStatCorrMethod
            SumeryFigMetric={'ShapeCat','ColorCat','Rule','CPI','ShapeAvgDist','ColorAvgDist','ResponseLoc'};
            nSumeryFigMetric=[  1           2        3      4         5               6          1];
            nThisMetric=nSumeryFigMetric(strcmp(SumeryFigMetric,MetricName));
            SumerySourceSpSetName={'STD','ucorpv','TJB15','TJB24','TJB30'};
            NRowSum=length(SumerySourceSpSetName);NColSum=6;
            SumeryFigSPs    ={[NRowSum NColSum  1],[NRowSum NColSum 2],[NRowSum NColSum 3],[NRowSum NColSum 4],[NRowSum NColSum 5],[NRowSum NColSum 6],[NRowSum NColSum 1]};
            SumeryTargetSpSet={[0*NColSum+nThisMetric],[1*NColSum+nThisMetric],[2*NColSum+nThisMetric],[3*NColSum+nThisMetric],[4*NColSum+nThisMetric]}; % referint to subplots inside the figure
            SumerySourceSpSet={[2,4,6],[2,4,4],[2,4,7],[2,4,8],[2,4,2]}; % referint to subplots inside the figure

            if nThisCorrMethod==1 & (strcmp(Metric,'Accuracy') & Dim==1)
                hsum=fp.RenderFigure(nMethods,[]);
                [Hsum,TargetSps_sum]=cellfun(@(x) fp.RenderSubplots(NRowSum,NColSum,x,[]),hsum,'UniformOutput',0);
            end

            % loop on figures and methods
            for tsp=1:length(SumerySourceSpSet)
                figure(SourceFig);
                SourceSp=subplot(SumerySourceSpSet{tsp}(1),SumerySourceSpSet{tsp}(2),SumerySourceSpSet{tsp}(3)); % choose the relevant Sp
                fp.CopySubplot(SourceSp,TargetSps_sum{nThisCorrMethod}(SumeryTargetSpSet{tsp})); % copy into target
                figure(Hsum{nThisCorrMethod})
                subplot(TargetSps_sum{nThisCorrMethod}(SumeryTargetSpSet{tsp}))
                title([MetricName ' ' SumerySourceSpSetName{tsp}])
            end
            figure(Hsum{nThisCorrMethod})
            sgtitle([RepSpace(ThisFileNameSpecs);{CorrMethod}])
        end
    end
end

    if AvgDim2==1 % then we take 3 which is our default
        Ext='AvgDim2';
    elseif AvgDim2==0
        Ext='';
    else
        Ext=['AvgDim2_' num2str(AvgDim2)];
    end
    if opts.PlotSummaryStatSpecs
        % save figures
        fprintf('\nSaving current figures ...')
        fp.SaveCurrentFigs(['AllStats_' NStgTxt Ext ThisFileNameSpecs],[AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep]);

        % save specific figures in one file
        for kk=[1 2 4 5]
            fp.SaveFigSeries(['AllStats_' NStgTxt Ext MetricName '_' SourceSpSetName{kk}],[AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep],Htot{kk});
        end
    end
end

