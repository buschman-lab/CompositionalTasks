function [a,p,X,Y,a_time]=CalCorrColShapRatioBhvPerf(ManData,NeuAna,ClassifierOpts,TargetFactorDim1,TargetFactorDim2,AllFactorVals_1ndD_avg,AccuracyDataDim1,AccuracyDataDim2,opts,PlotFlag,BhvFactor,NtrlsAvgBhv,NormType) % calculates correlation between belief and CPI

if ~(strcmp(TargetFactorDim1,'ShapeCat') & strcmp(TargetFactorDim2,'ColorCat'))
    warning('first dimension has to be shape and second dimension color');
    a=nan;p=nan;X=nan;Y=nan;return
end

% get sensory data
TimeIndStimAccuracy = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(1, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(1, 2);

 AvgAccuracyShape=(mean(AccuracyDataDim1(:,TimeIndStimAccuracy),2));
 AvgAccuracyColor=(mean(AccuracyDataDim2(:,TimeIndStimAccuracy),2));

 % get accuracy in time 
 AvgAccuracyShape_Time=ManData.SmoothData(AccuracyDataDim1,NeuAna.WidthSmoothing,'SmoothingMethod','movmean','DimSmoothing',2);
 AvgAccuracyColor_Time=ManData.SmoothData(AccuracyDataDim2,NeuAna.WidthSmoothing,'SmoothingMethod','movmean','DimSmoothing',2);

 % Compute the ratio for average and for time
 if contains(NormType,'Raw','IgnoreCase',true)
     RatioAvgAccuShpCol=AvgAccuracyColor./AvgAccuracyShape;
 elseif contains(NormType,'Norm','IgnoreCase',true)
     RatioAvgAccuShpCol=(AvgAccuracyColor-AvgAccuracyShape)./(AvgAccuracyColor+AvgAccuracyShape);
     RatioAvgAccuShpCol_Time=(AvgAccuracyColor_Time-AvgAccuracyShape_Time)./(AvgAccuracyColor_Time+AvgAccuracyShape_Time);
 elseif contains(NormType,'Chance','IgnoreCase',true)
     RatioAvgAccuShpCol=(AvgAccuracyColor-0.5)./(AvgAccuracyShape-0.5);
 end

   
% get the behavioral performance for each trial
nTrialRange=length(AllFactorVals_1ndD_avg);
Bhv_AllTrlRng=arrayfun(@(x) mean(AllFactorVals_1ndD_avg(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
Bhv_AllTrlRng=cell2mat(Bhv_AllTrlRng)';
BhvPerfTxt=['BhvPerf' BhvFactor num2str(NtrlsAvgBhv)];
BhvFactorInd=cellfun(@(x) find(strcmp(ClassifierOpts.AnalysisOpts.FactorInds2Keep,x)),{BhvPerfTxt});
BhvPerfTrl=Bhv_AllTrlRng(:,BhvFactorInd);

% Calculate the correlation between mean scores and compression index
X=RatioAvgAccuShpCol;
Y=BhvPerfTrl;
[a,p] = ManData.Correlation(X,Y,0,opts.CorrelationMetric);

% calculate correlation in time as well
[a_time] = ManData.Correlation(RatioAvgAccuShpCol_Time,Y,0,opts.CorrelationMetric);


if PlotFlag
    subplot(131)
    scatter(X,Y);
    xlabel('Norm Ratio of Shape/Color')
    ylabel(BhvPerfTxt)
    axis square
    subplot(132)
    plot(RatioAvgAccuShpCol)
    xlabel('Trials from switch')
    ylabel('Norm Ratio of Shape/Color')
    axis square
end
end