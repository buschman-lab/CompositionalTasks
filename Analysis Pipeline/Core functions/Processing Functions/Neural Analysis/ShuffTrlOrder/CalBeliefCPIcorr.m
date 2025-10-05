function [a,p,X,Y,a_time]=CalBeliefCPIcorr(ManData,NeuAna,ClassifierOpts,TargetFactorDim1,Scores_2ndD_avg_avg_Flip,CPI,opts,PlotFlag) % calculates correlation between belief and CPI

if contains(TargetFactorDim1,'Rule');warning('first dimension has to be shape or color');return
end
% get belief data 
TimeIndBelief = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(3, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(3, 2);

BeliefEncoding = ManData.ReshapeStruct2Mat(Scores_2ndD_avg_avg_Flip, 'TrialRange', 2);
MeanBelief= mean(BeliefEncoding(:,TimeIndBelief),2);

% get CPI data
CPI=log(CPI);
TimeIndComp = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(1, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(1, 2);
CompIndex = (mean(CPI(:, TimeIndComp), 2));


% Calculate the correlation between mean scores and compression index
X=MeanBelief;
Y=CompIndex;
[a,p] = ManData.Correlation(X,Y,0,opts.CorrelationMetric);

% calculate correlation in time 
CPI=ManData.SmoothData(CPI,NeuAna.WidthSmoothing,'SmoothingMethod','movmean','DimSmoothing',2);
[a_time] = arrayfun(@(x) ManData.Correlation(X,CPI(:,x),0,opts.CorrelationMetric),1:length(ClassifierOpts.AnalysisOpts.Time));


if PlotFlag
    subplot(121)
    scatter(X,Y);
    xlabel('Average belief Encoding')
    ylabel('Average CPI')
    axis square

    subplot(122)
    plot(ClassifierOpts.AnalysisOpts.Time,a_time)
    xlabel('Time')
    ylabel('Correlation')
    title('Correlation between CPI and belief encoding')
    axis square
end


end