function [a,p,X,Y,a_time]=CalBeliefBhvPerfcorr(ManData,NeuAna,ClassifierOpts,TargetFactorDim1,Scores_2ndD_avg_avg_Flip,AllFactorVals_1ndD_avg,opts,PlotFlag,BhvFactor,NtrlsAvgBhv) % calculates correlation between belief and CPI

if contains(TargetFactorDim1,'Rule');warning('first dimension has to be shape or color');return
end
% get belief data
TimeIndBelief = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(3, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(3, 2);

BeliefEncoding = ManData.ReshapeStruct2Mat(Scores_2ndD_avg_avg_Flip, 'TrialRange', 2);
MeanBelief= mean(BeliefEncoding(:,TimeIndBelief),2);

% take average in time
BeliefEncoding=ManData.SmoothData(BeliefEncoding,NeuAna.WidthSmoothing,'SmoothingMethod','movmean','DimSmoothing',2);


% get the behavioral performance for each trial
nTrialRange=length(AllFactorVals_1ndD_avg);
Bhv_AllTrlRng=arrayfun(@(x) mean(AllFactorVals_1ndD_avg(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
Bhv_AllTrlRng=cell2mat(Bhv_AllTrlRng)';
BhvPerfTxt=['BhvPerf' BhvFactor num2str(NtrlsAvgBhv)];
BhvFactorInd=cellfun(@(x) find(strcmp(ClassifierOpts.AnalysisOpts.FactorInds2Keep,x)),{BhvPerfTxt});
BhvPerfTrl=Bhv_AllTrlRng(:,BhvFactorInd);

% Calculate the correlation between mean scores and compression index
X=MeanBelief;
Y=BhvPerfTrl;
[a,p] = ManData.Correlation(X,Y,0,opts.CorrelationMetric);

if PlotFlag
    scatter(X,Y);
    xlabel('Average belief Encoding')
    ylabel(BhvPerfTxt)
    axis square
end
% calculate correlation in time as well
[a_time,p_time] = ManData.Correlation(BeliefEncoding,BhvPerfTrl,0,opts.CorrelationMetric);


end