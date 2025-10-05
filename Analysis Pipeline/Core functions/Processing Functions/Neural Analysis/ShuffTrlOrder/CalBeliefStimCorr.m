function [a,p,X,Y,a_avg,X_avg,Y_avg]=CalBeliefStimCorr(ManData,NeuAna,AvgDim2Val,ClassifierOpts,TargetFactorDim1,ScoresMag_1ndD_RS,Scores_2ndD_avg_avg_Flip,opts,PlotFlag) % calculates correlation between belief and stimulus encoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  'RulePreStimCompres'
%obj.PlotSVMProjectionsScores2DPEV_Learning(ClassifierResults_ndD{1},ClassifierResults_ndD{3},ClassifierOpts,...
%    Sp(5),Conds(c),'Line','RulePreStimCompres',[1 3],'ThisTargetFactor','Quadrants','LookatDim2',2,'ExtraData',CompIndexEncoding);

if contains(TargetFactorDim1,'Rule');warning('first dimension has to be shape or color');h=[];return
end

TimeIndBelief = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(3, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(3, 2);

% TimeIndSensory = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(1, 1) & ...
%     ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(1, 2);

TimeIndSensory = ClassifierOpts.AnalysisOpts.Time >= 0.2 & ...
    ClassifierOpts.AnalysisOpts.Time <= 0.3;

BeliefEncoding = ManData.ReshapeStruct2Mat(Scores_2ndD_avg_avg_Flip, 'TrialRange', 2);
MeanBelief= mean(BeliefEncoding(:,TimeIndBelief),2);

% calculate the mean of stimulus encoding
nTrialRange=size(ScoresMag_1ndD_RS,2);
MeanScoresMag_1ndD_RS=arrayfun(@(x) mean(ScoresMag_1ndD_RS(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
AvgStimEncoding=cell2mat(MeanScoresMag_1ndD_RS)';

AvgStimEncodingAvg=mean(AvgStimEncoding(:,TimeIndSensory),2);

AvgStimEncoding=ManData.SmoothData(AvgStimEncoding,NeuAna.WidthSmoothing,'SmoothingMethod','movmean','DimSmoothing',2);

X=MeanBelief;
Y=AvgStimEncoding;
[a,p] = ManData.Correlation(X,Y,0,opts.CorrelationMetric);

if PlotFlag
    plot(ClassifierOpts.AnalysisOpts.Time,a);
    xlabel('Time')
    ylabel([opts.CorrelationMetric ' Correlation'])
    axis square
end

% do it for average timing
X_avg=X;
Y_avg=AvgStimEncodingAvg;
[a_avg] = ManData.Correlation(X,Y_avg,0,opts.CorrelationMetric);

end