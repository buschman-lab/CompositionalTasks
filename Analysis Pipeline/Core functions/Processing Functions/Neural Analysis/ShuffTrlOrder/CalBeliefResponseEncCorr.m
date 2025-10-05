function [a,p,X,Y]=CalBeliefResponseEncCorr(ManData,NeuAna,AvgDim2Val,ClassifierOpts,TargetFactorDim1,ScoresMag_1ndD_RS,Scores_2ndD_avg_avg_Flip,opts,PlotFlag) % calculates correlation between belief and stimulus encoding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  'RulePreStimCompres'
%obj.PlotSVMProjectionsScores2DPEV_Learning(ClassifierResults_ndD{1},ClassifierResults_ndD{3},ClassifierOpts,...
%    Sp(5),Conds(c),'Line','RulePreStimCompres',[1 3],'ThisTargetFactor','Quadrants','LookatDim2',2,'ExtraData',CompIndexEncoding);

if contains(TargetFactorDim1,'Rule');warning('first dimension has to be shape or color');h=[];return
end

TimeIndBelief = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(3, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(3, 2);

BeliefEncoding = ManData.ReshapeStruct2Mat(Scores_2ndD_avg_avg_Flip, 'TrialRange', 2);
MeanBelief= mean(BeliefEncoding(:,TimeIndBelief),2);

% calculate the mean of response encoding
TimeIndResponse = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(1, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(1, 2);

nTrialRange=size(ScoresMag_1ndD_RS,2);
MeanScoresMag_1ndD_RS=arrayfun(@(x) mean(ScoresMag_1ndD_RS(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
AvgRespEncoding=cell2mat(MeanScoresMag_1ndD_RS)';
% Define time indices based on the specified time range
AvgRespEncoding=mean(AvgRespEncoding(:,TimeIndResponse),2);


X=MeanBelief;
Y=AvgRespEncoding;
[a,p] = ManData.Correlation(X,Y,0,opts.CorrelationMetric);

if PlotFlag
    scatter(X,Y);
    ylabel('Average Response Encoding')
    xlabel('Average Belief')
    axis square
end

end