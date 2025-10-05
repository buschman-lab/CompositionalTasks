function [a,p,X,Y]=CalCPIStimEncodingCorr(ManData,ClassifierOpts,TargetFactorDim1,ScoresMag_1ndD_RS,CPI,opts,PlotFlag) % calculates correlation between belief and stimulus encoding

% CPI=ManData.SmoothData(CPI,15,'SmoothingMethod','movmean','DimSmoothing',2);
% CPI=ManData.SmoothData(CPI,3,'SmoothingMethod','movmean','DimSmoothing',1);

% Calculate the mean of the compression index within a specific time period
 TimeIndComp = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(1, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(1, 2);
CPI=log(CPI);
CompIndex = (mean(CPI(:, TimeIndComp), 2));

% calculate the mean of stimulus encoding
nTrialRange=size(ScoresMag_1ndD_RS,2);
MeanScoresMag_1ndD_RS=arrayfun(@(x) mean(ScoresMag_1ndD_RS(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
AvgStimEncoding=cell2mat(MeanScoresMag_1ndD_RS)';
% Define time indices based on the specified time range
AvgStimEncoding=mean(AvgStimEncoding(:,TimeIndComp),2);

% Calculate the correlation between mean scores and compression index
X=AvgStimEncoding;
Y=CompIndex;
[a,p] = ManData.Correlation(X,Y,0,opts.CorrelationMetric);

if PlotFlag
    scatter(X,Y);
    xlabel('Average Stim Encoding')
    ylabel('Average CPI')
    axis square
end

end