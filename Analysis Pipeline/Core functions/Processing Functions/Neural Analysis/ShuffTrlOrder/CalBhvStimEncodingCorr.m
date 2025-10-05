function [a,p,X,Y]=CalBhvStimEncodingCorr(ManData,ClassifierOpts,TargetFactorDim1,ScoresMag_1ndD_RS,AllFactorVals_1ndD_avg,opts,PlotFlag,BhvFactor,NtrlsAvgBhv) % calculates correlation between belief and stimulus encoding


TimeIndStim = ClassifierOpts.AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(1, 1) & ...
    ClassifierOpts.AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(1, 2);

% get the behavioral performance for each trial
nTrialRange=length(AllFactorVals_1ndD_avg);
Bhv_AllTrlRng=arrayfun(@(x) mean(AllFactorVals_1ndD_avg(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
Bhv_AllTrlRng=cell2mat(Bhv_AllTrlRng)';
BhvPerfTxt=['BhvPerf' BhvFactor num2str(NtrlsAvgBhv)];
BhvFactorInd=cellfun(@(x) find(strcmp(ClassifierOpts.AnalysisOpts.FactorInds2Keep,x)),{BhvPerfTxt});
BhvPerfTrl=Bhv_AllTrlRng(:,BhvFactorInd);

% calculate the mean of stimulus encoding
nTrialRange=size(ScoresMag_1ndD_RS,2);
MeanScoresMag_1ndD_RS=arrayfun(@(x) mean(ScoresMag_1ndD_RS(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
AvgStimEncoding=cell2mat(MeanScoresMag_1ndD_RS)';

% AvgStimEncoding=ManData.SmoothData(AvgStimEncoding,15,'SmoothingMethod','movmean','DimSmoothing',2);
% AvgStimEncoding=ManData.SmoothData(AvgStimEncoding,3,'SmoothingMethod','movmean','DimSmoothing',1);

% Define time indices based on the specified time range
AvgStimEncoding=mean(AvgStimEncoding(:,TimeIndStim),2);

% Calculate the correlation between mean scores and compression index
X=AvgStimEncoding;
Y=BhvPerfTrl;
[a,p] = ManData.Correlation(X,Y,0,opts.CorrelationMetric);

if PlotFlag
    scatter(X,Y);
    xlabel('Average Stim Encoding')
    ylabel('Average behavioral perforemance')
    axis square
end

end