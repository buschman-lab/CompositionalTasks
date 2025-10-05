AnalysisOpts.Time=-0.55:0.01:0.6;
TimeInd=AnalysisOpts.Time>-0.15 & AnalysisOpts.Time<=0.6;
TimeTJB=AnalysisOpts.Time(TimeInd);
StatObsvTJB=squeeze(StatObsv);
StatShuffTJB=squeeze(StatShuff)';
StatObsvTJB=StatObsvTJB(TimeInd);
StatShuffTJB=StatShuffTJB(TimeInd,:);
ClusterCorrectionTJB(StatObsvTJB,StatShuffTJB,'Time',TimeTJB,'ThresholdPercentage',20,'DoubleSided',1)


[~, ~, cluster_p]=ClusterCorrectionTJB(smoothdata(rand(500,1),1,'movmean',15),smoothdata(rand(500,250),1,'movmean',15),'ThresholdPercentage',20,'DoubleSided',1)



AnalysisOpts.Time=TimeTJB;
StatShuffLim=StatShuff(:,:,TimeInd);
StatObsvLim=StatObsv(:,:,TimeInd);
[p_values, t_sums, clustIdx,permutation_distribution] = ...
                obj.ManData.ClusterMassCorrection_permutationTwoTail(StatShuffLim,StatObsvLim,AnalysisOpts.pvalClassifierAnalysisTrlShuff_ClusterCorrect,two_sided,'ShowClustCorrectionPlot',1);
     