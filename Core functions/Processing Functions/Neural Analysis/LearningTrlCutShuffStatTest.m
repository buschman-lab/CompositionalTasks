StatObsv=Observed(:,2,:)-Observed(:,1,:);
StatShuff=Shuffle(:,2,:)-Shuffle(:,1,:);

[p_values, t_sums, clustIdx,permutation_distribution] = ...
    obj.ManData.ClusterMassCorrection_permutationTwoTail(StatShuff,StatObsv,0.1,1,'ShowClustCorrectionPlot',1);


figure 
Col=copper(2);
subplot(221)
plot(AnalysisOpts.Time,squeeze(StatObsv),'r')
MeanShuffle=squeeze(mean(StatShuff,1));
hold on
plot(AnalysisOpts.Time,MeanShuffle,'b')
legend('avg obsv','avg shuff')
subplot(222)
hold on
plot(AnalysisOpts.Time,squeeze(StatShuff)')
 
