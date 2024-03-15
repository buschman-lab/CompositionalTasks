global AnalysisOpts
figure
subplot(231)
plot(AnalysisOpts.Time,squeeze(observed),'linewidth',5)
hold on
plot(AnalysisOpts.Time,squeeze(p_thresholdLower))
plot(AnalysisOpts.Time,squeeze(p_thresholdUpper))
yyaxis right
hold on
PVals=nan(size(clustIdx));for i=1:length(obs_clust_p);PVals(clustIdx==i)=obs_clust_p(i);end
plot(AnalysisOpts.Time,PVals,'k')
ylim([-0.1 1])
legend('Obsv',['th ' num2str(p_thresh)],'Obvs>th','Location','best')
yyaxis left
xlabel('Time(s)');ylabel('Classifier Accuracy')
title(sprintf('# of Shuffle Reps:%i',length(AnalysisOpts.NRep2Use4StatTest)));
subplot(232)
hold on
histogram(max_clustMass_dist,'Normalization','probability')
v=axis;
arrayfun(@(x) plot([x x],[v(3) v(4)],'r'),abs(clustMass))
xlabel('Cluster size');ylabel('Probability')
legend('clust siz shuff','clust siz obs','Location','best')
title({['distribution of cluster']; ['size for shuffle and observed']})

subplot(233)
histogram(observed,'Normalization','probability')
hold on
histogram(Shuffle,'Normalization','probability')
legend('Prob Observed','Prob Shuffle','Location','best')
title({['dist of values across all'];[ 'time for shuffle and observed']})


subplot(234)
plot(AnalysisOpts.Time,squeeze(observed),'r','LineWidth',5)
MeanShuffle=squeeze(mean(Shuffle,1));
hold on
plot(AnalysisOpts.Time,MeanShuffle,'b')
legend('Observed','avg Shuffle','Location','best')
xlabel('Time(s)');
title('Comparision of avg shuffle and observed')

subplot(235)
hold on
plot(AnalysisOpts.Time,squeeze(Shuffle)')
xlabel('Time(s)');
title('full distribution of shuffle with observed')
plot(AnalysisOpts.Time,squeeze(observed),'r','LineWidth',5)


