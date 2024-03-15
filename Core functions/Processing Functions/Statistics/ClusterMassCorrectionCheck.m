global AnalysisOpts
figure
subplot(221)
plot(AnalysisOpts.Time,squeeze(observed),'linewidth',5)
hold on
plot(AnalysisOpts.Time,squeeze(p_threshold))
yyaxis right
hold on
PVals=nan(size(clustIdx));for i=1:length(obs_clust_p);PVals(clustIdx==i)=obs_clust_p(i);end
plot(AnalysisOpts.Time,PVals,'k')
ylim([-0.1 1])
legend('Observed',['threshold alpha' num2str(p_thresh)],'Observed>threshold')
yyaxis left
xlabel('Time(s)');ylabel('Classifier Accuracy')
title(sprintf('# of Shuffle Reps:%i',length(AnalysisOpts.NRep2Use4StatTest)));
subplot(222)
hold on
histogram(max_clustMass_dist,'Normalization','probability')
v=axis;
arrayfun(@(x) plot([x x],[v(3) v(4)],'r'),abs(clustMass))
xlabel('Cluster size');ylabel('Probability')
subplot(223)
histogram(observed,'Normalization','probability')
hold on
histogram(Shuffle,'Normalization','probability')
legend('Prob Observed','Prob Shuffle')
