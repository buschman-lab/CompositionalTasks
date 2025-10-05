global AnalysisOpts
fp=fig_params;ManData=ManipulateData;
Row=2;Col=4;
figure
subplot(Row,Col,1)
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
ylabel('P val')
v=axis;
plot([v(1) v(2)],[0.05 0.05],'g--') % show 0.05 line 
% show observed
yyaxis left
xlabel('Time(s)');ylabel('Classifier Accuracy')
title(sprintf('# of Shuffle Reps:%i',length(AnalysisOpts.NRep2Use4StatTest)));
subplot(Row,Col,2)
hold on
histogram(max_clustMass_dist,'Normalization','probability')
v=axis;
arrayfun(@(x) plot([x x],[v(3) v(4)],'r'),abs(clustMass))
xlabel('Cluster size');ylabel('Probability')
legend('clust siz shuff','clust siz obs','Location','best')
title({['distribution of cluster']; ['size for shuffle and observed']})

subplot(Row,Col,3)
histogram(observed,'Normalization','probability')
hold on
histogram(Shuffle,'Normalization','probability')
legend('Prob Observed','Prob Shuffle','Location','best')
title({['dist of values across all'];[ 'time for shuffle and observed']})


subplot(Row,Col,4)% show p value shuffle 
yyaxis left
hold on
PvalShufflNonCorrected=ManData.CalpValShuffle(squeeze(Shuffle),squeeze(observed),'pvaltail','both');
PvalShufflNonCorrected_Left=ManData.CalpValShuffle(squeeze(Shuffle),squeeze(observed),'pvaltail','left');
PvalShufflNonCorrected_right=ManData.CalpValShuffle(squeeze(Shuffle),squeeze(observed),'pvaltail','right');

plot(AnalysisOpts.Time,mapminmax(squeeze(observed)',0,1),'r','LineWidth',1)
ylabel('Obs Val')
yyaxis right
plot(AnalysisOpts.Time,PvalShufflNonCorrected,'b','LineWidth',2)
plot(AnalysisOpts.Time,PvalShufflNonCorrected_Left,'k','LineWidth',2)
plot(AnalysisOpts.Time,PvalShufflNonCorrected_right,'c','LineWidth',2)
ylabel('P val')
legend({'Observed','pval twotail','pval left','pval right'},'Location','best','Box','off')
ylim([0 1]);
xlabel('Time(s)');
title('Uncorrected pvalue')

% show 0.05 line 
v=axis;
plot([v(1) v(2)],[0.05 0.05],'g--')

subplot(Row,Col,5)
hold on
plot(AnalysisOpts.Time,squeeze(Shuffle)')
xlabel('Time(s)');
title('full distribution of shuffle with observed')
plot(AnalysisOpts.Time,squeeze(observed),'r','LineWidth',5)

% subplot(Row,Col,6)
% hold on
% fp.PlotMeanStd(AnalysisOpts.Time,squeeze(Shuffle),[],'Time','',1,1,'MeanStd Shuffle','STD_method','bootstrap','LegendTxt','Shuf');
% fp.PlotMeanStd(AnalysisOpts.Time,squeeze(observed)',[],'Time','',2,1,'MeanStd Shuffle','STD_method','bootstrap','LegendTxt','Obsv');
% xlabel('Time(s)');
% title(obj.ThisTitle)

