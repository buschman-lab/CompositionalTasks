% ObservedBuff=Observed;
% ShuffleBuff=Shuffle;
% 
% Observed=obj.ManData.SmoothData(Observed,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',3);
% %Observed=obj.ManData.SmoothData(Observed,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
% Shuffle=obj.ManData.SmoothData(Shuffle,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',3);
% %Shuffle=obj.ManData.SmoothData(Shuffle,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);


ObservedCorr=squeeze(Observed);
ShuffleCorr=squeeze(Shuffle);

TrlRng=[1:size(Shuffle,2)]';
Nrep=size(ShuffleCorr,1);
NTim=size(ObservedCorr,2);

CorrType='spearman';
 %CorrObsv=arrayfun(@(x) corr(TrlRng,ObservedCorr(:,x),'type',CorrType),1:NTim);
 %CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(ShuffleCorr(rep,:,x))','type',CorrType),1:NTim)',1:Nrep,'UniformOutput',0))';
CorrObsv=arrayfun(@(x) LinearRegression(TrlRng,ObservedCorr(:,x)),1:NTim);
CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) LinearRegression(TrlRng,squeeze(ShuffleCorr(rep,:,x))),1:NTim)',1:Nrep,'UniformOutput',0))';

CorrObsv=permute(CorrObsv,[3 1 2]);
CorrShuff=permute(CorrShuff,[1 3 2]);

[p_values, t_sums, clustIdx,permutation_distribution] = ...
    obj.ManData.ClusterMassCorrection_permutationTwoTail(CorrShuff,CorrObsv,AnalysisOpts.pvalClassifierAnalysisTrlShuff_ClusterCorrect,two_sided,'ShowClustCorrectionPlot',1);
clusters=arrayfun(@(x) find(clustIdx==x),unique(clustIdx(clustIdx~=0)),'UniformOutput',0);


figure 
Col=copper(16);
subplot(421)
imagesc(ObservedCorr)
subplot(422)
MeanShuffle=squeeze(mean(ShuffleCorr,1));
imagesc(MeanShuffle)
subplot(423)
hold on
arrayfun(@(x) plot(AnalysisOpts.Time,ObservedCorr(x,:),'Color',Col(x,:)),1:16)

subplot(424)
hold on
arrayfun(@(x) errorbar(AnalysisOpts.Time,MeanShuffle(x,:),std(squeeze(ShuffleCorr(:,x,:)),1,1),'Color',Col(x,:)),1:16);

subplot(425)
hold on
 plot(AnalysisOpts.Time,squeeze(ShuffleCorr(:,1,:)))

shuffind=10; 
subplot(426)
cla
hold on
arrayfun(@(x) plot(AnalysisOpts.Time,squeeze(ShuffleCorr(shuffind,x,:)),'Color',Col(x,:)),1:16)


subplot(427)
cla
hold on
plot(AnalysisOpts.Time,squeeze(CorrShuff(shuffind,1,:)));

subplot(428)
cla
hold on
plot(AnalysisOpts.Time,squeeze(CorrShuff)');
plot(AnalysisOpts.Time,squeeze(CorrObsv),'LineWidth',5);


