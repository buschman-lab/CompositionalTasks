function SaveFileName=FitArima2CompressionIndex(Lag,Metric,AvgDim2,FileNameSpecs,Observed_Org,Shuffle_Org,Dim)
% AvgDim2 are we averaging across Dim 2 dimension?
% Lag Number of lags for ARIMA
% Metric we care about could be CPI ShapeDist or ColorDist
global AnalysisOpts

ManData=ManipulateData;
NeuAna=NeuralAnalysisFuncsTemp;
fp=fig_params;

PATH=AnalysisOpts.CodeTestDataPath;
Vars={[Metric '_Obsv'],[Metric '_Shuff']};

if isempty(Observed_Org)
    if ~exist([PATH FileNameSpecs 'ShData'],'file')
        File='CompressionIndexShuffData.mat';
        switch Metric
            case 'CPI'
                Vars={'ScoresMetricAll_obsv','ScoresMetricAll_Shuff'};
            case 'ShapeDist'
                Vars={'ShapeDistAvg_Obsv','ShapeDistAvg_Shuff'};
            case 'ColorDist'
                Vars={'ColorDistAvg_Obsv','ColorDistAvg_Shuff'};
        end
        Observed_Org=load([PATH File],Vars{1});Shuffle_Org=load([PATH File],Vars{2});
    else
        % try and load data from the file
        LoadFileName=[PATH FileNameSpecs 'ShData'];
        Observed_Org=load(LoadFileName,[Metric '_Obsv']);
        Shuffle_Org=load(LoadFileName, [Metric '_Shuff']);
    end
    
else
    Observed.([Metric '_Obsv'])=Observed_Org;Shuffle.([Metric '_Shuff'])=Shuffle_Org;
end

Observed=ManData.SmoothData(Observed.(Vars{1}),NeuAna.WidthSmoothing,'SmoothingMethod',NeuAna.SmoothingMethod,'DimSmoothing',2);
Shuffle=ManData.SmoothData(Shuffle.(Vars{2}),NeuAna.WidthSmoothing,'SmoothingMethod',NeuAna.SmoothingMethod,'DimSmoothing',2);

model = arima('ARLags', Lag, 'D', 0, 'Constant', 0); % AR(1) model
TrlRng=[1:size(Observed,1)]';
NTim=size(Observed,2);
% avergae across dim2?
if AvgDim2
    Observed=ManData.SmoothData(Observed,3,'SmoothingMethod','movmean','DimSmoothing',1);
    Shuffle=ManData.SmoothData(Shuffle,3,'SmoothingMethod','movmean','DimSmoothing',1);
end

[ArObserved,ArObserved_p]= (arrayfun(@(x) ARIMAtrendcorr(Observed(:,x),TrlRng,model,'Spearman'),1:NTim,'UniformOutput',0));

for i=1:size(Shuffle,3)
    fprintf('\n Cal ArShuff %i',i);
    tic
    [ArShuffle(:,:,i),ArShuffle_p(:,:,i)]= (arrayfun(@(x) ARIMAtrendcorr(squeeze(Shuffle(:,x,i)),TrlRng,model,'Spearman'),1:NTim,'UniformOutput',0));
    toc
end
if AvgDim2;Ext='AvgDim2';else;Ext='';end
SaveFileName=['ArimaTestShuffLag' num2str(Lag) Ext '_D' num2str(Dim) '_' Metric FileNameSpecs];
save([PATH SaveFileName],'ArShuffle','ArShuffle_p','ArObserved','ArObserved_p');

%% plot and save results
% model = arima('ARLags', 1, 'D', 0, 'Constant', 0); % AR(1) model
% [ArObserved,ArObserved_p]= (arrayfun(@(x) ARIMAtrendcorr(ObservedAvg(:,x),TrlRng,model,'Spearman'),40:NTim,'UniformOutput',0));
% plot tau and pval on observed
RepSpace=@(x) strrep(x,'_',' ');
col=copper(16);Time=-0.6:0.01:0.55;
figure;subplot(121);hold on;arrayfun(@(x) plot(Time,Observed(x,1:NTim),'color',col(x,:)),TrlRng');
yyaxis right;plot(Time,cell2mat(ArObserved_p),'b');ylabel pval;xlabel Time;axis square
title(['Pval Arima+Spearman on observed with lag ' num2str(Lag)])

subplot(122);hold on;arrayfun(@(x) plot(Time,Observed(x,1:NTim),'color',col(x,:)),TrlRng');
yyaxis right;plot(Time,cell2mat(ArObserved),'r');ylabel Tau;xlabel Time;axis square
title(['Tau Arima+Spearman on observed with lag ' num2str(Lag)])
sgtitle(RepSpace(FileNameSpecs))

% plot with cluster correction
AnalysisOpts.Time=Time; % this has to be modified
ManData.ClusterMassCorrection_permutationTwoTail(permute(cell2mat(ArShuffle),[3 1 2]),permute(cell2mat(ArObserved),[1 3 2]),0.2,1,'ShowClustCorrectionPlot',1);
sgtitle(RepSpace(FileNameSpecs))
end