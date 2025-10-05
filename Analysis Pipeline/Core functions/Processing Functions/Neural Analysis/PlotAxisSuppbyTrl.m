
TrlRange =obj.ManData.GenMovingInds(1,50,15,1);
NTargFactors=1;
TargetFactor={'Axis'};
TimInd=Time>=-0.4 & Time<=0;
clear CatPSTHPreStim
for i=1:length(TrlRange)
    FactorizedData=obj.LimitFactorizedDataTrialsBasedonFactor(AnalysisOpts.FactorizedData(ChNumSorted),'TrialNum',TrlRange{i});
    for f=1:NTargFactors % loop on the factors we want
        factind=f;%TargetFactorInds(f);
        FactorColor=AnalysisOpts.([TargetFactor{f} 'Colors']);
        for Neu=1:length(FactorizedData) % loop on the neurons
            % find the data for this factor
            [SortFactorizedData(Neu).data,SortFactorizedData(Neu).FactorLevels]=obj.SortFactorDatabyLevel(FactorizedData(Neu),TargetFactor{factind});
        end
        % now concatinate the PSTH for each neuron based on the value of the factor for each level
        [CatPSTH]=obj.CatRasterPSTH(SortFactorizedData,[],Time,[],obj.SubtractBaseLine);
        CatPSTHPreStim{1}(i)= arrayfun(@(x) mean(CatPSTH(x).ZSpkCountBin(:,TimInd),2),1,'UniformOutput',0);
        CatPSTHPreStim{2}(i)= arrayfun(@(x) mean(CatPSTH(x).ZSpkCountBin(:,TimInd),2),2,'UniformOutput',0);
    end
end
figure
arrayfun(@(x) obj.FigParams.PlotMeanStd(1:length(TrlRange),cell2mat(CatPSTHPreStim{x}),[],['Trial from Switch(15 trls window)'],'Norm PSTH',x,...
    obj.MeanStdPlotType,[],'AppendTitles',1,'NormalizebyMax',0,'WidthSmoothing',1,'SmoothingMethod','movmean','LegendTxt',['Axis' num2str(x)],'include_n',0),1:2);


