function h=PlotCompressionIndexData(FigNum,Area,DataPath)

FigParams=fig_params;
 
load([DataPath 'AnalysisOpts.mat']);
load([DataPath Area '_CompressionData.mat']);

[Sp,h]=GetFig;

for Rule=1:3

    Col=AnalysisOpts.RuleColors(Rule,:);
    ScoresMetric=CompressionEncoding{Rule}.TrlAvg.All;

    FigParams.PlotMeanStd(AnalysisOpts.Time,log(ScoresMetric),[],'Time relative to stimulus onset (s)',...
        [],Col,1,[],'Sp',Sp,'AppendTitles',1,'WidthSmoothing',15,'SmoothingMethod','movmean');  

    Ylbl=['Compression Index'];
    title([FigNum ': Compression Index']);ylabel(Ylbl)
   
end
axis tight
