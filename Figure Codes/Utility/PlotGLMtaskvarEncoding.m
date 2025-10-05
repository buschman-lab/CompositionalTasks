function varargout=PlotGLMtaskvarEncoding(FigName,DataPath) % plots results  for Fig.1h-m and Fig. S3
% FigName is the figure name in the manuscript
% DataPath is the path to Figure Data folder

load([DataPath 'GLMdata.mat']);

% load classes 
FigParams=fig_params;
NeuAna=NeuralAnalysisFuncsTemp;

% plot comparisions for each factor
[h1]=FigParams.RenderFigure(2,[]);
[h1{1},Sp1]=FigParams.RenderSubplots([],[],h1{1},NTargFactors);
[h1{2},Sp2]=FigParams.RenderSubplots([],[],h1{2},2);
k=1;

for Ar=AreaNum
    % plot each factor in a sperate subplot with its significance for all neurons
    arrayfun(@(x) NeuAna.PlotFactorizedData(CPDsubsh{k},AnalysisOpts.Time,1,Col(Ar,:),DataFieldNames{x},...
        [DataFieldNames{x} ' information'],AnalysisOpts.Xlabel,'CPD',h1(1),Sp1(x),'Neu2Use',...
        [],'MeanStdPlotType',3,'SubtractBaseLine',0,'ThisMarker','none',...
        'NormalizebyMax',0,'ThisLegendTxt',[AnalysisOpts.AreaNames{Ar}],'n_movavg',obj.WidthSmoothing),1:NTargFactors,'UniformOutput',1);%IndNeuSigFactor{Ar}(IndMatchCPDFieldNames(x),:)

   k=k+1;
end
% plot timing for each factor averaged across all areas
arrayfun(@(x) NeuAna.PlotFactorizedData(CPDsubshAll,AnalysisOpts.Time,1,ColFactors(GetFieldInd(DataFieldNamesNoRule{x}),:),DataFieldNamesNoRule{x},...
    ['Time course of task information for all areas'],AnalysisOpts.Xlabel,'CPD',h1(2),Sp2(1),'Neu2Use',...
    '','MeanStdPlotType',1,'SubtractBaseLine',1,'NPnts_SubtractBaseLine','auto',...
    'NormalizebyMax',1,'NormalizebyMean',0,'ThisLegendTxt',[DataFieldNamesNoRule{x}],'n_movavg',obj.WidthSmoothing,...
    'ThisMarker','none'),1:length(DataFieldNamesNoRule),'UniformOutput',1);

% plot rule information on its own
arrayfun(@(x) NeuAna.PlotFactorizedData(CPDsubshAll,AnalysisOpts.Time,1,ColFactors(GetFieldInd('Rule'),:),'Rule',...
    ['Time course of task information for all areas'],AnalysisOpts.Xlabel,'CPD',h1(2),Sp2(2),'Neu2Use',...
    '','MeanStdPlotType',1,'SubtractBaseLine',0,'NPnts_SubtractBaseLine',5,...
    'NormalizebyMax',1,'NormalizebyMean',0,'ThisLegendTxt',['Rule'],'n_movavg',obj.WidthSmoothing,...
    'ThisMarker','none'),1,'UniformOutput',1);
ylim([0 1.2])

varargout=h1;