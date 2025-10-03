function h=PlotClassifierEncodingAxisCompressionIndex(FigNum,DataPath)

global AnalysisOpts

ManData=ManipulateData;
FigParams=fig_params;
 
FigNum=strrep(FigNum,'Fig','');
FileNames={'Fig5_ShapeColorRuleS1C2C1_LearningData.mat'};

switch FigNum
    case '5a'
        FileNum=1;
        EncodingDistVar='ColorAvgDist';  
        WidthSmoothing=1;
        UseLogScale=0;
    case '5b'
        FileNum=1;
        EncodingDistVar='ShapeAvgDist';
        WidthSmoothing=1;
        UseLogScale=0;
    case '5d'
        FileNum=1;
        EncodingDistVar='CPI';
        WidthSmoothing=15;
        UseLogScale=1;
end

load([DataPath 'AnalysisOpts.mat']);
load([DataPath FileNames{FileNum}]);

[Sp,h]=GetFig;
Title=['Figure:' FigNum];
Time=ClassifierOpts.AnalysisOpts.Time;

ScoresMetric=EncodingDist{strcmp(EncodingDistVars,EncodingDistVar)};
if UseLogScale;ScoresMetric=ManData.SmoothData(ScoresMetric,3,'SmoothingMethod','movmean','DimSmoothing',1);end

nTrialRange=16;
colormap(Sp,copper(nTrialRange));
Col=copper(nTrialRange);

if ~iscell(ScoresMetric);ScoresMetric=arrayfun(@(x) ScoresMetric(x,:),1:nTrialRange,'UniformOutput',0);end

if UseLogScale;ScoresMetric=cellfun(@log, ScoresMetric,'UniformOutput',0);end

arrayfun(@(TrlRng) FigParams.PlotMeanStd(Time,ScoresMetric{TrlRng},[],ClassifierOpts.AnalysisOpts.Xlabel,...
    [],Col(TrlRng,:),1,[],'Sp',Sp,'AppendTitles',1,'WidthSmoothing',WidthSmoothing,...
    'SmoothingMethod','movmean'),1:nTrialRange,'UniformOutput',0);

% if we are performing trend stat test on this data as well
X=StatTestTrlShuff{1}.(EncodingDistVar).clusters;
P=StatTestTrlShuff{1}.(EncodingDistVar).statsummery;
AnalysisOpts.ShowStatPvalinPlot=0; % note that the all p values p<0.05 are converted to 0.049, p<0.01 to 0.0099 and p<0.001 to 0.0009
FigParams.plot_significance_level(X,...
    P,Time,'auto','k',[],[],'WidthSmoothing',WidthSmoothing,'SmoothingMethod','movmean');


Ylbl=[EncodingDistVar];
title(Title);ylabel(Ylbl)
if nTrialRange>1
    hCbar=colorbar;
    L=hCbar.Ticks(end);
    hCbar.Ticks=hCbar.Ticks(1):L/(nTrialRange-1):hCbar.Ticks(end);
    hCbar.TickLabels=arrayfun(@(x) num2str(x),45:5:120,'UniformOutput',0);
    hCbar.Label.String='Trials from switch';
    hCbar.FontSize=10;
end
axis tight
