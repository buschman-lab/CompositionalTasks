function h=PlotClassifierLearningResults(FigNum,DataPath)

global AnalysisOpts

ManData=ManipulateData;
FigParams=fig_params;
NeuAna=NeuralAnalysisFuncsTemp;

FigNum=strrep(FigNum,'Fig','');
FileNames={'Fig4_ShapeColorRuleS1C2C1_LearningData.mat','Fig4_ShapeColorRuleC1C2C1_LearningData.mat',...
    'Fig4_ResponseC1C2C1_LearningData.mat','Fig5_ShapeColorRuleS1C2C1LearningData.mat','Fig5_AxisDecoding_LearningData'};

AnalysisOpts.ThisTimeAxisStart=-0.2;AnalysisOpts.ThisTimeAxisEnd=0.55;

switch FigNum
    case '4c'
        FileNum=1;  
        TrnTstNum=3;   

    case '4d'
        FileNum=2;  
        TrnTstNum=3;

    case '4f'
        FileNum=1;  
        TrnTstNum=2;
    case '4g'
        FileNum=2;
        TrnTstNum=2;
    case '4i'
        FileNum=1;
        TrnTstNum=1;
    case '4j'
        FileNum=2;
        TrnTstNum=1;
    case '4l'
        FileNum=3;
        TrnTstNum=1;
    case '5h'
        FileNum=5;
        TrnTstNum=1;

end
%Area='PFC';

%load([DataPath 'AnalysisOpts.mat']);
load([DataPath FileNames{FileNum}]);

if strcmp(FigNum,'4c') | strcmp(FigNum,'4d') | strcmp(FigNum,'5h') 
    AnalysisOpts.ThisTimeAxisStart=-0.4;AnalysisOpts.ThisTimeAxisEnd=0;
end

Dim=TrnTstNum;Cond=1;

if TrnTstNum==1
    TargetFactorTxt=ClassifierOpts.TargetFactors{1};
    TrainCondTxt=sprintf('TrainCond');
    TestCondTxt=sprintf('TestCond');
elseif TrnTstNum==2
    TargetFactorTxt=ClassifierOpts.TargetFactors_2ndD{1};
    TrainCondTxt=sprintf('TrainCond%i',TrnTstNum);
    TestCondTxt=sprintf('TestCond%i',TrnTstNum);
elseif TrnTstNum==3
    TargetFactorTxt=ClassifierOpts.TargetFactors_3ndD{1};
    TrainCondTxt=sprintf('TrainCond%i',TrnTstNum);
    TestCondTxt=sprintf('TestCond%i',TrnTstNum);
end

[~,TrainTimInd,TestTimInd,nXtimePnt,TimeMatrixSize]=NeuAna.GetTimeRangeforThisCond(ClassifierOpts);

TrialRangeSet=ClassifierOpts.TestTrlNumRange{Cond}(2:4);
TrialRange=[TrialRangeSet(2):TrialRangeSet(3):abs(TrialRangeSet(1))];if isempty(TrialRange);TrialRange=1;end
TrainCond=ClassifierOpts.(TrainCondTxt){Cond};
TestCond=ClassifierOpts.(TestCondTxt){Cond};
nTrialRange=length(TrialRange);

MetricVals=squeeze(mean(MetricValsOrg{Dim},2))';
Time=TimeMetricValsDim{Dim}{1};

MetricVals=ManData.SmoothData(MetricVals,1,'SmoothingMethod','movmean','DimSmoothing',1);
CurrnTrialRange=nTrialRange;
if TrialRangeSet(3)==1;nTrialRange=length(1:5:CurrnTrialRange);MetricVals=MetricVals(1:5:CurrnTrialRange,:);end

Col=copper(nTrialRange);
if ~iscell(MetricVals);MetricVals=arrayfun(@(x) MetricVals(x,:),1:nTrialRange,'UniformOutput',0);end
[Sp,h]=GetFig;
arrayfun(@(TrlRng) FigParams.PlotMeanStd(Time,MetricVals{TrlRng},[],AnalysisOpts.Xlabel,...
    [],Col(TrlRng,:),1,[],'Sp',Sp,'AppendTitles',1,'WidthSmoothing',15,'performtrend_stattest',0,...
    'SmoothingMethod','movmean','STD_method','bootstrap'),1:nTrialRange,'UniformOutput',0);
Ylbl=[{'Accuracy'} ; {TargetFactorTxt}];
title(['Classifier Accuracy for Fig:' FigNum]);ylabel(Ylbl)
if nTrialRange>1
    hCbar=colorbar;
    hCbar.Ticks=1:nTrialRange;
    hCbar.Label.String='Trials from switch';hCbar.Label.FontSize=12;
    hCbar.FontSize=12;
end

%% add trend stat test for each time point
if ~isempty(StatTestTrlShuff{Cond})
    if ~isempty(StatTestTrlShuff{Cond}.(TargetFactorTxt).clusters)
        FigParams.plot_significance_level(StatTestTrlShuff{Cond}.(TargetFactorTxt).clusters,...
            StatTestTrlShuff{Cond}.(TargetFactorTxt).statsummery,Time,'auto','k',[],[],'WidthSmoothing',15,'SmoothingMethod','movmean');
    end
end

axis tight
v=axis;
xticks(v(1):0.2:v(2))
xtickangle(0)
if (v(4)-v(3))<0.5
    yticks((floor(v(3) * 10) / 10):0.05:v(4))
else
    yticks((floor(v(3) * 10) / 10):0.1:v(4))
end

AnalysisOpts.ThisTimeAxisStart=[];AnalysisOpts.ThisTimeAxisEnd=[];