function h=PlotClassifierResults(FigNum,DataPath)
% plot TDR anagle analysis results

global AnalysisOpts


ManData=ManipulateData;
FigParams=fig_params;
NeuAna=NeuralAnalysisFuncsTemp;

FigNum=strrep(FigNum,'Fig','');
switch FigNum
    case '2a'
        CondTrip=2;
    case '2b'
    case '2c'
        CondTrip=1;
    case '2e'

    case '2f' 
        CondTrip=6;
    case '2h'

    case '2i'
        CondTrip=3;

    case 'S4a'
        CondTrip=7;
    case 'S4e'

    case 'S4f'

    case 'S4g'

    case 'S4h'

    case 'S5c'
        CondTrip=10;
    case 'S5d'
        CondTrip=8;
    case 'S5e'
        CondTrip=9;
end
Area='PFC';

load([DataPath 'AnalysisOpts.mat']);
load([DataPath Area '_ClassifierData.mat']);


TriplePlots=ClassifierOpts.TriplePlots;
DimTxt=[0 2 3];
RunCrossTemporalClassifer=0;
RunClassifierStatTestAgainst50=0;
NormalizebyMax=0;
PerfMetric='Accuracy';
[Sp,h]=GetFig;

ThisTriplePlots=TriplePlots{CondTrip};
for x=1:size(ThisTriplePlots,1)
    LookatDim2=DimTxt(ThisTriplePlots(x,1));
    Cond=ThisTriplePlots(x,2);
    TrnTstNum=ThisTriplePlots(x,1);
    ThisColor=ThisTriplePlots(x,2);

    [~,TargetFactorTxt,~,~,~,~,~,~,~,TrainCondTxt,TestCondTxt,DimNum]=NeuAna.getClassifierDimInfo(ClassifierOpts,LookatDim2);
    LegTxt=[TargetFactorTxt(1:3) ' D' num2str(DimNum) ' Trn' ManData.ConvMat2Char(ClassifierOpts.(TrainCondTxt){Cond}),...
        'Tst' ManData.ConvMat2Char(ClassifierOpts.(TestCondTxt){Cond})];

    if RunCrossTemporalClassifer
        % plot accuracy as image
        FigParams.Image(Time,Time,mean(MetricVal,3),...
            {AnalysisOpts.Xlabel ;['Test ' ' R' ManData.ConvMat2Char(ClassifierOpts.(TestCondTxt){Cond})]},...
            {AnalysisOpts.Xlabel ;['Train ' ' R' ManData.ConvMat2Char(ClassifierOpts.(TrainCondTxt){Cond})]},...
            'Perf',ClassifierOpts.Name,Sp,'image_colormap',AnalysisOpts.ClassifierAccuracyColormap,'caxis_limits',ClassifierOpts.caxis_limits_XTemp,...
            'WidthSmoothing',15,'WidthSmoothingDim2',15,'SmoothingMethod','movmean');%ClassifierOpts.caxis_limits);
        v=axis;
        xticks(v(1):0.2:v(2));
        xtickangle(0);
        yticks(v(3):0.2:v(4));
        ytickangle(0);
    else
        Col=FigParams.getSingleColor(ThisColor);
        MetricVal=MetricValsOrg_SuperImposed{CondTrip}{x};
        Time=TimeMetricVals_SuperImposed{CondTrip}{x};

        FigParams.PlotMeanStd(Time,MetricVal,[],AnalysisOpts.Xlabel,'Accuracy',Col,...
            1,['Classifier Accuracy Fig:' FigNum],'Sp',Sp,'AppendTitles',0,'NormalizebyMax',0,...
            'WidthSmoothing',15,'SmoothingMethod','movmean','LegendTxt',LegTxt,...
            'STD_method','bootstrap','LegendLoc','best','p_line_style','-','p_MarkerNpnts',10,...
            'p_marker','none','p_marker_size',5);


        ThisStatTest=StatTest{DimNum};% grab statistical test

        if isfield(ThisStatTest,PerfMetric) % add a line for each of the significant clusters
            if length(ThisStatTest)==1;Cond=1;end % then this is the statistical test for this condition
            axis tight
            SigPlotDone=FigParams.plot_significance_level(ThisStatTest(Cond).Accuracy.clusters,...
                ThisStatTest(Cond).Accuracy.statsummery,Time,...
                ClassifierOpts.caxis_limits,Col,[],[],'WidthSmoothing',15,'SmoothingMethod','movmean');

        elseif RunClassifierStatTestAgainst50
            StatTest50_pval=ManData.CalpValShuffle(MetricVal,0.5*ones(1,size(MetricVal,2)),'pvaltail','left');
            AnalysisOpts.ShowStatPvalinPlot=0;
            SigPlotDone= FigParams.plot_significance_level([],StatTest50_pval,AnalysisOpts.Time,ClassifierOpts.caxis_limits,...
                Col,[],[],'skipmultiplecomp_stattest',1,'WidthSmoothing',15,'SmoothingMethod','movmean');
            AnalysisOpts.ShowStatPvalinPlot=1;
        else
            SigPlotDone=0;
        end
        % adjust ylim with regard to stat test results
        v=axis;
        % ylim([min(ClassifierOpts.caxis_limits(1),v(3)) max(ClassifierOpts.caxis_limits(2),v(4))]);
        axis tight
        xticks(v(1):0.2:v(2))
        xtickangle(0)
        if ~NormalizebyMax & (~SigPlotDone | ClassifierOpts.enforceAxisLimits)
            if strcmp(ClassifierOpts.caxis_limits,'auto')
                YLims=v(3:4);
            else
                YLims=ClassifierOpts.caxis_limits;
                FigParams.SetAxisTicks('y',ClassifierOpts.caxis_limits,0.1);
            end
            ylim([YLims(1) max(YLims(2),v(4))]);
        end
        if ~NormalizebyMax
            if isfield(ClassifierOpts,'ChanceLevel');ChanceLevel=ClassifierOpts.ChanceLevel;else;ChanceLevel=0.5;end
            FigParams.PlotHorizontalLine(ChanceLevel,Sp,[0.5 0.5 0.5],'p_line_style','--');
        end
        v=axis;
        if (v(4)-v(3))<0.5
            yticks((floor(v(3) * 10) / 10):0.05:v(4))
        else
            yticks((floor(v(3) * 10) / 10):0.1:v(4))
        end
    end
end
end