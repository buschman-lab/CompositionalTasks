function h=PlotCorrelationAnalysis(FigNum,DataPath)

global AnalysisOpts

ManData=ManipulateData;
FigParams=fig_params;

FigNum=strrep(FigNum,'Fig','');
FileNames={'Learning3D_Sh_Co_Ru_AltRule_45_120_CorrData.mat','Learning3D_Sh_Co_Ru_AltRule_40_110_CorrData','Learning3D_Sh_Co_Ru_SameRule_CorrData','Learning3D_Co_Res_Ru_AltRule_CorrData'};
MarkerSize=50;
CalShuffTrlOrderClassifier_ClusterTH=24;

switch FigNum
    case '4h'
        FileNum=2;
        PlotType=2; % plot in time
        SubField='BeliefStimCorr_Color';
        Xlabel='Time relative to stimulus onset (s)'; Ylabel='Zscorred Kendall corr';
        XLims=[-0.2 0.55];YLims=[-0.12 2];
    case '4e'
        FileNum=2;
        PlotType=2; % plot in time
        SubField='BeliefBhvPerfcorr_Color_Time';
        Xlabel='Time relative to stimulus onset (s)'; Ylabel='Zscorred Kendall corr';
        XLims=[-0.4 0];YLims=[-0.12 1.7];
    case '5e'
        FileNum=1;
        PlotType=1; % plot in scatter
        SubField='BeliefCPIcorr';
        Xlabel='Task Belief Encoding'; Ylabel='Compression Index';
        XLims=[0.4 0.55];YLims=[-0.2 0.8];
    case '5f'
        FileNum=1;
        PlotType=1; % plot in scatter
        SubField='CPIStimEncodingCorr';
        Xlabel='Color Category Encoding'; Ylabel='Compression Index';
        XLims=[0.56 0.67];YLims=[-0.15 0.75];
    case 'S8c'
        FileNum=2;
        PlotType=2; % plot in time
        SubField='BeliefBhvPerfcorr_Shape_Time';
        Xlabel='Time relative to stimulus onset (s)'; Ylabel='Zscorred Kendall corr';
        XLims=[-0.4 0];YLims=[-1.5 0];
    case 'S8d'
        FileNum=3;
        PlotType=2; % plot in time
        SubField='BeliefBhvPerfcorr_Color_Time';
        Xlabel='Time relative to stimulus onset (s)'; Ylabel='Zscorred Kendall corr';
        XLims=[-0.6 0.55];YLims=[-1.5 2];
    case 'S8e'
        FileNum=3;
        PlotType=2; % plot in time
        SubField='BeliefBhvPerfcorr_Shape_Time';
        Xlabel='Time relative to stimulus onset (s)'; Ylabel='Zscorred Kendall corr';
        XLims=[-0.6 0.55];YLims=[-1.7 1.5];
    case 'S8h'
        FileNum=2;
        PlotType=1; % plot in scatter
        SubField='BeliefStimCorr_Shape_avg';
        Xlabel='Task Belief Encoding'; Ylabel='Shape Category Encoding';
        XLims=[0.4 0.55];YLims=[0.45 0.62];
    case 'S8i'
        FileNum=4;
        PlotType=1; % plot in scatter
        SubField='BeliefRespCorr';
        Xlabel='Task Belief Encoding'; Ylabel='Response Direction Encoding';
        XLims=[0.45 0.52];YLims=[0.45 0.57];



end

load([DataPath 'AnalysisOpts.mat']);
load([DataPath FileNames{FileNum}],'CorrelationVar_AvgDim2');

[Sp,h]=GetFig;

if PlotType==1 % plot as scatter
    if iscell(CorrelationVar_AvgDim2.a_obsv.(SubField));CorrelationVar_AvgDim2.a_obsv.(SubField)=cell2mat(CorrelationVar_AvgDim2.a_obsv.(SubField));end
    if iscell(CorrelationVar_AvgDim2.a_sh.(SubField));CorrelationVar_AvgDim2.a_sh.(SubField)=cell2mat(CorrelationVar_AvgDim2.a_sh.(SubField));end

    aCorrected_obsv=ManData.CalZscoreShuffle(CorrelationVar_AvgDim2.a_obsv.(SubField),CorrelationVar_AvgDim2.a_sh.(SubField));
    PvalShuff=ManData.CalpValShuffle(CorrelationVar_AvgDim2.a_sh.(SubField)',CorrelationVar_AvgDim2.a_obsv.(SubField));
    X=cell2mat(CorrelationVar_AvgDim2.a_obsv.([SubField '_X'])); Y=cell2mat(CorrelationVar_AvgDim2.a_obsv.([SubField '_Y']));
    scatter(X,Y,MarkerSize,'Marker','o','MarkerFaceColor','k',MarkerEdgeColor='k');
    xlabel(Xlabel)
    ylabel(Ylabel)
    xlim(XLims);ylim(YLims);
    axis square
    title(sprintf('Figure %s: Correlation %s and %s aZ=%0.3f,p=%0.3f',FigNum,Xlabel,Ylabel,aCorrected_obsv,PvalShuff));
else

    if iscell(CorrelationVar_AvgDim2.a_obsv.(SubField));CorrelationVar_AvgDim2.a_obsv.(SubField)=cell2mat(CorrelationVar_AvgDim2.a_obsv.(SubField));end
    if iscell(CorrelationVar_AvgDim2.a_sh.(SubField));CorrelationVar_AvgDim2.a_sh.(SubField)=cell2mat(CorrelationVar_AvgDim2.a_sh.(SubField));end

    [a_obsv,a_sh]=ManData.CalZscoreShuffle(CorrelationVar_AvgDim2.a_obsv.(SubField),CorrelationVar_AvgDim2.a_sh.(SubField));
    ClusterCorrectionTJB(a_obsv,a_sh,{[],[1,1,1]},'Time',-0.6:0.01:0.55,'PlotSignificantClusts',1, ...
        'ThresholdPercentage',CalShuffTrlOrderClassifier_ClusterTH,'DoubleSided',1,'ZscoreData',0,'PlotSignificantClustsTH',0.05);
    xlim(XLims)
    axis square
    xlabel(Xlabel)
    ylabel(Ylabel)
    title(SubField)

end


