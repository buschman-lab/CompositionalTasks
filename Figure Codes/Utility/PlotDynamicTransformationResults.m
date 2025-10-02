function h=PlotDynamicTransformationResults(FigName,DataPath)
% FigName is the figure name in the manuscript
% DataPath is the path to Figure Data folder

% load the data that we need 
load([DataPath 'DynamicTransformationData.mat']);

% Define Classes we use
FigParams=fig_params;
ManData=ManipulateData;

% create figure and subplots
h=FigParams.RenderFigure(1,[]);
[h,Sp]=FigParams.RenderSubplots([1],[3],h{1},3);

if contains(FigName,'Fig3')
    ConditionTitle={'Corr ColorCat C1-X> Resp Dir C2','Corr ColorCat C1-> Resp Dir S1'};
    aCorr=[aCorrC1C2 aCorrC1S1];
    pCorr=[pCorrC1C2 pCorrC1S1];
elseif contains(FigName,'FigS7')
    ConditionTitle={'Corr ColorCat C2-X> Resp Dir S1','Corr ColorCat C2-> Resp Dir C2'};
    aCorr=[aCorrC2S1 aCorrC2C2];
    pCorr=[pCorrC2S1 pCorrC2C2];
end

% generate time axis 
opts.PopulationAna.PSTHbin=AnalysisOpts.PopulationAna.PSTHbin;
opts.SpkParams.PSTH_BinShift=AnalysisOpts.SpkParams.PSTH_BinShift;
opts.SpkParams.BaselineDelay=(AnalysisOpts.ThisTimeAxisEnd-AnalysisOpts.ThisTimeAxisStart);
opts.SpkParams.PeriodLength=2*opts.SpkParams.BaselineDelay;
Col={'r','b'};

for ec=1:2 % loop on conditions
    axes(Sp(ec));
    % plot image of correlation
    [~,TimeIndX,TimeIndY]=FigParams.Image(AnalysisOpts.Time,AnalysisOpts.Time,aCorr{ec},...
        {AnalysisOpts.Xlabel;'ColorCat'},{AnalysisOpts.Xlabel;'Response'},'Correlation',[ConditionTitle{ec}],Sp(ec),'OriginLine',5,'OriginLine_width',1,...
        'WidthSmoothing',15,'WidthSmoothingDim2',15,'SmoothingMethod','','imageplotfunc','pcolor','caxis_limits','auto');

    % plot significance image superimposed on it
    FigParams.PlotSignificanceImage(pCorr{ec}(TimeIndX,TimeIndY),AnalysisOpts.Time(TimeIndX),AnalysisOpts.Time(TimeIndY),aCorr{ec}(TimeIndX,TimeIndY),AnalysisOpts.pvalEntropyAnalysis);
    
    % Time axis for projection on diagonal
    ProjTime=-(sum(TimeIndY)-1)*opts.SpkParams.PSTH_BinShift:opts.SpkParams.PSTH_BinShift:(sum(TimeIndY)-1)*opts.SpkParams.PSTH_BinShift;%(-opts.SpkParams.BaselineDelay+opts.SpkParams.PSTH_BinShift):opts.SpkParams.PSTH_BinShift:(opts.SpkParams.BaselineDelay-opts.SpkParams.PSTH_BinShift);%obj.ManData.GenerateTimeAxis('leading',opts);

    % project the mean of the correct on the diagonal and add it to the plots
    SigCorr=aCorr{ec}(TimeIndX,TimeIndY);
    [~,~,MeanProjDiagCorrectPos]=ManData.ProjectontoMatDiag(SigCorr,0);
    
    FigParams.PlotMeanStd(ProjTime,MeanProjDiagCorrectPos,[],{AnalysisOpts.Xlabel;'ColorCat'},{AnalysisOpts.Xlabel;'Response'},...
        Col{ec},1,['Correlation Timing'],'WidthSmoothing',15,'SmoothingMethod','',...
        'Sp',Sp(3),'LegendTxt',ConditionTitle{ec},'IsthisAxisTime',0);
    % fit a gaussian to the curve
    gaussianModel = fit(ProjTime', MeanProjDiagCorrectPos', 'gauss1');

    % add vertical line
    MeanTime=gaussianModel.b1;%ProjTime(find(MeanProjDiagCorrectPos>=mean(MeanProjDiagCorrectPos),1,'first'));
    FigParams.PlotVerticalLine(MeanTime,Sp(3),Col{ec},'p_line_style','--');
    Ticks=[-0.2:0.05:0.2];
    v=axis;
    text(MeanTime,v(4),num2str(MeanTime,4));
    xticks(Ticks);
    xticklabels(ManData.CovertDouble2CellStr(Ticks));
    xtickangle(45);axis tight

end