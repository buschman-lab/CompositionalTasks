function h1=PlotBlockInfoTransferResults(FigName,DataPath) % plots results for information trafer across blocks for Fig. 4b and Fig. S11d
% FigName is the figure name in the manuscript
% DataPath is the path to Figure Data folder

% load the data that we need
switch FigName
    case {'Fig4b1'}
        load([DataPath 'Learning3D_RuleR2SeqAnalysisLearningCtrlR2.mat']);
    case {'Fig4b2'}
        load([DataPath 'Learning3D_RuleR2SeqAnalysisLearningCtrl.mat']);
    case {'FigS8f1'}
        load([DataPath 'Learning3D_RuleR2SeqAnalysisLearningCtrlR2FEF.mat']);
    case {'FigS8f2'}
        load([DataPath 'Learning3D_RuleR2SeqAnalysisLearningCtrlFEF.mat']);
end
% Define Classes we use
FigParams=fig_params;
ManData=ManipulateData;
NeuAna=NeuralAnalysisFuncsTemp;

% plot conditions from ClassifierOpts.TrainTriple
TriplePlots=ClassifierOpts.TriplePlots;
nPlts=length(TriplePlots);
nTrlRng=length(ClassifierResults_ndD{1}(1).TrialRange);
h=FigParams.RenderFigure(1,[]);
if nTrlRng>1;nPlts=size(TriplePlots{1},1);end % we can only have one condition in this case

% plot results for a average window
[h1,Sp1]=FigParams.RenderSubplots([],[],h{1},nPlts);
NeuAna.PutSGtitle4Figure(ClassifierOpts);
NeuAna.StatTest=StatTest; % load stat test

for CondTrip=1:nPlts
    ThisTriplePlots=TriplePlots{CondTrip};
    arrayfun(@(x) NeuAna.PlotSVMPerfMetric_Learning(ClassifierResults_ndD{ThisTriplePlots(x,1)},ClassifierOpts,'Accuracy',Sp1(CondTrip),...
        ThisTriplePlots(x,2),ThisTriplePlots(x,1),ThisTriplePlots(x,1),...
        'LookatDim2',DimTxt(ThisTriplePlots(x,1)),'ExtraData',x,'ThisColor',x,'MeanStdPlotType','violin'),1,'UniformOutput',0);
    xlabel('Condition')
    xlim([0 size(ThisTriplePlots,1)+1])
    ylim([0.2 0.8])
    xticks(1:size(ThisTriplePlots,1))
    if ~isempty(ClassifierOpts.ConditionLabel)
        xticklabels(ClassifierOpts.ConditionLabel);
        xtickangle(45);
    end
end

end