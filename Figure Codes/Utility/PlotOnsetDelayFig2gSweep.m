function h=PlotOnsetDelayFig2gSweep(DataPath)
% FigName is the figure name in the manuscript
% DataPath is the path to Figure Data folder

% load the data that we need 
load([DataPath 'Fig2gData.mat']);
load([DataPath 'AnalysisOpts.mat']);

% Define Classes we use
FigParams=fig_params;
ManData=ManipulateData;


%PlotCombs3D=Tasks2Compare.Cond;
AreaNum=[1 2 3 4 5];
nCol=4;
nb=1;
NBootStrapSet=[25:25:250];
for NBootStrap=NBootStrapSet
    for Ar=[1:5] % loop on all areas to compute the latency and stat tests
        % compare onset latency between conditions
        [LatencyAnalysis.Area(Ar,nb).OnSetLatencySec1,LatencyAnalysis.Area(Ar,nb).OnSetLatencySec2,...
            LatencyAnalysis.Area(Ar,nb).pval]=arrayfun(@(x)  ManData.CompareRiseTimesPaired(MetricVals{Ar}{1}{x},...
            MetricVals{Ar}{2}{x},TimeMetricVals{Ar}{1}{x},1,'auto',1,'NBootStrap',NBootStrap),1:2,'UniformOutput',0);

        [~,~,LatencyAnalysis.Area(Ar,nb).pvalPermute]=arrayfun(@(x)  ManData.CompareRiseTimesPairedWithDiffPermutationStatTest(MetricVals{Ar}{1}{x},...
             MetricVals{Ar}{2}{x},TimeMetricVals{Ar}{1}{x},1,'auto',1),1:2,'UniformOutput',0);
    end
    fprintf('\n%i',NBootStrap)
    nb=nb+1;
end

%% plot results 
% create figure and subplots
h=FigParams.RenderFigure(1,[]);
[h,Sp]=FigParams.RenderSubplots([2],[3],h{1},6);

t2c=2;
for Ar=AreaNum % loop on area
    for nb=1:length(NBootStrapSet)

        subplot(Sp(Ar));hold on

        LatencyDiffBetConds=LatencyAnalysis.Area(Ar,nb).OnSetLatencySec1{t2c}'-LatencyAnalysis.Area(Ar,nb).OnSetLatencySec2{t2c}';
        boxplot(LatencyDiffBetConds,nb*ones(size(LatencyDiffBetConds,1),1),'Whisker', Inf,'positions', ...
            nb,'colors',AnalysisOpts.AreaColors(Ar,:),'Widths', 0.5)

        ylim([-0.2 0.5]);
        axis square
        % pval bootstrap
        arrayfun(@(x) FigParams.AddDetailedSignificanceStar(nb,LatencyAnalysis.Area(Ar,nb).pval{x},...
            'r',Sp((Ar)),'SigStar_fontsize',5,'AddDetailedSigTxt',1,'SigStaeLocBias',-0.05),t2c,'UniformOutput',0);
        % pval permute
        arrayfun(@(x) FigParams.AddDetailedSignificanceStar(nb,LatencyAnalysis.Area(Ar,nb).pvalPermute{x},...
            'b',Sp((Ar)),'SigStar_fontsize',5,'AddDetailedSigTxt',1,'SigStaeLocBias',0),t2c,'UniformOutput',0);
    end
    xticks(1:length(NBootStrapSet))
    xticklabels([NBootStrapSet])
    xtickangle(45)
    title(AnalysisOpts.AreaNames{Ar})
    xlabel('Bootstrap size')
    ylabel('onset delay')
end
