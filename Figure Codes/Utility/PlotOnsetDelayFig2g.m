function h=PlotOnsetDelayFig2g(DataPath)
% FigName is the figure name in the manuscript
% DataPath is the path to Figure Data folder

% load the data that we need 
load([DataPath 'Fig2gData.mat']);
load([DataPath 'AnalysisOpts.mat']);

% Define Classes we use
FigParams=fig_params;
ManData=ManipulateData;

% create figure and subplots
h=FigParams.RenderFigure(1,[]);
[h,Sp]=FigParams.RenderSubplots([1],[1],h{1},1);

%PlotCombs3D=Tasks2Compare.Cond;
AreaNum=[1 3 4 5];
nCol=4;

for Ar=[1:5] % loop on all areas to compute the latency and stat tests
   % IndClassResults=(Ar-1)*Tasks2Compare.nMainTasks;
   
    % compare onset latency between conditions
    [LatencyAnalysis.Area(Ar).OnSetLatencySec1,LatencyAnalysis.Area(Ar).OnSetLatencySec2,...
        LatencyAnalysis.Area(Ar).pval]=arrayfun(@(x)  ManData.CompareRiseTimesPaired(MetricVals{Ar}{1}{x},...
        MetricVals{Ar}{2}{x},TimeMetricVals{Ar}{1}{x},1,'auto',1,'NBootStrap',250),1:2,'UniformOutput',0);

     [~,~,LatencyAnalysis.Area(Ar).pvalPermute]=arrayfun(@(x)  ManData.CompareRiseTimesPairedWithDiffPermutationStatTest(MetricVals{Ar}{1}{x},...
         MetricVals{Ar}{2}{x},TimeMetricVals{Ar}{1}{x},1,'auto',1),1:2,'UniformOutput',0);
end

%% plot results 
k=1;t2c=2;
for Ar=AreaNum % loop on area
    
    subplot(Sp(1));hold on
    
    LatencyDiffBetConds=LatencyAnalysis.Area(Ar).OnSetLatencySec1{t2c}'-LatencyAnalysis.Area(Ar).OnSetLatencySec2{t2c}';
    boxplot(LatencyDiffBetConds,k*ones(size(LatencyDiffBetConds,1),1),'Whisker', Inf,'positions', ...
        k,'colors',AnalysisOpts.AreaColors(Ar,:),'Widths', 0.5)
   
    ylim([-0.15 0.25]);
    axis square


    arrayfun(@(x) FigParams.AddDetailedSignificanceStar(k,LatencyAnalysis.Area(Ar).pval{x},...
        'r',Sp(1),'SigStar_fontsize',15,'AddDetailedSigTxt',1),t2c,'UniformOutput',0);

    arrayfun(@(x) FigParams.AddDetailedSignificanceStar(k,LatencyAnalysis.Area(Ar).pvalPermute{x},...
        'b',Sp(1),'SigStar_fontsize',15,'AddDetailedSigTxt',1,'SigStaeLocBias',0),t2c,'UniformOutput',0);

    k=k+1;
end
xticks(1:4)
xticklabels(AnalysisOpts.AreaNames(AreaNum))
xtickangle(45)