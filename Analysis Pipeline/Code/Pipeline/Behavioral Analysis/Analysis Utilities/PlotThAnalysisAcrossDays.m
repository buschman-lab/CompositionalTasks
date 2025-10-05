function [varargout] = PlotThAnalysisAcrossDays(Perf)
%PLOTAVGPERFACROSSDAYS plots average performance across days with respect
%to the simulated threshold for performance

% Perf: comoated in pipeline BehaviorAnalisis pipelines

% define functions we need
close all
fp=fig_params;
bhvAna=BhvAnalysisFuncs;
Thvals=[150:50:300];
nThvals=length(Thvals);
ColSet=distinguishable_colors(nThvals);
fh= fp.RenderFigure(4+nThvals,[]);hold all
h1=cell(1,3);h2=cell(1,3);
for i=Thvals
    n=find(i==Thvals);
%     if i==75
%         load(['Trl103_IntTrl' num2str(i) '_Data.mat']);
%     elseif i==200
        load(['Trl' num2str(i) '_IntTrl75_Data.mat']);
%     end
    %% start with concatinating all of the imporant info across days
    [AllPSMPerf,AllTrlPerf,IndSamp,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards] = bhvAna.CancatinateInfoDays(Perf);
    Col=ColSet(i==Thvals,:);
    %% plot psychometric average psychometric curve across days
    % creat a plot with labels and color codes 
    figure(fh{1});hold on
    v=axis;
    plot([v(1) v(2)],[n/2 n/2],'color',Col,'LineWidth',15)
    text((v(1)+v(2))/2,n/2+0.2,['Number of Integrated Trials ' num2str(i)],'FontSize',15);axis off
    AllPerfFig{1}=fh{1};
    ylim([0 length(Thvals)/2+1])
    
    [AllPerfFig{2},h1,h2]=bhvAna.PlotAvgPSM_ThAnalysis(i,AllPSMPerf,Col,fh(2),h1,h2); % 1figs
    [AllPerfFig{3:4}]=bhvAna.PlotTrlPerf_ThAnalysis(Col,n,fh(3:4),AllTrlPerf,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards); %2 figs
    [AllPerfFig{4+n}]=bhvAna.PlotAvgPSMCongInCong_ThAnalysis(fh(4+n),i,Col,IndSamp,AllPSMPerf);
end
FigSaveFileName=['TrialThresholdAnalysis75Trl'];
fp.SaveFigSeries(FigSaveFileName,'',[AllPerfFig],'SaveEachFrame',1)
varargout=AllPerfFig;
end

