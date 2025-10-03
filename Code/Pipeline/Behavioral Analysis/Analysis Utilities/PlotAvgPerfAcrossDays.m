function [varargout] = PlotAvgPerfAcrossDays(Perf)
%PLOTAVGPERFACROSSDAYS plots average performance across days
% Perf: comoated in pipeline BehaviorAnalisis pipelines

% define functions we need
fp=fig_params;
bhvAna=BhvAnalysisFuncs;


% start with concatinating all of the imporant info across days

[AllPSMPerf,AllTrlPerf,IndSamp,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist] =...
    bhvAna.CancatinateInfoDays(Perf);

% plot psychometric average psychometric curve across days 
 AllPerfFig=cell(1,6);
 %[AllPerfFig{1:2}]=bhvAna.PlotAvgPSM(AllPSMPerf); % 2figs
 [AllPerfFig{3:4}]=bhvAna.PlotTrlPerf(AllTrlPerf,AllTrlCount,AllTrlCountDay,NBlocksDay,RewardPulse,NCorrectTrl,NumRewards,AllSeqHist); %2 figs
 %[AllPerfFig{5}]=bhvAna.PlotSampPerfInfo(IndSamp);%bhvAna.PlotSampInfo(IndSamp); %1 fig
 %[AllPerfFig{6}]=bhvAna.PlotAvgPSMCongInCong(IndSamp,AllPSMPerf);


varargout=AllPerfFig;
end

