function  [fh,I,AllBlkTrls]=CalculateOverallPerformance(i,NSwitch,Trials,Nback,TrialInterval,I_WSLS,fh)
global Ncol N Itot LastTrialMat LastTrial PerfFit PerfFitTot TrialInfoBlk
if isempty(fh); fh = figure(10000); end
figure(fh)
set(0, 'currentfigure', fh);
hold on;

subplot(3,1,i)
sm_kern = ones(1, Nback); sm_kern = sm_kern./sum(sm_kern);
INDOBJ_AllTrials=(([Trials(TrialInterval).StopCondition] == 1 | [ Trials(TrialInterval).StopCondition] == -1));
INDOBJ_CorrectTrials=([Trials(TrialInterval).StopCondition] == 1);
CorretMat=INDOBJ_CorrectTrials(INDOBJ_AllTrials);
AllBlkTrls=Trials(TrialInterval); % output all of the trials of this block
LastTrial(i)=LastTrial(i)+length(CorretMat);
LastTrialMat{i}=[LastTrialMat{i} LastTrial(i)];
PlotInterval=LastTrialMat{i}(N(i)):(LastTrialMat{i}(N(i)+1)-Nback+1);
if ~isempty(PlotInterval) && length(PlotInterval)>3
    I=convn(CorretMat,sm_kern,'valid');
    plot(PlotInterval(1:length(I)),I,'linewidth',5); hold all;
    plot(PlotInterval(1:length(I_WSLS)),I_WSLS,':','linewidth',3); hold all;
    text(PlotInterval(1),1,['S:' num2str(NSwitch)],'FontSize',15);
    %  Perc_nonRes=PlotExtraTrials(Trials,I,TrialInterval,PlotInterval,INDOBJ_AllTrials,sm_kern) ;
    %%% now fit the logistic curve
    x=(0:length(I)-1);
    FitResult=FitLogistic(x',I');
    PerfFit.L_Sh{i}(N(i))=FitResult.L;PerfFit.ks_Sh{i}(N(i))=FitResult.k;PerfFit.x0_Sh{i}(N(i))=FitResult.x0;
    PerfFitTot.L=[PerfFitTot.L FitResult.L];
    PerfFitTot.k=[PerfFitTot.k FitResult.k];
    PerfFitTot.x0=[PerfFitTot.x0 FitResult.x0];
    PerfFitTot.Ntrl=[PerfFitTot.Ntrl length(CorretMat)];
    TrialInfoBlk=[TrialInfoBlk {TrialInterval}];
    
    plot(x+PlotInterval(1),PerfFit.L_Sh{i}(N(i))./(1+exp(-PerfFit.ks_Sh{i}(N(i))*(x-PerfFit.x0_Sh{i}(N(i))))),'--g')
    drawnow
    
    %  [Is1,Is2,Is3,Is4]=PlotIndSamples(bhv,Trials,i,TrialInterval);
    %   if Perc_nonRes<Anaopts.Perc_nonResTH
    Ncol(i)=Ncol(i)+1;
    Itot{i,Ncol(i)}=I;
    %                      Istot1{i,Ncol(i)}=Is1;Istot2{i,Ncol(i)}=Is2;Istot3{i,Ncol(i)}=Is3;Istot4{i,Ncol(i)}=Is4;
    %       end
    title(['Rule:' num2str(i)])
    ylim([0 1.2])
else
    
    PerfFitTot.L=[PerfFitTot.L nan];
    PerfFitTot.k=[PerfFitTot.k nan];
    PerfFitTot.x0=[PerfFitTot.x0 nan];
    PerfFitTot.Ntrl=[PerfFitTot.Ntrl nan];
    PerfFit.L_Sh{i}(N(i))=nan;PerfFit.ks_Sh{i}(N(i))=nan;PerfFit.x0_Sh{i}(N(i))=nan;
    TrialInfoBlk=[TrialInfoBlk {nan}];
    
end
%   InferenceAnalysisFunc(bhv,1,'Shape1')


N(i)=N(i)+1;

end