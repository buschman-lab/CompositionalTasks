function Perc_NoRES=PlotExtraTrials(Trials,I,TrialInterval,PlotInterval,INDOBJ_AllTrials,sm_kern)
if 0
     A=([Trials(TrialInterval).ColorIdentity] == 1 & [Trials(TrialInterval).Response] == 2 & [Trials(TrialInterval).Condition] == 1) | ([Trials(TrialInterval).ColorIdentity] == 2 & [Trials(TrialInterval).Response] == 1 & [Trials(TrialInterval).Condition] == 1);
     B=([Trials(TrialInterval).ColorIdentity] == 1 & [Trials(TrialInterval).Response] == 1 & [Trials(TrialInterval).Condition] == 1) | ([Trials(TrialInterval).ColorIdentity] == 2 & [Trials(TrialInterval).Response] == 2 & [Trials(TrialInterval).Condition] == 1);
     plot(PlotInterval(1:length(I)),convn(A(INDOBJ_AllTrials), sm_kern, 'valid'),'r-'); hold all; 
     plot(PlotInterval(1:length(I)),convn(B(INDOBJ_AllTrials), sm_kern, 'valid'),'k-'); hold all; 
     C=([Trials(TrialInterval).ColorIdentity] == 1 & [Trials(TrialInterval).Response] == 3 & [Trials(TrialInterval).Condition] == 1) | ([Trials(TrialInterval).ColorIdentity] == 2 & [Trials(TrialInterval).Response] == 4 & [Trials(TrialInterval).Condition] == 1);
     D=([Trials(TrialInterval).ColorIdentity] == 1 & [Trials(TrialInterval).Response] == 4 & [Trials(TrialInterval).Condition] == 1) | ([Trials(TrialInterval).ColorIdentity] == 2 & [Trials(TrialInterval).Response] == 3 & [Trials(TrialInterval).Condition] == 1);
     plot(PlotInterval(1:length(I)),convn(C(INDOBJ_AllTrials), sm_kern, 'valid'),'ro'); hold all; 
     plot(PlotInterval(1:length(I)),convn(D(INDOBJ_AllTrials), sm_kern, 'valid'),'ko'); hold all;
    %    plot(PlotInterval(1:length(I)),1.1*ones(1,length(I)),'LineWidth',10,'color',Trials(TrialInterval(1)).BackgroundColor)
     legend('Perf','Col12D','Col12R','Col34D','Col34R')
end
     INDOBJ_NoRES=([Trials(TrialInterval).StopCondition] == -2 | [Trials(TrialInterval).StopCondition] == -3 | [Trials(TrialInterval).StopCondition] == -5 | [Trials(TrialInterval).StopCondition] == -9);
   %  Ic=convn(INDOBJ_NoRES, sm_kern, 'valid');
    % plot(PlotInterval(1):PlotInterval(1)+length(Ic)-1,Ic,'linewidth',2); hold all;  
     Perc_NoRES=sum(INDOBJ_NoRES)/length(TrialInterval);      
     plot(PlotInterval(1),Perc_NoRES,'*','MarkerSize',10)
     
end
 