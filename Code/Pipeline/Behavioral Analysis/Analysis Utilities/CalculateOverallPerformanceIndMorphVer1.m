function  [fh]=CalculateOverallPerformanceIndMorphVer1(i,NSwitch,Trials,Nback,TrialInterval,fh)
global   N_Ind   LastTrialMat_Ind LastTrial_Ind
if isempty(fh); fh = figure(20000); end
figure(fh)
set(0, 'currentfigure', fh);
hold on;
LookMorph=[0 30 70 100 130 170];
NMorph=length(LookMorph);
Col=distinguishable_colors(NMorph);
subplot(3,1,i);hold on

if i==1
    Feature='ObjectMorphLevel';
else
    Feature='ColorMorphLevel';
end

k=1;
AttemptedTrials=[Trials(TrialInterval).StopCondition] >=-1;
for ML=LookMorph
    INDOBJ_AllTrials=double(([Trials(TrialInterval).StopCondition] == 1 | [Trials(TrialInterval).StopCondition] == -1) & [Trials(TrialInterval).(Feature)] == ML);
    INDOBJ_AllTrials(INDOBJ_AllTrials==0)=NaN;
    INDOBJ_CorrectTrials=([Trials(TrialInterval).StopCondition] == 1);
    CorrectMat=INDOBJ_AllTrials.*INDOBJ_CorrectTrials;
    CorrectMat=CorrectMat(AttemptedTrials);
    MeanCorrectMat{k}=movmean(CorrectMat,Nback,2,'omitnan','Endpoints','discard');
    k=k+1;
end

LastTrial_Ind(i)=LastTrial_Ind(i)+length(CorrectMat);
LastTrialMat_Ind{i}=[LastTrialMat_Ind{i} LastTrial_Ind(i)];
PlotInterval=LastTrialMat_Ind{i}(N_Ind(i)):(LastTrialMat_Ind{i}(N_Ind(i)+1)-Nback+1);
if ~isempty(PlotInterval) && length(PlotInterval)>0
    arrayfun(@(x) plot(PlotInterval(1:length(MeanCorrectMat{x})),MeanCorrectMat{x},'linewidth',2,'color',Col(x,:)),1:NMorph);
    LegendTxt=arrayfun(@(x) [num2str(x)],LookMorph,'UniformOutput',0);
    text(PlotInterval(1),1,['S:' num2str(NSwitch)],'FontSize',15);
    title(['Rule:' num2str(i)])
    ylim([0 1.2])
    legend(LegendTxt,'Location','northeastoutside')
end
v=axis;
plot([v(1) v(2)],[0.7 0.7],'--k')
plot([PlotInterval(end) PlotInterval(end)],[v(3) v(4)],'--k');
N_Ind(i)=N_Ind(i)+1;
end