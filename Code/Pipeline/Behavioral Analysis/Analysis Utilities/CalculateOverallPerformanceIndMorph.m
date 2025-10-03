% calculates overall performance of individual morph levels during
% learning. we should indicate which morphlevels we are intresting to look

function  [fh,I]=CalculateOverallPerformanceIndMorph(bhv,i,NSwitch,Trials,TrialInterval,NTrls2Look,fh)
global Ncol N Itot Anaopts LastTrialMat LastTrial PerfFit PerfFitTot TrialInfoBlk
if isempty(fh); fh = figure(10000); end
figure(fh)
set(0, 'currentfigure', fh);
hold on;
subplot(3,1,i)

StimIden=bhv.opts.MetaConditionList(Cond).S;
RespIden=bhv.opts.MetaConditionList(Cond).L;
Stim1Resp=RespIden((StimIden==1));
Stim2Resp=RespIden((StimIden==2));
Col={'r','b'};
CondTrans=[3 2 2 3]; %% whichi column to look for corresponding Morph (Prop ColorMorph ShapeMorph ShapeIden ColIden)
ii=1;
TrialInterval=[Trials.Condition] == Cond;
for M=0:MorphSteps:200
    SamplePool=find(bhv.cuemapshape(:,CondTrans(MCond))==M);
    INDOBJ_AllTrials    =(([Trials(TrialInterval).StopCondition] == 1 | [ Trials(TrialInterval).StopCondition] == -1));
  % INDOBJ_CorrectTrials= ([Trials(TrialInterval).StopCondition] == 1 );
    INDOBJ_AllTrials_RT= [Trials(TrialInterval).ReactionTime];

    AllTrialsMorphLevel=zeros(1,length(INDOBJ_AllTrials));CorrTrialsMorphLevel=zeros(1,length(INDOBJ_AllTrials));RT_ThisMorphLevel=[];

    for S=SamplePool'                         
          AllTrialsWithThisSample=(([Trials(TrialInterval).StopCondition] == 1 | [ Trials(TrialInterval).StopCondition] == -1) & [Trials(TrialInterval).Sample] == S)  ;
          CorrectTrialsWithThisSample=(([Trials(TrialInterval).StopCondition] == 1) & [Trials(TrialInterval).Sample] == S)  ;                          
          AllTrialsMorphLevel=AllTrialsWithThisSample|AllTrialsMorphLevel;
          CorrTrialsMorphLevel=CorrectTrialsWithThisSample | CorrTrialsMorphLevel;
    end
    CorrTrialsMorphLevel=CorrTrialsMorphLevel(AllTrialsMorphLevel);

end









sm_kern = ones(1, Nback); sm_kern = sm_kern./sum(sm_kern);
INDOBJ_AllTrials=(([Trials(TrialInterval).StopCondition] == 1 | [ Trials(TrialInterval).StopCondition] == -1));
INDOBJ_CorrectTrials=([Trials(TrialInterval).StopCondition] == 1);
CorretMat=INDOBJ_CorrectTrials(INDOBJ_AllTrials);
LastTrial(i)=LastTrial(i)+length(CorretMat);
LastTrialMat{i}=[LastTrialMat{i} LastTrial(i)];
PlotInterval=LastTrialMat{i}(N(i)):(LastTrialMat{i}(N(i)+1)-Nback+1);
I=convn(CorretMat,sm_kern,'valid');
plot(PlotInterval(1:length(I)),I,'linewidth',5); hold all;  
text(PlotInterval(1),1,['S:' num2str(NSwitch)],'FontSize',15);

title(['Rule:' num2str(i)])
ylim([0 1.2])

N(i)=N(i)+1;

end