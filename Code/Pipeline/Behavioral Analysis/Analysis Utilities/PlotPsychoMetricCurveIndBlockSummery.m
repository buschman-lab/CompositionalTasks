%%% added the psychometric curve fo rthe inactive rule called Swap rue
function [AllPerf,AllPerfLevls,AllCount,fs]=PlotPsychoMetricCurveIndBlockSummery(bhv,Trials,Cond,UseRuleSwap,N,fs)

if isempty(fs{Cond}); fs{Cond} = figure(1000+Cond); end
figure(fs{Cond})
set(0, 'currentfigure', fs{Cond});
hold on;

MorphSteps=unique(diff(bhv.cuemapshape(:,2)));MorphSteps=max(MorphSteps);USEREAL50=1;RuleSwap=[3 1 1];
if UseRuleSwap==0
    StimIden=bhv.opts.MetaConditionList(Cond).S;
    RespIden=bhv.opts.MetaConditionList(Cond).L;
    Stim1Resp=RespIden((StimIden==1));
    Stim2Resp=RespIden((StimIden==2));
    Col={'r','b'};
else
    %%% see if he is responding to the other rule
    RuleSwapCond=(RuleSwap(Cond));
    StimIden =bhv.opts.MetaConditionList(RuleSwapCond).S;
    RespIden =bhv.opts.MetaConditionList(RuleSwapCond).L;
    Stim1Resp=RespIden((StimIden==1));
    Stim2Resp=RespIden((StimIden==2));
    Col={'g','c'};
end

CondTrans=[3 2 2 3]; %% whichi column to look for corresponding Morph (Prop ColorMorph ShapeMorph ShapeIden ColIden)
ii=1;
TrialInterval=[Trials.Condition] == Cond;
for MCond=[Cond]
    for M=0:MorphSteps:200
        SamplePool=find(bhv.cuemapshape(:,CondTrans(MCond))==M);
        INDOBJ_AllTrials    =(([Trials(TrialInterval).StopCondition] == 1 | [ Trials(TrialInterval).StopCondition] == -1));
        % INDOBJ_CorrectTrials= ([Trials(TrialInterval).StopCondition] == 1 );
        INDOBJ_AllTrials_RT= [Trials(TrialInterval).ReactionTime];

        AllTrialsMorphLevel=zeros(1,length(INDOBJ_AllTrials));CorrTrialsMorphLevel=zeros(1,length(INDOBJ_AllTrials));RT_ThisMorphLevel=[];
        if (M==50 || M==150) & USEREAL50
            AllTrialsStim1_MidMorph=zeros(1,length(INDOBJ_AllTrials));  AllTrialsStim2_MidMorph=zeros(1,length(INDOBJ_AllTrials));

            for S=SamplePool'
                AllTrialsWithThisSample=(([Trials(TrialInterval).StopCondition] == 1 | [ Trials(TrialInterval).StopCondition] == -1) & [Trials(TrialInterval).Sample] == S)  ;
                Stim1Trials=(([Trials(TrialInterval).Response] == Stim1Resp) & [Trials(TrialInterval).Sample] == S);
                Stim2Trials=(([Trials(TrialInterval).Response] == Stim2Resp) & [Trials(TrialInterval).Sample] == S);

                AllTrialsStim1_MidMorph=AllTrialsStim1_MidMorph | Stim1Trials;
                AllTrialsStim2_MidMorph=AllTrialsStim2_MidMorph | Stim2Trials;
                AllTrialsMorphLevel=AllTrialsWithThisSample | AllTrialsMorphLevel;
                RT_ThisMorphLevel=[RT_ThisMorphLevel INDOBJ_AllTrials_RT( AllTrialsWithThisSample)];
            end

            Stim1TrialsMorphLevel=AllTrialsStim1_MidMorph(AllTrialsMorphLevel);
            SumAll=sum(AllTrialsMorphLevel);SumStim1= sum(Stim1TrialsMorphLevel);
            Perf=SumStim1/SumAll;
            if M>50 & M<150
                Perf=1-Perf;
            end
            ALLinfo(ii,:)=[M SumAll SumStim1 Perf];
            AllSamplePool{ii}=SamplePool;
            AllRTTimes{ii}=RT_ThisMorphLevel;
        else
            for S=SamplePool'
                if UseRuleSwap
                    AllTrialsWithThisSample=(([Trials(TrialInterval).Response] == Stim1Resp | [ Trials(TrialInterval).Response] == Stim2Resp) & [Trials(TrialInterval).Sample] == S)  ;
                    if Cond==1
                        CorrectTrialsWithThisSample=(([Trials(TrialInterval).Response] == Stim1Resp) & [Trials(TrialInterval).Sample] == S & [Trials(TrialInterval).ColorIdentity] == 1)  ;
                    else
                        CorrectTrialsWithThisSample=(([Trials(TrialInterval).Response] == Stim1Resp) & [Trials(TrialInterval).Sample] == S & [Trials(TrialInterval).ObjectIdentity] == 1)  ;
                    end
                else
                    AllTrialsWithThisSample=(([Trials(TrialInterval).StopCondition] == 1 | [ Trials(TrialInterval).StopCondition] == -1) & [Trials(TrialInterval).Sample] == S)  ;
                    CorrectTrialsWithThisSample=(([Trials(TrialInterval).StopCondition] == 1) & [Trials(TrialInterval).Sample] == S)  ;
                end
                AllTrialsMorphLevel=AllTrialsWithThisSample|AllTrialsMorphLevel;
                CorrTrialsMorphLevel=CorrectTrialsWithThisSample | CorrTrialsMorphLevel;
                RT_ThisMorphLevel=[RT_ThisMorphLevel INDOBJ_AllTrials_RT(CorrectTrialsWithThisSample)];
            end
            CorrTrialsMorphLevel=CorrTrialsMorphLevel(AllTrialsMorphLevel);
            SumAll=sum(AllTrialsMorphLevel);SumCorr= sum(CorrTrialsMorphLevel);
            Perf=SumCorr/SumAll;
            if M>50 & M<150
                Perf=1-Perf;
            end

            ALLinfo(ii,:)=[M SumAll SumCorr Perf];
            AllSamplePool{ii}=SamplePool;
            AllRTTimes{ii}=RT_ThisMorphLevel;
        end

        ii=ii+1;
    end
    % plotcircularimages(AllRTTimes,ALLinfo,AllSamplePool,bhv,Cond)
    HighInd=find(ALLinfo(:,1)>=0 & ALLinfo(:,1)<=100 & ~isnan(ALLinfo(:,end)));LowInd=find(ALLinfo(:,1)>=100 & ALLinfo(:,1)<=200 & ~isnan(ALLinfo(:,end)));

    LowInd=[1;LowInd(end:-1:1)]; %% adding the 0 percent
    HighPerf=1-ALLinfo(HighInd,end);HighLabl=ALLinfo(HighInd,1);HighCnt=ALLinfo(HighInd,2);
    LowPerf=1-ALLinfo(LowInd,end);LowLabl=ALLinfo(LowInd,1);LowCnt=ALLinfo(LowInd,2);
    %%%% fit both psychometric curves
    %           ColPars_High=mle_fit_psycho([0.01*HighLabl' ;[ALLinfo(HighInd,2)]' ; HighPerf'],'erf_psycho' );
    %           ColPars_Low=mle_fit_psycho([0.01*HighLabl' ;[ALLinfo(LowInd,2)]' ; LowPerf'],'erf_psycho' );

    %figure(500+Cond);
    title(['Rule'  num2str(Cond)])
    subplot(2,N(1),N(2))
    hold on
    %  plot(1:length(HighPerf),HighPerf,'LineWidth',2,'Color',Col{1});
    XTICKS=0:1/(length(HighPerf)-1):1;
    plot(XTICKS,HighPerf,'LineWidth',2,'Color',Col{1});
    arrayfun(@(x) text(XTICKS(x),HighPerf(x)+0.01,num2str(HighCnt(x))),1:length(XTICKS))
    %           plot([0:0.01:1],erf_psycho(ColPars_High,[0:0.01:1]),'--','LineWidth',2)

    % plot(LowPerf,'LineWidth',2,'Color','r');a
    set(gca,'XTick',0:1/(length(HighPerf)-1):1);set(gca,'XTickLabel',{HighLabl});
    %     for i=1:length(HighPerf);text(XTICKS(i),HighPerf(i)+0.05,num2str(ALLinfo(HighInd(i),2)));end
    if isfield(bhv,'PerfTH_30')
        plot([0 1],[bhv.PerfTH_30 bhv.PerfTH_30],'k')
        plot([0 1],[1-bhv.PerfTH_30 1-bhv.PerfTH_30],'k')

        plot([0 1],[bhv.PerfTH_0 bhv.PerfTH_0],'g')
        plot([0 1],[1-bhv.PerfTH_0 1-bhv.PerfTH_0],'g')
    end
    %    title(['Psychometic Curve first Half Rule:' num2str(Cond)])
    ylim([0 1])
    subplot(2,N(1),N(2)+N(1))
    hold on
    XTICKS=0:1/(length(LowPerf)-1):1;
    plot(XTICKS,LowPerf,'LineWidth',2,'Color',Col{2});
    arrayfun(@(x) text(XTICKS(x),LowPerf(x)+0.01,num2str(LowCnt(x))),1:length(XTICKS))

    %            plot([0:0.01:1],erf_psycho(ColPars_Low,[0:0.01:1]),'--','LineWidth',2)

    set(gca,'XTick',0:1/(length(LowPerf)-1):1);set(gca,'XTickLabel',{LowLabl});
    %   for i=1:length(LowPerf);text(XTICKS(i),LowPerf(i)+0.05,num2str(ALLinfo(LowInd(i),2)));end
    if isfield(bhv,'PerfTH_30')

        plot([0 1],[bhv.PerfTH_30 bhv.PerfTH_30],'k')
        plot([0 1],[1-bhv.PerfTH_30 1-bhv.PerfTH_30],'k')

        plot([0 1],[bhv.PerfTH_0 bhv.PerfTH_0],'g')
        plot([0 1],[1-bhv.PerfTH_0 1-bhv.PerfTH_0],'g')
    end
    %    title(['Psychometic Curve Second Half Rule:' num2str(Cond)])
    ylim([0 1])

end
%               HighMorphPerf.Perf=HighPerf;HighMorphPerf.Labl=HighLabl;
%               LowMorphPerf.Perf=LowPerf;LowMorphPerf.Labl=LowLabl;
HighMorphPerf=HighPerf;
LowMorphPerf=LowPerf;
title(['Rule'  num2str(Cond)])
set(gcf,'Name',['Psychometric curve for Rule: ' num2str(Cond)])

%% organize performances for output
AllHighPerfLevls=[0 20   30  40 50   60  70  80 100];
AllLowPerfLevls= [0 180 170 160 150 140 130 120 100];
AllPerfLevls=[AllHighPerfLevls AllLowPerfLevls];
AllPerf=NaN*ones(1,length(AllPerfLevls));
AllCount=NaN*ones(1,length(AllPerfLevls));

AnimalPerf=[LowPerf;HighPerf];AnimalPerfLvl=[LowLabl;HighLabl];
AnimalCount=[LowCnt;HighCnt];
k=1;
for i=AllPerfLevls
    if sum(AnimalPerfLvl==i)
        AllPerf(k)=unique(AnimalPerf(AnimalPerfLvl==i));
        AllCount(k)=unique(AnimalCount(AnimalPerfLvl==i));
    end
    k=k+1;
end


end