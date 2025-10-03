function h1=PlotGLMmixedselectivityResults(FigName,DataPath) % plots results  for Fig.1h-m and Fig. S3
% FigName is the figure name in the manuscript
% DataPath is the path to Figure Data folder

load([DataPath 'GLMdata.mat']);

% load classes
FigParams=fig_params;
NeuAna=NeuralAnalysisFuncsTemp;

% plot comparisions for each factor
[h1]=FigParams.RenderFigure(2,[]);
[h1{1},Sp1]=FigParams.RenderSubplots([2],[5],h1{1},[]);


%BetaCPD,CPDpval,CPDcluster,IndNeuSigFactor,CPDFieldNames

BetaCPD=CPDsubshAll;
CPDpval=[CPDpval{:}];
CPDcluster=[CPDcluster{:}];
IndNeuSigFactor=[IndNeuSigFactor{:}];
CPDFieldNames=CPDFieldNames{1};

Factors2Look={'ResponseLoc','Reward','Rule','ColorCat','ShapeCat'};
Factors2LookTiming={[0.25,0.4],[0.4,0.6],[-0.2,0],[0.05,0.3],[0.05,0.3]}; % timing of averaging for each of these factors
nFactors=length(Factors2Look);
GetTimeInd=@(x) AnalysisOpts.Time>=x(1) & AnalysisOpts.Time<=x(2);
NNeu=size(BetaCPD.ResponseLoc,1);
% pick pairs of these factors and plot the scatter plot of selectivity
Pairs=nchoosek(1:nFactors,2); PairCol=distinguishable_colors(size(Pairs,1));%obj.FigParams.getColorPalet(size(Pairs,1));
for i=1:size(Pairs,1)
    Fact1Name=Factors2Look{Pairs(i,1)};Fact2Name=Factors2Look{Pairs(i,2)};
    FactNames=[Factors2Look(Pairs(i,1)),Factors2Look(Pairs(i,2))];
    Fact1NameInd=strcmp(CPDFieldNames,Fact1Name);
    Fact2NameInd=strcmp(CPDFieldNames,Fact2Name);
    % go through each neuron and average the CPD where there is significant encoding of that variable
    clear Fact
    for p=1:2
        for n=1:NNeu
            temp=[];
            for cl=1:size(CPDcluster(n).(FactNames{p}),1)
                SharedTimes=intersect(find(GetTimeInd(Factors2LookTiming{Pairs(i,p)})),CPDcluster(n).(FactNames{p}){cl}') ;
                if ~isempty(SharedTimes)
                    if CPDpval(n).(FactNames{p}){cl}(1)<0.05
                        temp=[temp BetaCPD.(FactNames{p})(n,SharedTimes)];
                    end
                end
            end
            if ~isempty(temp)
                Fact{p}(n)=nanmean(temp);
            else
                Fact{p}(n)=0;
            end
        end
    end

    Fact1=Fact{1}';
    Fact2=Fact{2}';
    FactConjct=IndNeuSigFactor(Fact1NameInd,:) | IndNeuSigFactor(Fact2NameInd,:);
    Fact1=Fact1(FactConjct);Fact2=Fact2(FactConjct);
    Fact1=Fact1/nanstd(Fact1);Fact2=Fact2/nanstd(Fact2);

    LegTxt{i}=[Fact1Name '-' Fact2Name];

    subplot(Sp1(i))
    %subplot(2,5,i);
    hold on

    Fact1Count=sum(IndNeuSigFactor(Fact1NameInd,:),2);
    Fact2Count=sum(IndNeuSigFactor(Fact2NameInd,:),2);
    FactConjCount=sum(IndNeuSigFactor(Fact1NameInd,:) & IndNeuSigFactor(Fact2NameInd,:),2);
    RatioJoinConj=(FactConjCount)/min(Fact1Count,Fact2Count);
    % do statistical test
    Pstat=PopOverlapPermutStatTest(Fact1Count,Fact2Count,FactConjCount,NNeu);
    bar(1,Fact1Count/NNeu,'Facecolor',PairCol(i,:))
    bar(2,Fact2Count/NNeu,'Facecolor',PairCol(i,:))
    bar(3,FactConjCount/NNeu,'Facecolor',PairCol(i,:))
    ylabel('Proportion of Neurons')
  
    xticks(1:3);
    xticklabels({Fact1Name,Fact2Name,[Fact1Name '-' Fact2Name]});
    xtickangle(45);
    title(sprintf('%s, p=%0.4f',LegTxt{i},Pstat));
    axis square
end