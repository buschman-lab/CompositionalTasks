% create a shuffle distribution for each neuron
function [ClassifierOpts]=CreateShuffleDist4ClassifierCorrectShuff(obj,FactorData,ClassifierOpts,FactorLevelComb,DimNum,Cond)
global AnalysisOpts
ManData=ManipulateData;
fprintf('\ngenerating shuffle distribution for dimension:%i',DimNum);
DimTxt=[0 2 3];

if ClassifierOpts.StimCongruency(DimNum)~=0
    TargetFields=obj.StimCongruencyFactorName{ClassifierOpts.StimCongruency(DimNum)};
else
    TargetFields={'none'};
end

[~,TargetFactor,ObservedFieldName,TrainCond,TestCond,DimTxt,TestFactorTxt,...
    TrainFactorTxt,LevelFieldName]=obj.getClassifierDimInfo(ClassifierOpts,DimTxt(DimNum));

TrainCond=TrainCond{Cond};
IndTrainCond=cell2mat(arrayfun(@(x) find(FactorLevelComb(:,2)==x)',TrainCond,'UniformOutput',0));


TargetFactorSet={'Color','Shape','ResponseLoc','Reward'}; % list of keywords for target factor
TargetFactorSetAlter={'ShapeCat','ColorCat','ResponseLoc','Reward'}; % list of alternate of the target factors

% now based on the situation create a shuffle distribution
if contains(TargetFactor,'Color') |  contains(TargetFactor,'Shape')

    AlternateFactor= TargetFactorSetAlter{contains(TargetFactorSet,TargetFactor(1:5))};

    if sum(contains(TargetFields,'Congruency','IgnoreCase',true))
        [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,AlternateFactor,FactorData,DimNum,IndTrainCond);
    elseif sum(contains(TargetFields,'ResponseLoc','IgnoreCase',true))
        [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,'ResponseLoc',FactorData,DimNum,IndTrainCond);
    elseif sum(contains(TargetFields,'Reward','IgnoreCase',true))
        [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,'Reward',FactorData,DimNum,IndTrainCond);
    end

elseif sum(contains(TargetFactor,'ResponseLoc','IgnoreCase',true))

    % then make it stable with respect to the feature in effect
    if TrainCond==1;AlternateFactor='ShapeCat';else;AlternateFactor='ColorCat';end
    [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,AlternateFactor,FactorData,DimNum,IndTrainCond);

elseif sum(contains(TargetFactor,'Rule','IgnoreCase',true)) % we don't balance for anything
   
    [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,'',FactorData,DimNum,IndTrainCond);

%     for ss=1:ClassifierOpts.NrepShufperFold
%         ShuffLbl(:,ss)=randsample(AllLabels,length(AllLabels));
%     end
end

% copy indexes for same recording neurons
NeuGroups=cellfun(@(x) find(strcmp(AnalysisOpts.Ch_2look_RecDate,x)),unique(AnalysisOpts.Ch_2look_RecDate),'UniformOutput',0);
for NeuG=1:length(NeuGroups) % loop through the groups and copy the shuffle for each one
    ThisNe=NeuGroups{NeuG};
    if length(ThisNe)>1
        ShuffLbl(ThisNe(2:end))=ShuffLbl(ThisNe(1));
    end
end
ClassifierOpts.ClassifierShuffleLabel{DimNum}=ShuffLbl;
ClassifierOpts.ClassifierShuffleTrialIndex{DimNum}=TrialIndex;
ClassifierOpts.ClassifierShuffleTrainCondIndex{DimNum}=TrainCondIndex;
end


function [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,AlternateFactor,FactorData,DimNum,IndTrainCond)

global AnalysisOpts

% if this a shape or color category or shape or color morphlevel change the associated level
if contains(TargetFactor,'ML')
    TargetFactor=[TargetFactor(1:5) 'Cat'];
end

AlternateFactorInd=strcmp(AnalysisOpts.factornames,AlternateFactor);
TargetFactorInd=strcmp(AnalysisOpts.factornames,TargetFactor);

NNeu=length(FactorData);
for NN=1:NNeu %loop on neurons
    FactorDataNeu=FactorData(NN);
    ThisTrainFactors=ManData.ReshapeCell2Mat(FactorDataNeu.AllFactors(IndTrainCond),62); % now we have the train factors to be shuffled
    % get trial index for each Index
    TrialIndex{NN}=cell2mat(cellfun(@(x) 1:size(x,1),FactorDataNeu.AllFactors(IndTrainCond),'uniformoutput',0));
    TrlSiz=cellfun(@(x) size(x,1),FactorDataNeu.AllFactors(IndTrainCond));
    TrainCondIndex{NN}=cell2mat(arrayfun(@(x) IndTrainCond(x)*ones(1,TrlSiz(x)),1:length(IndTrainCond),'uniformoutput',0));
    ConcatIndex{NN}=1:sum(TrlSiz);

    AllLabels=ThisTrainFactors(:,TargetFactorInd);
    AltFactor=ThisTrainFactors(:,AlternateFactorInd);

    ShuffLbl{NN}=nan*ones(length(AllLabels),ClassifierOpts.NrepShufperFold);
    % get the leve of the factors
    UniqAltFactorVals=unique(AltFactor)';
    UniqueLabels=unique(AllLabels)';

    if isempty(AlternateFactor)
        if length(UniqueLabels)~=2;error('we have more than two labels');end

        for ss=1:ClassifierOpts.NrepShufperFold
            LabelInd=randsample(length(AllLabels),length(AllLabels));
            ShuffLbl{NN}(:,ss)=AllLabels(LabelInd);
        end
    else
        if length(UniqueLabels)~=2 | length(UniqAltFactorVals)~=2;error('we have more than two labels');end
        for ss=1:ClassifierOpts.NrepShufperFold
            Ncat1=sum(AltFactor==UniqAltFactorVals(1));
            Ncat2=sum(AltFactor==UniqAltFactorVals(2));

            LabelInd1=randsample(find(AltFactor==UniqAltFactorVals(1)),Ncat1);
            LabelInd2=randsample(find(AltFactor==UniqAltFactorVals(2)),Ncat2);

            ShuffLbl{NN}(AltFactor==UniqAltFactorVals(1),ss)=AllLabels(LabelInd1);
            ShuffLbl{NN}(AltFactor==UniqAltFactorVals(2),ss)=AllLabels(LabelInd2);
        end
    end

end
end

