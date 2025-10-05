% create a shuffle distribution for each neuron
function [ClassifierOpts]=CreateShuffleDist4ClassifierCorrectShuffV3(obj,FactorData,ClassifierOpts,FactorLevelComb,DimNum,Cond,TrlRng)
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

if sum(contains(TargetFields,'SeqHist','IgnoreCase',true))
    SeqHistVal=ClassifierOpts.SeqHistCond{Cond};
else
    SeqHistVal=nan;
end

if obj.CalShuffTrlOrder % if we are shuffleing the trial order then we will need to do it within the target factor and all of the balances we are doing
    TargetFields=[TargetFields {TargetFactor}];
    [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,'TrialNum',TargetFields,FactorData,DimNum,IndTrainCond,SeqHistVal);
else
    [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,TargetFields,FactorData,DimNum,IndTrainCond,SeqHistVal);
end
% copy indexes for same recording neurons
NeuGroups=cellfun(@(x) find(strcmp(AnalysisOpts.Ch_2look_RecDate,x)),unique(AnalysisOpts.Ch_2look_RecDate),'UniformOutput',0);
for NeuG=1:length(NeuGroups) % loop through the groups and copy the shuffle for each one
    ThisNe=NeuGroups{NeuG};
    if length(ThisNe)>1
        ShuffLbl(ThisNe(2:end))=ShuffLbl(ThisNe(1));
    end
end
ClassifierOpts.ClassifierShuffleLabel{Cond}{TrlRng}{DimNum}=ShuffLbl;
ClassifierOpts.ClassifierShuffleTrialIndex{Cond}{TrlRng}{DimNum}=TrialIndex;
ClassifierOpts.ClassifierShuffleTrainCondIndex{Cond}{TrlRng}{DimNum}=TrainCondIndex;
end


function [ShuffLbl,TrialIndex,TrainCondIndex]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,TargetFields,FactorData,DimNum,IndTrainCond,SeqHistVal)

global AnalysisOpts

% if this a shape or color category or shape or color morphlevel change the associated level
% if contains(TargetFactor,'ML')
%     TargetFactor=[TargetFactor(1:5) 'Cat'];
% end
if ischar(TargetFields);TargetFields={TargetFields};end
% remove Seqhist from TargetFields
TargetFields=TargetFields(~contains(TargetFields,'SeqHist'));
nTargetFields=length(TargetFields);
TargetFieldsInd=cellfun(@(x) find(strcmp(AnalysisOpts.factornames,x)),TargetFields);
TargetFactorInd=strcmp(AnalysisOpts.factornames,TargetFactor);
SeqHistInd=strcmp(AnalysisOpts.factornames,'SeqHist');

NNeu=length(FactorData);
for NN=1:NNeu %loop on neurons
    FactorDataNeu=FactorData(NN);
    ThisTrainFactors=ManData.ReshapeCell2Mat(FactorDataNeu.AllFactors(IndTrainCond),62); % now we have the train factors to be shuffled
    % get trial index for each Index
    TrialIndex{NN}=cell2mat(cellfun(@(x) 1:size(x,1),FactorDataNeu.AllFactors(IndTrainCond),'uniformoutput',0));
    TrlSiz=cellfun(@(x) size(x,1),FactorDataNeu.AllFactors(IndTrainCond));
    TrainCondIndex{NN}=cell2mat(arrayfun(@(x) IndTrainCond(x)*ones(1,TrlSiz(x)),1:length(IndTrainCond),'uniformoutput',0));
    ConcatIndex{NN}=1:sum(TrlSiz);
   
    % if we don have seqhist then all of the trials are included
    if isnan(SeqHistVal)
        ThisSeqHistVal=-1;
        SeqhistValThisNeu=-1*ones(sum(TrlSiz),1);
    else
        SeqhistValThisNeu=ThisTrainFactors(:,SeqHistInd);
    end

    AllLabels=ThisTrainFactors(:,TargetFactorInd);
    AltFactor=arrayfun(@(x) ThisTrainFactors(:,x),TargetFieldsInd,'UniformOutput',0);

    ShuffLbl{NN}=nan*ones(length(AllLabels),ClassifierOpts.NrepShufperFold);

    % get the leve of the factors
    UniqueLabels=unique(AllLabels)';
    if nTargetFields>1
        UniqAltFactorVals=cellfun(@(x) unique(x)',AltFactor,'UniformOutput',0);
        UniqAltFactorVals=ManData.CreatPairCombsCell(UniqAltFactorVals{1},UniqAltFactorVals{2});
    else
        UniqAltFactorVals=unique(AltFactor{1});
    end

    if isempty(TargetFields)

        for ss=1:ClassifierOpts.NrepShufperFold
            LabelInd=randsample(length(AllLabels),length(AllLabels));
            ShuffLbl{NN}(:,ss)=AllLabels(LabelInd);
        end
    else
        % sample form the labels but consider the value of seqhist if it is present 
        if nTargetFields>1

            for ss=1:ClassifierOpts.NrepShufperFold
                for ll=1:size(UniqAltFactorVals,1)
                    Ncat=sum(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu==ThisSeqHistVal);

                    LabelInd=randsample(find(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu==ThisSeqHistVal),Ncat);

                    % change the labels thaat are affected
                    ShuffLbl{NN}(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu==ThisSeqHistVal,ss)=AllLabels(LabelInd);
                end
                % keep the labels that are not affected 
                ShuffLbl{NN}(isnan(ShuffLbl{NN}(:,ss)),ss)=AllLabels(isnan(ShuffLbl{NN}(:,ss)));
            end

        else

            for ss=1:ClassifierOpts.NrepShufperFold
                for ll=1:size(UniqAltFactorVals,1)
                    Ncat=sum(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal);

                    LabelInd=randsample(find(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal),Ncat);

                    % change the labels thaat are affected
                    ShuffLbl{NN}(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal,ss)=AllLabels(LabelInd);
                end
                % keep the labels that are not affected 
                ShuffLbl{NN}(isnan(ShuffLbl{NN}(:,ss)),ss)=AllLabels(isnan(ShuffLbl{NN}(:,ss)));
            end
        end
        
    end

end
end

