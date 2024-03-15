for trlind=1:18
IndTest=find(sum(cell2mat(arrayfun(@(i) (FactorDataTest(Neu).AllFactors{CondInds(x)}(:,TrialNumInd)==FactorData(Neu).AllFactors{CondInds(x)}(i,TrialNumInd) & ...
FactorDataTest(Neu).AllFactors{CondInds(x)}(:,BlkOrderInd)==FactorData(Neu).AllFactors{CondInds(x)}(i,BlkOrderInd)),TrainStimInds{Neu}{x}(trlind),'UniformOutput',0)),2))

out(trlind,:)=FactorDataTest(Neu).AllFactors{CondInds(x)}(IndTest,:)==FactorData(Neu).AllFactors{CondInds(x)}(TrainStimInds{Neu}{x}(trlind),:);
end


for i=1:151
   IndCommon{i}=find(FactorData(1).AllFactors{1}(:,39)==FactorDataTest(1).AllFactors{1}(i,39) & FactorData(1).AllFactors{1}(:,41)==FactorDataTest(1).AllFactors{1}(i,41));
end

hold on
for i=1:25;Acc(i,:)=cell2mat(cellfun(@(x) x.Accurary,ClassifierResults_Shuffled(1).TrialRange(1).Rep(i).Observed_3ndD,'UniformOutput',0));
end
plot(mean(Acc,1))



hold on
for i=1:16;Acc(i,:)=cell2mat(cellfun(@(x) x.Accurary,ClassifierResults_Shuffled(1).TrialRange(1).Rep(i).Observed,'UniformOutput',0));
end
plot(mean(Acc,1))


for ll=1:size(UniqAltFactorVals,1)
    for A=1:length(UniqueLabels)
       Labelcount(A,ll)=sum(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal & AllLabels==UniqueLabels(A));
       LabelShuff(A,ll)=sum(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal & ShuffLbl{NN}==UniqueLabels(A));
    end
end

