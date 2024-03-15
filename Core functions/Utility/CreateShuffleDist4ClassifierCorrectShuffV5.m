% create a shuffle distribution for each neuron
% (V4 generate one shuffle sample for each neuron, this is to avoid savign a shuffle file)
function [ClassifierOpts]=CreateShuffleDist4ClassifierCorrectShuffV5(obj,FactorData,ClassifierOpts,FactorLevelComb,DimNum,Cond,TrlRng)
global AnalysisOpts

ManData=ManipulateData;
%fprintf('\ngenerating shuffle distribution for dimension:%i',DimNum);
DimTxt=[0 2 3];
CurrentNrepShufperFold=ClassifierOpts.NrepShufperFold;
ClassifierOpts.NrepShufperFold=1; % we just need one for this instance


%% detetimene the target fields
if ClassifierOpts.StimCongruency(DimNum)~=0
    TargetFields=obj.StimCongruencyFactorName{ClassifierOpts.StimCongruency(DimNum)};
else
    TargetFields={'none'};
end

[~,TargetFactor,ObservedFieldName,TrainCond,TestCond,DimTxt,TestFactorTxt,...
    TrainFactorTxt,LevelFieldName]=obj.getClassifierDimInfo(ClassifierOpts,DimTxt(DimNum));

%% if we have rule condition we need to add the first factor as target field
if ((strcmp(TargetFactor,'Rule') & sum(contains(TargetFields,'Congruency')) & ...
        (strcmp(ClassifierOpts.TargetFactors{1},'ShapeCat') | ...
        strcmp(ClassifierOpts.TargetFactors{1},'ColorCat')) & ...
        ClassifierOpts.StimCongruency(DimNum)==10) | ...
        (strcmp(TargetFactor,'Rule') & strcmp(ClassifierOpts.TargetFactors{1},'ResponseLoc') & ClassifierOpts.StimCongruency(DimNum)==13)) & ...
        obj.CalShuffTrlOrder

    TargetFields=[TargetFields, ClassifierOpts.TargetFactors{1}];
end

%% limit trials based on SeqHist
if sum(contains(TargetFields,'SeqHist','IgnoreCase',true))
    SeqHistVal=ClassifierOpts.SeqHistCond{Cond};
else
    SeqHistVal=nan;
end
% limit based on the trial number of the train or test 
[TrainTrlRange,TestTrlRange,~,~]=obj.GetTrialRangeforThisCond(ClassifierOpts,Cond);
if obj.CalShuffTrlOrder
    TrialRngInc=unique(cell2mat(TestTrlRange));
else
    TrialRngInc=unique(cell2mat(TrainTrlRange));
end

%% Limit trials based on reward
[TrlTypeVal_Train,TrlTypeVal_Test]=AdjustLimitFactorCondition(ClassifierOpts,DimNum,Cond);

if obj.CalShuffTrlOrder % if we are shuffleing the trial order then we will need to do it within the target factor and all of the balances we are doing
    TargetFields=[TargetFields {TargetFactor}];
    TestCond=TestCond{Cond};
    IndTestCond=cell2mat(arrayfun(@(x) find(FactorLevelComb(:,2)==x)',TestCond,'UniformOutput',0));
    ShuffleTrlSet='Test'; % we are shuffling the test set
    [ShuffLbl,TrialIndex,TrainCondIndex,BlkOrderShuff,RevTrlNumShuff]=SamplefromFactor(ManData,ClassifierOpts,'TrialNum',TargetFields,FactorData,DimNum,IndTestCond,SeqHistVal,Cond,TrlTypeVal_Test,TrialRngInc,ShuffleTrlSet);
else
    FactorData=SkipSeqHistTrain(FactorData,ClassifierOpts,Cond); % if ClassifierOpts.SkipSeqHistTrain=1; then we don't apply SeqHist to train data
    TrainCond=TrainCond{Cond};
    IndTrainCond=cell2mat(arrayfun(@(x) find(FactorLevelComb(:,2)==x)',TrainCond,'UniformOutput',0));
    % if this is category then adjust it based on the rule
    if sum(strcmp(TargetFields,'Cat'))
        ThisFactorRuleBased=AnalysisOpts.RuleMainFeature{TrainCond};
        TargetFields{strcmp(TargetFields,'Cat')}=ThisFactorRuleBased;
    end
    ShuffleTrlSet='Train'; % we are shuffling the train set
    [ShuffLbl,TrialIndex,TrainCondIndex,BlkOrderShuff,RevTrlNumShuff]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,TargetFields,FactorData,DimNum,IndTrainCond,SeqHistVal,Cond,TrlTypeVal_Train,TrialRngInc,ShuffleTrlSet);
end

%% copy indexes for same recording neurons
NeuGroups=cellfun(@(x) find(strcmp(AnalysisOpts.Ch_2look_RecDate,x)),unique(AnalysisOpts.Ch_2look_RecDate),'UniformOutput',0);
for NeuG=1:length(NeuGroups) % loop through the groups and copy the shuffle for each one
    ThisNe=NeuGroups{NeuG};
    if length(ThisNe)>1
        ShuffLbl(ThisNe(2:end))=ShuffLbl(ThisNe(1));
        BlkOrderShuff(ThisNe(2:end))=BlkOrderShuff(ThisNe(1));
        RevTrlNumShuff(ThisNe(2:end))=RevTrlNumShuff(ThisNe(1));
    end
end

%% save the varibales in ClassifierOpts
ClassifierOpts.NrepShufperFold=CurrentNrepShufperFold;
ClassifierOpts.ClassifierShuffleLabel{DimNum}=ShuffLbl;
ClassifierOpts.ClassifierShuffleTrialIndex{DimNum}=TrialIndex;
ClassifierOpts.ClassifierShuffleTrainCondIndex{DimNum}=TrainCondIndex;
ClassifierOpts.ClassifierShuffleTrainCondBlkOrder{DimNum}=BlkOrderShuff;
ClassifierOpts.ClassifierShuffleRevTrlNum{DimNum}=RevTrlNumShuff;
end

function [ShuffLbl,TrialIndex,TrainCondIndex,BlkOrderShuff,RevTrlNumShuff]=SamplefromFactor(ManData,ClassifierOpts,TargetFactor,TargetFields,FactorData,DimNum,IndTrainCond,SeqHistVal,Cond,TrlTypeVal,TrialRngInc,ShuffleTrlSet)

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
FromSwitchBhvPerfInd=find(strcmp(AnalysisOpts.factornames,'FromSwitchBhvPerf'));
BlkOrderInd=find(strcmp(AnalysisOpts.factornames,'BlkOrder'));
RewardInd=strcmp(AnalysisOpts.factornames,'Reward');
RevTrlNumInd=strcmp(AnalysisOpts.factornames,'TrialNumReverse');
TrlNumInd=strcmp(AnalysisOpts.factornames,'TrialNum');

if sum(TrialRngInc<0);FactorIndTrlNum=RevTrlNumInd;else;FactorIndTrlNum=TrlNumInd;end


NNeu=length(FactorData);
for NN=1:NNeu %loop on neurons
    FactorDataNeu=FactorData(NN);
    if isfield(ClassifierOpts,'LimitFromSwitchPerf')
        if strcmp(ShuffleTrlSet,'Test')
            IncludedTrials=LimitTrialsBasedonFactor(FactorDataNeu,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(2),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{2});
        elseif strcmp(ShuffleTrlSet,'Train')
            IncludedTrials=LimitTrialsBasedonFactor(FactorDataNeu,'FromSwitchBhvPerf',ClassifierOpts.LimitFromSwitchPerf{Cond}(1),ClassifierOpts.LimitFromSwitchPerf_Operation{Cond}{1});
        end
    else
        IncludedTrials=[];
    end

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
        ThisSeqHistVal=SeqHistVal;
        SeqhistValThisNeu=ThisTrainFactors(:,SeqHistInd);
    end

    % are we limiting the trial type?
    if isnan(TrlTypeVal)
        ThisTrlTypeVal=-1;
        TrlTypeValThisNeu=-1*ones(sum(TrlSiz),1);
    else
        ThisTrlTypeVal=TrlTypeVal;
        TrlTypeValThisNeu=ThisTrainFactors(:,RewardInd);
    end

    % if we are limiting any factor
    if ~isempty(IncludedTrials)
        ThisIncTrlsLimitFactor=ManData.ReshapeCell2Mat(IncludedTrials(IndTrainCond),64)'; % now we have the train factors to be shuffled
    else
        ThisIncTrlsLimitFactor=ones(sum(TrlSiz),1);
    end

    % limit based on the trial number
    ThisTrainTrlNum=ThisTrainFactors(:,FactorIndTrlNum);
    ThisTrainTrlNumInc=arrayfun(@(x) logical(sum(TrialRngInc==x)),ThisTrainTrlNum);
   
    ThisIncTrlsLimitFactor=ThisIncTrlsLimitFactor & ThisTrainTrlNumInc;

    AllLabels=ThisTrainFactors(:,TargetFactorInd);
    AltFactor=arrayfun(@(x) ThisTrainFactors(:,x),TargetFieldsInd,'UniformOutput',0);
    BlkOrder=ThisTrainFactors(:,BlkOrderInd); % get BlkOrderInfo
    RevTrlNum=ThisTrainFactors(:,RevTrlNumInd); %get RevTrlNum

    ShuffLbl{NN}=nan*ones(length(AllLabels),ClassifierOpts.NrepShufperFold);
    BlkOrderShuff{NN}=ShuffLbl{NN};
    RevTrlNumShuff{NN}=ShuffLbl{NN};

    % get the leve of the factors
    UniqueLabels=unique(AllLabels)';
    if nTargetFields==2
        UniqAltFactorVals=cellfun(@(x) unique(x)',AltFactor,'UniformOutput',0);
        if ClassifierOpts.StimCongruency(DimNum)==12 | ClassifierOpts.StimCongruency(DimNum)==16
            UniqAltFactorVals{strcmp(TargetFields,'Congruency')}=0;
        end
        UniqAltFactorVals=ManData.CreatPairCombsCell(UniqAltFactorVals{1},UniqAltFactorVals{2});
    elseif nTargetFields==3
        UniqAltFactorVals=cellfun(@(x) unique(x)',AltFactor,'UniformOutput',0);
        UniqAltFactorVals=ManData.CreatTripleCombsCell(UniqAltFactorVals{1},UniqAltFactorVals{2},UniqAltFactorVals{3});
    else
        UniqAltFactorVals=unique(AltFactor{1});
    end

    if isempty(TargetFields) | ClassifierOpts.SkipShuffleBalancing(DimNum) % if we don't have any balancing condition for shuffling
        % LabelInd=randsample(length(AllLabels),length(AllLabels));
        % ShuffLbl{NN}(:,ss)=AllLabels(LabelInd);
        for ss=1:ClassifierOpts.NrepShufperFold
            Ncat=sum(SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor);
            LabelInd=randsample(find(SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor),Ncat);
            
            % change the labels thaat are affected
            ShuffLbl{NN}(SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=AllLabels(LabelInd);
            
            % change the block order with the same random sample as well
             BlkOrderShuff{NN}(SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=BlkOrder(LabelInd);
            RevTrlNumShuff{NN}(SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=RevTrlNum(LabelInd);
        end
        % keep the labels that are not affected
        ShuffLbl{NN}(isnan(ShuffLbl{NN}(:,ss)),ss)=AllLabels(isnan(ShuffLbl{NN}(:,ss)));
        
        % keep block order of the trials that are not affected
        BlkOrderShuff{NN}(isnan(BlkOrderShuff{NN}(:,ss)),ss)=BlkOrder(isnan(BlkOrderShuff{NN}(:,ss)));
        RevTrlNumShuff{NN}(isnan(RevTrlNumShuff{NN}(:,ss)),ss)=RevTrlNum(isnan(RevTrlNumShuff{NN}(:,ss)));        
    else
        %% sample form the labels but consider the value of seqhist if it is present
        if nTargetFields==2            
            for ss=1:ClassifierOpts.NrepShufperFold
                for ll=1:size(UniqAltFactorVals,1)
                    Ncat=sum(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor);
                    
                    LabelInd=randsample(find(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor),Ncat);
                    
                    % change the labels thaat are affected
                    ShuffLbl{NN}(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=AllLabels(LabelInd);
                    
                    % change the block order with the same random sample as well
                    BlkOrderShuff{NN}(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=BlkOrder(LabelInd);
                    RevTrlNumShuff{NN}(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=RevTrlNum(LabelInd);
                end
                % keep the labels that are not affected
                ShuffLbl{NN}(isnan(ShuffLbl{NN}(:,ss)),ss)=AllLabels(isnan(ShuffLbl{NN}(:,ss)));
                
                % keep block order of the trials that are not affected
                BlkOrderShuff{NN}(isnan(BlkOrderShuff{NN}(:,ss)),ss)=BlkOrder(isnan(BlkOrderShuff{NN}(:,ss)));
                RevTrlNumShuff{NN}(isnan(RevTrlNumShuff{NN}(:,ss)),ss)=RevTrlNum(isnan(RevTrlNumShuff{NN}(:,ss)));
                
            end            
            %% check if the number of trials per category remains the same before and after the shuffleing for any trial trange
            if strcmp(TargetFactor,'TrialNum')
                tt=1;
                for trl=51:150
                    for ss=1:ClassifierOpts.NrepShufperFold
                        for ll=1:size(UniqAltFactorVals,1)
                            ThisTrials=AllLabels<trl & AllLabels>(trl-50);
                            AltFactorCheck=arrayfun(@(x) ThisTrainFactors(ThisTrials,x),TargetFieldsInd,'UniformOutput',0);
                            NcatCheck(ss,ll)=sum(AltFactorCheck{1}==UniqAltFactorVals(ll,1) & AltFactorCheck{2}==UniqAltFactorVals(ll,2)  &  SeqhistValThisNeu(ThisTrials)==ThisSeqHistVal & TrlTypeValThisNeu(ThisTrials)==ThisTrlTypeVal & TrlTypeValThisNeu(ThisTrials)==ThisTrlTypeVal & ThisIncTrlsLimitFactor(ThisTrials));
                            
                            ThisTrialsShuff=ShuffLbl{NN}<trl & ShuffLbl{NN}>(trl-50);
                            AltFactorShuffCheck=arrayfun(@(x) ThisTrainFactors(ThisTrialsShuff,x),TargetFieldsInd,'UniformOutput',0);
                            NcatCheckShuff(ss,ll)=sum(AltFactorShuffCheck{1}==UniqAltFactorVals(ll,1) & AltFactorShuffCheck{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu(ThisTrialsShuff)==ThisSeqHistVal & TrlTypeValThisNeu(ThisTrialsShuff)==ThisTrlTypeVal & ThisIncTrlsLimitFactor(ThisTrialsShuff));
                            
                            if sum(NcatCheckShuff==NcatCheck)~=length(NcatCheckShuff);
                                error('in the num of trials');
                            end
                        end
                    end
                    tt=tt+1;
                end
            else % check if the same number of categories are present
                for ss=1:ClassifierOpts.NrepShufperFold
                    ThisTrials=logical(ones(length(AllLabels),1));
                    ThisTrialsShuff=ThisTrials;
                    
                    for ll=1:size(UniqAltFactorVals,1)
                        AltFactorCheck=arrayfun(@(x) ThisTrainFactors(ThisTrials,x),TargetFieldsInd,'UniformOutput',0);
                        NcatCheck(ss,ll)=sum(AltFactorCheck{1}==UniqAltFactorVals(ll,1) & AltFactorCheck{2}==UniqAltFactorVals(ll,2)  &  SeqhistValThisNeu(ThisTrials)==ThisSeqHistVal & TrlTypeValThisNeu(ThisTrials)==ThisTrlTypeVal & TrlTypeValThisNeu(ThisTrials)==ThisTrlTypeVal & ThisIncTrlsLimitFactor(ThisTrials));
                        
                        AltFactorShuffCheck=arrayfun(@(x) ThisTrainFactors(ThisTrialsShuff,x),TargetFieldsInd,'UniformOutput',0);
                        NcatCheckShuff(ss,ll)=sum(AltFactorShuffCheck{1}==UniqAltFactorVals(ll,1) & AltFactorShuffCheck{2}==UniqAltFactorVals(ll,2) & SeqhistValThisNeu(ThisTrialsShuff)==ThisSeqHistVal & TrlTypeValThisNeu(ThisTrialsShuff)==ThisTrlTypeVal & ThisIncTrlsLimitFactor(ThisTrialsShuff));
                        
                        if sum(NcatCheckShuff(:)==NcatCheck(:))~=length(NcatCheckShuff(:))
                            error('in the num of trials');
                        end
                    end
                end
            end
        elseif nTargetFields==3            
            for ss=1:ClassifierOpts.NrepShufperFold
                for ll=1:size(UniqAltFactorVals,1)
                    Ncat=sum(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & AltFactor{3}==UniqAltFactorVals(ll,3) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor);
                    
                    LabelInd=randsample(find(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & AltFactor{3}==UniqAltFactorVals(ll,3) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor),Ncat);
                    
                    % change the labels thaat are affected
                    ShuffLbl{NN}(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & AltFactor{3}==UniqAltFactorVals(ll,3) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=AllLabels(LabelInd);
                    
                    % change the block order with the same random sample as well
                    BlkOrderShuff{NN}(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & AltFactor{3}==UniqAltFactorVals(ll,3) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=BlkOrder(LabelInd);
                    RevTrlNumShuff{NN}(AltFactor{1}==UniqAltFactorVals(ll,1) & AltFactor{2}==UniqAltFactorVals(ll,2) & AltFactor{3}==UniqAltFactorVals(ll,3) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=RevTrlNum(LabelInd);
                    
                end
                % keep the labels that are not affected
                ShuffLbl{NN}(isnan(ShuffLbl{NN}(:,ss)),ss)=AllLabels(isnan(ShuffLbl{NN}(:,ss)));
                
                % keep block order of the trials that are not affected
                BlkOrderShuff{NN}(isnan(BlkOrderShuff{NN}(:,ss)),ss)=BlkOrder(isnan(BlkOrderShuff{NN}(:,ss)));
                RevTrlNumShuff{NN}(isnan(RevTrlNumShuff{NN}(:,ss)),ss)=RevTrlNum(isnan(RevTrlNumShuff{NN}(:,ss)));
            end            
            %% check if the number of trials per category remains the same before and after the shuffleing for any trial trange
            if strcmp(TargetFactor,'TrialNum')
                tt=1;
                for trl=51:150
                    for ss=1:ClassifierOpts.NrepShufperFold
                        for ll=1:size(UniqAltFactorVals,1)
                            ThisTrials=AllLabels<trl & AllLabels>(trl-50);
                            AltFactorCheck=arrayfun(@(x) ThisTrainFactors(ThisTrials,x),TargetFieldsInd,'UniformOutput',0);
                            NcatCheck(ss,ll)=sum(AltFactorCheck{1}==UniqAltFactorVals(ll,1) & AltFactorCheck{2}==UniqAltFactorVals(ll,2) & AltFactorCheck{3}==UniqAltFactorVals(ll,3) &  SeqhistValThisNeu(ThisTrials)==ThisSeqHistVal & TrlTypeValThisNeu(ThisTrials)==ThisTrlTypeVal & ThisIncTrlsLimitFactor(ThisTrials));
                            
                            ThisTrialsShuff=ShuffLbl{NN}<trl & ShuffLbl{NN}>(trl-50);
                            AltFactorShuffCheck=arrayfun(@(x) ThisTrainFactors(ThisTrialsShuff,x),TargetFieldsInd,'UniformOutput',0);
                            NcatCheckShuff(ss,ll)=sum(AltFactorShuffCheck{1}==UniqAltFactorVals(ll,1) & AltFactorShuffCheck{2}==UniqAltFactorVals(ll,2) & AltFactorShuffCheck{3}==UniqAltFactorVals(ll,3) & SeqhistValThisNeu(ThisTrialsShuff)==ThisSeqHistVal & TrlTypeValThisNeu(ThisTrialsShuff)==ThisTrlTypeVal & ThisIncTrlsLimitFactor(ThisTrialsShuff));
                            
                            if sum(NcatCheckShuff==NcatCheck)~=length(NcatCheckShuff);
                                error('in the num of trials');
                            end
                        end
                    end
                    tt=tt+1;
                end
            end
        else           
            for ss=1:ClassifierOpts.NrepShufperFold
                for ll=1:size(UniqAltFactorVals,1)
                    Ncat=sum(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor);
                    
                    LabelInd=randsample(find(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor),Ncat);
                    
                    % change the labels thaat are affected
                    ShuffLbl{NN}(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=AllLabels(LabelInd);
                    
                    % change the block order with the same random sample as well
                    BlkOrderShuff{NN}(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=BlkOrder(LabelInd);
                    RevTrlNumShuff{NN}(AltFactor{1}==UniqAltFactorVals(ll) & SeqhistValThisNeu==ThisSeqHistVal & TrlTypeValThisNeu==ThisTrlTypeVal & ThisIncTrlsLimitFactor,ss)=RevTrlNum(LabelInd);
                    
                end
                % keep the labels that are not affected
                ShuffLbl{NN}(isnan(ShuffLbl{NN}(:,ss)),ss)=AllLabels(isnan(ShuffLbl{NN}(:,ss)));
                
                % keep block order of the trials that are not affected
                BlkOrderShuff{NN}(isnan(BlkOrderShuff{NN}(:,ss)),ss)=BlkOrder(isnan(BlkOrderShuff{NN}(:,ss)));
                
                RevTrlNumShuff{NN}(isnan(RevTrlNumShuff{NN}(:,ss)),ss)=RevTrlNum(isnan(RevTrlNumShuff{NN}(:,ss)));
                
            end            
            if 1% check if the same number of categories are present
                for ss=1:ClassifierOpts.NrepShufperFold
                    ThisTrials=logical(ones(length(AllLabels),1));
                    ThisTrialsShuff=ThisTrials;
                    
                    for ll=1:size(UniqAltFactorVals,1)
                        AltFactorCheck=arrayfun(@(x) ThisTrainFactors(ThisTrials,x),TargetFieldsInd,'UniformOutput',0);
                        NcatCheck(ss,ll)=sum(AltFactorCheck{1}==UniqAltFactorVals(ll,1)    &  SeqhistValThisNeu(ThisTrials)==ThisSeqHistVal & TrlTypeValThisNeu(ThisTrials)==ThisTrlTypeVal & TrlTypeValThisNeu(ThisTrials)==ThisTrlTypeVal & ThisIncTrlsLimitFactor(ThisTrials));
                        
                        AltFactorShuffCheck=arrayfun(@(x) ThisTrainFactors(ThisTrialsShuff,x),TargetFieldsInd,'UniformOutput',0);
                        NcatCheckShuff(ss,ll)=sum(AltFactorShuffCheck{1}==UniqAltFactorVals(ll,1) & SeqhistValThisNeu(ThisTrialsShuff)==ThisSeqHistVal & TrlTypeValThisNeu(ThisTrialsShuff)==ThisTrlTypeVal & ThisIncTrlsLimitFactor(ThisTrialsShuff));
                        
                        if sum(NcatCheckShuff(:)==NcatCheck(:))~=length(NcatCheckShuff(:))
                            error('in the num of trials');
                        end
                    end
                end
            end
        end
    end
end
end

function IncludedTrls=LimitTrialsBasedonFactor(FactorData,FactorName,FactorValues,FactorOperation) % limits trials of a selection of based on limits imposed by an specific factor
%@FactorData is the output of PrepareData4ClassifierAnalysis
%@FactorName name of the factor we want to change
%@FactorValues accepted values of this factor
%@FactorOperation define what operation we are doign with
%factor can be 'equal','bigger','smaller','interval(provide values in pairs of columns)'
global AnalysisOpts

if isnan(FactorValues);IncludedTrls=[];return;end % we are not changing anything
% if it is trial number but they are negative then reverse the order
if strcmp('TrialNum',FactorName) & sum(FactorValues<0); FactorName='TrialNumReverse';end
if ~exist('FactorOperation','var');FactorOperation='equal';end
if isnan(FactorOperation);return;end % we are not changing anything

FactorInd=find(strcmp(AnalysisOpts.factornames,FactorName));
% loop through all of the structures and change the fields (remove data mean because it is not correct anymore

switch FactorOperation
    case 'equal'
        IncludedTrls=cellfun(@(x) cell2mat(arrayfun(@(y) (x(:,FactorInd)==y)',FactorValues,'Uniformoutput',0)),FactorData.AllFactors,'uniformoutput',0);
    case 'bigger'
        IncludedTrls=cellfun(@(x) cell2mat(arrayfun(@(y) (x(:,FactorInd)>=y)',FactorValues,'Uniformoutput',0)),FactorData.AllFactors,'uniformoutput',0);
    case 'smaller'
        IncludedTrls=cellfun(@(x) cell2mat(arrayfun(@(y) (x(:,FactorInd)<=y)',FactorValues,'Uniformoutput',0)),FactorData.AllFactors,'uniformoutput',0);
    case 'interval'
        IncludedTrls=cellfun(@(x) cell2mat(arrayfun(@(y) (x(:,FactorInd)>=FactorValues(y,1) & x(:,FactorInd)<=FactorValues(y,2))',1:size(FactorValues,1),'Uniformoutput',0)),FactorData(i).AllFactors,'uniformoutput',0);
end

end

function [TrlTypeVal_Train,TrlTypeVal_Test]=AdjustLimitFactorCondition(ClassifierOpts,DimNum,Cond)

DimTxt={'','2','3'};

% if we have a specific trial type then we need to take that into account
if ClassifierOpts.TrialType(DimNum)~=3 % if we are taking correct or incorrect trials

    if ClassifierOpts.TrialType(DimNum)==2;ClassifierOpts.TrialType(DimNum)=0;end % this is incorrect trials

    TrlTypeVal_Train=ClassifierOpts.TrialType(DimNum);
    TrlTypeVal_Test=TrlTypeVal_Train;
    

else  % if we are taking all trials check if we are limiting the distribution of train or test

    if isfield(ClassifierOpts,['LimFactorTrialVal' num2str(DimTxt{DimNum})])
        eval(sprintf('LimiFactorVal=ClassifierOpts.LimFactorTrialVal%s{Cond};',DimTxt{DimNum}));
        eval(sprintf('ClassifierOpts.LimFactorTrialName=ClassifierOpts.LimFactorTrialName%s;',DimTxt{DimNum}));

    elseif isfield(ClassifierOpts,'LimFactorTrialVal')
        LimiFactorVal=ClassifierOpts.LimFactorTrialVal{Cond};
        ClassifierOpts.LimFactorTrialName=ClassifierOpts.LimFactorTrialName;

    else
         TrlTypeVal_Train=nan;TrlTypeVal_Test=nan;
        return;
    end
   
    
    if isnan(LimiFactorVal) || DimNum~=LimiFactorVal
        TrlTypeVal_Train=nan;TrlTypeVal_Test=nan;
        return;
    end

    if strcmp(ClassifierOpts.LimFactorTrialName,'CorrectTrls')
        TrlTypeVal_Train=1; TrlTypeVal_Test=1;
    elseif strcmp(ClassifierOpts.LimFactorTrialName,'CorrectTrlsTrainIncorTrlsTest') % keep correct trials at the train but change the test
        TrlTypeVal_Train=1;TrlTypeVal_Test=0; % only take correct trials for train and keep incorrect for test
    elseif strcmp(ClassifierOpts.LimFactorTrialName,'AllTrlsTrainIncorTrlsTest') % keep all trials at the train but change the test
        TrlTypeVal_Train=nan;TrlTypeVal_Test=0; % only take all trials for train and keep incorrect for test
    elseif strcmp(ClassifierOpts.LimFactorTrialName,'AllTrlsTrainCorrTrlsTest') % keep all trials at the train but change the test
        TrlTypeVal_Train=nan;TrlTypeVal_Test=1; % only take all trials for train and keep correct for test
    else
        TrlTypeVal_Train=nan;TrlTypeVal_Test=nan;
    end

end
end
function FactorData=SkipSeqHistTrain(FactorData,ClassifierOpts,Cond)
% if ClassifierOpts.SkipSeqHistTrain=1; then we don't apply SeqHist to train data

if ClassifierOpts.SkipSeqHistTrain~=1;return;end
SeqHistInd=ClassifierOpts.SeqHist.Ind;
for Neu=1:length(FactorData)
    for i=1:length(FactorData(Neu).AllFactors)
        FactorData(Neu).AllFactors{i}(:,SeqHistInd)=ClassifierOpts.SeqHistCond{Cond};
    end
end
end

