function ShuffLbl=CreateShuffleDist4Classifier(obj,TrainDataAllFactors,TestDataAllFactors,response,ClassifierOpts,DimNum,Cond)


fprintf('\ngenerating shuffle distribution for dimension:%i',DimNum);
DimTxt=[0 2 3];

% first take the train label and test labels
if ClassifierOpts.CVspecs(DimNum,Cond)==0 % then use only training data
    AllLabels=[response{1}];
elseif ClassifierOpts.CVspecs(DimNum,Cond)==1 % then use test and train data
    AllLabels=[response{1} ;response{2}];
end


if ClassifierOpts.StimCongruency(DimNum)~=0
    TargetFields=obj.StimCongruencyFactorName{ClassifierOpts.StimCongruency(DimNum)};
else
    TargetFields={'none'};
end

[~,TargetFactor,ObservedFieldName,TrainCond,TestCond,DimTxt,TestFactorTxt,...
    TrainFactorTxt,LevelFieldName]=obj.getClassifierDimInfo(ClassifierOpts,DimTxt(DimNum));

TrainCond=TrainCond{Cond};

TargetFactorSet={'Color','Shape','ResponseLoc','Reward'}; % list of keywords for target factor
TargetFactorSetAlter={'ShapeCat','ColorCat','ResponseLoc','Reward'}; % list of alternate of the target factors

% now based on the situation create a shuffle distribution
if contains(TargetFactor,'Color') |  contains(TargetFactor,'Shape')

    AlternateFactor= TargetFactorSetAlter{contains(TargetFactorSet,TargetFactor(1:5))};

    if sum(contains(TargetFields,'Congruency','IgnoreCase',true))
        ShuffLbl=SamplefromFactor(ClassifierOpts,AlternateFactor,TrainDataAllFactors,TestDataAllFactors,AllLabels,DimNum,Cond);
    elseif sum(contains(TargetFields,'ResponseLoc','IgnoreCase',true))
        ShuffLbl=SamplefromFactor(ClassifierOpts,'ResponseLoc',TrainDataAllFactors,TestDataAllFactors,AllLabels,DimNum,Cond);
    elseif sum(contains(TargetFields,'Reward','IgnoreCase',true))
        ShuffLbl=SamplefromFactor(ClassifierOpts,'Reward',TrainDataAllFactors,TestDataAllFactors,AllLabels,DimNum,Cond);
    end

elseif sum(contains(TargetFactor,'ResponseLoc','IgnoreCase',true))

    % then make it stable with respect to the feature in effect
    if TrainCond==1;AlternateFactor='ShapeCat';else;AlternateFactor='ColorCat';end
    ShuffLbl=SamplefromFactor(ClassifierOpts,AlternateFactor,TrainDataAllFactors,TestDataAllFactors,AllLabels,DimNum,Cond);

elseif sum(contains(TargetFactor,'Rule','IgnoreCase',true)) % we don't balance for anything
    for ss=1:ClassifierOpts.NrepShufperFold
        ShuffLbl(:,ss)=randsample(AllLabels,length(AllLabels));
    end
end
end


function ShuffLbl=SamplefromFactor(ClassifierOpts,AlternateFactor,TrainDataAllFactors,TestDataAllFactors,AllLabels,DimNum,Cond)

try
    if ClassifierOpts.CVspecs(DimNum,Cond)==1 % then use train and test  data
        ThisFactor=[TrainDataAllFactors(:,ClassifierOpts.(AlternateFactor).Ind);...
            TestDataAllFactors(:,ClassifierOpts.(AlternateFactor).Ind)];
    elseif  ClassifierOpts.CVspecs(DimNum,Cond)==0 % then use only training data
        ThisFactor=[TrainDataAllFactors(:,ClassifierOpts.(AlternateFactor).Ind)];
    end
    % what loop on the unique values and unique labels
    UniqFactorVals=unique(ThisFactor)';
    UniqueLabels=unique(AllLabels)';
    if length(UniqueLabels)~=2 | length(UniqFactorVals)~=2 ;error('we have more than two labels ');end

    % check if the number of shape category within each color category is the same
    EqualNum1=sum(AllLabels(ThisFactor==UniqFactorVals(1))==UniqueLabels(1))==sum(AllLabels(ThisFactor==UniqFactorVals(1))==UniqueLabels(2));
    EqualNum2=sum(AllLabels(ThisFactor==UniqFactorVals(2))==UniqueLabels(1))==sum(AllLabels(ThisFactor==UniqFactorVals(2))==UniqueLabels(2));
    ShuffLbl=nan*ones(length(AllLabels),ClassifierOpts.NrepShufperFold);

    for ss=1:ClassifierOpts.NrepShufperFold
        if (EqualNum1==1 & EqualNum2==1) & ClassifierOpts.UseBalancedShuffling==1
            Ncat1=sum(ThisFactor==1);Ncat2=sum(ThisFactor==2);
            ShuffLbl(ThisFactor==1,ss)=randsample(AllLabels(ThisFactor==1),Ncat1);
            ShuffLbl(ThisFactor==2,ss)=randsample(AllLabels(ThisFactor==2),Ncat2);
        else
            ShuffLbl(:,ss)=randsample(AllLabels,length(AllLabels));
        end
    end
catch me
    disp(me.message)
end
end

