function [Scores_1ndD_Out,Scores_2ndD_Out,TargetFactorDim1,TargetFactorDim2,Nrep,NrepShuff,AllFactorVals_1ndD,AllFactorVals_2ndD,ThisFactorResp,ThisFactorResp2]=ExtractScores4ShuffleData(ClassifierResults_1ndD,ClassifierResults_2ndD,ClassifierOpts,Cond,Dims,TargetFactor) % plots PEV of two distibutions of scores
% extracts scores for shuffle data to process
global AnalysisOpts
NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
opts.SkipProcessingScores=1;

IsFoldRep=isfield(ClassifierResults_1ndD(1).TrialRange(1).Rep(1),'FoldRep'); % is there a fold rep for the suffle data
if IsFoldRep
    Nrep=length(ClassifierResults_1ndD(Cond).TrialRange(1).Rep(1).FoldRep);
    NrepShuff=length(ClassifierResults_1ndD(Cond).TrialRange(1).Rep);
else
    Nrep=length(ClassifierResults_1ndD(Cond).TrialRange(1).Rep);
    NrepShuff=1;
end

% look if test for congruent or incongruent trials is included here?
if contains(ClassifierOpts.Name,'Congruency') 
    LimFactorTrialVal=ClassifierOpts.LimFactorTrialVal{Cond};
    if isnan(LimFactorTrialVal);LimFactorTrialVal=[2 2];end
    if (strcmp(TargetFactor,'QuadrantsInCong') & LimFactorTrialVal(2)==1) | (strcmp(TargetFactor,'QuadrantsCong') & LimFactorTrialVal(2)==0)
        return
    end
end
nTrialRange=length(ClassifierResults_1ndD(Cond).TrialRange);
TrialRangeSet=ClassifierOpts.TestTrlNumRange{Cond}(2:4);
TrialRange=[TrialRangeSet(2):TrialRangeSet(3):TrialRangeSet(1)];
[LvlsTestFactor1,TargetFactorDim1,~,TrainCond1,TestCond1,Dim1Txt,TestFactorTxt1,TrainFactorTxt1,LevelFieldName1]=...
    NeuAna.getClassifierDimInfo(ClassifierOpts,Dims(1));
[LvlsTestFactor2,TargetFactorDim2,~,TrainCond2,TestCond2,Dim2Txt,TestFactorTxt2,TrainFactorTxt2,LevelFieldName2]=...
    NeuAna.getClassifierDimInfo(ClassifierOpts,Dims(2));

if ~ClassifierOpts.EqualizeTrialsXConds & strcmp(TargetFactor,'Quadrants')
    error('We are assuming exalized sampling of trials X dimensions...')
end

% get classifier perfromance
Accuracy_1ndD=ManData.ReshapeStruct2Mat(ClassifierResults_1ndD(Cond).TrialRange,'Accuracy',4);
Accuracy_2ndD=ManData.ReshapeStruct2Mat(ClassifierResults_2ndD(Cond).TrialRange,'Accuracy',4);

%get scores
for rep_sh=1:NrepShuff
    fprintf('\n processing shuffle repetition %i',rep_sh)
    for TrlRng=1:nTrialRange
        for rep=1:Nrep
            % grab trials and tile them based on different values
            % of factor we care about % this is tiling test factors
            % based on the TargetFactor
            if ~IsFoldRep
                [FactorDataInd1{rep},FactorLevels]=NeuAna.TileTrialsbasedonFactor(ClassifierResults_1ndD(Cond).TrialRange(TrlRng).Rep(rep).(TestFactorTxt1),TargetFactor,ClassifierOpts.(LevelFieldName1){1});
                [FactorDataInd2{rep}  ~          ]=NeuAna.TileTrialsbasedonFactor(ClassifierResults_2ndD(Cond).TrialRange(TrlRng).Rep(rep).(TestFactorTxt1),TargetFactor,ClassifierOpts.(LevelFieldName2){1});
               
                % grab score values now based on number of the trial in a specific repetition
                Scores_1ndD(TrlRng).TrialRange(rep,:)=cellfun(@(x) ClassifierResults_1ndD(Cond).TrialRange(TrlRng).Scores{rep}(x,:,:),FactorDataInd1{rep},'UniformOutput',0);
                Scores_2ndD(TrlRng).TrialRange(rep,:)=cellfun(@(x) ClassifierResults_2ndD(Cond).TrialRange(TrlRng).Scores{rep}(x,:,:),FactorDataInd2{rep},'UniformOutput',0);

                % get all of the factor values for these trials
                AllFactorVals_1ndD(TrlRng).TrialRange(rep,:)=cellfun(@(x) ClassifierResults_1ndD(Cond).TrialRange(TrlRng).Rep(rep).(TestFactorTxt1)(x,:),FactorDataInd1{rep},'UniformOutput',0);
                AllFactorVals_2ndD(TrlRng).TrialRange(rep,:)=cellfun(@(x) ClassifierResults_2ndD(Cond).TrialRange(TrlRng).Rep(rep).(TestFactorTxt1)(x,:),FactorDataInd2{rep},'UniformOutput',0);

            else
                [FactorDataInd1{rep},FactorLevels]=NeuAna.TileTrialsbasedonFactor(ClassifierResults_1ndD(Cond).TrialRange(TrlRng).Rep(rep_sh).FoldRep(rep).(TestFactorTxt1),TargetFactor,ClassifierOpts.(LevelFieldName1){1});
                [FactorDataInd2{rep},  ~          ]=NeuAna.TileTrialsbasedonFactor(ClassifierResults_2ndD(Cond).TrialRange(TrlRng).Rep(rep_sh).FoldRep(rep).(TestFactorTxt1),TargetFactor,ClassifierOpts.(LevelFieldName2){1});
              
                % grab score values now based on number of the trial in a specific repetition
                Scores_1ndD(TrlRng).TrialRange(rep,:)=cellfun(@(x) ClassifierResults_1ndD(Cond).TrialRange(TrlRng).Scores{rep_sh}{rep}(x,:,:),FactorDataInd1{rep},'UniformOutput',0);
                Scores_2ndD(TrlRng).TrialRange(rep,:)=cellfun(@(x) ClassifierResults_2ndD(Cond).TrialRange(TrlRng).Scores{rep_sh}{rep}(x,:,:),FactorDataInd2{rep},'UniformOutput',0);
              
                % get all of the factor values for these trials
                AllFactorVals_1ndD(TrlRng).TrialRange(rep,:)=cellfun(@(x) ClassifierResults_1ndD(Cond).TrialRange(TrlRng).Rep(rep_sh).FoldRep(rep).(TestFactorTxt1)(x,:),FactorDataInd1{rep},'UniformOutput',0);
                AllFactorVals_2ndD(TrlRng).TrialRange(rep,:)=cellfun(@(x) ClassifierResults_2ndD(Cond).TrialRange(TrlRng).Rep(rep_sh).FoldRep(rep).(TestFactorTxt1)(x,:),FactorDataInd2{rep},'UniformOutput',0);
            end

            % get the response category of the morph level based on the desired factor only for the first dimension
            ThisFactorResp=GetTestFactorLevels(ManData,ClassifierResults_1ndD,TargetFactorDim1,TestFactorTxt1,FactorDataInd1,Cond,TrlRng,rep,rep_sh);
            ThisFactorResp2=GetTestFactorLevels(ManData,ClassifierResults_2ndD,TargetFactorDim2,TestFactorTxt1,FactorDataInd2,Cond,TrlRng,rep,rep_sh);
        end
        if ~opts.SkipProcessingScores
            % reorganize factor values
            AllFactorVals_1ndD(TrlRng).TrialRange=arrayfun(@(x) cell2mat(AllFactorVals_1ndD(TrlRng).TrialRange(:,x)),1:length(FactorLevels),'UniformOutput',0);
            AllFactorVals_2ndD(TrlRng).TrialRange=arrayfun(@(x) cell2mat(AllFactorVals_2ndD(TrlRng).TrialRange(:,x)),1:length(FactorLevels),'UniformOutput',0);

            AllFactorVals_1ndD_RS(TrlRng).TrialRange=cell2mat(AllFactorVals_1ndD(TrlRng).TrialRange');
            AllFactorVals_2ndD_RS(TrlRng).TrialRange=cell2mat(AllFactorVals_2ndD(TrlRng).TrialRange');

            % take average of the factors in case we need them in the future
            AllFactorVals_1ndD_avg(TrlRng).TrialRange=mean(ManData.ReshapeCell2Mat(AllFactorVals_1ndD(TrlRng).TrialRange,3),3);
            AllFactorVals_2ndD_avg(TrlRng).TrialRange=mean(ManData.ReshapeCell2Mat(AllFactorVals_2ndD(TrlRng).TrialRange,3),3);


            %% check if factors across dimensions match with each other
            if size(AllFactorVals_1ndD_RS(TrlRng).TrialRange,1)==size(AllFactorVals_2ndD_RS(TrlRng).TrialRange,1)
                TestFactors=NeuAna.AdjustObjCategory4Factor(AllFactorVals_1ndD_RS(TrlRng).TrialRange)==NeuAna.AdjustObjCategory4Factor(AllFactorVals_2ndD_RS(TrlRng).TrialRange);
                if sum(TestFactors(:)==0) & ClassifierOpts.EqualizeTrialsXConds
                    error('Trial types are not corresponding in this analysis');
                elseif sum(TestFactors(:)==0) & ~ClassifierOpts.EqualizeTrialsXConds
                    warning('Trial types are not corresponding in this analysis');
                end
            end
            %% warning this magnitude flip is only correct to use when factor of 1ndD matches Target factor
            % get the magnitude for first dimension as well
            % to perfrom this either we shoud have Target factor as
            % qudrant and the first dimension factor as shape or clor
            % Or the factor of the first dimension should match the TargetFactor
            if contains(TargetFactor(1:5),TargetFactorDim1(1:5)) | (contains(TargetFactor,'Quadrants') & (contains(TargetFactorDim1,'Shape') | contains(TargetFactorDim1,'Color')))
                MagnitudeFlip=ThisFactorResp';MagnitudeFlip=sort(unique(MagnitudeFlip),2,'ascend');
                ScoresMag_1ndD(TrlRng).TrialRange=arrayfun(@(x) cell2mat(Scores_1ndD(TrlRng).TrialRange(:,x)),MagnitudeFlip,'UniformOutput',0);
                ScoresMag_1ndD(TrlRng).TrialRange=arrayfun(@(x) squeeze(ScoresMag_1ndD(TrlRng).TrialRange{x}(:,x,:)),MagnitudeFlip,'UniformOutput',0);
                ScoresMag_1ndD_RS(TrlRng).TrialRange=cell2mat(ScoresMag_1ndD(TrlRng).TrialRange');
            end

            % get the scores for both dimensions
            Scores_1ndD_Org(TrlRng).TrialRange=arrayfun(@(x) cell2mat(Scores_1ndD(TrlRng).TrialRange(:,x)),1:length(FactorLevels),'UniformOutput',0);
            Scores_1ndD(TrlRng).TrialRange=cellfun(@(x) squeeze(x(:,1,:)),Scores_1ndD_Org(TrlRng).TrialRange,'UniformOutput',0);
            Scores_1ndD_Flip(TrlRng).TrialRange=cellfun(@(x) squeeze(x(:,2,:)),Scores_1ndD_Org(TrlRng).TrialRange,'UniformOutput',0);


            Scores_2ndD_Org(TrlRng).TrialRange=arrayfun(@(x) cell2mat(Scores_2ndD(TrlRng).TrialRange(:,x)),1:length(FactorLevels),'UniformOutput',0);
            Scores_2ndD(TrlRng).TrialRange=cellfun(@(x) squeeze(x(:,1,:)),Scores_2ndD_Org(TrlRng).TrialRange,'UniformOutput',0);
            Scores_2ndD_Flip(TrlRng).TrialRange=cellfun(@(x) squeeze(x(:,2,:)),Scores_2ndD_Org(TrlRng).TrialRange,'UniformOutput',0);

            % concatinate the scores
            Scores_1ndD_RS(TrlRng).TrialRange=cell2mat(Scores_1ndD(TrlRng).TrialRange');
            Scores_2ndD_RS(TrlRng).TrialRange=cell2mat(Scores_2ndD(TrlRng).TrialRange');

            % take the mean of the scores for the second dimension as  we might need it for compression
            Scores_2ndD_avg(TrlRng).TrialRange=mean(ManData.ReshapeCell2Mat(Scores_2ndD(TrlRng).TrialRange,3),3);
            Scores_2ndD_avg_avg(TrlRng).TrialRange=mean(mean(ManData.ReshapeCell2Mat(Scores_2ndD(TrlRng).TrialRange,3),3),1);

            Scores_2ndD_avg_Flip(TrlRng).TrialRange=mean(ManData.ReshapeCell2Mat(Scores_2ndD_Flip(TrlRng).TrialRange,3),3);
            Scores_2ndD_avg_avg_Flip(TrlRng).TrialRange=mean(mean(ManData.ReshapeCell2Mat(Scores_2ndD_Flip(TrlRng).TrialRange,3),3),1);

            % get the magnitude based on each dimension category
            TargetFactorDim1Ch=TargetFactorDim1;TargetFactorDim2Ch=TargetFactorDim2;
            if strcmp(TargetFactorDim1,'ColorCat') & length(AnalysisOpts.FactorInds2Keep)==9;TargetFactorDim1Ch='ColorML';elseif strcmp(TargetFactorDim1,'ShapeCat') & length(AnalysisOpts.FactorInds2Keep)==9;TargetFactorDim1Ch='ShapeML';end
            if strcmp(TargetFactorDim2,'ColorCat') & length(AnalysisOpts.FactorInds2Keep)==9;TargetFactorDim2Ch='ColorML';elseif strcmp(TargetFactorDim2,'ShapeCat') & length(AnalysisOpts.FactorInds2Keep)==9;TargetFactorDim2Ch='ShapeML';end
            AdjFactors_1nD=NeuAna.AdjustObjCategory4Factor(AllFactorVals_1ndD_RS(TrlRng).TrialRange);
            AdjFactors_2nD=NeuAna.AdjustObjCategory4Factor(AllFactorVals_2ndD_RS(TrlRng).TrialRange);

            TargFactorIndDim_1nD=contains(AnalysisOpts.FactorInds2Keep,TargetFactorDim1Ch);
            TargFactorIndDim_2nD=contains(AnalysisOpts.FactorInds2Keep,TargetFactorDim2Ch);
            Levels_1nD=sort(unique(AdjFactors_1nD(:,TargFactorIndDim_1nD)));
            Levels_2nD=sort(unique(AdjFactors_2nD(:,TargFactorIndDim_2nD)));

            Scores_1ndD_RS_MagCat(TrlRng).TrialRange=nan*ones(size(Scores_1ndD_RS(TrlRng).TrialRange));
            Scores_2ndD_RS_MagCat(TrlRng).TrialRange=nan*ones(size(Scores_2ndD_RS(TrlRng).TrialRange));
            Scores_1ndD_Org_RS(TrlRng).TrialRange=ManData.ReshapeCell2Mat(Scores_1ndD_Org(TrlRng).TrialRange,62);
            Scores_2ndD_Org_RS(TrlRng).TrialRange=ManData.ReshapeCell2Mat(Scores_2ndD_Org(TrlRng).TrialRange,62);

            LevelsTrls_1nD=arrayfun(@(x) AdjFactors_1nD(:,TargFactorIndDim_1nD)==x,Levels_1nD,'uniformoutput',0);
            LevelsTrls_2nD=arrayfun(@(x) AdjFactors_2nD(:,TargFactorIndDim_2nD)==x,Levels_2nD,'uniformoutput',0);

            Scores_1ndD_RS_MagCat(TrlRng).TrialRange(LevelsTrls_1nD{1},:)=squeeze(Scores_1ndD_Org_RS(TrlRng).TrialRange(LevelsTrls_1nD{1},1,:));
            Scores_2ndD_RS_MagCat(TrlRng).TrialRange(LevelsTrls_2nD{1},:)=squeeze(Scores_2ndD_Org_RS(TrlRng).TrialRange(LevelsTrls_2nD{1},1,:));

            if length(Levels_1nD)>1
                Scores_1ndD_RS_MagCat(TrlRng).TrialRange(LevelsTrls_1nD{2},:)=squeeze(Scores_1ndD_Org_RS(TrlRng).TrialRange(LevelsTrls_1nD{2},2,:));
            end
            if length(Levels_2nD)>1
                Scores_2ndD_RS_MagCat(TrlRng).TrialRange(LevelsTrls_2nD{2},:)=squeeze(Scores_2ndD_Org_RS(TrlRng).TrialRange(LevelsTrls_2nD{2},2,:));
            end
        end
        % report the final result
        Scores_1ndD_Out(TrlRng,rep_sh).TrialRange=Scores_1ndD(TrlRng).TrialRange;
        Scores_2ndD_Out(TrlRng,rep_sh).TrialRange=Scores_2ndD(TrlRng).TrialRange;
    end
end
end
function ThisFactorResp=GetTestFactorLevels(ManData,ClassifierResults,TargetFactorDim,TestFactorTxt,FactorDataInd,Cond,TrlRng,rep,rep_sh)
global AnalysisOpts

IsFoldRep=isfield(ClassifierResults(Cond).TrialRange(1).Rep(1),'FoldRep'); % is there a fold rep for the suffle data

% get the labels for the test samples
if strcmp(TargetFactorDim,'ColorCat') | strcmp(TargetFactorDim,'ShapeCat');TargetFactorDim=[TargetFactorDim(1:5) 'ML'];end
FactorInd=contains(AnalysisOpts.FactorInds2Keep,TargetFactorDim);%(1:4));
if IsFoldRep
    ThisFactorResp=ClassifierResults(Cond).TrialRange(TrlRng).Rep(rep_sh).FoldRep(rep).(TestFactorTxt)(cell2mat(FactorDataInd{rep}),FactorInd);
else
    ThisFactorResp=ClassifierResults(Cond).TrialRange(TrlRng).Rep(rep).(TestFactorTxt)(cell2mat(FactorDataInd{rep}),FactorInd);
end
if contains(TargetFactorDim,'Shape') | contains(TargetFactorDim,'Color')
    ThisFactorResp=ManData.CategorizeMorphlevel(ThisFactorResp);
end
end