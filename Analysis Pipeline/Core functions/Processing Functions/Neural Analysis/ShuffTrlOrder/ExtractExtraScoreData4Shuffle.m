function outdata=ExtractExtraScoreData4Shuffle(ClassifierOpts,AllFactorVals_1ndD,AllFactorVals_2ndD,Scores_1ndD,Scores_2ndD,ThisFactorResp,TargetFactorDim1,TargetFactorDim2,TargetFactor)

ManData=ManipulateData;
NeuAna=NeuralAnalysisFuncsTemp;
nTrialRange=length(Scores_1ndD);

if strcmp(TargetFactor,'ResponseLoc')
    FactorLevels=1:2;
elseif strcmp(TargetFactor,'Quadrants')
    FactorLevels=1:4;
end

for TrlRng=1:nTrialRange

    AllFactorVals_1ndD_RS(TrlRng).TrialRange=cell2mat(AllFactorVals_1ndD(TrlRng).TrialRange');
    AllFactorVals_2ndD_RS(TrlRng).TrialRange=cell2mat(AllFactorVals_2ndD(TrlRng).TrialRange');

    % take average of the factors in case we need them in the future
    % AllFactorVals_1ndD_avg(TrlRng).TrialRange=mean(ManData.ReshapeCell2Mat(AllFactorVals_1ndD(TrlRng).TrialRange,3),3);
    % AllFactorVals_2ndD_avg(TrlRng).TrialRange=mean(ManData.ReshapeCell2Mat(AllFactorVals_2ndD(TrlRng).TrialRange,3),3);
    A1=ManData.ReshapeCell2Mat(AllFactorVals_1ndD(TrlRng).TrialRange,42);
    A2=ManData.ReshapeCell2Mat(AllFactorVals_2ndD(TrlRng).TrialRange,42);
    AllFactorVals_1ndD_avg(TrlRng).TrialRange=mean(cell2mat(arrayfun(@(x) mean(A1{x},4),FactorLevels,'uniformoutput',0)),2)';
    AllFactorVals_2ndD_avg(TrlRng).TrialRange=mean(cell2mat(arrayfun(@(x) mean(A2{x},4),FactorLevels,'uniformoutput',0)),2)';
    
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
        MagnitudeFlip=ThisFactorResp';%MagnitudeFlip=sort(unique(MagnitudeFlip),2,'ascend');
        ScoresMag_1ndD(TrlRng).TrialRange=arrayfun(@(x) cell2mat(Scores_1ndD(TrlRng).TrialRange(:,x)),1:length(FactorLevels),'UniformOutput',0);
        ScoresMag_1ndD(TrlRng).TrialRange=arrayfun(@(x) squeeze(ScoresMag_1ndD(TrlRng).TrialRange{x}(:,MagnitudeFlip(x),:)),1:length(FactorLevels),'UniformOutput',0);
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
    if strcmp(TargetFactorDim1,'ColorCat') & length(ClassifierOpts.AnalysisOpts.FactorInds2Keep)==9;TargetFactorDim1Ch='ColorML';elseif strcmp(TargetFactorDim1,'ShapeCat') & length(ClassifierOpts.AnalysisOpts.FactorInds2Keep)==9;TargetFactorDim1Ch='ShapeML';end
    if strcmp(TargetFactorDim2,'ColorCat') & length(ClassifierOpts.AnalysisOpts.FactorInds2Keep)==9;TargetFactorDim2Ch='ColorML';elseif strcmp(TargetFactorDim2,'ShapeCat') & length(ClassifierOpts.AnalysisOpts.FactorInds2Keep)==9;TargetFactorDim2Ch='ShapeML';end
    AdjFactors_1nD=NeuAna.AdjustObjCategory4Factor(AllFactorVals_1ndD_RS(TrlRng).TrialRange);
    AdjFactors_2nD=NeuAna.AdjustObjCategory4Factor(AllFactorVals_2ndD_RS(TrlRng).TrialRange);

    TargFactorIndDim_1nD=contains(ClassifierOpts.AnalysisOpts.FactorInds2Keep,TargetFactorDim1Ch);
    TargFactorIndDim_2nD=contains(ClassifierOpts.AnalysisOpts.FactorInds2Keep,TargetFactorDim2Ch);
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
outdata=struct('ScoresMag_1ndD',ScoresMag_1ndD,'ScoresMag_1ndD_RS',ScoresMag_1ndD_RS,'Scores_1ndD_RS',Scores_1ndD_RS,...
    'Scores_2ndD_RS',Scores_2ndD_RS,'Scores_2ndD_avg',Scores_2ndD_avg,'Scores_2ndD_avg_avg',Scores_2ndD_avg_avg,...
    'Scores_2ndD_avg_Flip',Scores_2ndD_avg_Flip,'Scores_2ndD_avg_avg_Flip',Scores_2ndD_avg_avg_Flip, ...
    'AllFactorVals_1ndD_avg',AllFactorVals_1ndD_avg,'AllFactorVals_2ndD_avg',AllFactorVals_2ndD_avg);


end

% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   StimEncodingCompresScatterAvg   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % obj.PlotSVMProjectionsScores2DPEV_Learning(ClassifierResults_ndD{2},ClassifierResults_ndD{3},ClassifierOpts,...
%     Sp(8),Conds(c),'Scatter','StimEncodingCompresScatterAvg',[2 3],'ThisTargetFactor','Quadrants','LookatDim2',2,'ExtraData',CompIndexEncoding);
% 
% % Check if the target factor in the first dimension is 'Rule'
% if contains(TargetFactorDim1, 'Rule')
%     % Display a warning and exit the function if the first dimension is not 'shape' or 'color'
%     warning('first dimension has to be shape or color');
%     h = []; % Clear the handle (not shown in the provided code)
%     return; % Exit the function
% end
% 
% % Calculate the compression index for the second dimension (rule) for each trial
% CompIndex=obj.ExtraData.TrlAvg.All;
% CompIndex=log(CompIndex);
% CompIndex=obj.ManData.SmoothData(CompIndex,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
% 
% % Calculate the mean of the compression index within a specific time period
% TimRange=ClassifierOpts.SpkCountPeriod(1,:);
% TimeIndComp = AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(1, 1) & ...
%     AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(1, 2);
% CompIndex = mean(CompIndex(:, TimeIndComp), 2);
% 
% % calculate the mean of stimulus encoding
% MeanScoresMag_1ndD_RS=arrayfun(@(x) mean(ScoresMag_1ndD_RS(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
% AvgStimEncoding=cell2mat(MeanScoresMag_1ndD_RS)';
% 
% % Define time indices based on the specified time range
% AvgStimEncoding=mean(AvgStimEncoding(:,TimeIndComp),2);
% 
% % Generate line legend text based on time range
% LineLegTxt = sprintf('%0.1f->%0.1f', TimRange(1), TimRange(2));
% 
% % Calculate the correlation between mean scores and compression index
% X=AvgStimEncoding;
% Y=CompIndex;
% [a,p] = obj.ManData.Correlation(X,Y,0,'Modified_MannKendall');
% 
% % Correct p-value for multiple comparisons
% p = p * 1;
% 
% % Perform significance test based on adjusted p-value threshold
% pa = p < (0.05);
% Xlbl=['Stimulus' TargetFactorDim1 'encoding'];Ylbl='Compression Index';
% PlotType='Scatter';
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CorrBeliefBhvPerfScatterAvg
% % plot correlation of belief with behavioral performance in color and shape in scatter plot
% obj.PlotSVMProjectionsScores2DPEV_Learning(ClassifierResults_ndD{1},ClassifierResults_ndD{3},ClassifierOpts,...
%     Sp(11),Conds(c),'Scatter','CorrBeliefBhvPerfScatterAvg',[0 3],'ThisTargetFactor','Quadrants','LookatDim2',2,...
%     'ThisTitle','Corr of Color Bhv Performance 10 trls and Task Belief');
% if contains(TargetFactorDim1,'Rule');warning('first dimension has to be shape or color');h=[];return
% end
% % find the factors we want to correlate the data with
% if contains(obj.ThisTitle,'color','IgnoreCase',true)
%     BhvFactor='Color';
% elseif contains(obj.ThisTitle,'shape','IgnoreCase',true)
%     BhvFactor='Shape';
% end
% 
% % adjust the level of avegraing
% if contains(obj.ThisTitle,'10','IgnoreCase',true)
%     NtrlsAvgBhv='10';
% elseif contains(obj.ThisTitle,'50','IgnoreCase',true)
%     NtrlsAvgBhv='50';
% else
%     NtrlsAvgBhv='10';
% end
% 
% TimRange= sort( ClassifierOpts.SpkCountPeriod(3,:),2,'descend'); % Define time points of interest
% TimeInd = AnalysisOpts.Time >= ClassifierOpts.SpkCountPeriod(1, 1) & ...
%     AnalysisOpts.Time <= ClassifierOpts.SpkCountPeriod(1, 2);
% 
% % get average belief encoding
% BeliefEncoding = obj.ManData.ReshapeStruct2Mat(Scores_2ndD_avg_avg_Flip, 'TrialRange', 2);
% MeanBelief= mean(BeliefEncoding(:,TimeInd),2);
% 
% % get the behavioral performance for each trial
% Bhv_AllTrlRng=arrayfun(@(x) mean(AllFactorVals_1ndD_avg(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
% Bhv_AllTrlRng=cell2mat(Bhv_AllTrlRng)';
% 
% BhvFactorInd=cellfun(@(x) find(strcmp(AnalysisOpts.FactorInds2Keep,x)),{['BhvPerf' BhvFactor NtrlsAvgBhv]});
% BhvPerfTrl=Bhv_AllTrlRng(:,BhvFactorInd);
% 
% % Generate line legend text based on time range
% LineLegTxt = sprintf('%0.1f->%0.1f', TimRange(1), TimRange(2));
% 
% % Calculate the correlation between mean scores and compression index
% [a,p] = obj.ManData.Correlation(MeanBelief,BhvPerfTrl,0,'Modified_MannKendall');
% X=MeanBelief;Y=BhvPerfTrl;
% 
% % Correct p-value for multiple comparisons
% p = p * 1;
% 
% % Perform significance test based on adjusted p-value threshold
% pa = p < (0.05);
% Xlbl='Task belief encoding';Ylbl='Avg Behavioral Performance';
% PlotType='Scatter';
% %% do the permutation statistical test for differences
% Bhv_AllTrlRng_AllTrls=cell2mat(arrayfun(@(x) AllFactorVals_1ndD_avg(x).TrialRange(:,BhvFactorInd)',1:nTrialRange,'uniformoutput',0));
% BeliefEncoding_AllTrls =cell2mat(arrayfun(@(x) mean(Scores_2ndD_avg_Flip(x).TrialRange(:,TimeInd),2)',1:nTrialRange,'uniformoutput',0));
% Ppercentile=15;
% Lower=prctile(BeliefEncoding_AllTrls,Ppercentile);
% Upper=prctile(BeliefEncoding_AllTrls,100-Ppercentile);
% UpperTrials=(BeliefEncoding_AllTrls>=Upper);
% LowerTrials=(BeliefEncoding_AllTrls<=Lower);
% BhvUpperTrials=Bhv_AllTrlRng_AllTrls(UpperTrials);
% BhvLowerTrials=Bhv_AllTrlRng_AllTrls(LowerTrials);
% [p]=obj.ManData.DiffPermutationStatTestNonPaired(BhvLowerTrials,BhvUpperTrials);
% pa = p < (0.05);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% obj.PlotSVMProjectionsScores2DPEV_Learning(ClassifierResults_ndD{1},ClassifierResults_ndD{3},ClassifierOpts,...
%     Sp(6),Conds(c),'Line','PreStimScoreCorr',[0 3],'ThisTargetFactor','Quadrants','performtrend_stattest',1);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% CorrBeliefMotorRep
% obj.PlotSVMProjectionsScores2DPEV_Learning(ClassifierResults_ndD{1},ClassifierResults_ndD{3},ClassifierOpts,...
%     Sp(1),Conds(c),'Line','BeleifRespScoreCorr',[0 3],'ThisTargetFactor','ResponseLoc','performtrend_stattest',1);
% 
% obj.PlotSVMProjectionsScores2DPEV_Learning(ClassifierResults_ndD{2},ClassifierResults_ndD{3},ClassifierOpts,...
%     Sp(2),Conds(c),'Line','BeleifRespScoreCorr',[2 3],'ThisTargetFactor','ResponseLoc','performtrend_stattest',1);
% 
% if contains(TargetFactorDim1,'Shape') | contains(TargetFactorDim1,'Color')
%     warning('first dimension has to be response location');h=[];return
% end
% TimRange=sort( ClassifierOpts.SpkCountPeriod(3,:),2,'descend'); % look at different time points 0 -0.1;-0.1 -0.2;-0.2 -0.3;-0.3 -0.4;-0.4 -0.6;
% nTimRange=size(TimRange,1);
% TimRangeResp=sort( ClassifierOpts.SpkCountPeriod(1,:),2,'descend'); % look at time points for response
% 
% %% below code is when we correlate with each individual trial
% if ~AnalysisOpts.PopulationAna.UseAvg_BeleifRespScoreCorr
%     % correlate magnetude of representation
%     Scores_1ndD_AllTrlRng=obj.ManData.ReshapeStruct2Mat(ScoresMag_1ndD_RS,'TrialRange',62);
%     %   if contains(TargetFactor,'InCong');Scores_1ndD_AllTrlRng=-1*Scores_1ndD_AllTrlRng;end
%     %    Scores_1ndD_AllTrlRng=obj.ManData.ReshapeStruct2Mat(Scores_1ndD_RS,'TrialRange',62);
%     for k=1:nTimRange
%         LineLegTxt{k}=sprintf('%0.1f->%0.1f',TimRange(k,1),TimRange(k,2));
%         TimeInd=AnalysisOpts.Time<=TimRange(k,1) & AnalysisOpts.Time>=TimRange(k,2);
%         Scores_2ndD_AllTrlRng=obj.ManData.ReshapeStruct2Mat(Scores_2ndD_RS,'TrialRange',62);
%         [ScoresMetric{k},p]=obj.ManData.Correlation(mean(Scores_2ndD_AllTrlRng(:,TimeInd),2),...
%             mean(Scores_1ndD_AllTrlRng,1),0,'Modified_MannKendall');
%         p=p*AnalysisOpts.NTim; % correct for multiple comparisions
%         pa{k}=p<(0.05);
%     end
% else
%     %% below code is when we correlate with average
%     % Retrieve scores for rule encoding from the averaged 2nd dimension scores
%     MeanScoresMag_1ndD_RS=arrayfun(@(x) mean(ScoresMag_1ndD_RS(x).TrialRange,1)',1:nTrialRange,'uniformoutput',0);
%     AvgRespEncoding=cell2mat(MeanScoresMag_1ndD_RS)';
% 
%     BeliefEncoding = obj.ManData.ReshapeStruct2Mat(Scores_2ndD_avg_avg_Flip, 'TrialRange', 2);
% 
%     TimeInd=AnalysisOpts.Time<=TimRange(1) & AnalysisOpts.Time>=TimRange(2);
%     TimeIndResp=AnalysisOpts.Time<=TimRangeResp(1) & AnalysisOpts.Time>=TimRangeResp(2);
% 
%     MeanBelief= mean(BeliefEncoding(:,TimeInd),2);
%     MeanResponseEncoding=mean(AvgRespEncoding(:,TimeIndResp),2);
%     [a,p] = obj.ManData.Correlation(MeanBelief,MeanResponseEncoding,0,'Modified_MannKendall');
%     X=MeanBelief;Y=MeanResponseEncoding;
% 
%     % Correct p-value for multiple comparisons
%     p = p * 1;
%     % Perform significance test based on adjusted p-value threshold
%     pa = p < (0.05);
%     Xlbl='Task belief encoding';Ylbl='Response Location Encoding';
%     PlotType='Scatter';
% 
% end