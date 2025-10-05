function [BeliefStimCorr,CPIStimEncodingCorr,BeliefCPIcorr,BeliefBhvPerfcorr]=ProcessCorrelationStatTests(ClassifierOpts,AvgDim2)
% computes correlations and statistical tests for figures 4 and 5

global AnalysisOpts
SetAnalysisOptions_RuleRepresentation

fprintf('\nComputing correlation stats data from classifier results')

AnalysisOpts.AreaNum=1;
AnalysisOpts.CalShuffTrlOrderClassifier_ClusterTH=24; % for all correlation analysis 
ManData=ManipulateData;
NeuAna=NeuralAnalysisFuncs;
FigParams=fig_params;

if strcmp(NeuAna.PSTHTimRef,'leading')
    ClassifierOpts.AnalysisOpts.Time=ClassifierOpts.AnalysisOpts.Time-ClassifierOpts.AnalysisOpts.PopulationAna.PSTHbin*0.001/2;
    warning('We have subtracted half of the PSTH bin from time vector')
end

if AvgDim2==1 % then we take 3 which is our default
    Ext='AvgDim2';
    AvgDim2Val=3;
elseif AvgDim2==0
    Ext='';
    AvgDim2Val=1;
else
    Ext=['AvgDim2_' num2str(AvgDim2)];
    AvgDim2Val=AvgDim2;
end

FileNameSyntax=['_CPI' ClassifierOpts.Name '_PFC_' AnalysisOpts.SpkCntStartFieldName '_' AnalysisOpts.TrlSpkTimeFieldName '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
[Scores_1ndD_obs,FileExists]=ManData.LoadVar('Classifier','Scores_1ndD_obs',FileNameSyntax,0,'WantedDate','ALL');
[Scores_2ndD_obs]     =ManData.LoadVar('Classifier','Scores_2ndD_obs',FileNameSyntax,0,'WantedDate','ALL');
[Scores_3ndD_obs]     =ManData.LoadVar('Classifier','Scores_3ndD_obs',FileNameSyntax,0,'WantedDate','ALL');
[Scores_1ndD_sh]      =ManData.LoadVar('Classifier','Scores_1ndD_sh',FileNameSyntax,0,'WantedDate','ALL');
[Scores_2ndD_sh]      =ManData.LoadVar('Classifier','Scores_2ndD_sh',FileNameSyntax,0,'WantedDate','ALL');
[Scores_3ndD_sh]      =ManData.LoadVar('Classifier','Scores_3ndD_sh',FileNameSyntax,0,'WantedDate','ALL');
[TargetFactorDim1]    =ManData.LoadVar('Classifier','TargetFactorDim1',FileNameSyntax,0,'WantedDate','ALL');
[TargetFactorDim2]    =ManData.LoadVar('Classifier','TargetFactorDim2',FileNameSyntax,0,'WantedDate','ALL');
[TargetFactorDim3]    =ManData.LoadVar('Classifier','TargetFactorDim3',FileNameSyntax,0,'WantedDate','ALL');
AllFactorVals_1ndD_obs=ManData.LoadVar('Classifier','AllFactorVals_1ndD_obs',FileNameSyntax,0,'WantedDate','ALL');
AllFactorVals_2ndD_obs=ManData.LoadVar('Classifier','AllFactorVals_2ndD_obs',FileNameSyntax,0,'WantedDate','ALL');
AllFactorVals_3ndD_obs=ManData.LoadVar('Classifier','AllFactorVals_3ndD_obs',FileNameSyntax,0,'WantedDate','ALL');
AllFactorVals_1ndD_sh =ManData.LoadVar('Classifier','AllFactorVals_1ndD_sh',FileNameSyntax,0,'WantedDate','ALL');
AllFactorVals_2ndD_sh =ManData.LoadVar('Classifier','AllFactorVals_2ndD_sh',FileNameSyntax,0,'WantedDate','ALL');
AllFactorVals_3ndD_sh =ManData.LoadVar('Classifier','AllFactorVals_3ndD_sh',FileNameSyntax,0,'WantedDate','ALL');
ThisFactorResp1_obs=ManData.LoadVar('Classifier','ThisFactorResp1_obs',FileNameSyntax,0,'WantedDate','ALL');
ThisFactorResp2_obs=ManData.LoadVar('Classifier','ThisFactorResp2_obs',FileNameSyntax,0,'WantedDate','ALL');
ThisFactorResp3_obs=ManData.LoadVar('Classifier','ThisFactorResp3_obs',FileNameSyntax,0,'WantedDate','ALL');
ThisFactorResp1_sh =ManData.LoadVar('Classifier','ThisFactorResp1_sh',FileNameSyntax,0,'WantedDate','ALL');
ThisFactorResp2_sh =ManData.LoadVar('Classifier','ThisFactorResp2_sh',FileNameSyntax,0,'WantedDate','ALL');
ThisFactorResp3_sh =ManData.LoadVar('Classifier','ThisFactorResp3_sh',FileNameSyntax,0,'WantedDate','ALL');

%% load also compression index data
[~,TargetFactorDim]=arrayfun(@(dim) NeuAna.getClassifierDimInfo(ClassifierOpts,dim),[0 1 3],'UniformOutput',0);
FileNameSpecs=@(x) [ClassifierOpts.Name '_' TargetFactorDim{x}];
ThisFileNameSpecs=FileNameSpecs(1);
PATH=AnalysisOpts.CodeTestDataPath;
LoadFileName=[ThisFileNameSpecs '_D1ShData'];
Metric='CPI';
CPIData=load([PATH LoadFileName],[Metric '_Obsv'],[Metric '_Shuff']);
if ~(isfield(CPIData,[Metric '_Obsv']) & isfield(CPIData,[Metric '_Shuff']))
    fprintf('\nWe dont have the data to process this condition')
    return
end
% check if we have CPI data 
ExistCPIdata=~isempty(CPIData.CPI_Obsv);
%% load accuracy for shffule and observed
for ThisDim=1:3
    ThisFileNameSpecs=FileNameSpecs(ThisDim);
    LoadFileName=[ThisFileNameSpecs '_D' num2str(ThisDim) 'ShData'];
    AccuracyData{ThisDim}=load([PATH LoadFileName],'Accuracy_Obsv','Accuracy_Shuff');
end

%% now extract all of the extra information we need to process
NrepShufperFold=length(Scores_1ndD_sh);
% get extra scores var for observed and shuffle
%ExtraScore12_Obsv=arrayfun(@(x) ExtractExtraScoreData4Shuffle(ClassifierOpts,AllFactorVals_1ndD_obs{end},AllFactorVals_2ndD_obs{end},...
%    Scores_1ndD_obs{end},Scores_2ndD_obs{end},ThisFactorResp1_obs{end},TargetFactorDim1,TargetFactorDim2,'Quadrants'),1,'uniformoutput',0); % for observed
if strcmp(ClassifierOpts.TargetFactors{1},'ResponseLoc')
    TargetFactor='ResponseLoc';
else
    TargetFactor='Quadrants';
end

ExtraScore13_Obsv=arrayfun(@(x) ExtractExtraScoreData4Shuffle(ClassifierOpts,AllFactorVals_1ndD_obs{end},AllFactorVals_3ndD_obs{end},...
    Scores_1ndD_obs{end},Scores_3ndD_obs{end},ThisFactorResp1_obs{end},TargetFactorDim1,TargetFactorDim3,TargetFactor),1,'uniformoutput',0); % for observed

ExtraScore23_Obsv=arrayfun(@(x) ExtractExtraScoreData4Shuffle(ClassifierOpts,AllFactorVals_2ndD_obs{end},AllFactorVals_3ndD_obs{end},...
    Scores_2ndD_obs{end},Scores_3ndD_obs{end},ThisFactorResp2_obs{end},TargetFactorDim2,TargetFactorDim3,TargetFactor),1,'uniformoutput',0); % for observed

%%% CHANGE THE REP
%ExtraScore12_Shuff=arrayfun(@(x) ExtractExtraScoreData4Shuffle(ClassifierOpts,AllFactorVals_1ndD_sh{x},AllFactorVals_2ndD_sh{x},...
%    Scores_1ndD_sh{x},Scores_2ndD_sh{x},ThisFactorResp1_sh{x},TargetFactorDim1,TargetFactorDim2,'Quadrants'),1:NrepShufperFold,'uniformoutput',0); % for shuffle

ExtraScore13_Shuff=arrayfun(@(x) ExtractExtraScoreData4Shuffle(ClassifierOpts,AllFactorVals_1ndD_sh{x},AllFactorVals_3ndD_sh{x},...
    Scores_1ndD_sh{x},Scores_3ndD_sh{x},ThisFactorResp1_sh{x},TargetFactorDim1,TargetFactorDim3,TargetFactor),1:NrepShufperFold,'uniformoutput',0); % for shuffle

ExtraScore23_Shuff=arrayfun(@(x) ExtractExtraScoreData4Shuffle(ClassifierOpts,AllFactorVals_2ndD_sh{x},AllFactorVals_3ndD_sh{x},...
    Scores_2ndD_sh{x},Scores_3ndD_sh{x},ThisFactorResp2_sh{x},TargetFactorDim2,TargetFactorDim3,TargetFactor),1:NrepShufperFold,'uniformoutput',0); % for shuffle

%% process each condition of correlation for observed and shuffle
opts.CorrelationMetric='Kendall';
ZscoreData=0;

%% if this is for response correlation we only need one correlation of response representation and belief 
if strcmp(TargetFactor,'ResponseLoc')
    [a_obsv.BeliefRespCorr,~,a_obsv.BeliefRespCorr_X,a_obsv.BeliefRespCorr_Y]= (arrayfun(@(x) CalBeliefResponseEncCorr(ManData,NeuAna,AvgDim2Val,ClassifierOpts,TargetFactorDim1,ExtraScore13_Obsv{x}.ScoresMag_1ndD_RS,ExtraScore13_Obsv{x}.Scores_2ndD_avg_avg_Flip,opts,1),1,'uniformoutput',0)); % calculates correlation between belief and CPI
    [a_sh.BeliefRespCorr  ,~,a_sh.BeliefRespCorr_X,a_sh.BeliefRespCorr_Y]= (arrayfun(@(x) (CalBeliefResponseEncCorr(ManData,NeuAna,AvgDim2Val,ClassifierOpts,TargetFactorDim1,ExtraScore13_Shuff{x}.ScoresMag_1ndD_RS,ExtraScore13_Shuff{x}.Scores_2ndD_avg_avg_Flip,opts,0)),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
    a_obsv.BeliefRespCorr=cell2mat(a_obsv.BeliefRespCorr)';
    a_sh.BeliefRespCorr=ManData.ReshapeCell2Mat(a_sh.BeliefRespCorr,1);
    aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.BeliefRespCorr,a_sh.BeliefRespCorr);
    PvalShuff=ManData.CalpValShuffle(a_sh.BeliefRespCorr',a_obsv.BeliefRespCorr);
    title(sprintf('Corr Resp Encoding and Belief a=%0.3f,aZ=%0.3f,p=%0.3f',a_obsv.BeliefRespCorr,aCorrected_obsv,PvalShuff) )
    BeliefStimCorr=[];CPIStimEncodingCorr=[];BeliefCPIcorr=[];BeliefBhvPerfcorr=[];
else
    % plot correlation of belief and color category in time
    BeliefStimCorr=figure;
    subplot(241)
    [a_obsv.BeliefStimCorr_Color,~,a_obsv.BeliefStimCorr_Color_X,a_obsv.BeliefStimCorr_Color_Y,a_obsv.BeliefStimCorr_Color_avg]= (arrayfun(@(x) CalBeliefStimCorr(ManData,NeuAna,AvgDim2Val,ClassifierOpts,TargetFactorDim2,ExtraScore23_Obsv{x}.ScoresMag_1ndD_RS,ExtraScore23_Obsv{x}.Scores_2ndD_avg_avg_Flip,opts,1),1,'uniformoutput',0)); % calculates correlation between belief and CPI
    [a_sh.BeliefStimCorr_Color,~,a_sh.BeliefStimCorr_Color_X,a_sh.BeliefStimCorr_Color_Y,a_sh.BeliefStimCorr_Color_avg]= (arrayfun(@(x)  (CalBeliefStimCorr(ManData,NeuAna,AvgDim2Val,ClassifierOpts,TargetFactorDim2,ExtraScore23_Shuff{x}.ScoresMag_1ndD_RS,ExtraScore23_Shuff{x}.Scores_2ndD_avg_avg_Flip,opts,0)),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
    a_obsv.BeliefStimCorr_Color=cell2mat(a_obsv.BeliefStimCorr_Color)';
    a_sh.BeliefStimCorr_Color=ManData.ReshapeCell2Mat(a_sh.BeliefStimCorr_Color,1);
    [a_obsv.BeliefStimCorr_Color,a_sh.BeliefStimCorr_Color]=ManData.CalZscoreShuffle(a_obsv.BeliefStimCorr_Color,a_sh.BeliefStimCorr_Color);
    ClusterCorrectionTJB(a_obsv.BeliefStimCorr_Color,a_sh.BeliefStimCorr_Color,{[],[2,4,2]},'Time',ClassifierOpts.AnalysisOpts.Time,'ThresholdPercentage',AnalysisOpts.CalShuffTrlOrderClassifier_ClusterTH,'DoubleSided',1,'ZscoreData',ZscoreData);
    axis square

    aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.BeliefStimCorr_Color,a_sh.BeliefStimCorr_Color);
    PvalShuff=ManData.CalpValShuffle(a_sh.BeliefStimCorr_Color',a_obsv.BeliefStimCorr_Color);
    PvalShuff_avg_Col=ManData.CalpValShuffle(cell2mat(a_sh.BeliefStimCorr_Color_avg)',cell2mat(a_obsv.BeliefStimCorr_Color_avg));
    subplot(243)
    plot(ClassifierOpts.AnalysisOpts.Time,PvalShuff);
    xlabel('Time');ylabel('pval')
    hold on
    plot([ClassifierOpts.AnalysisOpts.Time(1) ClassifierOpts.AnalysisOpts.Time(end)],[0.05 0.05],'--g')
    yyaxis right
    plot(ClassifierOpts.AnalysisOpts.Time,aCorrected_obsv)
    ylabel('tau corrected')
    axis square
    title('Correlation between belief and Color stimulus encoding')
    % plot zscore of the correlation in time
    subplot(244)
    ZBeliefStimCorr_Color=ManData.CalZscoreShuffle(a_obsv.BeliefStimCorr_Color,a_sh.BeliefStimCorr_Color);
    plot(ClassifierOpts.AnalysisOpts.Time,ZBeliefStimCorr_Color)
    ylabel('tau zscored')
    axis square
    title('Z Corr between belief and Color stimulus encoding')

    % plot correlation of belief and shape category in time
    subplot(245)
    [a_obsv.BeliefStimCorr_Shape,~,a_obsv.BeliefStimCorr_Shape_X,a_obsv.BeliefStimCorr_Shape_Y,a_obsv.BeliefStimCorr_Shape_avg,a_obsv.BeliefStimCorr_Shape_avg_X,a_obsv.BeliefStimCorr_Shape_avg_Y]= (arrayfun(@(x) CalBeliefStimCorr(ManData,NeuAna,AvgDim2Val,ClassifierOpts,TargetFactorDim1,ExtraScore13_Obsv{x}.ScoresMag_1ndD_RS,ExtraScore13_Obsv{x}.Scores_2ndD_avg_avg_Flip,opts,1),1,'uniformoutput',0)); % calculates correlation between belief and CPI
    [a_sh.BeliefStimCorr_Shape  ,~,a_sh.BeliefStimCorr_Shape_X,a_sh.BeliefStimCorr_Shape_Y,a_sh.BeliefStimCorr_Shape_avg,a_sh.BeliefStimCorr_Shape_avg_X,a_sh.BeliefStimCorr_Shape_avg_Y]= (arrayfun(@(x) (CalBeliefStimCorr(ManData,NeuAna,AvgDim2Val,ClassifierOpts,TargetFactorDim1,ExtraScore13_Shuff{x}.ScoresMag_1ndD_RS,ExtraScore13_Shuff{x}.Scores_2ndD_avg_avg_Flip,opts,0)),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
    a_obsv.BeliefStimCorr_Shape=cell2mat(a_obsv.BeliefStimCorr_Shape)';
    a_sh.BeliefStimCorr_Shape=ManData.ReshapeCell2Mat(a_sh.BeliefStimCorr_Shape,1);
    [a_obsv.BeliefStimCorr_Shape,a_sh.BeliefStimCorr_Shape]=ManData.CalZscoreShuffle(a_obsv.BeliefStimCorr_Shape,a_sh.BeliefStimCorr_Shape);
    ClusterCorrectionTJB(a_obsv.BeliefStimCorr_Shape,a_sh.BeliefStimCorr_Shape,{[],[2,4,6]},'Time',ClassifierOpts.AnalysisOpts.Time,'ThresholdPercentage',AnalysisOpts.CalShuffTrlOrderClassifier_ClusterTH,'DoubleSided',1,'ZscoreData',ZscoreData);
    axis square

    aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.BeliefStimCorr_Shape,a_sh.BeliefStimCorr_Shape);
    PvalShuff=ManData.CalpValShuffle(a_sh.BeliefStimCorr_Shape',a_obsv.BeliefStimCorr_Shape);

    subplot(247)
    plot(ClassifierOpts.AnalysisOpts.Time,PvalShuff);
    xlabel('Time');ylabel('pval')
    hold on
    plot([ClassifierOpts.AnalysisOpts.Time(1) ClassifierOpts.AnalysisOpts.Time(end)],[0.05 0.05],'--g')
    yyaxis right
    plot(ClassifierOpts.AnalysisOpts.Time,aCorrected_obsv)
    ylabel('tau corrected')
    axis square
    title('Correlation between belief and Shape stimulus encoding')

    subplot(248)
    ZBeliefStimCorr_Shape=ManData.CalZscoreShuffle(a_obsv.BeliefStimCorr_Shape,a_sh.BeliefStimCorr_Shape);
    plot(ClassifierOpts.AnalysisOpts.Time,ZBeliefStimCorr_Shape(:,1))
    ylabel('tau zscored')
    axis square
    title('Z Corr between belief and Shape stimulus encoding')
     
    % plot average correlation with shape
    AvgCorrBeliefShape=figure;
    PvalShuff_avg_Shp=ManData.CalpValShuffle(cell2mat(a_sh.BeliefStimCorr_Shape_avg)',cell2mat(a_obsv.BeliefStimCorr_Shape_avg));
    aCorrected_obsv_avg_Shp=ManData.CalZscoreShuffle(cell2mat(a_obsv.BeliefStimCorr_Shape_avg),cell2mat(a_sh.BeliefStimCorr_Shape_avg));
    
    scatter(cell2mat(a_obsv.BeliefStimCorr_Shape_avg_X),cell2mat(a_obsv.BeliefStimCorr_Shape_avg_Y));
    xlabel('Average belief Encoding')
    ylabel('Average stim encoding shape')
    axis square
    title(sprintf('Corr shape Stim Encoding and belief a=%0.3f,aZ=%0.3f,p=%0.3f',cell2mat(a_obsv.BeliefStimCorr_Shape_avg), ...
        aCorrected_obsv_avg_Shp,PvalShuff_avg_Shp) )

    %% correlation of CPI with stim encoding
    if ExistCPIdata
        %[1 2],[2 1] %[red bunny, green tee, green bunny, red tee]
        CPIStimEncodingCorr=figure;
        [a_obsv.CPIStimEncodingCorr,~,a_obsv.CPIStimEncodingCorr_X,a_obsv.CPIStimEncodingCorr_Y]=(arrayfun(@(x) CalCPIStimEncodingCorr(ManData,ClassifierOpts,TargetFactorDim2,ExtraScore23_Obsv{x}.ScoresMag_1ndD_RS,CPIData.CPI_Obsv{1},opts,1),1,'uniformoutput',0)); % calculates correlation between belief and CPI
        [a_sh.CPIStimEncodingCorr  ,~,a_sh.CPIStimEncodingCorr_X  ,a_sh.CPIStimEncodingCorr_Y]=(arrayfun(@(x) (CalCPIStimEncodingCorr(ManData,ClassifierOpts,TargetFactorDim2,ExtraScore23_Shuff{x}.ScoresMag_1ndD_RS,CPIData.CPI_Shuff{1}(:,:,x),opts,0)),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
        a_obsv.CPIStimEncodingCorr=cell2mat(a_obsv.CPIStimEncodingCorr)';
        a_sh.CPIStimEncodingCorr=cell2mat(a_sh.CPIStimEncodingCorr);
        aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.CPIStimEncodingCorr,a_sh.CPIStimEncodingCorr);
        PvalShuff=ManData.CalpValShuffle(a_sh.CPIStimEncodingCorr',a_obsv.CPIStimEncodingCorr);
        title(sprintf('Corr Stim Encoding and CPI a=%0.3f,aZ=%0.3f,p=%0.3f',a_obsv.CPIStimEncodingCorr,aCorrected_obsv,PvalShuff) )
    end
    %% correlation of color behavior and color encoding
    BhvStimEncodingCorr=figure;
    [a_obsv.BhvStimEncodingCorr,~,a_obsv.BhvStimEncodingCorr_X,a_obsv.BhvStimEncodingCorr_Y]=(arrayfun(@(x) CalBhvStimEncodingCorr(ManData,ClassifierOpts,TargetFactorDim2,ExtraScore23_Obsv{x}.ScoresMag_1ndD_RS,ExtraScore23_Obsv{x}.AllFactorVals_1ndD_avg,opts,1,'Color',10),1,'uniformoutput',0)); % calculates correlation between belief and CPI
    [a_sh.BhvStimEncodingCorr  ,~,a_sh.BhvStimEncodingCorr_X  ,a_sh.BhvStimEncodingCorr_Y]=(arrayfun(@(x) (CalBhvStimEncodingCorr(ManData,ClassifierOpts,TargetFactorDim2,ExtraScore23_Shuff{x}.ScoresMag_1ndD_RS,ExtraScore23_Shuff{x}.AllFactorVals_1ndD_avg,opts,0,'Color',10)),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
    a_obsv.BhvStimEncodingCorr=cell2mat(a_obsv.BhvStimEncodingCorr)';
    a_sh.BhvStimEncodingCorr=cell2mat(a_sh.BhvStimEncodingCorr);
    aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.BhvStimEncodingCorr,a_sh.BhvStimEncodingCorr);
    PvalShuff=ManData.CalpValShuffle(a_sh.BhvStimEncodingCorr',a_obsv.BhvStimEncodingCorr);
    title(sprintf('Corr Stim Encoding and Color behavior a=%0.3f,aZ=%0.3f,p=%0.3f',a_obsv.BhvStimEncodingCorr,aCorrected_obsv,PvalShuff) )


    %% calculate correlation of Belief with CPI
    if ExistCPIdata
        BeliefCPIcorr=figure;
        [a_obsv.BeliefCPIcorr,~,a_obsv.BeliefCPIcorr_X,a_obsv.BeliefCPIcorr_Y,a_obsv.BeliefCPIcorr_Time]=(arrayfun(@(x) CalBeliefCPIcorr(ManData,NeuAna,ClassifierOpts,TargetFactorDim2,ExtraScore23_Obsv{x}.Scores_2ndD_avg_avg_Flip,CPIData.CPI_Obsv{1},opts,1),1,'uniformoutput',0)); % calculates correlation between belief and CPI
        [a_sh.BeliefCPIcorr  ,~,a_sh.BeliefCPIcorr_X  ,a_sh.BeliefCPIcorr_Y,a_sh.BeliefCPIcorr_Time]=(arrayfun(@(x) (CalBeliefCPIcorr(ManData,NeuAna,ClassifierOpts,TargetFactorDim2,ExtraScore23_Shuff{x}.Scores_2ndD_avg_avg_Flip,CPIData.CPI_Shuff{1}(:,:,x),opts,0)),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
        a_obsv.BeliefCPIcorr=cell2mat(a_obsv.BeliefCPIcorr)';a_sh.BeliefCPIcorr=cell2mat(a_sh.BeliefCPIcorr);
        aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.BeliefCPIcorr,a_sh.BeliefCPIcorr);
        PvalShuff=ManData.CalpValShuffle(a_sh.BeliefCPIcorr',a_obsv.BeliefCPIcorr);%,'pvaltail','right');
        subplot(121)
        title(sprintf('Corr Belief and CPI a=%0.3f,aZ=%0.3f,p=%0.3f',a_obsv.BeliefCPIcorr,aCorrected_obsv,PvalShuff))

        % plot correaltion of belief and CPI in time
        a_obsv.BeliefCPIcorr_Time=ManData.ReshapeCell2Mat(a_obsv.BeliefCPIcorr_Time,1);
        a_sh.BeliefCPIcorr_Time=ManData.ReshapeCell2Mat(a_sh.BeliefCPIcorr_Time,1);
        [a_obsv.BeliefCPIcorr_Time,a_sh.BeliefCPIcorr_Time]=ManData.CalZscoreShuffle(a_obsv.BeliefCPIcorr_Time,a_sh.BeliefCPIcorr_Time);
        ClusterCorrectionTJB(a_obsv.BeliefCPIcorr_Time,a_sh.BeliefCPIcorr_Time,{[],[1,2,2]},'Time',ClassifierOpts.AnalysisOpts.Time,'ThresholdPercentage',AnalysisOpts.CalShuffTrlOrderClassifier_ClusterTH,'DoubleSided',1,'ZscoreData',ZscoreData);
    end
    %% calculate correlation of Belief with behavioral performance in color and shape
    BeliefBhvPerfcorr=figure;
    subplot(221)
    [a_obsv.BeliefBhvPerfcorr_Color,~,a_obsv.BeliefBhvPerfcorr_Color_X,a_obsv.BeliefBhvPerfcorr_Color_Y,a_obsv.BeliefBhvPerfcorr_Color_Time]=(arrayfun(@(x) CalBeliefBhvPerfcorr(ManData,NeuAna,ClassifierOpts,TargetFactorDim2,ExtraScore23_Obsv{x}.Scores_2ndD_avg_avg_Flip,ExtraScore23_Obsv{x}.AllFactorVals_1ndD_avg,opts,1,'Color',10),1,'uniformoutput',0)); % calculates correlation between belief and CPI
    [a_sh.BeliefBhvPerfcorr_Color  ,~,a_sh.BeliefBhvPerfcorr_Color_X  ,a_sh.BeliefBhvPerfcorr_Color_Y,a_sh.BeliefBhvPerfcorr_Color_Time]  =(arrayfun(@(x)  (CalBeliefBhvPerfcorr(ManData,NeuAna,ClassifierOpts,TargetFactorDim2,ExtraScore23_Shuff{x}.Scores_2ndD_avg_avg_Flip,ExtraScore23_Shuff{x}.AllFactorVals_1ndD_avg,opts,0,'Color',10)),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
    a_obsv.BeliefBhvPerfcorr_Color=cell2mat(a_obsv.BeliefBhvPerfcorr_Color)';
    a_sh.BeliefBhvPerfcorr_Color=cell2mat(a_sh.BeliefBhvPerfcorr_Color);
    aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.BeliefBhvPerfcorr_Color,a_sh.BeliefBhvPerfcorr_Color);
    PvalShuff=ManData.CalpValShuffle(a_sh.BeliefBhvPerfcorr_Color',a_obsv.BeliefBhvPerfcorr_Color);%,'pvaltail','right');
    title(sprintf('Corr Belief and Bhv perf a=%0.3f,aZ=%0.3f,p=%0.3f',a_obsv.BeliefBhvPerfcorr_Color,aCorrected_obsv,PvalShuff))

    subplot(222)
    [a_obsv.BeliefBhvPerfcorr_Shape,~,a_obsv.BeliefBhvPerfcorr_Shape_X,a_obsv.BeliefBhvPerfcorr_Shape_Y,a_obsv.BeliefBhvPerfcorr_Shape_Time]= (arrayfun(@(x) CalBeliefBhvPerfcorr(ManData,NeuAna,ClassifierOpts,TargetFactorDim1,ExtraScore13_Obsv{x}.Scores_2ndD_avg_avg_Flip,ExtraScore13_Obsv{x}.AllFactorVals_1ndD_avg,opts,1,'Shape',10),1,'uniformoutput',0)); % calculates correlation between belief and CPI
    [a_sh.BeliefBhvPerfcorr_Shape  ,~,a_sh.BeliefBhvPerfcorr_Shape_X  ,a_sh.BeliefBhvPerfcorr_Shape_Y,a_sh.BeliefBhvPerfcorr_Shape_Time]    = (arrayfun(@(x) (CalBeliefBhvPerfcorr(ManData,NeuAna,ClassifierOpts,TargetFactorDim1,ExtraScore13_Shuff{x}.Scores_2ndD_avg_avg_Flip,ExtraScore13_Shuff{x}.AllFactorVals_1ndD_avg,opts,0,'Shape',10)),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
    a_obsv.BeliefBhvPerfcorr_Shape=cell2mat(a_obsv.BeliefBhvPerfcorr_Shape)';
    a_sh.BeliefBhvPerfcorr_Shape=cell2mat(a_sh.BeliefBhvPerfcorr_Shape);
    aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.BeliefBhvPerfcorr_Shape,a_sh.BeliefBhvPerfcorr_Shape);
    PvalShuff=ManData.CalpValShuffle(a_sh.BeliefBhvPerfcorr_Shape',a_obsv.BeliefBhvPerfcorr_Shape);%,'pvaltail','left');
    title(sprintf('Corr Blf and Bhv perf a=%0.3f,aZ=%0.3f,p=%0.3f, %0.3f-%0.3f',a_obsv.BeliefBhvPerfcorr_Shape,aCorrected_obsv,PvalShuff,ClassifierOpts.SpkCountPeriod(3,1),ClassifierOpts.SpkCountPeriod(3,2)))

    %% plot correlation of belief and bhv performance in time
    BeliefTime=logical(ones(1,116));%ClassifierOpts.AnalysisOpts.Time<=0.001 & ClassifierOpts.AnalysisOpts.Time>=-0.4;%
    a_obsv.BeliefBhvPerfcorr_Color_Time=cell2mat(a_obsv.BeliefBhvPerfcorr_Color_Time);
    a_sh.BeliefBhvPerfcorr_Color_Time=cell2mat(a_sh.BeliefBhvPerfcorr_Color_Time);
    [a_obsv.BeliefBhvPerfcorr_Color_Time,a_sh.BeliefBhvPerfcorr_Color_Time]=ManData.CalZscoreShuffle(a_obsv.BeliefBhvPerfcorr_Color_Time,a_sh.BeliefBhvPerfcorr_Color_Time);
    ClusterCorrectionTJB(a_obsv.BeliefBhvPerfcorr_Color_Time(BeliefTime),a_sh.BeliefBhvPerfcorr_Color_Time(BeliefTime,:),{[],[2,2,3]},'Time',ClassifierOpts.AnalysisOpts.Time(BeliefTime),'ThresholdPercentage',AnalysisOpts.CalShuffTrlOrderClassifier_ClusterTH,'DoubleSided',1);
    axis tight
    %PvalShuff=ManData.CalpValShuffle(a_sh.BeliefBhvPerfcorr_Shape_Time',a_obsv.BeliefBhvPerfcorr_Shape_Time);

    a_obsv.BeliefBhvPerfcorr_Shape_Time=cell2mat(a_obsv.BeliefBhvPerfcorr_Shape_Time);
    a_sh.BeliefBhvPerfcorr_Shape_Time=cell2mat(a_sh.BeliefBhvPerfcorr_Shape_Time);
    [a_obsv.BeliefBhvPerfcorr_Shape_Time,a_sh.BeliefBhvPerfcorr_Shape_Time]=ManData.CalZscoreShuffle(a_obsv.BeliefBhvPerfcorr_Shape_Time,a_sh.BeliefBhvPerfcorr_Shape_Time);
    ClusterCorrectionTJB(a_obsv.BeliefBhvPerfcorr_Shape_Time(BeliefTime),a_sh.BeliefBhvPerfcorr_Shape_Time(BeliefTime,:),{[],[2,2,4]},'Time',ClassifierOpts.AnalysisOpts.Time(BeliefTime),'ThresholdPercentage',AnalysisOpts.CalShuffTrlOrderClassifier_ClusterTH,'DoubleSided',1);
    axis tight

    %% calculate correlation of normalized shape and color accuracy ratio and behavioral performance
    AccuRatioBhvPerfcorr=figure;
    [a_obsv.AccuRatioBhvPerfcorr,~,a_obsv.AccuRatioBhvPerfcorr_X,a_obsv.AccuRatioBhvPerfcorr_Y,a_obsv.AccuRatioBhvPerfcorr_Time]=(arrayfun(@(x) CalCorrColShapRatioBhvPerf(ManData,NeuAna,ClassifierOpts,TargetFactorDim1,TargetFactorDim2,ExtraScore23_Obsv{x}.AllFactorVals_1ndD_avg,AccuracyData{1}.Accuracy_Obsv{1},AccuracyData{2}.Accuracy_Obsv{1},opts,1,'Color',10,'Norm'),1,'uniformoutput',0)); % calculates correlation between belief and CPI
    [a_sh.AccuRatioBhvPerfcorr  ,~,a_sh.AccuRatioBhvPerfcorr_X  ,a_sh.AccuRatioBhvPerfcorr_Y,a_sh.AccuRatioBhvPerfcorr_Time]=(arrayfun(@(x)  (CalCorrColShapRatioBhvPerf(ManData,NeuAna,ClassifierOpts,TargetFactorDim1,TargetFactorDim2,ExtraScore23_Shuff{x}.AllFactorVals_1ndD_avg,AccuracyData{1}.Accuracy_Shuff{1}(:,:,x),AccuracyData{2}.Accuracy_Shuff{1}(:,:,x),opts,0,'Color',10,'Norm')),1:NrepShufperFold,'uniformoutput',0)); % calculates correlation between belief and CPI
    a_obsv.AccuRatioBhvPerfcorr=cell2mat(a_obsv.AccuRatioBhvPerfcorr)';
    a_sh.AccuRatioBhvPerfcorr=cell2mat(a_sh.AccuRatioBhvPerfcorr);
    aCorrected_obsv=ManData.CalZscoreShuffle(a_obsv.AccuRatioBhvPerfcorr,a_sh.AccuRatioBhvPerfcorr);
    PvalShuff=ManData.CalpValShuffle(a_sh.AccuRatioBhvPerfcorr',a_obsv.AccuRatioBhvPerfcorr);
    title(sprintf('Corr ratio Shp/Col and Bhv perf a=%0.3f,aZ=%0.3f,p=%0.3f',a_obsv.AccuRatioBhvPerfcorr,aCorrected_obsv,PvalShuff))

    % correlation of ratio in time 
    RatioTime=ClassifierOpts.AnalysisOpts.Time<=0.6 & ClassifierOpts.AnalysisOpts.Time>=-0.2;
    a_obsv.AccuRatioBhvPerfcorr_Time=cell2mat(a_obsv.AccuRatioBhvPerfcorr_Time);
    a_sh.AccuRatioBhvPerfcorr_Time=cell2mat(a_sh.AccuRatioBhvPerfcorr_Time);
    [a_obsv.AccuRatioBhvPerfcorr_Time,a_sh.AccuRatioBhvPerfcorr_Time]=ManData.CalZscoreShuffle(a_obsv.AccuRatioBhvPerfcorr_Time,a_sh.AccuRatioBhvPerfcorr_Time);
    ClusterCorrectionTJB(a_obsv.AccuRatioBhvPerfcorr_Time(RatioTime),a_sh.AccuRatioBhvPerfcorr_Time(RatioTime,:),{[],[1,3,3]},'Time',ClassifierOpts.AnalysisOpts.Time(RatioTime),'ThresholdPercentage',AnalysisOpts.CalShuffTrlOrderClassifier_ClusterTH,'DoubleSided',1);
    axis tight
    axis square
end
%% save off the data into the file again
SaveVarName=['CorrelationVar_'  Ext];

eval(sprintf('%s=struct(''%s'',%s,''%s'',%s);',SaveVarName,'a_obsv','a_obsv','a_sh','a_sh'));

ThisFileNameSpecs=FileNameSpecs(1);
LoadFileName=[ThisFileNameSpecs '_D' num2str(1) 'ShData'];
fprintf('\nSaving correlation metrics into file %s',LoadFileName)
save([PATH LoadFileName],SaveVarName,'-append');
save([PATH LoadFileName],'opts','-append');


%% save off figures 
mkdir([AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep ClassifierOpts.Name 'CorrStatTests'] )
FigParams.SaveCurrentFigs([ClassifierOpts.Name 'CorrStatTests'],[AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep  ClassifierOpts.Name 'CorrStatTests' filesep],'SaveEachFrame',1)

end