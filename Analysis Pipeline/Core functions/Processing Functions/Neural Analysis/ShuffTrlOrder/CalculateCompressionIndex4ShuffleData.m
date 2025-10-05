% calculates compression index and stat test for shuffle data
% you can load data from Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\Code Test Data\CompressionIndexShuffData100Reps.mat
%CalculateCompressionIndex4ShuffleData(ClassifierOpts,Scores_1ndD_sh,Scores_2ndD_sh,Scores_1ndD_obsv,Scores_2ndD_obsv,TargetFactorDim1,TargetFactorDim2,NrepShuff)
function [CPI_Obsv,CPI_Shuff,ShapeAvgDist_Obsv,ShapeAvgDist_Shuff,ColorAvgDist_Obsv,ColorAvgDist_Shuff,CompressionEncoding_Obsv,CompressionEncoding_Shuff]=CalculateCompressionIndex4ShuffleData(ClassifierOpts,PlotFlag,UseSavedData,varargin)

global AnalysisOpts
SetAnalysisOptions_RuleRepresentation

fprintf('\nCalculating Compression index data from classifier results')

AnalysisOpts.AreaNum=1;

ManData=ManipulateData;
obj=NeuralAnalysisFuncsTemp;
obj.ThisExtraTargetFactor='';
UseSavedData=1;

if ~isempty(varargin) % if we are aleading passing the data along
    Scores_1ndD_sh=varargin{1};Scores_2ndD_sh=varargin{2};
    Scores_1ndD_obsv=varargin{3};Scores_2ndD_obsv=varargin{4};
    TargetFactorDim1=varargin{5};
    TargetFactorDim2=varargin{6};
    NrepShuff=varargin{7};

elseif UseSavedData
     FileNameSyntax=['_CPI' ClassifierOpts.Name '_PFC_' AnalysisOpts.SpkCntStartFieldName '_' AnalysisOpts.TrlSpkTimeFieldName '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];    
     [Scores_1ndD_obs,FileExists]=ManData.LoadVar('Classifier','Scores_1ndD_obs',FileNameSyntax,0,'WantedDate','ALL');
     [Scores_2ndD_obs]=ManData.LoadVar('Classifier','Scores_2ndD_obs',FileNameSyntax,0,'WantedDate','ALL');
     [Scores_1ndD_sh]=ManData.LoadVar('Classifier','Scores_1ndD_sh',FileNameSyntax,0,'WantedDate','ALL');
     [Scores_2ndD_sh]=ManData.LoadVar('Classifier','Scores_2ndD_sh',FileNameSyntax,0,'WantedDate','ALL');
     [TargetFactorDim1]=ManData.LoadVar('Classifier','TargetFactorDim1',FileNameSyntax,0,'WantedDate','ALL');
     [TargetFactorDim2]=ManData.LoadVar('Classifier','TargetFactorDim2',FileNameSyntax,0,'WantedDate','ALL');
     if ~FileExists
         CPI_Obsv=[];CPI_Shuff=[];
         ShapeAvgDist_Obsv=[];ShapeAvgDist_Shuff=[];
         ColorAvgDist_Obsv=[];ColorAvgDist_Shuff=[];
         return
     end
     if isempty(Scores_1ndD_obs{1});Scores_1ndD_obs=Scores_1ndD_obs(end);end
     if isempty(Scores_2ndD_obs{1});Scores_2ndD_obs=Scores_2ndD_obs(end);end

     %% IF WE ARE USING PREVIOUS DATA 
%     FileNameSyntax=['_Shuff_'  'C1_' ClassifierOpts.Name '_PFC_' AnalysisOpts.SpkCntStartFieldName '_' AnalysisOpts.TrlSpkTimeFieldName '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
%     [ClassifierOpts]=ManData.LoadVar('Classifier','ClassifierOpts',FileNameSyntax,0,'WantedDate','ALL');
% 
%     for C=1:501
%         fprintf('\nLoading CPI data %i',C);
%         ExtraStrSaveCPI=['_CPI' num2str(C) '_Learning3D_Sh_Co_Ru_AltRule_RB_F_SkTrSq_MS_FEA_FL-75_50_4_50AT_' AnalysisOpts.Area_2look{1} '_' AnalysisOpts.SpkCntStartFieldName '_' AnalysisOpts.TrlSpkTimeFieldName '_' num2str(AnalysisOpts.PopulationAna.PSTHbin)];
%         Scores_1ndD=ManData.LoadVar('Classifier','Scores_1ndD',ExtraStrSaveCPI,0,'WantedDate','ALL');
%         Scores_2ndD=ManData.LoadVar('Classifier','Scores_2ndD',ExtraStrSaveCPI,0,'WantedDate','ALL');
%         if C==1
%             TargetFactorDim1=ManData.LoadVar('Classifier','TargetFactorDim1',ExtraStrSaveCPI,0,'WantedDate','ALL');
%             TargetFactorDim2=ManData.LoadVar('Classifier','TargetFactorDim2',ExtraStrSaveCPI,0,'WantedDate','ALL');
%         end
%         if C~=501
%             Scores_1ndD_sh(:,C)=Scores_1ndD(C); Scores_2ndD_sh(:,C)=Scores_2ndD(C);
%         else
%             Scores_1ndD_obs(:,1)=Scores_1ndD(C); Scores_2ndD_obs(:,1)=Scores_2ndD(C);
%         end
%     end
    NrepShuff=size(Scores_1ndD_sh,2);
else
    % process shuffle data first
    [Scores_1ndD_sh,Scores_2ndD_sh,TargetFactorDim1,TargetFactorDim2,Nrep,NrepShuff]=ExtractScores4ShuffleData(obj.ClassifierResults_Shuff{1},obj.ClassifierResults_Shuff{2},ClassifierOpts,1,[0 2],'Quadrants'); % plots PEV of two distibutions of scores
    % process observed data now
    [Scores_1ndD_obsv,Scores_2ndD_obsv,TargetFactorDim1,TargetFactorDim2]=ExtractScores4ShuffleData(obj.ClassifierResults_Observed{1},obj.ClassifierResults_Observed{2},ClassifierOpts,1,[0 2],'Quadrants'); % plots PEV of two distibutions of scores
end

% calculate compression index and shape and color encoding for shuffle
[ScoresMetricAll,EncodingDist_Shuff,CompressionEncoding_Shuff]=arrayfun(@(x) CalculateIndex('ALL',obj,ManData,'',Scores_1ndD_sh{x},Scores_2ndD_sh{x},TargetFactorDim1,TargetFactorDim2,x),1:NrepShuff,'UniformOutput',0);
ScoresMetricAll_Shuff=ManData.ReshapeCell2Mat(ScoresMetricAll,3);
CPI_Shuff=ScoresMetricAll_Shuff;

% calculate compression index and shape and color encoding for shuffle
[ScoresMetricAll_obsv,EncodingDist_Obsv,CompressionEncoding_Obsv]=arrayfun(@(x) CalculateIndex('ALL',obj,ManData,'',Scores_1ndD_obs{x},Scores_2ndD_obs{x},TargetFactorDim1,TargetFactorDim2,x),1,'UniformOutput',0);
ScoresMetricAll_obsv=ManData.ReshapeCell2Mat(ScoresMetricAll_obsv,3);
CPI_Obsv=ScoresMetricAll_obsv;

% get shape and color distance as well
ShapeAvgDist_Obsv=EncodingDist_Obsv{1}.ShapeDistAvg;
ColorAvgDist_Obsv=EncodingDist_Obsv{1}.ColorDistAvg;
ShapeAvgDist_Shuff=arrayfun(@(x) EncodingDist_Shuff{x}.ShapeDistAvg,1:NrepShuff,'UniformOutput',0);
ShapeAvgDist_Shuff=ManData.ReshapeCell2Mat(ShapeAvgDist_Shuff,3);
ColorAvgDist_Shuff=arrayfun(@(x) EncodingDist_Shuff{x}.ColorDistAvg,1:NrepShuff,'UniformOutput',0);
ColorAvgDist_Shuff=ManData.ReshapeCell2Mat(ColorAvgDist_Shuff,3);

if ~PlotFlag;return;end

Prmt=@(x) permute(x,[1 3 2]); % permute
Time=-0.55:0.01:0.6;%ClassifierOpts.AnalysisOpts.Time;
NTim=length(Time);

% smooth data in time
Observed=ManData.SmoothData(ScoresMetricAll_obsv,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
Shuffle=ManData.SmoothData(ScoresMetricAll_Shuff,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);

% smooth data in trials?
%Observed=ManData.SmoothData(Observed,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
%Shuffle=ManData.SmoothData(Shuffle,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);

NStg=16;
Shuffle=permute(Shuffle,[3 1 2]);
TrlRng=[1:NStg]';
col=copper(16);
figure;hold on;arrayfun(@(x) plot(Observed(x,:),'color',col(x,:)),1:16);
MeanShuff=squeeze(mean(Shuffle,1));
figure;hold on;arrayfun(@(x) plot(MeanShuff(x,:),'color',col(x,:)),1:16);

figure;hold on;
for i=1:500
    arrayfun(@(x) plot(squeeze(Shuffle(i,x,:)),'color',col(x,:)),1:16);
    pause
    cla
end

%% remove the trend with arima
model = arima('ARLags', 10, 'D', 0, 'Constant', 0); % AR(1) model
[ArObserved,ArObserved_p]= (arrayfun(@(x) ARIMAtrendcorr(Observed(:,x),TrlRng,model,'Spearman'),40:NTim,'UniformOutput',0));
clear ArShuffle ArShuffle_p
for i=1:NrepShuff
    %   fprintf('\n Cal ArShuff %i');
    tic
    [ArShuffle(:,:,i),ArShuffle_p(:,:,i)]= (arrayfun(@(x) ARIMAtrendcorr(squeeze(Shuffle(i,:,x))',TrlRng,model,'Spearman'),40:NTim,'UniformOutput',0));
    toc
end
ManData.ClusterMassCorrection_permutationTwoTail(permute(cell2mat(ArShuffle),[3 1 2]),permute(cell2mat(ArObserved),[1 3 2]),0.2,1,'ShowClustCorrectionPlot',1);

figure;hold on;arrayfun(@(y) plot(arrayfun(@(x) corr(TrlRng,squeeze(ArShuffle(:,x,y)),'type','Spearman'),1:NTim)),1:10)
plot(arrayfun(@(x) corr(TrlRng,ArObserved(:,x),'type','Spearman'),1:NTim),'LineWidth',20)
%%
ArObserved=cell2mat(arrayfun(@(x) GetArimaResiduals(Observed(:,x),2),85,'UniformOutput',0));
for i=1:500%NrepShuff
    fprintf('\n Cal ArShuff %i');
    ArShuffle2(:,:,i)=cell2mat(arrayfun(@(x) GetArimaResiduals(squeeze(Shuffle(i,:,x)),2)',85,'UniformOutput',0));
end
CorrShuffAr=arrayfun(@(x) corr(TrlRng,squeeze(ArShuffle2(1,:,x))','type','Spearman'),1:500);
CorrObsvAr=corr(TrlRng,ArObserved,'type','Spearman');

sum(CorrShuffAr>=CorrObsvAr)/501

pval=ManData.CalpValShuffle(CorrShuff_PrdAvg,CorrObsv_PrdAvg,'pvaltail','right');

AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod='';
if strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'Kendall')
    % calculate correlation metric
    CorrObsv=arrayfun(@(x) corr(TrlRng,Observed(1:NStg,x),'type','Kendall'),1:NTim);
    CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(Shuffle(rep,1:NStg,x))','type','Kendall'),1:NTim)',1:NrepShuff,'UniformOutput',0))';
elseif strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'GPTKendall')
    CorrObsv=arrayfun(@(x) mannKendallTestAC(Observed(1:NStg,x),0.05),1:NTim);
    CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) mannKendallTestAC(squeeze(Shuffle(rep,1:NStg,x)),0.05),1:NTim)',1:NrepShuff,'UniformOutput',0))';
 
elseif strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'TheilSen')
    CorrObsv=arrayfun(@(x) TheilSen([TrlRng,Observed(1:NStg,x)]),1:NTim);
    CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) TheilSen([TrlRng,squeeze(Shuffle(rep,1:NStg,x))']),1:NTim)',1:NrepShuff,'UniformOutput',0))';
elseif strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'Spearman')
    % calculate correlation metric
    CorrObsv=arrayfun(@(x) corr(TrlRng,Observed(1:NStg,x),'type','Spearman'),1:NTim);
    CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(Shuffle(rep,1:NStg,x))','type','Spearman'),1:NTim)',1:NrepShuff,'UniformOutput',0))';

elseif strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'Regression')
    CorrObsv=arrayfun(@(x) LinearRegression(TrlRng,Observed(1:NStg,x)),1:NTim);
    CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) LinearRegression(TrlRng,squeeze(Shuffle(rep,1:NStg,x))),1:NTim)',1:NrepShuff,'UniformOutput',0))';
elseif strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'MMK')
    significance_value_tau = AnalysisOpts.Modified_MannKendall_significance_value_tau;
    significance_value_ac = AnalysisOpts.Modified_MannKendall_significance_value_ac;

    CorrObsv=arrayfun(@(x) Modified_MannKendall_test(TrlRng, Observed(1:NStg,x)', significance_value_tau, significance_value_ac),1:NTim);
    CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) Modified_MannKendall_test(TrlRng, squeeze(Shuffle(rep,1:NStg,x))', significance_value_tau, significance_value_ac),1:NTim)',1:NrepShuff,'UniformOutput',0))';
elseif strcmp(AnalysisOpts.Classifier_TrlShuff_TrendCorrMethod,'Sigmoid')
    [SlopeObsv,~,~,DifMinMaxObsv,NormDifMinMaxObsv,NormDifMinObsv]=arrayfun(@(x) fitAndAdjustSigmoidLsq(TrlRng,Observed(1:NStg,x),0,num2str(Time(x))),1:NTim,'UniformOutput',0);
    for rep=1:NrepShuff
        [SlopeShuff(rep,:),~,~,DifMinMaxShuff(rep,:),NormDifMinMaxShuff(rep,:),NormDifMinShuff(rep,:)]=arrayfun(@(x) fitAndAdjustSigmoidLsq(TrlRng,squeeze(Shuffle(rep,1:NStg,x))',0,num2str(Time(x))),1:NTim,'UniformOutput',0);
    end
    CorrObsv=cell2mat(NormDifMinMaxObsv);
    CorrShuff=cell2mat(NormDifMinMaxShuff);

    % CorrObsv(isoutlier(CorrObsv,1))=nan;
    % CorrShuff(isoutlier(CorrShuff,1))=nan;
end

%% look at stat tests
TimeInd=Time>=-0.15 & Time<=0.6;
for i=20
    ClusterCorrectionTJB(CorrObsv(1,TimeInd)',CorrShuff(:,TimeInd)','Time',Time(TimeInd),'ThresholdPercentage',i,'DoubleSided',1);
end
ManData.ClusterMassCorrection_permutationTwoTail(Prmt(CorrShuff(:,TimeInd)),Prmt(CorrObsv(1,TimeInd)),0.1,1,'ShowClustCorrectionPlot',1);

%% look at stat test with only avegrage of 100-200ms
TimPrdAvg=Time>=0.1 & Time<=0.3;
Observed_PrdAvg=mean(Observed(:,TimPrdAvg),2);
Shuffle_PrdAvg=mean(Shuffle(:,:,TimPrdAvg),3);

CorrObsv_PrdAvg=arrayfun(@(x) corr(TrlRng,Observed_PrdAvg(1:NStg),'type','Spearman'),1);
CorrShuff_PrdAvg=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,Shuffle_PrdAvg(rep,1:NStg)','type','Spearman'),1)',1:NrepShuff,'UniformOutput',0))';

pval=ManData.CalpValShuffle(CorrShuff_PrdAvg,CorrObsv_PrdAvg,'pvaltail','right');

%% calculate stats by calculating the difference between the first and the last
LogObserved=  (Observed);LogShuffle=  (Shuffle);

ObservedDiff=Prmt((LogObserved(NStg,:)-LogObserved(1,:))./(LogObserved(1,:)+LogObserved(NStg,:)));
ShuffleDiff=Prmt(squeeze((LogShuffle(:,NStg,:)-LogShuffle(:,1,:))./(LogShuffle(:,1,:)+LogShuffle(:,NStg,:))));

ObservedDiff=Prmt((LogObserved(NStg,:)-LogObserved(1,:)));
ShuffleDiff=Prmt(squeeze((LogShuffle(:,NStg,:)-LogShuffle(:,1,:))));

%[ObservedDiff,ShuffleDiff]=ManipulateMeticVal(ObservedDiff,ShuffleDiff,'BaseLine');

TimeInd=Time>=-0.15 & Time<=0.6;
ObservedDiffTJB=squeeze(ObservedDiff);ShuffleDiffTJB=squeeze(ShuffleDiff);
for i=20
    ClusterCorrectionTJB(ObservedDiffTJB(TimeInd),ShuffleDiffTJB(:,TimeInd)','Time',Time(TimeInd),'ThresholdPercentage',i,'DoubleSided',0);
end
ManData.ClusterMassCorrection_permutationTwoTail(ShuffleDiff(:,:,TimeInd),ObservedDiff(:,:,TimeInd),0.2,1,'ShowClustCorrectionPlot',1);
%% look at stat test of difference with only avegrage of 100-200ms
ObservedDiff_PrdAvg=mean(ObservedDiff(:,TimPrdAvg),2);
ShuffleDiff_PrdAvg=mean(ShuffleDiff(:,:,TimPrdAvg),3);

end
function [ScoresMetric,EncodingDist,CompressionEncoding]=CalculateIndex(MaxTrls,obj,ManData,IndexName,Scores_1ndD,Scores_2ndD,TargetFactorDim1,TargetFactorDim2,rep_sh)
global AnalysisOpts

% 4 objects correspond to [Shape Color] category [1 1],[2 2],
% [1 2],[2 1] %[red bunny, green tee, green bunny, red tee]
% calculate compression index across time for each learning stage
if strcmp(TargetFactorDim1,'ColorCat') & strcmp(TargetFactorDim2,'ShapeCat')
    % then the first dimension is color and second dimension is shape
    ScoresColor=Scores_1ndD;ScoresShape=Scores_2ndD;
elseif strcmp(TargetFactorDim2,'ColorCat') & strcmp(TargetFactorDim1,'ShapeCat')
    % then the first dimension is shape and second dimension is color
    ScoresColor=Scores_2ndD;ScoresShape=Scores_1ndD;
else
    warning('This Analysis is only possible when dim1 is color and dim2 is shape or vice versa')
    h=[];return
end
nTrialRange=size(ScoresColor,1);

if size(ScoresColor(1).TrialRange{1},3)>1
    for TrlRng=1:nTrialRange
        ScoresColor_Org(TrlRng).TrialRange=arrayfun(@(x) cell2mat(ScoresColor(TrlRng).TrialRange(:,x)),1:4,'UniformOutput',0);
        ScoresColor(TrlRng).TrialRange=cellfun(@(x) squeeze(x(:,1,:)),ScoresColor_Org(TrlRng).TrialRange,'UniformOutput',0);
        ScoresShape_Org(TrlRng).TrialRange=arrayfun(@(x) cell2mat(ScoresShape(TrlRng).TrialRange(:,x)),1:4,'UniformOutput',0);
        ScoresShape(TrlRng).TrialRange=cellfun(@(x) squeeze(x(:,1,:)),ScoresShape_Org(TrlRng).TrialRange,'UniformOutput',0); 
    end
end

if strcmp(MaxTrls,'ALL')
    ColorAxisTime=arrayfun(@(TrlRng) arrayfun(@(obj) ScoresColor(TrlRng).TrialRange{obj},...
        1:4,'UniformOutput',0),1:nTrialRange,'UniformOutput',0);
    ShapeAxisTime=arrayfun(@(TrlRng)  arrayfun(@(obj) ScoresShape(TrlRng).TrialRange{obj},...
        1:4,'UniformOutput',0),1:nTrialRange,'UniformOutput',0);
else
    ColorAxisTime=arrayfun(@(TrlRng) arrayfun(@(obj) ScoresColor(TrlRng).TrialRange{obj}(1:MaxTrls,:),...
        1:4,'UniformOutput',0),1:nTrialRange,'UniformOutput',0);
    ShapeAxisTime=arrayfun(@(TrlRng)  arrayfun(@(obj) ScoresShape(TrlRng).TrialRange{obj}(1:MaxTrls,:),...
        1:4,'UniformOutput',0),1:nTrialRange,'UniformOutput',0);
end
[CompressionEncoding,EncodingDist]=obj.CalQuadrantCompressionEncodingAxis(ShapeAxisTime,ColorAxisTime);
% choose which factor we are looking at
if contains(obj.ThisExtraTargetFactor,'EncodingDist')% then we are looking at encoding distances
    EncodingField=erase(obj.ThisExtraTargetFactor,'EncodingDist_');
    if strcmp(AnalysisOpts.TrendTestMethod,'Permutation') % if we are using permutation to find trends
        ScoresMetricResamp=EncodingDist.([EncodingField 'ReSampAvg']);
        ScoresMetricResamp=ManData.SmoothData(ScoresMetricResamp,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
        ScoresMetric=squeeze(mean(ScoresMetricResamp,1));
    else
        ScoresMetric=EncodingDist.([EncodingField 'Avg']);
        ScoresMetric=ManData.SmoothData(ScoresMetric,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
    end
elseif contains(IndexName,'EncodingAxis') % then we are looking at encoidng axis for four objects
    PlotType='QuadrantObjsTime';
else% then we are looking at compression index
    if strcmp(IndexName,'QuadrantsCong')
        ScoresMetric=CompressionEncoding.TrlAvg.Cong;
        ScoresMetric=log(ScoresMetric);
        ScoresMetric=ManData.SmoothData(ScoresMetric,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
    elseif strcmp(IndexName,'QuadrantsInCong')
        ScoresMetric=CompressionEncoding.TrlAvg.InCong;
        ScoresMetric=log(ScoresMetric);
        ScoresMetric=ManData.SmoothData(ScoresMetric,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
    else
        %         if strcmp(AnalysisOpts.TrendTestMethod,'Permutation') % if we are using permutation to find trends
        %             ScoresMetric=CompressionEncoding.TrlAvg.All;
        %             ScoresMetric=log(ScoresMetric);
        %             ScoresMetric=ManData.SmoothData(ScoresMetric,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
        %             ScoresMetricResamp{1}=ManData.SmoothData(EncodingDist.ColorDistReSampAvg,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
        %             ScoresMetricResamp{2}=ManData.SmoothData(EncodingDist.ShapeDistReSampAvg,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
        %         else
        ScoresMetric=CompressionEncoding.TrlAvg.All;
        %     ScoresMetric=log(ScoresMetric);
        %     ScoresMetric=ManData.SmoothData(ScoresMetric,obj.WidthSmoothingDim2,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',1);
        %         end
    end
end

end