function StatTest=ReportTrialOrderShuffleStatTestResultsTemp(ClassifierOpts,Cond,varargin) % perform statistical test on classifier results
% reports the trial shuffle statistical test results that are already saved and is generated before
% Saved in Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\Code Test Data
% use the Syntax provided to load the data for that specific metric
global AnalysisOpts

NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;

[ClassifierOpts]=NeuAna.DefineClassifierTestOptions(ClassifierOpts.TestName,'CalShuff',1,'CalShuffTrlOrder',1);

%prepare variables
StatCompMethod=AnalysisOpts.CalShuffTrlOrderClassifier_StatCompMethod;
ClusterTH=AnalysisOpts.CalShuffTrlOrderClassifier_ClusterTH;

[~,~,~,TrlRng]=NeuAna.GetTrialRangeforThisCond(ClassifierOpts,Cond);

TrlRng=TrlRng';NStg=size(TrlRng,1);NStgTxt='';
if NeuAna.WidthSmoothingDim2==3;AvgDim2=1;else;AvgDim2=NeuAna.WidthSmoothingDim2;end
[~,TargetFactorDim]=arrayfun(@(dim) NeuAna.getClassifierDimInfo(ClassifierOpts,dim),[0 1 3],'UniformOutput',0);
FileNameSpecs=@(x) [ClassifierOpts.Name '_' TargetFactorDim{x}];
MetricNameSet=[TargetFactorDim,'CPI','ShapeAvgDist','ColorAvgDist'];
MetricDimSet =[1 2 3 1 1 1];
TimeIndSet=repmat({logical(ones(1,length(AnalysisOpts.Time)))},[1 3]);%FigParams.LimitTimeAxis4Classifier(AnalysisOpts.Time,ClassifierOpts);
% load data that was already processed
try
    for MetricName=MetricNameSet
        Dim=MetricDimSet(strcmp(MetricNameSet,MetricName));
        LoadFileName=[FileNameSpecs(Dim) '_D' num2str(Dim) 'ShData'];

        % load fits in time
        [OutVar,~,SaveVarName]=FitStatMetric2Data(2,StatCompMethod,AvgDim2,NStg,NStgTxt,LoadFileName,[],[],[MetricName{1}]);

        % load fits in average time
        [OutVar_TimAvg,~,SaveVarName_TimAvg]=FitStatMetric2Data(2,StatCompMethod,AvgDim2,NStg,NStgTxt,LoadFileName,[],[],[MetricName{1} '_TimAvg']);

        if ~isempty(OutVar)
            CorrObsv=OutVar.(SaveVarName).Observed_a;
            CorrShuff=OutVar.(SaveVarName).Shuffle_a;
            CorrObsv_TimAvg=OutVar_TimAvg.(SaveVarName_TimAvg).Observed_a;
            CorrShuff_TimAvg=OutVar_TimAvg.(SaveVarName_TimAvg).Shuffle_a;

            TimeInd=TimeIndSet{Dim};
            % TimeInd=logical(ones(1,length(AnalysisOpts.Time)));

            [~,PvalTxtTimAvg]=ManData.PlotPvals(CorrObsv_TimAvg,CorrShuff_TimAvg,0);

            [start_ind, end_ind, cluster_p]=ClusterCorrectionTJB(CorrObsv(1,TimeInd)',CorrShuff(:,TimeInd)',[],'Time',AnalysisOpts.Time(TimeInd),...
                'ThresholdPercentage',ClusterTH,'DoubleSided',1,'PlotFigure',1,'Title',[SaveVarName;{PvalTxtTimAvg}]);
            % convert this so we can read it later.
            [X,P]=ManData.ConvertTJBclustercorrection(start_ind, end_ind, cluster_p,AnalysisOpts.Time);
            % report stat test results
            StatTest(Cond).(MetricName{1}).clusters=X; 
            StatTest(Cond).(MetricName{1}).statsummery=P;
            StatTest(Cond).(MetricName{1}).p_values=[];
            StatTest(Cond).(MetricName{1}).t_sums=[];
            StatTest(Cond).(MetricName{1}).clustIdx=[];
            StatTest(Cond).(MetricName{1}).permutation_distribution=[];
        end
    end
catch me
    ManData.ThrowErrorTxt(me)
    StatTest=[];
end
end
