% processes trial shuffle statistical test based on the distribution 
function ProcessTrlShuffStatTests(ClassifierResults_Observed,ClassifierResults_Shuff,ClassifierOpts,Lag,AvgDim2,Dim,FileNameSpecs,SaveFileName)
%AvgDim2 % are we averaging in second dimension?

global AnalysisOpts

NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
fp=fig_params;
SaveArimaPlots=0; % don't save them for now because we plot them later

%% load data
% the data should be already loaded and passed along to this 
%% preprocess the data
DimNum=[0 1 3];

ClassifierResults_Shuff=arrayfun(@(x) NeuAna.PreprocessClassifierResults(ClassifierResults_Shuff,ClassifierOpts,x,1),DimNum(Dim),'UniformOutput',0);
ClassifierResults_Observed=arrayfun(@(x) NeuAna.PreprocessClassifierResults(ClassifierResults_Observed,ClassifierOpts,x,0),DimNum(Dim),'UniformOutput',0);

%% extract important information
[~,Accuracy_Obsv] =arrayfun(@(x) NeuAna.GetPerfMetricData4Classifier(ClassifierResults_Observed{x},'Accuracy',1,0),1,'UniformOutput',0);
[~,Accuracy_Shuff]=arrayfun(@(x) NeuAna.GetPerfMetricData4Classifier(ClassifierResults_Shuff{x},   'Accuracy',1,1),1,'UniformOutput',0);

% structue of MetricVal and MEtricValShuff is Rep*TrlRng*Time
Accuracy_Obsv=cellfun(@(x) ManData.ReshapeCell2Mat(x,3),Accuracy_Obsv,'UniformOutput',0);
Accuracy_Obsv=cellfun(@(x) permute(x,[3 2 1]),Accuracy_Obsv,'UniformOutput',0);
Accuracy_Shuff=cellfun(@(x) ManData.ReshapeCell2Mat(x,3),Accuracy_Shuff,'UniformOutput',0);
Accuracy_Shuff=cellfun(@(x) permute(x,[3 1 2]),Accuracy_Shuff,'UniformOutput',0);

%% run stat test for CPI as well
if strcmp(ClassifierOpts.TestName,'Learning3D_Shape_Color_Rule_Xgen_AltRule_RB') & ClassifierOpts.IncludeAllClassifierInfo & Dim==1
    UseSavedData=0;PlotFlag=0;
    [CPI_Obsv,CPI_Shuff,ShapeAvgDist_Obsv,ShapeAvgDist_Shuff,ColorAvgDist_Obsv,ColorAvgDist_Shuff,CompressionEncoding_Obsv,CompressionEncoding_Shuff]=...
         arrayfun(@(x) CalculateCompressionIndex4ShuffleData(ClassifierOpts,PlotFlag,UseSavedData),1,'UniformOutput',false);

else
    CPI_Obsv=[];CPI_Shuff=[];ShapeAvgDist_Obsv=[];ShapeAvgDist_Shuff=[];ColorAvgDist_Obsv=[];ColorAvgDist_Shuff=[];FitArima2CPI=0;CompressionEncoding_Obsv=[];CompressionEncoding_Shuff=[];
end

%% Save this data so we can use it in the future 
fprintf('\nSaving data of classifier results to %s',SaveFileName)
save(SaveFileName,'Accuracy_Obsv','Accuracy_Shuff','CPI_Obsv','CPI_Shuff','ShapeAvgDist_Obsv','ShapeAvgDist_Shuff', ...
    'ColorAvgDist_Obsv','ColorAvgDist_Shuff','CompressionEncoding_Obsv','CompressionEncoding_Shuff');

%% Fit Arima Model + plot and save results

%Lag=2; % ARIMA lag
% SaveFileName=arrayfun(@(x) FitArima2CompressionIndex(Lag,'Accuracy',AvgDim2,FileNameSpecs(Dim(x)),Accuracy_Obsv{x},Accuracy_Shuff{x},Dim),1,'Uniformoutput',0);
% 
% if FitArima2CPI
%     FitArima2CompressionIndex(Lag,'CPI',AvgDim2,ClassifierOpts.Name,CPI_obsv,CPI_Shuff);
%     FitArima2CompressionIndex(Lag,'ShapeDist',AvgDim2,ClassifierOpts.Name,EncodingDistShape_Obsv,EncodingDistShape_Shuff);
%     FitArima2CompressionIndex(Lag,'ColoDist',AvgDim2,ClassifierOpts.Name,EncodingDistColor_Obsv,EncodingDistColor_Shuff);
% end
% 
if SaveArimaPlots
    fp.SaveCurrentFigs(SaveFileName{1},[AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep]);
end

end