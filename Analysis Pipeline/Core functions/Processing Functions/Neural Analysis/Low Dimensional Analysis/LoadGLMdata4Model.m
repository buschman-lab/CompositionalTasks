function [GLMWeightsAll,GLMCPDAll,Mdl,Time,TimeIndGLM,GLMModelComp,WantedFieldName,BestLambdaInd]=LoadGLMdata4Model(TimeRef,Area,MdlName,BestLambda)

global AnalysisOpts
ManData=ManipulateData;
NeuAna=NeuralAnalysisFuncsTemp;

%% load varibales
% TimeRef='SAMPLE_ON';%'SACCADE_START';
% Area='PFC';
GLMpath=[AnalysisOpts.DataSavePath 'GLM' filesep 'ALL' filesep];

%AnalsysOptsOrg=load([GLMpath 'GLM_ALL_ALL_2SingCellAna_' Area '_' TimeRef '_100_Summery.mat'], 'AnalysisOpts');
NNeu=1:AnalysisOpts.NNeuAnalysis;
 
GLMModelTxt=['GLMfitSingFact_' MdlName];
for i=NNeu
    tic
    FileName=sprintf(['GLM_ALL_ALL_%iSingCellAna_' Area '_' TimeRef '_100_Summery.mat'],i);
    Mdl(i)=load([GLMpath FileName],GLMModelTxt);
    fprintf('\nLoading file %i in %0.2f secs...',i,toc);
end
% load time data
Mdl4Time=load([GLMpath 'GLM_ALL_ALL_1SingCellAna_' Area '_' TimeRef '_100_ShufSingFact.mat'],['GLMfit_' MdlName]);
Time=Mdl4Time.(['GLMfit_' MdlName]).Time;
% generate time axis
if strcmp(TimeRef,'SAMPLE_ON')
    TimeIndGLM=ManData.GetExatTimeAxis(Time,[-0.2 0.6]);%Time>=TimePeriod(1) & Time<=TimePeriod(end); % take wanted time points
elseif strcmp(TimeRef,'SACCADE_START')
    TimeIndGLM=ManData.GetExatTimeAxis(Time,[-0.3 0.2]);%Time>=TimePeriod(1) & Time<=TimePeriod(end); % take wanted time points
end

FieldNames=Mdl(1).(GLMModelTxt).(MdlName).GLMfactorsnames;
WantedFieldName=FieldNames(~contains(FieldNames,'Time') & ~contains(FieldNames,'Reward') & ...
    ~contains(FieldNames,'Bias'));
FieldNamesCPD=fieldnames(Mdl(1).(GLMModelTxt).CPD)';
WantedFieldNameCPD=FieldNamesCPD(~contains(FieldNamesCPD,'Time') & ~contains(FieldNamesCPD,'Reward') & ...
    ~contains(FieldNamesCPD,'Bias'));

%% load model comparision results to find the best lambda
ModelCompName=['Compare' MdlName];
[GLMModelComp,GLMModelCompFileExist]=NeuAna.LoadGLMfitfile(NNeu,ModelCompName,'ModelComp'); % load model comparision summery file
[IndLambdaMSE,GLMBestLambdaInd,IndLambdaR2TimePoint]=arrayfun(@(x) NeuAna.FindBestGLMlambda(GLMModelComp(x),MdlName),1:length(GLMModelComp));
if ~isempty(BestLambda)
    BestLambdaInd=BestLambda*ones(1,NNeu(end));
else
    BestLambdaInd=GLMBestLambdaInd;
end

%% organize the GLM for factors
for Factor=WantedFieldName
    IndFactor=find(contains(FieldNames,Factor{1},'IgnoreCase',true));
    % get Beta values
    for i=NNeu
        GLMWeights=squeeze(Mdl(i).(GLMModelTxt).(MdlName).GLM_weights_full(:,BestLambdaInd(i),:));
        GLMWeightsAll.(Factor{1})(:,i)=GLMWeights(IndFactor,TimeIndGLM)';
    end
end
for Factor=WantedFieldNameCPD
    IndFactor=find(contains(FieldNamesCPD,Factor{1},'IgnoreCase',true));
    % get CPD values
    GLMCPDAll.(Factor{1})=cell2mat(arrayfun(@(x) Mdl(x).(GLMModelTxt).CPD.(Factor{1})(BestLambdaInd(x),TimeIndGLM)',NNeu,'UniformOutput',0));
end
end