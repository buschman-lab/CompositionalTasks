function [GLMWeightsAll,Mdl,Time,TimeIndGLM,WantedFieldName]=LoadGLMdata4ModelFitSamp(TimeRef,Area,MdlName,Lambda)

global AnalysisOpts
ManData=ManipulateData;
NeuAna=NeuralAnalysisFuncsTemp;

%% load varibales
% TimeRef='SAMPLE_ON';%'SACCADE_START';
% Area='PFC';
GLMpath=[AnalysisOpts.DataSavePath 'GLM' filesep 'ALL' filesep];

%AnalsysOptsOrg=load([GLMpath 'GLM_ALL_ALL_2SingCellAna_' Area '_' TimeRef '_100_Summery.mat'], 'AnalysisOpts');
NNeu=1:AnalysisOpts.NNeuAnalysis;
 
GLMModelTxt=['GLMfit_' MdlName];
GLMModelObsTxt=['GLMfitObsrv_' MdlName];
GLMModelShTxt=['GLMfitSh_' MdlName];

for i=NNeu  
    tic
    FileName=sprintf(['GLM_ALL_ALL_%iSingCellAna_' Area '_' TimeRef '_100_GLMfitSamp.mat'],i);
    Mdl(i)=load([GLMpath FileName],GLMModelTxt,GLMModelObsTxt,GLMModelShTxt);    % load resampling, observed and shuffle data  
    fprintf('\nLoading file %i in %0.2f secs...',i,toc);
end

% load time data
Time=Mdl(1).(['GLMfit_' MdlName]).Time;
% generate time axis
if strcmp(TimeRef,'SAMPLE_ON')
    TimeIndGLM=ManData.GetExatTimeAxis(Time,[-0.2 0.6]);%Time>=TimePeriod(1) & Time<=TimePeriod(end); % take wanted time points
elseif strcmp(TimeRef,'SACCADE_START')
    TimeIndGLM=ManData.GetExatTimeAxis(Time,[-0.3 0.2]);%Time>=TimePeriod(1) & Time<=TimePeriod(end); % take wanted time points
end

FieldNames=Mdl(1).(GLMModelTxt)(1).(MdlName).GLMfactorsnames;
WantedFieldName=FieldNames(~contains(FieldNames,'Time')  & ...
    ~contains(FieldNames,'Bias'));

%% organize the GLM for factors
% for Factor=WantedFieldName
%     IndFactor=find(strcmp(FieldNames,Factor{1}));
%     % get Beta values
%     for i=NNeu
%         GLMWeights=squeeze(Mdl(i).(GLMModelTxt).(MdlName).GLM_weights_full(:,Lambda,:));
%         GLMWeightsAll.(Factor{1})(:,i)=GLMWeights(IndFactor,TimeIndGLM)';
%     end
% end
GLMWeightsAll=[];
end