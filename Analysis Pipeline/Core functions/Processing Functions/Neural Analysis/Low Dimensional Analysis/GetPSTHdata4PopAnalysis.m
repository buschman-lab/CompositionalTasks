% gets PSTH data average and their z-score for population analysis 
function [CatPSTH,zCatPSTH,ConcatzPSTHall,CatInds,CatIndsTot,NewTime]=GetPSTHdata4PopAnalysis(Area,TimePeriod,Reward,OnsetRef,TargetFactor,WidthSmoothing,equlCongIncongWeight)
% Area is the name of the Area we want the data from
% Reward=1 or 0 for correct and incorrect trials 2 for all trials
% OnsetRef is the time reference where we are looking at the data SAMLE_ON or SACCADE_START
% TargetFactor what are the target factors we are using. include xRuleX to have that factor during specific rule 
% Outputs
% CatPSTH= average PSTH for each factor
%zCatPSTH= zscored(by mean and std of all timepoints and trials)
%ConcatzPSTHall =concatinated PSTH for all trials 
%CatInds identity in the form of [factor, stimind] for each time point
%CatIndsTot total combinations of identites for trials
global AnalysisOpts

NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
%equlCongIncongWeight=1; % are we giving equal weight to congruent and incongruent
 
%% load PFC data
if ~isfield(AnalysisOpts,[Area 'data_' OnsetRef])
    fprintf(2,'\nLoading PSTH data ....')
    PATHData=[AnalysisOpts.DataSavePath  'Core Data' filesep 'SpikingData' filesep];
    AnalysisOpts.([Area 'data_' OnsetRef])=load([PATHData 'ALL_' AnalysisOpts.Ch_2look_AreaName{1} '_AllTrials_' OnsetRef '.mat']);
    AnalysisOpts.NNeu=length(AnalysisOpts.([Area 'data_' OnsetRef]).PSTHRaster_100ms);
    [AnalysisOpts.FactorizedData]=NeuAna.GetFactorizedData(AnalysisOpts.([Area 'data_' OnsetRef]).PSTHRaster_100ms, ...
        AnalysisOpts.([Area 'data_' OnsetRef]).TrialSpikeData{2},0,'MeanSubtractByRule',0, ...
        'SpkCountPeriod',[]); % get factorized data
end

%% create population response for conditions 
[~,Time]=ManData.GenerateTimeAxis(NeuAna.PSTHTimRef);

TrlRange =ManData.GenMovingInds(1,800,800,0);
if isempty(TargetFactor)
    TargetFactor={'ColorMLxRule1','ColorMLxRule2','ColorMLxRule3',...
        'ShapeMLxRule1','ShapeMLxRule2','ShapeMLxRule3',...
        'ResponseLocxRule1','ResponseLocxRule2','ResponseLocxRule3'};
end
NTargFactors=length(TargetFactor);
%GetTime=@(x,y) floor(x*1000)>=floor(y(1)*1000) &  floor(x*1000)<=floor(y(end)*1000);
TimInd=ManData.GetExatTimeAxis(Time,TimePeriod);%Time>=TimePeriod(1) & Time<=TimePeriod(end); % take wanted time points
NTim=sum(TimInd);
NewTime=Time(TimInd);
CatPSTHtot=[];
for i=1:length(TrlRange)
    FactorizedData=NeuAna.LimitFactorizedDataTrialsBasedonFactor(AnalysisOpts.FactorizedData,'TrialNum',TrlRange{i});
    if Reward~=2 % limit trials based on reward
        FactorizedData=NeuAna.LimitFactorizedDataTrialsBasedonFactor(FactorizedData,'Reward',Reward);
    end

    for f=1:NTargFactors % loop on the factors we want
        factind=f;%TargetFactorInds(f);
        for Neu=1:length(FactorizedData) % loop on the neurons
            % find the data for this factor
            [SortFactorizedData(Neu).data,SortFactorizedData(Neu).FactorLevels,SortFactorizedData(Neu).CongruencyData]=...
                NeuAna.SortFactorDatabyLevel(FactorizedData(Neu),TargetFactor{factind});
            if equlCongIncongWeight & (contains(TargetFactor{factind},'shape','IgnoreCase',1) | contains(TargetFactor{factind},'color','IgnoreCase',1))
                % take average of the trails based on congruent and incongruent
                SortFactorizedData(Neu).data=arrayfun(@(x) ...
                    [mean(SortFactorizedData(Neu).data{x}(SortFactorizedData(Neu).CongruencyData{x}==0,:),1);...
                    mean(SortFactorizedData(Neu).data{x}(SortFactorizedData(Neu).CongruencyData{x}==1,:),1)], ...
                    1:length(SortFactorizedData(Neu).data),'UniformOutput',0);

                SortFactorizedData(Neu).(TargetFactor{factind}).CountCongTrl=cell2mat(arrayfun(@(x) [sum(SortFactorizedData(Neu).CongruencyData{x}==0);...
                    sum(SortFactorizedData(Neu).CongruencyData{x}==1)], ...
                    1:length(SortFactorizedData(Neu).data),'UniformOutput',0));
            end
        end
        nFactLvl(f)=length(SortFactorizedData(Neu).FactorLevels);
                           
        % now concatinate the PSTH for each neuron based on the value of the factor for each level
        [CatPSTH.(TargetFactor{f})]=NeuAna.CatRasterPSTH(SortFactorizedData,[],Time,[],0);
        % smooth the data for each factor now 
        for fl=1:nFactLvl(f)
            CatPSTH.(TargetFactor{f})(fl).SpkCountBin=ManData.SmoothData(CatPSTH.(TargetFactor{f})(fl).SpkCountBin(:,TimInd), ...
                WidthSmoothing,'SmoothingMethod','movmean');
        end
        % concatinate across trials so we can do the z-score
        CatPSTHtot=[CatPSTHtot cell2mat(arrayfun(@(x) CatPSTH.(TargetFactor{f})(x).SpkCountBin,1:nFactLvl(f),'UniformOutput',0))];    
    end

    % after all of this now do the z-score
    NeuMean=repmat(mean(CatPSTHtot,2),[1 NTim]);
    NeuStd=repmat(std(CatPSTHtot,0,2),[1 NTim]);

    % now create z-score by subtracting the mean and dividing by std for each neuron and factor 
    for f=1:NTargFactors
        for fl=1:nFactLvl(f)
            zCatPSTH.(TargetFactor{f})(fl).SpkCountBin=(CatPSTH.(TargetFactor{f})(fl).SpkCountBin-NeuMean)./NeuStd;
        end
    end
end

% concatinate all factors together 
ConcatzPSTHall=[];CatInds=[];CatIndsTot=[];
 
for f=1:NTargFactors
    for fl=1:nFactLvl(f)
        ConcatzPSTHall=[ConcatzPSTHall;zCatPSTH.(TargetFactor{f})(fl).SpkCountBin'] ;
        CatInds=[CatInds [fl*ones(1,NTim);f*ones(1,NTim)]];
        CatIndsTot=[CatIndsTot [fl;f]];
    end
end
 
end