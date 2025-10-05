global AnalysisOpts

% response to reviewer 2 for linear and non-linear embedding 
NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
fp=fig_params;
SM=@(x,dim,k) smoothdata(x,dim,"movmean",k);

%% load PFC data
if ~isfield(AnalysisOpts,'PFCdata')
    fprintf(2,'\nLoading PSTH data ....')
    PATHData=[AnalysisOpts.DataSavePath  'Core Data' filesep 'SpikingData' filesep];
    AnalysisOpts.PFCdata=load([PATHData 'ALL_' AnalysisOpts.Ch_2look_AreaName{1} '_AllTrials_SAMPLE_ON.mat']);
    [AnalysisOpts.FactorizedData]=obj.GetFactorizedData(AnalysisOpts.PFCdata.PSTHRaster_100ms,AnalysisOpts.PFCdata.TrialSpikeData{2},0,'MeanSubtractByRule',0,'SpkCountPeriod',ClassifierOpts.SpkCountPeriod); % get factorized data
end

%% create population response for conditions 
[~,TimeRef]=ManData.GenerateTimeAxis(NeuAna.PSTHTimRef);
Time=TimeRef-AnalysisOpts.PopulationAna.PSTHbin*0.001/2;

TrlRange =ManData.GenMovingInds(1,800,800,0);
TargetFactor={'ColorMLxRule1','ColorMLxRule2','ColorMLxRule3'};
%TargetFactor={'ResponseLocxRule1','ResponseLocxRule2','ResponseLocxRule3'};

NTargFactors=length(TargetFactor);

TimInd=Time>=-0.2 & Time<=0.6;
NTim=sum(TimInd);
clear  CatPSTH
CatPSTHtot=[];
for i=1:length(TrlRange)
    FactorizedData=NeuAna.LimitFactorizedDataTrialsBasedonFactor(AnalysisOpts.FactorizedData,'TrialNum',TrlRange{i});
    FactorizedData=NeuAna.LimitFactorizedDataTrialsBasedonFactor(FactorizedData,'Reward',1);

    for f=1:NTargFactors % loop on the factors we want
        factind=f;%TargetFactorInds(f);
        for Neu=1:length(FactorizedData) % loop on the neurons
            % find the data for this factor
            [SortFactorizedData(Neu).data,SortFactorizedData(Neu).FactorLevels]=...
                NeuAna.SortFactorDatabyLevel(FactorizedData(Neu),TargetFactor{factind});
        end
        nFactLvl(f)=length(SortFactorizedData(Neu).FactorLevels);
        % now concatinate the PSTH for each neuron based on the value of the factor for each level
        [CatPSTH.(TargetFactor{f})]=NeuAna.CatRasterPSTH(SortFactorizedData,[],Time,[],0);
        % smooth the data for each factor now 
        for fl=1:nFactLvl(f)
            CatPSTH.(TargetFactor{f})(fl).SpkCountBin=ManData.SmoothData(CatPSTH.(TargetFactor{f})(fl).SpkCountBin(:,TimInd), ...
                1,'SmoothingMethod','movmean');
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
ColorCatColors=AnalysisOpts.MorphlevelsColRGB;

for f=1:NTargFactors
    for fl=1:nFactLvl(f)
        ConcatzPSTHall=[ConcatzPSTHall;zCatPSTH.(TargetFactor{f})(fl).SpkCountBin'] ;
        CatInds=[CatInds [fl*ones(1,NTim);f*ones(1,NTim)]];
        CatIndsTot=[CatIndsTot [fl;f]];
    end
end
data=ConcatzPSTHall;
wSM=1;


%% IF WE ARE LOOKING AT BETAS
if 0
    GLMWeightsConcatTot=ManData.ReshapeCell2Mat(GLMWeightsConcat,1)';
    ColorCatColors=AnalysisOpts.MorphlevelsColRGBInc50([1 5 2 3 4],:);
    data= GLMWeightsConcatTot;
    NTim=length(TimInd);
    Np=size(GLMWeightsConcatTot,1);
    NTargFactors=3;
    nFactLvl=ones(1,NTargFactors)*Np/(NTim*NTargFactors);
    CatInds=[];CatIndsTot=[];
    for f=1:NTargFactors
        for fl=1:nFactLvl(f)
            %  ConcatzPSTHall=[ConcatzPSTHall;zCatPSTH.(TargetFactor{f})(fl).SpkCountBin'] ;
            CatInds=[CatInds [fl*ones(1,NTim);f*ones(1,NTim)]];
            CatIndsTot=[CatIndsTot [fl;f]];
        end
    end
    wSM=15;
end

%% project into low dimensional space
% take  pca
[coeff,score,latent,tsquared,explained,mu]=pca(data);
[coefftsne,losstsne]=tsne(data,'Algorithm','exact','NumDimensions',3);
[coeffuma]=run_umap(data,'n_components',3,'verbose','none');
% do multidimensional scaling
CDMD=mdscale(squareform(pdist(data,'euclidean')),3);

% smooth data of we need to 
score=SM(score,1,wSM);coefftsne=SM(coefftsne,1,wSM);coeffuma=SM(coeffuma,1,wSM);CDMD=SM(CDMD,1,wSM);

%% Plot results
% plot data for each target factor 
LineStyle={'-','--',':'};
PlotType=0;% sperate figure 0 subplot 1

% plot using pca
figure
if PlotType
    subplot(221)
end
hold on
arrayfun(@(x) ManData.PlotDatainTime(Time(TimInd), ...
    score(CatInds(1,:)==CatIndsTot(1,x) & CatInds(2,:)==CatIndsTot(2,x),1:3),[1 2 3],...
    ColorCatColors(CatIndsTot(1,x),:),LineStyle{CatIndsTot(2,x)},[0:0.15:0.6]),1:sum(nFactLvl));

fp.ForcePutLegends({'R1','R2','R3'},ColorPallet(1:3,:),LineStyle);
title(['PCA, explained var PC1-PC3=' num2str(sum(explained(1:3))) '%'])

% plot using tsne
if PlotType
    subplot(222)
else
    figure
end
hold on
arrayfun(@(x) ManData.PlotDatainTime(Time(TimInd), ...
    coefftsne(CatInds(1,:)==CatIndsTot(1,x) & CatInds(2,:)==CatIndsTot(2,x),1:3),[1 2 3],...
    ColorCatColors(CatIndsTot(1,x),:),LineStyle{CatIndsTot(2,x)},[0:0.15:0.6]),1:sum(nFactLvl));

fp.ForcePutLegends({'R1','R2','R3'},ColorPallet(1:3,:),LineStyle);
xlabel('tsne1');ylabel('tsne2');zlabel('tsne3');
title(['TSNE, loss=' num2str(losstsne) ]);

% plot using umap
if PlotType
    subplot(223)
else
    figure
end
hold on
arrayfun(@(x) ManData.PlotDatainTime(Time(TimInd), ...
    coeffuma(CatInds(1,:)==CatIndsTot(1,x) & CatInds(2,:)==CatIndsTot(2,x),1:3),[1 2 3],...
    ColorCatColors(CatIndsTot(1,x),:),LineStyle{CatIndsTot(2,x)},[0:0.15:0.6]),1:sum(nFactLvl));

fp.ForcePutLegends({'R1','R2','R3'},ColorPallet(1:3,:),LineStyle);
xlabel('umap1');ylabel('umap2');zlabel('umap3');
title('UMAP')
  view(0,90)
if PlotType
    subplot(224)
else
    figure;
end
hold on
arrayfun(@(x) ManData.PlotDatainTime(Time(TimInd), ...
    CDMD(CatInds(1,:)==CatIndsTot(1,x) & CatInds(2,:)==CatIndsTot(2,x),1:3),[1 2 3],...
    ColorCatColors(CatIndsTot(1,x),:),LineStyle{CatIndsTot(2,x)},[0:0.15:0.6]),1:sum(nFactLvl));

fp.ForcePutLegends({'R1','R2','R3'},ColorPallet(1:3,:),LineStyle);
xlabel('MDS1');ylabel('MDS2');xlabel('MDS3');
title('MDS')

