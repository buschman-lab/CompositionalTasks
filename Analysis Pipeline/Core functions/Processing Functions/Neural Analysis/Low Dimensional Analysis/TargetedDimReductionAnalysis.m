% does the analysis of Targeted Dimensionality Reduction of Mante et al 2013
function TargetedDimReductionAnalysis(TimeRef,Area,BestLambdaGLM,equlCongIncongWeight,PCexPercent)
%% load varibales
global AnalysisOpts

NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
fp=fig_params;

%TimeRef='SAMPLE_ON';%'SACCADE_START';%
%Area='PFC';
%GLMpath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\GLM\ALL\';
AreaNum=find(strcmp(AnalysisOpts.AreaNames,Area));
AnalysisOpts.NNeuAnalysis=AnalysisOpts.AreaTotalNumNeurons(strcmp(AnalysisOpts.AreaNames,Area));% get the number of neurons based on the area
ColorCatColors=AnalysisOpts.MorphlevelsColRGB;
ShapeCatColors=AnalysisOpts.MorphlevelsShpRGB;
ColorPallet=AnalysisOpts.ColorPalettCamden;
LineStyle={'-','--',':','-','--',':','-','--',':'};
WsmoothingPCA=1;WidthSmoothingPSTH=15;WidthSmoothingGLM=5; % how much we are smoothing the PSTH and how much we are smoothing the final results
SM=@(x,dim,k) smoothdata(x,dim,"movmean",k);
MeanCenter=@(x) x-mean(x,1);

% define timing variabels (what is the time limit we want for compression index)
if strcmp(TimeRef,'SAMPLE_ON') 
    TimeAvgCPI=[0.1 0.4];
    TimeLook4MaxArea={[0.1 0.2],[0.1 0.2],[0.1 0.2],[0.2 0.3],[0.1 0.2]};
    TimeLook4Max=TimeLook4MaxArea{AreaNum};
    Time2ShowInPlots=[-0.1 0.4];
elseif strcmp(TimeRef,'SACCADE_START')
    TimeAvgCPI=[-0.2 0.1];TimeLook4Max=[-0.2 0.1];Time2ShowInPlots=[-0.2 0.2];
end
 
%% load GLM fit data 
AT=@(x) [x '_' TimeRef]; %concatinate Varname with time ref
GLMmdlName='SensoryMotorIntegerMdl';
[GLMWeightsAll,GLMCPDAll,Mdl,Time,TimeIndGLM,GLMModelComp,WantedFieldName,BestLambdaInd]=...
    LoadGLMdata4Model(TimeRef,Area,GLMmdlName,BestLambdaGLM);

ManData.CopyVars2AnalysisOpts(AT('GLMWeightsAll'),GLMWeightsAll,AT('GLMCPDAll'),GLMCPDAll,AT('Mdl'),Mdl,AT('Time'),Time, ...
    AT('TimeIndGLM'),TimeIndGLM,AT('GLMModelComp'),GLMModelComp,AT('BestLambdaInd'),BestLambdaInd,AT('WantedFieldName'),WantedFieldName); % copy data to 

[GLMWeightsAll,GLMCPDAll,Mdl,Time,TimeIndGLM,GLMModelComp,WantedFieldName,BestLambdaInd]=ManData.ExtractVarsFromAnalysisOpts(AT('GLMWeightsAll'),AT('GLMCPDAll'),AT('Mdl'),AT('Time'), ...
    AT('TimeIndGLM'),AT('GLMModelComp'),AT('WantedFieldName'),AT('BestLambdaInd'));

NWantedFactors=length(WantedFieldName);
TimeGLM=Time(TimeIndGLM);

%if we are changing the data based on the Lambda
% BestLambdaIndNew=4*ones(1,AnalysisOpts.NNeuAnalysis);
% FieldNames=Mdl(1).(['GLMfitSingFact_' MdlName]).(MdlName).GLMfactorsnames;
% for Factor=WantedFieldName
%     IndFactor=find(contains(FieldNames,Factor{1},'IgnoreCase',true));
%     % get Beta values
%     for i=1:AnalysisOpts.NNeuAnalysis
%         GLMWeights=squeeze(Mdl(i).(['GLMfitSingFact_' MdlName]).(MdlName).GLM_weights_full(:,BestLambdaIndNew(i),:));
%         GLMWeightsAll.(Factor{1})(:,i)=GLMWeights(IndFactor,TimeIndGLM)';
%     end
% end

% smooth GLM fit data 
GLMWeightsAll=structfun(@(x) SM(x,1,WidthSmoothingGLM),GLMWeightsAll,'UniformOutput',0);

%% Get PSTH data for factors
[CatPSTH,zCatPSTH,ConcatzPSTHall,CatInds,CatIndsTot,PSTHTime]=GetPSTHdata4PopAnalysis(Area,TimeGLM,1,TimeRef,[],WidthSmoothingPSTH,equlCongIncongWeight);
AnalysisOpts.Time=PSTHTime;

%% Do PCA on Neural data to smooth the data 
[D,score,latent,tsquared,explained]=pca(ConcatzPSTHall);

if strcmp(PCexPercent,'all')
    UsePCdenoising=0;
    nPC=AnalysisOpts.NNeuAnalysis;
else
    UsePCdenoising=1;
    nPC=find(cumsum(explained)>=PCexPercent,1,"first");
end

% Ddenoise=D(:,1:nPC);
% GLMWeightsPCA=cellfun(@(x) [Ddenoise'*GLMWeightsAll.(x)']',WantedFieldName,'UniformOutput',0);

Ddenoise=D;Ddenoise(:,nPC+1:end)=0;
GLMWeightsPCA=cellfun(@(x) [Ddenoise'*GLMWeightsAll.(x)']',WantedFieldName,'UniformOutput',0);

%% find the Max norm for betas
TimeLimMaxColShape=find(ManData.GetExatTimeAxis(TimeGLM,TimeLook4Max));
TimeLimMax=repmat({1:length(TimeGLM)},[1 NWantedFactors]);
TimeLimMax(contains(WantedFieldName,'Shape') | contains(WantedFieldName,'Color'))={TimeLimMaxColShape};

if UsePCdenoising
    GLMWeightsNorm=cellfun(@(x) ManData.CalNormVectMatrix(x,2),GLMWeightsPCA,'UniformOutput',0);
   % [~,MaxTimeNorm]=cellfun(@max,GLMWeightsNorm);
    [~,MaxTimeNorm]=arrayfun(@(x) max(GLMWeightsNorm{x}(TimeLimMax{x})),1:NWantedFactors);
    MaxTimeNorm=arrayfun(@(x) MaxTimeNorm(x)+(TimeLimMax{x}(1)-1),1:NWantedFactors);
    WeightsMaxNorm=cell2mat(arrayfun(@(x) GLMWeightsPCA{x}(MaxTimeNorm(x),:)',1:NWantedFactors,'UniformOutput',0));
else
    GLMWeightsNorm=cellfun(@(x) ManData.CalNormVectMatrix(GLMWeightsAll.(x),2),WantedFieldName,'UniformOutput',0);
    %[~,MaxTimeNorm]=cellfun(@max,GLMWeightsNorm);
    [~,MaxTimeNorm]=arrayfun(@(x) max(GLMWeightsNorm{x}(TimeLimMax{x})),1:NWantedFactors);
    MaxTimeNorm=arrayfun(@(x) MaxTimeNorm(x)+(TimeLimMax{x}(1)-1),1:NWantedFactors);
    WeightsMaxNorm=cell2mat(arrayfun(@(x) GLMWeightsAll.(WantedFieldName{x})(MaxTimeNorm(x),:)',1:NWantedFactors,'UniformOutput',0));
end

%% Do the QR decomposition of 
[Q,R]=qr(WeightsMaxNorm);
TaskAxis=Q(:,1:NWantedFactors); % first colums that have orthogonizlied task axis

%% Project average population response into 
FieldPSTH=fieldnames(CatPSTH)';
if UsePCdenoising
    for f=FieldPSTH
        ProjPSTHonTaskAxis.(f{1})=arrayfun(@(x) TaskAxis'*([zCatPSTH.(f{1})(x).SpkCountBin'*Ddenoise]'),1:size(CatPSTH.(f{1}),2),'UniformOutput',0);
    end
else
    for f=FieldPSTH
        ProjPSTHonTaskAxis.(f{1})=arrayfun(@(x) TaskAxis'*zCatPSTH.(f{1})(x).SpkCountBin,1:size(CatPSTH.(f{1}),2),'UniformOutput',0);
    end
end

%% plot the co-efficients for each axis
h1=fp.RenderFigure(1,[]);
[h1,Sp]=fp.RenderSubplots([],[],h1{1},4);
GroupsAxisY={[1 2 3],[4],[5],[6 7]}; % groups of factors to plot together
GroupsAxisYY={[1 3],[5]}; % groups of factors to plot together

hold on
arrayfun(@(s) arrayfun(@(y) fp.PlotMeanStd(TimeGLM,GLMWeightsNorm{y},[],['Time from ' TimeRef],'Value',ColorPallet(y,:), ...
    1,WantedFieldName{y},'LegendTxt',WantedFieldName{y},'Sp',Sp(s)),GroupsAxisY{s},'UniformOutput',0),1:4,'UniformOutput',0);

% add a star for the maximum value 
arrayfun(@(s) arrayfun(@(y) fp.Text(TimeGLM(MaxTimeNorm(y)),GLMWeightsNorm{y}(MaxTimeNorm(y)),'*',ColorPallet(y,:),'font_size',15,'Sp', ...
    Sp(s)),GroupsAxisY{s},'UniformOutput',0),1:4,'UniformOutput',0);

% % plot axis yy conditions
% arrayfun(@(s) arrayfun(@(y) fp.PlotMeanStd(TimeGLM,GLMWeightsNorm{y},[],['Time from ' TimeRef],'Value',ColorPallet(y,:), ...
%     1,WantedFieldName{y},'LegendTxt',WantedFieldName{y},'Sp',Sp(s),'RightYYAxis',1),GroupsAxisYY{s},'UniformOutput',0),1:2,'UniformOutput',0);
% 
% arrayfun(@(s) arrayfun(@(y) fp.Text(TimeGLM(MaxTimeNorm(y)),GLMWeightsNorm{y}(MaxTimeNorm(y)),'*',ColorPallet(y,:),'font_size',15,'Sp', ...
%     Sp(s),'RightYYAxis',1),GroupsAxisY{s},'UniformOutput',0),1:3,'UniformOutput',0);
% 

%% plot task variables in different spaces
% show color information goes to axis 1 and axis 2 in R3 and R2
ProjCondNameColor={'ColorMLxRule1','ColorMLxRule2','ColorMLxRule3'};
ProjCondNameShape={'ShapeMLxRule1','ShapeMLxRule2','ShapeMLxRule3'};
ProjCondNameSet=[{ProjCondNameColor} {ProjCondNameShape}];

AxisNamesColor={'ColorIntML','Axis1Int','Axis2Int'};
AxisNamesShape={'ShapeIntML','Axis1Int','Axis2Int'};
AxisNamesSet=[{AxisNamesColor},{AxisNamesShape}];
CatColorsSet=[{ColorCatColors},{ShapeCatColors}];
Time2PlotInd=PSTHTime>=Time2ShowInPlots(1) & PSTHTime<=Time2ShowInPlots(2);
AxisCombinations=[1 2;1 3];%nchoosek(1:3,2);

clear AxisVall
for ii=1:2
    ni=1;
    ProjCondName=ProjCondNameSet{ii};
    CondNum=cellfun(@(x) find(strcmp(FieldPSTH,x)),ProjCondName);
    AxisNames=AxisNamesSet{ii};
    CatColors=CatColorsSet{ii};
    h2(ii)=figure;
    sgtitle(sprintf('Projection of %s into encoding axes',AxisNamesSet{ii}{1}))

    for x=CondNum    
        for ax=1:size(AxisCombinations,1)
            subplot(3,2,ax+2*(find(x==CondNum)-1))
            nf=cellfun(@(x) find(strcmp(WantedFieldName,x)),AxisNames(AxisCombinations(ax,:))); % find the index of the axes

            arrayfun(@(n) ManData.PlotDatainTime(PSTHTime(Time2PlotInd),smoothdata(ProjPSTHonTaskAxis.(FieldPSTH{x}){n}(nf,Time2PlotInd)', ...
                1,'movmean',WsmoothingPCA),[1 2],...
                CatColors(n,:),LineStyle{x},-1:0.2:1,1,'PlotWithConstantCol',1),1:6);
            title(ProjCondNameShape{(x==CondNum)});
            xlabel(WantedFieldName{nf(1)});ylabel(WantedFieldName{nf(2)});
            % xlim([-15 15]);ylim([-20 20]);
            %   fp.ForcePutLegends(FieldPSTH((x)),ColorCatColors((x),:),LineStyle((x)))
            AxisVall{ii}(ni,:)=axis;            
            ni=ni+1;
        end
    end
    % now set the best limits for the subplots
    for x=CondNum    
        for ax=1:size(AxisCombinations,1)
            subplot(3,2,ax+2*(find(x==CondNum)-1));
            xlim([min(AxisVall{ii}(:,1)) max(AxisVall{ii}(:,2))]);
            ylim([min(AxisVall{ii}(:,3)) max(AxisVall{ii}(:,4))]);
        end
    end
end
arrayfun(@(x) fp.FormatAllAxesFig(x),h2) % format axis

% plot each of the projections across time for color
ProjCondName={'ColorMLxRule1','ColorMLxRule2','ColorMLxRule3'};
CondNum=cellfun(@(x) find(strcmp(FieldPSTH,x)),ProjCondName);

AxisNames={'ColorIntML','Axis1Int','Axis2Int'};

h3=figure;
sgtitle(sprintf('Projection of encoding axes in time'));

for x=CondNum    
    for ax=1:3
        subplot(3,3,ax+3*(find(x==CondNum)-1))
        nf=cellfun(@(x) find(strcmp(WantedFieldName,x)),AxisNames(ax)); % find the index of the axes
      
        arrayfun(@(y) fp.PlotMeanStd(PSTHTime(Time2PlotInd),ProjPSTHonTaskAxis.(FieldPSTH{x}){y}(nf,Time2PlotInd),[],['Time from ' TimeRef],AxisNames{ax}, ...
            ColorCatColors(y,:),1,FieldPSTH{x}),1:6,'UniformOutput',0)
    end   
end
fp.FormatAllAxesFig(h3)

% and shape
ProjCondName={'ShapeMLxRule1','ShapeMLxRule2','ShapeMLxRule3'};
CondNum=cellfun(@(x) find(strcmp(FieldPSTH,x)),ProjCondName);

AxisNames={'ShapeIntML','Axis1Int','Axis2Int'};

h4=figure;
sgtitle(sprintf('Projection of encoding axes in time'));

for x=CondNum    
    for ax=1:3
        subplot(3,3,ax+3*(find(x==CondNum)-1))
        nf=cellfun(@(x) find(strcmp(WantedFieldName,x)),AxisNames(ax)); % find the index of the axes
      
        arrayfun(@(y) fp.PlotMeanStd(PSTHTime(Time2PlotInd),ProjPSTHonTaskAxis.(FieldPSTH{x}){y}(nf,Time2PlotInd),[],['Time from ' TimeRef],AxisNames{ax}, ...
            ShapeCatColors(y,:),1,FieldPSTH{x}),1:6,'UniformOutput',0)
    end   
end
fp.FormatAllAxesFig(h4)

%% plot color Axis 1 and Axis 2 in a three dimensional plot for rule 2 and 3
ProjCondName={'ColorMLxRule1','ColorMLxRule2','ColorMLxRule3'};
AxisNames={'ColorIntML','Axis1Int','Axis2Int'};
AxisCombinations3D=[1 2 3];
h5=figure;
for x=2:3   
    sgtitle(sprintf('Projection of %s into encoding axes',ProjCondName{x}))
    for ax=1:size(AxisCombinations3D,1)
         nf=cellfun(@(x) find(strcmp(WantedFieldName,x)),AxisNames(AxisCombinations3D(ax,:))); % find the index of the axes

        arrayfun(@(n) ManData.PlotDatainTime(PSTHTime(Time2PlotInd),smoothdata(ProjPSTHonTaskAxis.(FieldPSTH{x}){n}(nf,Time2PlotInd)', ...
            1,'movmean',WsmoothingPCA),[1 2 3],...
            ColorCatColors(n,:),LineStyle{x},-1:0.2:1,1,'PlotWithConstantCol',1),1:6);                      
    end
end
xlabel(WantedFieldName{nf(1)});ylabel(WantedFieldName{nf(2)}); zlabel(WantedFieldName{nf(3)});

% add an axis in the center 
AxisV=axis;
Org=[0 0 0];ar=1;
x_dir=[1 0 0]*AxisV(2)/ar;y_dir=[0 1 0]*AxisV(4)/ar;z_dir=[0 0 1]*AxisV(6)/ar;
quiver3(Org(1),Org(2),Org(3),x_dir(1),x_dir(2),x_dir(3),'r','filled','MaxHeadSize',0.5)
quiver3(Org(1),Org(2),Org(3),y_dir(1),y_dir(2),y_dir(3),'r','filled','MaxHeadSize',0.5)
quiver3(Org(1),Org(2),Org(3),z_dir(1),z_dir(2),z_dir(3),'r','filled','MaxHeadSize',0.5)

fp.FormatAllAxesFig(h5)
fp.ForcePutLegends(FieldPSTH(CondNum(x)),ColorCatColors(CondNum(x),:),LineStyle(CondNum(x)))

%% now measure angle between subspaces for color in each rule and subspaces of response 
if 0
    [GLMWeightsAll_InteractMdl,GLMCPDAll_InteractMdl,Mdl_InteractMdl,Time_InteractMdl,TimeIndGLM_InteractMdl,GLMModelComp_InteractMdl]=...
        LoadGLMdata4Model(TimeRef,Area,'SensoryMotorIntegerInteractMdl');

    % take PCA from the color weights to calculate the angle
    FieldsGLMInteract=fieldnames(GLMWeightsAll_InteractMdl);
    PCAField={'ColorIntMLxRule1','ColorIntMLxRule2','ColorIntMLxRule3'};
    [coeff,score,latent,tsquared,explained,mu]=cellfun(@(x) pca(GLMWeightsAll_InteractMdl.(x)),PCAField,'UniformOutput',0);

    % measure angles
    COefs={coeff};
    Explaineds={explained};
    RuleCombs=[1 2;1 3;2 3];
    for f=1:length(COefs)
        NPC=length(Explaineds{f}{1});
        clear PrincipalAngle SumPrincipalAngle
        for npc=1:NPC
            for n=1:3
                VaVb=COefs{f}{RuleCombs(n,1)}(:,1:npc)'*COefs{f}{RuleCombs(n,2)}(:,1:npc);
                [U,S,V]=svd(VaVb);
                PrincipalAngle(n,npc)=acosd(S(1,1));
                SumPrincipalAngle(n,npc)=sum(acosd(diag(S)));
            end
        end
        subplot(2,2,2+f)
        yyaxis left
        hold on
        arrayfun(@(x) plot(PrincipalAngle(x,:),'color',ColorPallet(x,:),'LineWidth',5,'LineStyle','-'),1:3)
        xlabel('included PC')
        ylabel('Principal Angle')
        % title(Txt{f});
        axis square
        % legend({'1-2','1-3','2-3'});

        yyaxis right
        plot(cumsum(Explaineds{f}{1}))
        ylabel('Explained Variance')
        axis tight
        fp.ForcePutLegends({'1-2','1-3','2-3','ev'},ColorPallet(1:4,:),[{'-','-','-','--'}])
    end
end
%% measure compression index for each rule by removing the choice axis
ProjCondName={'QuadrantObjsxRule1','QuadrantObjsxRule2','QuadrantObjsxRule3'};
AxisNames={'ColorIntML','ShapeIntML'};
AxisCombinations=[1 2];ax=1;
nf=cellfun(@(x) find(strcmp(WantedFieldName,x)),AxisNames(AxisCombinations(ax,:))); % find the index of the axes

[CatPSTH_CPI,zCatPSTH_CPI,ConcatzPSTHall_CPI,CatInds_CPI,CatIndsTot_CPI,PSTHTime_CPI]=...
GetPSTHdata4PopAnalysis(Area,Time(TimeIndGLM),1,TimeRef,ProjCondName,WidthSmoothingPSTH,0);

FieldPSTH_CPI=fieldnames(CatPSTH_CPI)';
CondNum=cellfun(@(x) find(strcmp(FieldPSTH_CPI,x)),ProjCondName);

if UsePCdenoising
    for f=FieldPSTH_CPI
        ProjPSTHonTaskAxis_CPI.(f{1})=arrayfun(@(x) TaskAxis'*([zCatPSTH_CPI.(f{1})(x).SpkCountBin'*Ddenoise]'),1:size(CatPSTH_CPI.(f{1}),2),'UniformOutput',0);
    end
else
    for f=FieldPSTH_CPI
        ProjPSTHonTaskAxis_CPI.(f{1})=arrayfun(@(x) TaskAxis'*zCatPSTH_CPI.(f{1})(x).SpkCountBin,1:size(CatPSTH_CPI.(f{1}),2),'UniformOutput',0);
    end
end

% now calculate compression index for each rule
Data=arrayfun(@(x) arrayfun(@(n) smoothdata(ProjPSTHonTaskAxis_CPI.(FieldPSTH_CPI{x}){n}(nf,:)', ...
    1,'movmean',WsmoothingPCA),1:4,'uniformoutput',0),CondNum,'UniformOutput',0);
[CompressionIndex,ColorDistData,ShapeDistData]=ManData.CalCPIOrthEncodingAxis(Data);
ColorDistDataMean=cellfun(@(x) mean(x(PSTHTime_CPI>=TimeAvgCPI(1) & PSTHTime_CPI<=TimeAvgCPI(2))),ColorDistData,'uniformoutput',0);
ShapeDistDataMean=cellfun(@(x) mean(x(PSTHTime_CPI>=TimeAvgCPI(1) & PSTHTime_CPI<=TimeAvgCPI(2))),ShapeDistData,'uniformoutput',0);
CompressionIndexMean=arrayfun(@(x) log(ColorDistDataMean{x}/ShapeDistDataMean{x}),1:3,'UniformOutput',0);

h6=figure;
subplot(221)
hold on
arrayfun(@(y) fp.PlotMeanStd(PSTHTime_CPI,CompressionIndex{y}',[],['Time from ' TimeRef],'CPI',AnalysisOpts.RuleColors(y,:), ...
    1,'Compression Index for Rules','LegendTxt',ProjCondName{y}),1:3,'UniformOutput',0);

subplot(222)
hold on
% arrayfun(@(y) fp.BarPlot(y,mean(CompressionIndex{y}(PSTHTime_CPI>=TimeAvgCPI(1) & PSTHTime_CPI<=TimeAvgCPI(2))),AnalysisOpts.RuleColors(y,:),['Rule'],'CPI', ...
%     'Compression Index for Rules','LegendTxt',ProjCondName{y}),1:3,'UniformOutput',0);
arrayfun(@(y) fp.BarPlot(y,CompressionIndexMean{y},AnalysisOpts.RuleColors(y,:),['Rule'],'CPI', ...
    'Compression Index for Rules','LegendTxt',ProjCondName{y}),1:3,'UniformOutput',0);

subplot(223)

hold on
arrayfun(@(y) fp.PlotMeanStd(ColorDistDataMean{y},ShapeDistDataMean{y},[],['Distance on Color'],'Distance on Shape', ...
    AnalysisOpts.RuleColors(y,:),1,'Compression Index for Rules','LegendTxt',ProjCondName{y},'p_marker','d',...
    'p_marker_size',20),1:3,'UniformOutput',0);

v=axis;
Amin=min(v);Amax=max(v);
xlim([Amin Amax]);ylim([Amin Amax])
fp.PlotXYLine([],1)

% show normalized distances as well
subplot(224)
hold on
arrayfun(@(y) fp.PlotMeanStd(ColorDistDataMean{y}/(ColorDistDataMean{y}+ShapeDistDataMean{y}), ...
    ShapeDistDataMean{y}/(ColorDistDataMean{y}+ShapeDistDataMean{y}),[],['Distance on Color'],'Distance on Shape', ...
    AnalysisOpts.RuleColors(y,:),1,'Compression Index for Rules','LegendTxt',ProjCondName{y},'p_marker','d',...
    'p_marker_size',20),1:3,'UniformOutput',0);
v=axis;
Amin=min(v)-0.1;Amax=max(v)+0.1;
xlim([Amin Amax]);ylim([Amin Amax])
fp.PlotXYLine([],1)

% format all of the axes
fp.FormatAllAxesFig(h6);

%% Measure Compression index during learning
AnalysisOpts.FactorizedDataBuff=AnalysisOpts.FactorizedData;
ProjCondName={'QuadrantObjsxRule1','QuadrantObjsxRule3'};
AxisNames={'ColorIntML','ShapeIntML'};
AxisCombinations=[1 2];ax=1;
nf=cellfun(@(x) find(strcmp(WantedFieldName,x)),AxisNames(AxisCombinations(ax,:))); % find the index of the axes

% generate compression index during learning and total 
TrlRange =ManData.GenMovingInds(1,125,50,10);
clear CompressionIndexR1 CompressionIndexR3
for trl=1:length(TrlRange)
    fprintf('\nCalculating CPI for Trial Range %i',trl)
    AnalysisOpts.FactorizedData=AnalysisOpts.FactorizedDataBuff;
    AnalysisOpts.FactorizedData=NeuAna.LimitFactorizedDataTrialsBasedonFactor(AnalysisOpts.FactorizedData,'TrialNum',TrlRange{trl});

    [CatPSTH_CPI,zCatPSTH_CPI,~,~,~,PSTHTime_CPI]=...
        GetPSTHdata4PopAnalysis(Area,Time(TimeIndGLM),1,TimeRef,ProjCondName,WidthSmoothingPSTH,0);

    FieldPSTH_CPI=fieldnames(CatPSTH_CPI)';
    CondNum=cellfun(@(x) find(strcmp(FieldPSTH_CPI,x)),ProjCondName);

    if UsePCdenoising
        for f=FieldPSTH_CPI
            ProjPSTHonTaskAxis_CPI.(f{1})=arrayfun(@(x) TaskAxis'*([zCatPSTH_CPI.(f{1})(x).SpkCountBin'*Ddenoise]'),1:size(CatPSTH_CPI.(f{1}),2),'UniformOutput',0);
        end
    else
        for f=FieldPSTH_CPI
            ProjPSTHonTaskAxis_CPI.(f{1})=arrayfun(@(x) TaskAxis'*zCatPSTH_CPI.(f{1})(x).SpkCountBin,1:size(CatPSTH_CPI.(f{1}),2),'UniformOutput',0);
        end
    end

    % now calculate compression index for each rule
    Data=arrayfun(@(x) arrayfun(@(n) smoothdata(ProjPSTHonTaskAxis_CPI.(FieldPSTH_CPI{x}){n}(nf,:)', ...
        1,'movmean',WsmoothingPCA),1:4,'uniformoutput',0),CondNum,'UniformOutput',0);
    CompressionIndexR1{trl}=ManData.CalCPIOrthEncodingAxis(Data(1));
    CompressionIndexR3{trl}=ManData.CalCPIOrthEncodingAxis(Data(2));
end
ColLearning=copper(length(TrlRange));

h7=figure;
subplot(221)
arrayfun(@(y) fp.PlotMeanStd(PSTHTime_CPI,CompressionIndexR1{y}{1}',[],['Time from ' TimeRef],'CPI',ColLearning(y,:), ...
    1,'Compression Index for Rule 1'),1:length(TrlRange),'UniformOutput',0) 

subplot(222)
arrayfun(@(y) fp.BarPlot(y,mean(CompressionIndexR1{y}{1}(PSTHTime_CPI>=TimeAvgCPI(1) & PSTHTime_CPI<=TimeAvgCPI(2))),ColLearning(y,:),['Trl'],'CPI', ...
    'Compression Index for Rule 1'),1:length(TrlRange),'UniformOutput',0);

subplot(223)
arrayfun(@(y) fp.PlotMeanStd(PSTHTime_CPI,CompressionIndexR3{y}{1}',[],['Time from ' TimeRef],'CPI',ColLearning(y,:), ...
    1,'Compression Index for Rule 3'),1:length(TrlRange),'UniformOutput',0)

subplot(224)
arrayfun(@(y) fp.BarPlot(y,mean(CompressionIndexR3{y}{1}(PSTHTime_CPI>=TimeAvgCPI(1) & PSTHTime_CPI<=TimeAvgCPI(2))),ColLearning(y,:),['Trl'],'CPI', ...
    'Compression Index for Rule 3'),1:length(TrlRange),'UniformOutput',0);

fp.FormatAllAxesFig(h7);
AnalysisOpts.FactorizedData=AnalysisOpts.FactorizedDataBuff;

h=[h1 h2 h3 h4 h5 h6 h7];
%% save what we have  
if equlCongIncongWeight;EqTxt='_EqCongInCong';else;EqTxt='';end

TDRAnaFigsFileName=[AnalysisOpts.ResultsSavePath filesep 'Low Dim Analysis' filesep,...
    'TargetedDimRed_' TimeRef '_'  Area  '_BestLambda',num2str(BestLambdaGLM) EqTxt '_PC' num2str(PCexPercent)];

fp.SaveFigSeries([],TDRAnaFigsFileName,h,'SaveEachFrame',1)

% fp.SaveCurrentFigs(['TargetedDimRed_' TimeRef '_'  Area  '_BestLambda',num2str(BestLambdaGLM) EqTxt ], ...
%    [AnalysisOpts.ResultsSavePath filesep 'Low Dim Analysis' filesep])
close all
%% look at the CPDs
%first at R2
% Model='SensoryMotorIntegerInteractMdl';
% R2Time=cell2mat(arrayfun(@(x) Mdl(x).(['GLMfitSingFact_' Model]).(Model).R2time',1:480,'UniformOutput',0));
% Fields=fieldnames(Mdl_InteractMdl(1).(['GLMfitSingFact_' Model]).CPD);
% Lam=4;
% figure;hold on
% Facts=7:9;%1:length(Fields);
% for f=Fields(Facts)'
%     Val=cell2mat(arrayfun(@(x) Mdl_InteractMdl(x).(['GLMfitSingFact_' Model]).CPD.(f{1})(Lam,:)',1:480,'UniformOutput',0));
%     plot(Time,nanmean(SM(Val,1,15),2))
% end
% legend(Fields(Facts))
end
