 %% load varibales
global AnalysisOpts

NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
fp=fig_params;

TimeRef='SAMPLE_ON';%'SACCADE_START';%
Area='PFC';
BestLambdaGLM=4;
GLMpath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\GLM\ALL\';
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
% AT=@(x) [x '_' TimeRef]; %concatinate Varname with time ref
% GLMmdlName='SensoryMotorIntegerInteractMdl';
% [GLMWeightsAll,GLMCPDAll,Mdl,Time,TimeIndGLM,GLMModelComp,WantedFieldName,BestLambdaInd]=...
%     LoadGLMdata4Model(TimeRef,Area,GLMmdlName,BestLambdaGLM);
% 
% ManData.CopyVars2AnalysisOpts(AT('GLMWeightsAll'),GLMWeightsAll,AT('GLMCPDAll'),GLMCPDAll,AT('Mdl'),Mdl,AT('Time'),Time, ...
%     AT('TimeIndGLM'),TimeIndGLM,AT('GLMModelComp'),GLMModelComp,AT('BestLambdaInd'),BestLambdaInd,AT('WantedFieldName'),WantedFieldName); % copy data to 
% 
% [GLMWeightsAll,GLMCPDAll,Mdl,Time,TimeIndGLM,GLMModelComp,WantedFieldName,BestLambdaInd]=ManData.ExtractVarsFromAnalysisOpts(AT('GLMWeightsAll'),AT('GLMCPDAll'),AT('Mdl'),AT('Time'), ...
%     AT('TimeIndGLM'),AT('GLMModelComp'),AT('WantedFieldName'),AT('BestLambdaInd'));
% 
% NWantedFactors=length(WantedFieldName);
% TimeGLM=Time(TimeIndGLM);
% 
% % smooth GLM fit data 
% GLMWeightsAll=structfun(@(x) SM(x,1,WidthSmoothingGLM),GLMWeightsAll,'UniformOutput',0);

 
ColorCatColors=AnalysisOpts.MorphlevelsColRGBInc50([1 5 2 3 4],:);
Path='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\GLM\ALL\';
%[~,Time]=ManData.GenerateTimeAxis(NeuAna.PSTHTimRef);
load([Path 'GLM_ALL_ALL_480SingCellAna_PFC_SAMPLE_ON_100_Summery.mat'], 'AnalysisOpts');
MdlName='SensoryMotorIntegerInteractMdl';
SingFactTxt=['GLMfitSingFact_' MdlName];
for i=1:480
    tic
    FileName=sprintf('GLM_ALL_ALL_%iSingCellAna_PFC_SAMPLE_ON_100_Summery.mat',i);
    Mdl(i)=load([Path FileName],['GLMfitSingFact_' MdlName]);
    fprintf('\nLoading file %i in %0.2f secs...',i,toc);
end
FieldNames=Mdl(1).(SingFactTxt).(MdlName).GLMfactorsnames;
% Define if we are using CPD as  Weights or the Beta weights
WeightDef='Beta';NNeu=1:480;
Time=AnalysisOpts.Time;

% clear GLMWeightsConcat GLMWeightAll
% ColorPallet=AnalysisOpts.ColorPalettCamden;
% ColorPallet2=AnalysisOpts.ColorPalett1COSYNE23;
% LineStyle={'-','--',':'};
% MeanCenter=@(x) x-mean(x,1);

Lam=6;% Lambda
Txt={'Per category','Per Condition'};
for EndTime=0.2
for Factor={'color'};%'color'};%,'rule'}
    IndFactor=find(contains(FieldNames,Factor{1},'IgnoreCase',true));
    clear GLMWeightsConcat GLMWeightAll
    %% prepare data
    switch WeightDef
        case 'Beta'
            for i=NNeu
                GLMWeights=squeeze(Mdl(i).(SingFactTxt).(MdlName).GLM_weights_full(:,Lam,:));
                switch Factor{1}
                    case {'color','shape','axis'}
                        Inds=[1 2 3];
                        Ncomps=2;
                        TimInd=find(Time>=0.1 & Time<=EndTime);
                        npoints=length(TimInd);
                        TimeAxis=Time(TimInd);

                        GLMWeightsConcat{1}(:,i)=[GLMWeights(IndFactor(1),TimInd)]';
                        GLMWeightsConcat{2}(:,i)=[GLMWeights(IndFactor(2),TimInd)]';
                        GLMWeightsConcat{3}(:,i)=[GLMWeights(IndFactor(3),TimInd)]';

                        for k=1:3
                            GLMWeightAll(:,i,k)=[mean(GLMWeights(IndFactor(k),TimInd),2)]';
                        end
                    case 'responseloc'
                        Inds=[1 2 3];
                        Ncomps=4;
                        TimInd=find(Time>=0.3 & Time<=0.5);
                        npoints=length(TimInd);
                        TimeAxis=Time(TimInd);

                        GLMWeightsConcat{1}(:,i)=[GLMWeights(IndFactor(1),TimInd)]';
                        GLMWeightsConcat{2}(:,i)=[GLMWeights(IndFactor(2),TimInd)]';
                        GLMWeightsConcat{3}(:,i)=[GLMWeights(IndFactor(3),TimInd)]';

                        for k=1:3
                            GLMWeightAll(:,i,k)=[GLMWeights(IndFactor(k),TimInd)]';
                        end
                    case 'TDR'  % use targeted dimensionality reduction 
                        RuleFactors=[{'Rule1','Reward','ColorIntMLxRule1','ShapeIntMLxRule1','Axis1IntxRule1'},...
                            {'Rule2','Reward','ColorIntMLxRule2','ShapeIntMLxRule2','Axis2IntxRule2'},...
                            {'Rule3','Reward','ColorIntMLxRule3','ShapeIntMLxRule3','Axis2IntxRule3'}];
                        for Rule=1:3
                            % find the indexes for this factors
                            TDRfactorInds{Rule}=cellfun(@(x) find(strcmp(FieldNames,x)),RuleFactors{Rule})
                        end
                end
            end
    end

    %% take PCA of color category per condition now
    figure
    sgtitle(['plotting subspaces for ' Factor{1}])
    [coeffCond,score,latent,tsquared,explainedCond,mu]=arrayfun(@(x) pca(GLMWeightsConcat{x}),1:3,'UniformOutput',0);
% 
%     % project all to the space of color in R2
%     ProjScore=arrayfun(@(x) MeanCenter(GLMWeightsConcat{x})*coeffCond{3},1:3,'UniformOutput',0);    
%     subplot(221)
%     hold on
%     Rule=1:3;Ncomps=1;
%     arrayfun(@(n) arrayfun(@(x) ManData.PlotDatainTime(TimeAxis,smoothdata(ProjScore{x}((1:npoints)+(n-1)*npoints,1:3), ...
%         1,'movmean',15),[1 2 3],...
%         ColorCatColors(n,:),LineStyle{x},[0 0.2 0.4 0.6],0),Rule),1:Ncomps);
% 
%     fp.ForcePutLegends({'R1','R2','R3'},ColorPallet(1:3,:),LineStyle)
%     title(Txt{1});

    %% Take PCA of each condition now
    [coeff,score,latent,tsquared,explained,mu]=arrayfun(@(x) pca(GLMWeightAll(:,:,x)),1:3,'UniformOutput',0);


    %% measure angles
    COefs=[{coeffCond},{coeff}];
    Explaineds=[{explainedCond},{explained}];
    RuleCombs=[1 2;1 3;2 3];
    for f=1
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
        title(Txt{f});
        axis square
       % legend({'1-2','1-3','2-3'});

        yyaxis right
        plot(cumsum(Explaineds{f}{1}))
        ylabel('Explained Variance')
        axis tight
        fp.ForcePutLegends({'1-2','1-3','2-3','ev'},ColorPallet(1:4,:),[{'-','-','-','--'}])
    end
    %% add the angle between average values 
    [angd]=arrayfun(@(x) ManData.GetAngleBetVectors(GLMWeightAll(1,:,RuleCombs(x,1)),GLMWeightAll(1,:,RuleCombs(x,2))),1:3);
    subplot(224)
    bar(angd);
    xticklabels({'1-2','1-3','2-3'})
end
title(num2str(EndTime,2))
end
