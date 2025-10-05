% does the analysis of Targeted Dimensionality Reduction of Mante et al 2013
function AngleAnalysisWithTDR(GLMmdlName,TimeRef,Area,BestLambdaGLM,equlCongIncongWeight,PCexPercent)
%% load varibales
global AnalysisOpts

NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
fp=fig_params;


if isscalar(GLMmdlName)
    GLMmdlName=AnalysisOpts.SingCellAna.GLMMdlNameSet{GLMmdlName};
end
%TimeRef='SAMPLE_ON';%'SACCADE_START';%
%Area='PFC';
%GLMpath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\GLM\ALL\';
AreaNum=find(strcmp(AnalysisOpts.AreaNames,Area));
AnalysisOpts.NNeuAnalysis=AnalysisOpts.AreaTotalNumNeurons(strcmp(AnalysisOpts.AreaNames,Area));% get the number of neurons based on the area
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
%GLMmdlName='SensoryMotorIntegerInteractMdl';
[GLMWeightsAll,Mdl,Time,TimeIndGLM,WantedFieldName]=LoadGLMdata4ModelFitSamp(TimeRef,Area,GLMmdlName,BestLambdaGLM);
 
ManData.CopyVars2AnalysisOpts(AT('GLMWeightsAll'),GLMWeightsAll,AT('Mdl'),Mdl,AT('Time'),Time, ...
    AT('TimeIndGLM'),TimeIndGLM,AT('WantedFieldName'),WantedFieldName); % copy data to 

[GLMWeightsAll,Mdl,Time,TimeIndGLM,WantedFieldName]=ManData.ExtractVarsFromAnalysisOpts(AT('GLMWeightsAll'),AT('Mdl'),AT('Time'), ...
    AT('TimeIndGLM'),AT('WantedFieldName'));

NWantedFactors=length(WantedFieldName);
TimeGLM=Time(TimeIndGLM);

%% Get PSTH data for factors
AnalysisOpts.Ch_2look_AreaName{1}=Area;
[CatPSTH,zCatPSTH,ConcatzPSTHall,CatInds,CatIndsTot,PSTHTime]=GetPSTHdata4PopAnalysis(Area,TimeGLM,1,TimeRef,[],WidthSmoothingPSTH,equlCongIncongWeight);
AnalysisOpts.Time=PSTHTime;
[D,score,latent,tsquared,explained]=pca(ConcatzPSTHall);


for Lambda=[1:6]
    %% change Lambda
    %if we are changing the data based on the Lambda
    BestLambdaIndNew=Lambda*ones(1,AnalysisOpts.NNeuAnalysis);
    FitTypeSet={'GLMfit','GLMfitObsrv','GLMfitSh'};
    FieldNames=Mdl(1).(['GLMfit_' GLMmdlName])(1).(GLMmdlName).GLMfactorsnames;
    clear GLMWeightsAll GLMWeights
    for FitType=FitTypeSet
        fprintf('\nReplacing Lambda for %s',FitType{1})
        for Factor=WantedFieldName
            IndFactor=find(strcmp(FieldNames,Factor{1}));
            % get Beta values
            for i=1:AnalysisOpts.NNeuAnalysis
                for rep=1:length(Mdl(i).([FitType{1} '_' GLMmdlName]))
                    GLMWeights=squeeze(Mdl(i).([FitType{1} '_' GLMmdlName])(rep).(GLMmdlName).GLM_weights_full(:,BestLambdaIndNew(i),:));
                    GLMWeightsAll.(FitType{1}).(Factor{1})(:,i,rep)=GLMWeights(IndFactor,TimeIndGLM)';
                end
            end
        end
    end

    %% Do the QR decomposition for each rule seperately
    RuleFactors={{'ColorIntMLxRule1','ShapeIntMLxRule1','Axis1IntxRule1','Rule1','Reward'},...
        {'ColorIntMLxRule2','ShapeIntMLxRule2','Axis2IntxRule2','Rule2','Reward'},...
        {'ColorIntMLxRule3','ShapeIntMLxRule3','Axis1IntxRule3','Rule3','Reward'}};

    % cal angle for resampling
    [angd_GLMfit,angd_WithinRule_GLMfit]=CalAngleWithTDR(GLMWeightsAll.GLMfit,WidthSmoothingGLM,PCexPercent,WantedFieldName,TimeGLM,TimeLook4Max,RuleFactors,D,explained);

    % cal angle for observed
    [angd_GLMfitObsrv,angd_WithinRule_GLMfitObsrv]=CalAngleWithTDR(GLMWeightsAll.GLMfitObsrv,WidthSmoothingGLM,PCexPercent,WantedFieldName,TimeGLM,TimeLook4Max,RuleFactors,D,explained);

    % cal angle for shuffle
    [angd_GLMfitSh,angd_WithinRule_GLMfitSh]=CalAngleWithTDR(GLMWeightsAll.GLMfitSh,WidthSmoothingGLM,PCexPercent,WantedFieldName,TimeGLM,TimeLook4Max,RuleFactors,D,explained);

    AllAnglesStruct=[];
    AllAnglesStruct=ManData.CopyVars2Struct(AllAnglesStruct,'angd_GLMfit',angd_GLMfit,'angd_WithinRule_GLMfit',angd_WithinRule_GLMfit,'angd_GLMfitObsrv',angd_GLMfitObsrv,...
        'angd_WithinRule_GLMfitObsrv',angd_WithinRule_GLMfitObsrv,'angd_GLMfitSh',angd_GLMfitSh,'angd_WithinRule_GLMfitSh',angd_WithinRule_GLMfitSh);

    AllAnglesStruct=structfun(@(x) ManData.degreesto0to90(x),AllAnglesStruct,'UniformOutput',0);

    AllAnglesStruct_Mean=structfun(@(x) rad2deg(circ_mean(deg2rad(x),[],3)),AllAnglesStruct,'UniformOutput',0);
    AllAnglesStruct_Std=structfun(@(x) rad2deg(circ_std(deg2rad(x),[],[],3)),AllAnglesStruct,'UniformOutput',0);

    % do the stat test for each condition
    [row,col] = ind2sub([4 3],1:12);
    pval_acrossrule=arrayfun(@(x) ManData.CalpValShuffle(squeeze(AllAnglesStruct.angd_GLMfitSh(row(x),col(x),:)),...
        AllAnglesStruct.angd_GLMfitObsrv(row(x),col(x)),'pvaltail','left'),1:12);
    pval_acrossrule=reshape(pval_acrossrule,[4,3]);

    pval_withinrule=arrayfun(@(x) ManData.CalpValShuffle(squeeze(AllAnglesStruct.angd_WithinRule_GLMfitSh(row(x),col(x),:)),...
        AllAnglesStruct.angd_WithinRule_GLMfitObsrv(row(x),col(x)),'pvaltail','left'),1:12);
    pval_withinrule=reshape(pval_withinrule,[4,3]);

    pval_all=[pval_acrossrule pval_withinrule(:,2:3)];

    % now plot the data with stat test and everything
    FactorNames={'Color','Shape','Axis','Rule'};

    figure
    n=1;
    for f=[1 3]
        subplot(1,2,n)
        hold on
        fp.PlotMeanStd(1:5,[AllAnglesStruct_Mean.angd_GLMfit(f,:) AllAnglesStruct_Mean.angd_WithinRule_GLMfit(f,2:3)],...
            [AllAnglesStruct_Std.angd_GLMfit(f,:)  AllAnglesStruct_Std.angd_WithinRule_GLMfit(f,2:3)],[],[],1,2,[],'STD_method','bootstrap');

        fp.PlotMeanStd(1.2:5.2,[AllAnglesStruct_Mean.angd_GLMfitSh(f,:) AllAnglesStruct_Mean.angd_WithinRule_GLMfitSh(f,2:3)],...
            [AllAnglesStruct_Std.angd_GLMfitSh(f,:)  AllAnglesStruct_Std.angd_WithinRule_GLMfitSh(f,2:3)],[],[],[0.5 0.5 0.5],0,[],...
            'p_line_style','none','p_marker_size',15,'p_marker','.','STD_method','bootstrap');

        arrayfun(@(x) fp.AddDetailedSignificanceStar(x,pval_all(f,x),'k',[],'SigStar_fontsize',20),1:5)

        xlabel('Task Pair')
        ylabel('Angle')
        xticks([1:5]);
        xticklabels({'1-2','1-3','2-3','within 2','within 3'});
        xtickangle(45)
        title(['Angle bet tasks for ' FactorNames{f}]);
        axis square
        n=n+1;
    end
    AllAnglesStructLambda{Lambda}=AllAnglesStruct;
end
fp.SaveCurrentFigs(['GLM_AngleAnalysis_' GLMmdlName],[AnalysisOpts.ResultsSavePath 'Low Dim Analysis' filesep],'SaveEachFrame',1)

% save the angle data so we can plot it in the later time 
ManData.SaveVar('GLM',AllAnglesStructLambda,['AllAnglesStructLambda_' GLMmdlName],[GLMmdlName '_AngleStats'],'WantedDate','ALL','SaveAnalysisOpts',0); 
           
%% plot the co-efficients for each axis
% h1=fp.RenderFigure(1,[]);
% [h1,Sp]=fp.RenderSubplots([],[],h1{1},4);
% GroupsAxisY={[1 2 3],[5,6,7],[8:10],[11:13]}; % groups of factors to plot together
% GroupsAxisYY={[1 3],[5]}; % groups of factors to plot together
% 
% hold on
% arrayfun(@(s) arrayfun(@(y) fp.PlotMeanStd(TimeGLM,GLMWeightsNorm{y},[],['Time from ' TimeRef],'Value',ColorPallet(y,:), ...
%     1,WantedFieldName{y},'LegendTxt',WantedFieldName{y},'Sp',Sp(s)),GroupsAxisY{s},'UniformOutput',0),1:4,'UniformOutput',0);
% 
% % add a star for the maximum value 
% arrayfun(@(s) arrayfun(@(y) fp.Text(TimeGLM(MaxTimeNorm(y)),GLMWeightsNorm{y}(MaxTimeNorm(y)),'*',ColorPallet(y,:),'font_size',15,'Sp', ...
%     Sp(s)),GroupsAxisY{s},'UniformOutput',0),1:4,'UniformOutput',0);
end
