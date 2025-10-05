% uses GLM Sampling for measuring the angles
global AnalysisOpts

NeuAna=NeuralAnalysisFuncs;
ManData=ManipulateData;
fp=fig_params;
ColorCatColors=AnalysisOpts.MorphlevelsColRGBInc50([1 5 2 3 4],:);
Path='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\GLM\ALL\';
%[~,Time]=ManData.GenerateTimeAxis(NeuAna.PSTHTimRef);
load([Path 'GLM_ALL_ALL_480SingCellAna_PFC_SAMPLE_ON_100_Summery.mat'], 'AnalysisOpts');
MdlName='GLMfit_SensoryCatMotorRulePureInteracMdl';
for i=1:480
    tic
    FileName=sprintf('%sGLM_ALL_ALL_%iSingCellAna_PFC_SAMPLE_ON_100_GLMfitSamp.mat',Path,i);
    Mdl{i}=load(FileName,MdlName);
    fprintf('\nLoading file %i in %0.2f secs...',i,toc);
end

FieldNames=Mdl{1}.(MdlName)(1).SensoryCatMotorRulePureInteracMdl.GLMfactorsnames;

% Define if we are using CPD as  Weights or the Beta weights
WeightDef='Beta';NNeu=1:480;
Time=AnalysisOpts.Time;
Lam=6;% Lambda
Tcol=[0.2 0.25];
Tresp=[0.3 0.35];
 
ColorPallet=AnalysisOpts.ColorPalettCamden;
ColorPallet2=AnalysisOpts.ColorPalett1COSYNE23;
LineStyle={'-','--',':'};
MeanCenter=@(x) x-mean(x,1);
Txt={'Per category','Per Condition'};
for Factor={'color'};%,'responseloc'}
    clear GLMWeightsConcat GLMWeightAll PrincipalAngle SumPrincipalAngle
    for nrep=1:100

        IndFactor=find(contains(FieldNames,Factor{1},'IgnoreCase',true));
        %% prepare data
        for i=NNeu
            GLMWeights=squeeze(Mdl{i}.(MdlName)(nrep).SensoryCatMotorRulePureInteracMdl.GLM_weights_full(:,Lam,:));
            switch Factor{1}
                case {'color','shape'}
                    Inds=[1 2;3 4;5 6];
                    Ncomps=2;
                    TimInd=find(Time>=Tcol(1) & Time<=Tcol(2));
                    npoints=length(TimInd);
                    TimeAxis=Time(TimInd);
                    YLIM=[30 100];

                    GLMWeightsConcat{1}(:,i)=[GLMWeights(IndFactor(1),TimInd),GLMWeights(IndFactor(2),TimInd)]';
                    GLMWeightsConcat{2}(:,i)=[GLMWeights(IndFactor(3),TimInd),GLMWeights(IndFactor(4),TimInd)]';
                    GLMWeightsConcat{3}(:,i)=[GLMWeights(IndFactor(5),TimInd),GLMWeights(IndFactor(6),TimInd)]';

                    for k=1:3
                        GLMWeightAll(:,i,k)=[arrayfun(@(x) norm(GLMWeights(IndFactor(Inds(k,1)):IndFactor(Inds(k,2)),x)),TimInd)]';
                        % GLMWeightAll(:,i,k)=[arrayfun(@(x) GLMWeights(IndFactor(Inds(k,1)),x)-GLMWeights(IndFactor(Inds(k,2)),x),TimInd)]';
                    end
                case 'rule'
                    Inds=[1 1;2 2;3 3];
                    Ncomps=1;
                    TimInd=find(Time>=-0.4 & Time<=0);
                    npoints=length(TimInd);
                    TimeAxis=Time(TimInd);

                    GLMWeightsConcat{1}(:,i)=[GLMWeights(IndFactor(1),TimInd)]';
                    GLMWeightsConcat{2}(:,i)=[GLMWeights(IndFactor(2),TimInd)]';
                    GLMWeightsConcat{3}(:,i)=[GLMWeights(IndFactor(3),TimInd)]';

                    for k=1:3
                        GLMWeightAll(:,i,k)=[arrayfun(@(x) norm(GLMWeights(IndFactor(Inds(k,1)):IndFactor(Inds(k,2)),x)),TimInd)]';
                    end
                case 'responseloc'
                    Inds=[1 4;5 8;9 12];
                    Ncomps=4;
                    TimInd=find(Time>=Tresp(1) & Time<=Tresp(2));
                    npoints=length(TimInd);
                    TimeAxis=Time(TimInd);
                    YLIM=[0 100];

                    GLMWeightsConcat{1}(:,i)=[GLMWeights(IndFactor(1),TimInd),GLMWeights(IndFactor(2),TimInd),GLMWeights(IndFactor(3),TimInd),GLMWeights(IndFactor(4),TimInd)]';
                    GLMWeightsConcat{2}(:,i)=[GLMWeights(IndFactor(5),TimInd),GLMWeights(IndFactor(6),TimInd),GLMWeights(IndFactor(7),TimInd),GLMWeights(IndFactor(8),TimInd)]';
                    GLMWeightsConcat{3}(:,i)=[GLMWeights(IndFactor(9),TimInd),GLMWeights(IndFactor(10),TimInd),GLMWeights(IndFactor(11),TimInd),GLMWeights(IndFactor(12),TimInd)]';

                    for k=1:3
                        %   GLMWeightAll(:,i,k)=[arrayfun(@(x) norm(GLMWeights(IndFactor(Inds(k,1)):IndFactor(Inds(k,2)),x)),TimInd)]';
                        GLMWeightAll(:,i,k)=[arrayfun(@(x) GLMWeights(IndFactor(Inds(k,1)),x)-GLMWeights(IndFactor(Inds(k,2)),x),TimInd)]';
                    end
            end
        end

        %% take PCA of color category per condition now
        [coeffCond,score,latent,tsquared,explainedCond,mu]=arrayfun(@(x) pca(GLMWeightsConcat{x}),1:3,'UniformOutput',0);
        [coeff,score,latent,tsquared,explained,mu]=arrayfun(@(x) pca(GLMWeightAll(:,:,x)),1:3,'UniformOutput',0);

        %% measure angles
        COefs=[{coeffCond},{coeff}];
        Explaineds=[{explainedCond},{explained}];
        RuleCombs=[1 2;1 3;2 3];
        %  clear PrincipalAngle SumPrincipalAngle
        for f=1:2
            NPC=length(Explaineds{f}{1});
            for npc=1:NPC
                for n=1:3
                    VaVb=COefs{f}{RuleCombs(n,1)}(:,1:npc)'*COefs{f}{RuleCombs(n,2)}(:,1:npc);
                    [U,S,V]=svd(VaVb);
                    PrincipalAngle{f}(n,npc,nrep)=acosd(S(1,1));
                    SumPrincipalAngle{f}(n,npc,nrep)=sum(acosd(diag(S)));
                end
            end
        end
    end

    %% plot results
    figure
    for f=1:2
        subplot(2,2,f)
        hold on
        arrayfun(@(x) fp.PlotMeanStd(1:size(PrincipalAngle{f},2),squeeze(PrincipalAngle{f}(x,:,:))',[],[],[],ColorPallet(x,:),1,[],'STD_method','bootstrap'),1:3)

        xlabel('included PC')
        ylabel('Principal Angle')
        title([Txt{f} Factor{1}]);
        axis square

        subplot(2,2,2+f)
        hold on
        %X,Y,YSTD,Xlbl,Ylbl,Col,Shaded,Title
        arrayfun(@(x) fp.PlotMeanStd(x,squeeze(PrincipalAngle{f}(x,end,:)),[],[],[],ColorPallet(x,:),2,[],...
            'STD_method','bootstrap','ForceShaded',1),1:3)
        xlabel('included PC')
        ylabel('Principal Angle')
        title(Txt{f});
        axis square
        xticks(1:3)
        xticklabels({'1-2','1-3','2-3'})
        xtickangle(45)
        ylim(YLIM)
        fp.PlotHorizontalLine(90,gca,[0.5 0.5 0.5],'p_line_style','--')
        % do stat test on pairs of conditions
        pairs=nchoosek(1:3,2);
        [~,p_pairs]=arrayfun(@(x) ttest2(squeeze(PrincipalAngle{f}(pairs(x,1),end,:)),squeeze(PrincipalAngle{f}(pairs(x,2),end,:))),1:3);
        p_pairs=p_pairs*3;
        groups={pairs(1,:),pairs(2,:),pairs(3,:)};
        sigstar(groups,p_pairs)
    end
end


%% put a figure with details of the angle analysis
fp.SaveCurrentFigs('GLM_AngleAnalysis','Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Results\Low Dim Analysis\','SaveEachFrame',1)

