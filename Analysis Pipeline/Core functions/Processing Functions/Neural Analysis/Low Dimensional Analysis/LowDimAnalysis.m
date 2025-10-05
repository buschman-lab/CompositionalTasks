global AnalysisOpts

NeuAna=NeuralAnalysisFuncs;
ManData=ManipulateData;
fp=fig_params;
ColorCatColors=AnalysisOpts.MorphlevelsColRGBInc50([1 5 2 3 4],:);
Path='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Input Output Data\GLM\ALL\';
%[~,Time]=ManData.GenerateTimeAxis(NeuAna.PSTHTimRef);
load([Path 'GLM_ALL_ALL_480SingCellAna_PFC_SAMPLE_ON_100_Summery.mat'], 'AnalysisOpts');

for i=1:480
    tic
    FileName=sprintf('%sGLM_ALL_ALL_%iSingCellAna_PFC_SAMPLE_ON_100_Summery.mat',Path,i);
    Mdl(i)=load(FileName,'GLMfitSingFact_SensoryCatMotorRulePureInteracMdl');
    fprintf('\nLoading file %i in %0.2f secs...',i,toc);
end

FieldNames=Mdl(1).GLMfitSingFact_SensoryCatMotorRulePureInteracMdl.SensoryCatMotorRulePureInteracMdl.GLMfactorsnames;

% Define if we are using CPD as  Weights or the Beta weights
WeightDef='Beta';NNeu=1:480;
Time=AnalysisOpts.Time;
Lam=5;% Lambda
clear GLMWeightsConcat GLMWeightAll

ColorPallet=AnalysisOpts.ColorPalettCamden;
ColorPallet2=AnalysisOpts.ColorPalett1COSYNE23;
LineStyle={'-','--',':'};
MeanCenter=@(x) x-mean(x,1);
Txt={'Per category','Per Condition'};
for Factor={'color','shape','responseloc'} 
    IndFactor=find(contains(FieldNames,Factor{1},'IgnoreCase',true));
    clear GLMWeightsConcat GLMWeightAll
    %% prepare data
    switch WeightDef
        case 'Beta'
            for i=NNeu
                GLMWeights=squeeze(Mdl(i).GLMfitSingFact_SensoryCatMotorRulePureInteracMdl.SensoryCatMotorRulePureInteracMdl.GLM_weights_full(:,Lam,:));
                %           GLMWeightsConcat= cell2mat(arrayfun(@(x) [norm(GLMWeights(1:2,x));norm(GLMWeights(3:4,x));...
                %               norm(GLMWeights(5:6,x));norm(GLMWeights(7:8,x));norm(GLMWeights(9:10,x));norm(GLMWeights(11:12,x))],...
                %               TimInd,'UniformOutput',0));
                %           GLMWeightsConcatRS(:,i)=reshape(GLMWeightsConcat',[6*size(TimInd,2),1]);
                switch Factor{1}
                    case {'color','shape'}
                        Inds=[1 2;3 4;5 6];
                        Ncomps=2;
                        TimInd=find(Time>=0.2 & Time<=0.25);                      
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
                        TimInd=find(Time>=0.3 & Time<=0.35);
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

        case 'CPD'
            Fields=fieldnames(Mdl(1).GLMfitSingFact_SensoryCatMotorRulePureInteracMdl.CPD);
            Lam=5;
            nFeild=4:9;
            for i=NNeu
                for f=1:6
                    FieldName=Fields{nFeild(f)};
                    GLMWeightAll(:,:,f)=cell2mat(arrayfun(@(x) Mdl(x).GLMfitSingFact_SensoryCatMotorRulePureInteracMdl.CPD.(FieldName)(Lam,TimInd)',1:480,'UniformOutput',0));
                end
            end
    end

    %% take PCA of color category per condition now
    figure
    sgtitle(['plotting subspaces for ' Factor{1}])
    [coeffCond,score,latent,tsquared,explainedCond,mu]=arrayfun(@(x) pca(GLMWeightsConcat{x}),1:3,'UniformOutput',0);

    % project all to the space of color in R2
    ProjScore=arrayfun(@(x) MeanCenter(GLMWeightsConcat{x})*coeffCond{3},1:3,'UniformOutput',0);    
    subplot(221)
    hold on
    Rule=1:3;
    arrayfun(@(n) arrayfun(@(x) ManData.PlotDatainTime(TimeAxis,smoothdata(ProjScore{x}((1:npoints)+(n-1)*npoints,1:3), ...
        1,'movmean',15),[1 2 3],...
        ColorCatColors(n,:),LineStyle{x},[0 0.2 0.4 0.6],0),Rule),1:Ncomps);

    fp.ForcePutLegends({'R1','R2','R3'},ColorPallet(1:3,:),LineStyle)
    title(Txt{1});

    %% Take PCA of each condition now
    [coeff,score,latent,tsquared,explained,mu]=arrayfun(@(x) pca(GLMWeightAll(:,:,x)),1:3,'UniformOutput',0);

    % project all to the space of color in R2
    ProjScore=arrayfun(@(x) MeanCenter(GLMWeightAll(:,:,x))*coeff{3},1:3,'UniformOutput',0);
    subplot(222)
    hold on
    arrayfun(@(x) ManData.PlotDatainTime(TimeAxis,smoothdata(ProjScore{x}(:,1:3),1,'movmean',15),[1 2 3],...
        ColorPallet(x,:),LineStyle{x},[0 .1 0.2 .3 0.4 .5 0.6],0),1:3);

    fp.ForcePutLegends({'R1','R2','R3'},ColorPallet(1:3,:),LineStyle)
    title(Txt{2});

    %% measure angles
    COefs=[{coeffCond},{coeff}];
    Explaineds=[{explainedCond},{explained}];
    RuleCombs=[1 2;1 3;2 3];
    for f=1:2
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


    figure
    for f=1:2
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
        NPC=length(Explaineds{f}{1});

        subplot(2,2,f)
        hold on
        arrayfun(@(x) bar(x,PrincipalAngle(x,end),'facecolor',ColorPallet(x,:)),1:3)
        xlabel('included PC')
        ylabel('Principal Angle')
        title(Txt{f});
        axis square
        xticks(1:3)
        xticklabels({'1-2','1-3','2-3'})
        xtickangle(45)
        ylim(YLIM)
        fp.PlotHorizontalLine(90,gca,[0.5 0.5 0.5],'p_line_style','--')

        subplot(2,2,2+f)
        hold on
        arrayfun(@(x) bar(x,SumPrincipalAngle(x,end),'facecolor',ColorPallet(x,:)),1:3)
        xlabel('included PC')
        ylabel('Sum of Principal Angle')
        title(Txt{f});
        axis square
        xticks(1:3)
        xticklabels({'1-2','1-3','2-3'})
        xtickangle(45)

    end

end

%% put a figure with details of the angle analysis 
fp.SaveCurrentFigs('GLM_AngleAnalysis','Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Analysis Pipeline\Results\Low Dim Analysis\','SaveEachFrame',1)

%%
[PC,Score]=pca(GLMWeightsConcatRS(1:6*npoints,:));

figure
hold on
NPoints2Plot=40;
arrayfun(@(x) ManData.PlotDatainTime(TimeAxis(1:NPoints2Plot),smoothdata(Score((1:NPoints2Plot)+(x-1)*npoints,1:3),1,'movmean',15),[1 2 3],...
    ColorPallet(x,:),LineStyle{x},[0 0.2 0.4 0.6]),1:3);

fp.ForcePutLegends({'R1','R2','R3'},ColorPallet(1:3,:),LineStyle)
% for i=1:3
%     arrayfun(@(x) plot3(movmean(Score((1:20)+(i-1)*59,1),5,1),movmean(Score((1:20)+(i-1)*59,2),5,1),movmean(Score((1:20)+(i-1)*59,3),5,1),'Color',AnalysisOpts.ColorPalett3(i,:)*x,'LineWidth',15),1)
% end
% figure
%legend({'R1','R2','R3'})


%% look at the CPDs
%first at R2
R2Time=cell2mat(arrayfun(@(x) Mdl(x).GLMfitSingFact_SensoryCatMotorRulePureInteracMdl.SensoryCatMotorRulePureInteracMdl.R2time',1:480,'UniformOutput',0));
Fields=fieldnames(Mdl(1).GLMfitSingFact_SensoryCatMotorRulePureInteracMdl.CPD);
Lam=5;
figure;hold on
for f=Fields(4:6)'
    Val=cell2mat(arrayfun(@(x) Mdl(x).GLMfitSingFact_SensoryCatMotorRulePureInteracMdl.CPD.(f{1})(Lam,:)',1:480,'UniformOutput',0));
    plot(nanmean(Val,2))
end
legend(Fields(4:6))