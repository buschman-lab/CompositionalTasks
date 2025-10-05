function h=PlotAngleAnaTDRresults(DataPath)
% plot TDR anagle analysis results

global AnalysisOpts


ManData=ManipulateData;
fp=fig_params;


% prepare the data
Models={'SensoryMotorIntegerInteractMdl',...
    'SensoryMotorIntegerInteractMdlFromSwitch',...
    'SensoryMotorIntegerInteractMdlToSwitch'};

PreTxt='AllAnglesStructLambda_';
[row,col] = ind2sub([4 3],1:12);
Lam=4;
for Mdl=Models
    ModelPATH=[DataPath 'GLM_ALL_ALL_75_cl2' Mdl{1} '_AngleStats.mat'];
    ModelData.(Mdl{1})=load(ModelPATH);   

    FullMdlName=[PreTxt Mdl{1}];

    AllAnglesStruct=ModelData.(Mdl{1}).(FullMdlName);
    AllAnglesStruct=AllAnglesStruct{Lam};
    AllAnglesStruct=structfun(@(x) ManData.degreesto0to90(x),AllAnglesStruct,'UniformOutput',0);

    Data.(FullMdlName).AllAnglesStruct_Mean=structfun(@(x) rad2deg(circ_mean(deg2rad(x),[],3)),AllAnglesStruct,'UniformOutput',0);
    Data.(FullMdlName).AllAnglesStruct_Std=structfun(@(x) rad2deg(circ_std(deg2rad(x),[],[],3)),AllAnglesStruct,'UniformOutput',0);
    Data.(FullMdlName).AllAnglesStruct_All=structfun(@(x) x,AllAnglesStruct,'UniformOutput',0);

    % do the stat test for each condition
    Data.(FullMdlName).pval_acrossrule=arrayfun(@(x) ManData.CalpValShuffle(squeeze(AllAnglesStruct.angd_GLMfitSh(row(x),col(x),:)),...
        AllAnglesStruct.angd_GLMfitObsrv(row(x),col(x)),'pvaltail','left'),1:12);
    Data.(FullMdlName).pval_acrossrule=reshape(Data.(FullMdlName).pval_acrossrule,[4,3]);

    Data.(FullMdlName).pval_withinrule=arrayfun(@(x) ManData.CalpValShuffle(squeeze(AllAnglesStruct.angd_WithinRule_GLMfitSh(row(x),col(x),:)),...
        AllAnglesStruct.angd_WithinRule_GLMfitObsrv(row(x),col(x)),'pvaltail','left'),1:12);
    Data.(FullMdlName).pval_withinrule=reshape(Data.(FullMdlName).pval_withinrule,[4,3]);

    Data.(FullMdlName).pval_all=[Data.(FullMdlName).pval_acrossrule Data.(FullMdlName).pval_withinrule(:,2:3)];
end
FMdlNmAll=[PreTxt Models{1}];
FMdlNmFS=[PreTxt Models{2}];
FMdlNmTS=[PreTxt Models{3}];

% now plot the data with stat test and everything
FactorNames={'Color','Shape','Axis','Rule'};
figure
n=1;
DataOrder=[1 4 2 5 3 6 7 8];
for f=[1 3]
    subplot(1,2,n)
    hold on
    MeanData2Show=[Data.(FMdlNmFS).AllAnglesStruct_Mean.angd_GLMfit(f,:) Data.(FMdlNmTS).AllAnglesStruct_Mean.angd_GLMfit(f,:)  Data.(FMdlNmAll).AllAnglesStruct_Mean.angd_WithinRule_GLMfit(f,2:3)];
    StdData2Show= [Data.(FMdlNmFS).AllAnglesStruct_Std.angd_GLMfit(f,:) Data.(FMdlNmTS).AllAnglesStruct_Std.angd_GLMfit(f,:)  Data.(FMdlNmAll).AllAnglesStruct_Std.angd_WithinRule_GLMfit(f,2:3)];
    
    fp.PlotMeanStd(1:8,MeanData2Show(DataOrder),...
       StdData2Show(DataOrder),[],[],1,2,[],'STD_method','bootstrap');

    MeanDataSh2Show=[Data.(FMdlNmFS).AllAnglesStruct_Mean.angd_GLMfitSh(f,:) Data.(FMdlNmTS).AllAnglesStruct_Mean.angd_GLMfitSh(f,:)  Data.(FMdlNmAll).AllAnglesStruct_Mean.angd_WithinRule_GLMfitSh(f,2:3)];
    StdData2ShShow= [Data.(FMdlNmFS).AllAnglesStruct_Std.angd_GLMfitSh(f,:) Data.(FMdlNmTS).AllAnglesStruct_Std.angd_GLMfitSh(f,:)  Data.(FMdlNmAll).AllAnglesStruct_Std.angd_WithinRule_GLMfitSh(f,2:3)];
    
fp.PlotMeanStd((1:8)+0.35,MeanDataSh2Show(DataOrder),...
        StdData2ShShow(DataOrder),[],[],[0.5 0.5 0.5],0,[],...
        'p_line_style','none','p_marker_size',15,'p_marker','.','STD_method','bootstrap');

%     fp.PlotMeanStd(1.2:5.2,[AllAnglesStruct_Mean.angd_GLMfitSh(f,:) AllAnglesStruct_Mean.angd_WithinRule_GLMfitSh(f,2:3)],...
%         [AllAnglesStruct_Std.angd_GLMfitSh(f,:)  AllAnglesStruct_Std.angd_WithinRule_GLMfitSh(f,2:3)],[],[],[0.5 0.5 0.5],0,[],...
%         'p_line_style','none','p_marker_size',15,'p_marker','.','STD_method','bootstrap');
    pval_all=[Data.(FMdlNmFS).pval_acrossrule Data.(FMdlNmTS).pval_acrossrule Data.(FMdlNmAll).pval_withinrule];
    arrayfun(@(x) fp.AddDetailedSignificanceStar(x,pval_all(f,DataOrder(x)),'k',[],'SigStar_fontsize',10,'ThisYval',MeanData2Show(DataOrder(x))),1:8)

    xlabel('Task Pair')
    ylabel('Angle')
    xticks([1:8]);
    xticklabels({'1-2 From Swc', '1-2 To Swc','1-3 From Swc', '1-3 To Swc','3-2 From Swc', '3-2 To Swc','within 2','within 3'});
    xtickangle(45)
    title(['Angle bet tasks for ' FactorNames{f}]);
    axis square
    ylim([0 95])
    n=n+1;
end

% second version of the figure with a line for within rule angles 
figure
n=1;
DataOrder=[1 4 2 5 3 6];
for f=[1 3]
    subplot(1,2,n)
    hold on
    MeanData2Show=[Data.(FMdlNmFS).AllAnglesStruct_Mean.angd_GLMfit(f,:) Data.(FMdlNmTS).AllAnglesStruct_Mean.angd_GLMfit(f,:)  ];
    StdData2Show= [Data.(FMdlNmFS).AllAnglesStruct_Std.angd_GLMfit(f,:) Data.(FMdlNmTS).AllAnglesStruct_Std.angd_GLMfit(f,:) ];
     
    fp.PlotMeanStd(1:6,MeanData2Show(DataOrder),...
       StdData2Show(DataOrder),[],[],1,2,[],'STD_method','bootstrap');

    MeanDataSh2Show=[Data.(FMdlNmFS).AllAnglesStruct_Mean.angd_GLMfitSh(f,:) Data.(FMdlNmTS).AllAnglesStruct_Mean.angd_GLMfitSh(f,:) ];
    StdData2ShShow= [Data.(FMdlNmFS).AllAnglesStruct_Std.angd_GLMfitSh(f,:) Data.(FMdlNmTS).AllAnglesStruct_Std.angd_GLMfitSh(f,:) ];
    
    fp.PlotMeanStd((1:6)+0.35,MeanDataSh2Show(DataOrder),...
        StdData2ShShow(DataOrder),[],[],[0.5 0.5 0.5],0,[],...
        'p_line_style','none','p_marker_size',15,'p_marker','.','STD_method','bootstrap');

%     fp.PlotMeanStd(1.2:5.2,[AllAnglesStruct_Mean.angd_GLMfitSh(f,:) AllAnglesStruct_Mean.angd_WithinRule_GLMfitSh(f,2:3)],...
%         [AllAnglesStruct_Std.angd_GLMfitSh(f,:)  AllAnglesStruct_Std.angd_WithinRule_GLMfitSh(f,2:3)],[],[],[0.5 0.5 0.5],0,[],...
%         'p_line_style','none','p_marker_size',15,'p_marker','.','STD_method','bootstrap');
    pval_all=[Data.(FMdlNmFS).pval_acrossrule Data.(FMdlNmTS).pval_acrossrule];
    arrayfun(@(x) fp.AddDetailedSignificanceStar(x,pval_all(f,DataOrder(x)),'k',[],'SigStar_fontsize',10,'ThisYval',MeanData2Show(DataOrder(x))),1:6)

    % plot within rule as horizontal line
    fp.PlotHorizontalLine(Data.(FMdlNmAll).AllAnglesStruct_Mean.angd_WithinRule_GLMfit(f,2),[],'g','p_line_style','--')
    fp.PlotHorizontalLine(Data.(FMdlNmAll).AllAnglesStruct_Mean.angd_WithinRule_GLMfit(f,3),[],'b','p_line_style','--')
   
    xlabel('Task Pair')
    ylabel('Angle')
    xticks([1:6]);
    xticklabels({'1-2 From Swc', '1-2 To Swc','1-3 From Swc', '1-3 To Swc','3-2 From Swc', '3-2 To Swc'});
    xtickangle(45)
    title(['Angle bet tasks for ' FactorNames{f}]);
    axis square
    ylim([0 95])
    n=n+1;
end

close all
% Third version of the figure with a line for within rule angles and whisker and box plot for others 
figure
n=1;
DataOrder=[1 4 2 5 3 6];
for f=[1 3]
    subplot(1,2,n)
    hold on
    % plot main data distribution using whisker and box plot
    AllData2Show=[squeeze(Data.(FMdlNmFS).AllAnglesStruct_All.angd_GLMfit(f,:,:))' squeeze(Data.(FMdlNmTS).AllAnglesStruct_All.angd_GLMfit(f,:,:))'  ];
    groups=cell2mat(arrayfun(@(x) x*ones(size(AllData2Show,1),1),1:6,'uniformoutput',0));
   
    % rearrange the order of the data
    AllData2Show=AllData2Show(:,DataOrder);
    % groups=groups(:,DataOrder);
    % do a box plot instead and get rid of outlier identification 
    boxplot(AllData2Show(:),groups(:),'Whisker', Inf,'Widths', 0.2)%,'PlotStyle','compact')
    
    
    % plot the shuffle distribution
    AllData2ShowSh=[squeeze(Data.(FMdlNmFS).AllAnglesStruct_All.angd_GLMfitSh(f,:,:))' squeeze(Data.(FMdlNmTS).AllAnglesStruct_All.angd_GLMfitSh(f,:,:))'  ];
    groupsSh=cell2mat(arrayfun(@(x) x*ones(size(AllData2ShowSh,1),1),1:6,'uniformoutput',0));
   
    boxplot(AllData2ShowSh(:),groupsSh(:),'Whisker', Inf,'positions',(1:6)+0.35,'colors',[0.5 0.5 0.5],'Widths', 0.2)%,'PlotStyle','compact')

    pval_all=[Data.(FMdlNmFS).pval_acrossrule Data.(FMdlNmTS).pval_acrossrule];
    arrayfun(@(x) fp.AddDetailedSignificanceStar(x,pval_all(f,DataOrder(x)),'k',[],'SigStar_fontsize',10,'ThisYval',MeanData2Show(DataOrder(x))),1:6)

    % plot within rule as horizontal line
    fp.PlotHorizontalLine(Data.(FMdlNmAll).AllAnglesStruct_Mean.angd_WithinRule_GLMfit(f,2),[],'g','p_line_style','--')
    fp.PlotHorizontalLine(Data.(FMdlNmAll).AllAnglesStruct_Mean.angd_WithinRule_GLMfit(f,3),[],'b','p_line_style','--')
   
    xlabel('Task Pair')
    ylabel('Angle')
    xticks([1:6]);
    xticklabels({'1-2 From Swc', '1-2 To Swc','1-3 From Swc', '1-3 To Swc','3-2 From Swc', '3-2 To Swc'});
    xtickangle(45)
    title(['Angle bet tasks for ' FactorNames{f}]);
    axis square
    ylim([0 95])
    n=n+1;
end

h=gcf;
