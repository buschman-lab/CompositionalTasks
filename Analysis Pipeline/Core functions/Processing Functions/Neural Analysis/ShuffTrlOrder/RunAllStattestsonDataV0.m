% process all of the data of the classifier with different methods
function RunAllStattestsonData(ClassifierOpts,Metric,AvgDim2,Dim,WantedNStg)
global AnalysisOpts

warning off
NeuAna=NeuralAnalysisFuncsTemp;
ManData=ManipulateData;
fp=fig_params;
Prmt=@(x) permute(x,[1 3 2]); % permute

%ClassifierOpts=NeuAna.DefineClassifierTestOptions(TaskName,'CalShuffTrlOrder',1);
%% load data
% the data should be already loaded and passed along to this
%% preprocess the data
DimNum=[0 1 3];
[~,TargetFactorDim]=arrayfun(@(dim) NeuAna.getClassifierDimInfo(ClassifierOpts,dim),[0 1 3],'UniformOutput',0);
FileNameSpecs=@(x) [ClassifierOpts.Name '_' TargetFactorDim{x}];
ThisFileNameSpecs=FileNameSpecs(Dim);

PATH=AnalysisOpts.CodeTestDataPath;
if strcmp(Metric,'Accuracy')
    LoadFileName=[ThisFileNameSpecs '_D' num2str(Dim(1)) 'ShData'];
    % prepare data
    Data=load([PATH LoadFileName],[Metric '_Obsv'],[Metric '_Shuff']);
    Observed=ManData.SmoothData(Data.([Metric '_Obsv']){1},NeuAna.WidthSmoothing,'SmoothingMethod',NeuAna.SmoothingMethod,'DimSmoothing',2);
    Shuffle=ManData.SmoothData(Data.([Metric '_Shuff']){1},NeuAna.WidthSmoothing,'SmoothingMethod',NeuAna.SmoothingMethod,'DimSmoothing',2);
    MetricName=TargetFactorDim{Dim};
else
    LoadFileName='CompressionIndexShuffData.mat';   
    % prepare data
    Data=load([PATH LoadFileName],[Metric '_Obsv'],[Metric '_Shuff']);
    Observed=ManData.SmoothData(Data.([Metric '_Obsv']),NeuAna.WidthSmoothing,'SmoothingMethod',NeuAna.SmoothingMethod,'DimSmoothing',2);
    Shuffle=ManData.SmoothData(Data.([Metric '_Shuff']),NeuAna.WidthSmoothing,'SmoothingMethod',NeuAna.SmoothingMethod,'DimSmoothing',2);
  %  Shuffle=permute(Shuffle,[3 2 1]);
    ThisFileNameSpecs=  [ClassifierOpts.Name '_' Metric];
    MetricName=Metric;
end
fprintf('\nLoading data of classifier results from %s',LoadFileName)


if isempty(WantedNStg) | size(Observed,1)<WantedNStg
    TrlRng=[1:size(Observed,1)]';NStg=size(Observed,1);NStgTxt='';
else
    TrlRng=[1:WantedNStg]';NStg=WantedNStg;NStgTxt=['Nstg' num2str(NStg)];
end

NTim=size(Observed,2);
% avergae across dim2?
if AvgDim2
    Observed=ManData.SmoothData(Observed,NeuAna.WidthSmoothingDim2,'SmoothingMethod','movmean','DimSmoothing',1);
    Shuffle=ManData.SmoothData(Shuffle,NeuAna.WidthSmoothingDim2,'SmoothingMethod','movmean','DimSmoothing',1);
end
Shuffle=permute(Shuffle,[3 1 2]);
NrepShuff=size(Shuffle,1);
%% start plotting the data

% plot observed data and shuffle
RepSpace=@(x) strrep(x,'_',' ');
col=copper(NStg);Time=-0.6:0.01:0.55;
figure;
subplot(121);hold on;arrayfun(@(x) plot(Time,Observed(x,1:NTim),'color',col(x,:)),TrlRng');
xlabel('Time');title('Observed')
MeanShuff=squeeze(mean(Shuffle,1));
subplot(122);hold on;arrayfun(@(x) plot(Time,MeanShuff(x,1:NTim),'color',col(x,:)),TrlRng');
xlabel('Time');title('Mean Shuffle')
sgtitle(['Shuffle and observed for ' RepSpace(ThisFileNameSpecs)])


%CorrMethodSet={'Kendall','Spearman','3PW_p','3PW_slop','TheilSen','Regression','Sigmoid','ARIMAlag2','ARIMAlag4','ARIMAlag6'};

CorrMethodSet={'Kendall','Spearman','TheilSen','Regression','Sigmoid','ARIMAlag2','ARIMAlag4','ARIMAlag6'};

nMethods=length(CorrMethodSet);

for nThisCorrMethod=1:nMethods
    CorrMethod=CorrMethodSet{nThisCorrMethod};
    PlotData=1;
    fprintf('\n calculating trend correlation using %s ...',CorrMethod)
    %% loop on different methods and calculate statistics
    if strcmp(CorrMethod,'Kendall')
        CorrObsv=arrayfun(@(x) corr(TrlRng,Observed(1:NStg,x),'type','Kendall'),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(Shuffle(rep,1:NStg,x))','type','Kendall'),1:NTim)',1:NrepShuff,'UniformOutput',0))';
    elseif strcmp(CorrMethod,'GPTKendall')
        CorrObsv=arrayfun(@(x) mannKendallTestAC(Observed(1:NStg,x),0.05),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) mannKendallTestAC(squeeze(Shuffle(rep,1:NStg,x)),0.05),1:NTim)',1:NrepShuff,'UniformOutput',0))';
    elseif strcmp(CorrMethod,'TheilSen')
        CorrObsv=arrayfun(@(x) TheilSen([TrlRng,Observed(1:NStg,x)]),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) TheilSen([TrlRng,squeeze(Shuffle(rep,1:NStg,x))']),1:NTim)',1:NrepShuff,'UniformOutput',0))';
    elseif strcmp(CorrMethod,'Spearman')
        CorrObsv=arrayfun(@(x) corr(TrlRng,Observed(1:NStg,x),'type','Spearman'),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) corr(TrlRng,squeeze(Shuffle(rep,1:NStg,x))','type','Spearman'),1:NTim)',1:NrepShuff,'UniformOutput',0))';
    elseif strcmp(CorrMethod,'Regression')
        CorrObsv=arrayfun(@(x) LinearRegression(TrlRng,Observed(1:NStg,x)),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) LinearRegression(TrlRng,squeeze(Shuffle(rep,1:NStg,x))),1:NTim)',1:NrepShuff,'UniformOutput',0))';
    elseif strcmp(CorrMethod,'MMK')
        significance_value_tau = AnalysisOpts.Modified_MannKendall_significance_value_tau;
        significance_value_ac = AnalysisOpts.Modified_MannKendall_significance_value_ac;
        CorrObsv=arrayfun(@(x) Modified_MannKendall_test(TrlRng, Observed(1:NStg,x)', significance_value_tau, significance_value_ac),1:NTim);
        CorrShuff=cell2mat(arrayfun(@(rep) arrayfun(@(x) Modified_MannKendall_test(TrlRng, squeeze(Shuffle(rep,1:NStg,x))', significance_value_tau, significance_value_ac),1:NTim)',1:NrepShuff,'UniformOutput',0))';
    elseif strcmp(CorrMethod,'Sigmoid')
        [SlopeObsv,~,~,DifMinMaxObsv,NormDifMinMaxObsv,NormDifMinObsv]=arrayfun(@(x) fitAndAdjustSigmoidLsq(TrlRng,Observed(1:NStg,x),0,num2str(Time(x))),1:NTim,'UniformOutput',0);
        for rep=1:NrepShuff
            [SlopeShuff(rep,:),~,~,DifMinMaxShuff(rep,:),NormDifMinMaxShuff(rep,:),NormDifMinShuff(rep,:)]=arrayfun(@(x) fitAndAdjustSigmoidLsq(TrlRng,squeeze(Shuffle(rep,1:NStg,x))',0,num2str(Time(x))),1:NTim,'UniformOutput',0);
        end
        CorrObsv=cell2mat(DifMinMaxObsv);
        CorrShuff=cell2mat(DifMinMaxShuff);
        % CorrObsv(isoutlier(CorrObsv,1))=nan;
        % CorrShuff(isoutlier(CorrShuff,1))=nan;
    elseif contains(CorrMethod,'ARIMA')
        % then load ARIMA file
        if AvgDim2;Ext='AvgDim2';else;Ext='';end
        Lag=str2double(CorrMethod(end));
        ARIMALoadFileName=['ArimaTestShuffLag' num2str(Lag) Ext '_D' num2str(Dim) '_' Metric ThisFileNameSpecs];
        if exist([PATH ARIMALoadFileName '.mat'],'file')
            load([PATH ARIMALoadFileName],'ArShuffle','ArShuffle_p','ArObserved','ArObserved_p');
            CorrObsv=cell2mat(ArObserved);
            CorrShuff=permute(cell2mat(ArShuffle),[3 2 1]);

        else
            fprintf('\n ARIMA file %s does not exist',[ARIMALoadFileName]);
            PlotData=0;            
        end
    elseif strcmp(CorrMethod,'3PW_p') | strcmp(CorrMethod,'PW_p') |  strcmp(CorrMethod,'TFPW_Y_p')
        if strcmp(CorrMethod,'3PW_p')
            PW_method='3PW';
        elseif strcmp(CorrMethod,'PW_p')
            PW_method='PW';
        elseif strcmp(CorrMethod,'TFPW_Y_p')
            PW_method='TFPW_Y';
        end
        clear p3PW_Obsv Slop3PW_Obvs p3PW_Shuff Slop3PW_Shuff
        test_data=timetable(datetime(datevec([datenum([2001 01 01]):365:datenum([2000+NStg 01 01])])));
    %    test_data=timetable(datetime(datevec([datenum([2001 01 02]):1:datenum([2001 01+NStg 01])])));
        % calculate observed
        for t=1:NTim
            test_data.param=Observed(1:NStg,t);
            results=MK_tempAggr(test_data, 0.001, 'PW_method',PW_method,'alpha_MK',99,'alpha_ak',95, 'alpha_CL',95, 'alpha_Xhomo',99);
            p3PW_Obsv(1,t)=results.P;
            Slop3PW_Obvs(1,t)=results.slope;
        end
      %  figure;plot(Time,p3PW_Obsv)
        % calculate observed shuffle
        for rep=1:NrepShuff
            fprintf('\n Cl 3PW for rep %i',rep)
            for t=1:NTim
                test_data.param=squeeze(Shuffle(rep,1:NStg,t))';
                results=MK_tempAggr(test_data, 0.001, 'PW_method',PW_method,'alpha_MK',99,'alpha_ak',95, 'alpha_CL',95, 'alpha_Xhomo',90);
                p3PW_Shuff(rep,t)=results.P;
                Slop3PW_Shuff(rep,t)=results.slope;
            end
        end
        CorrObsv=p3PW_Obsv;
        CorrShuff=p3PW_Shuff;
    elseif strcmp(CorrMethod,'3PW_slop')  | strcmp(CorrMethod,'PW_slop') |  strcmp(CorrMethod,'TFPW_Y_slop')
        CorrObsv=Slop3PW_Obvs;
        CorrShuff=Slop3PW_Shuff;
    end

    if PlotData
        % plot with cluster correction
        if strcmp(TargetFactorDim{Dim},'Rule')
            TimeInd=Time>=-0.5 & Time<=0;
        else
            TimeInd=Time>=-0.15 & Time<=0.6;
        end

        AnalysisOpts.Time=Time(TimeInd);
        ManData.ClusterMassCorrection_permutationTwoTail(Prmt(CorrShuff(:,TimeInd)),Prmt(CorrObsv(1,TimeInd)),0.1,1,'ShowClustCorrectionPlot',1);
        ClusterCorrectionTJB(CorrObsv(1,TimeInd)',CorrShuff(:,TimeInd)',{[],[2,4,7]},'Time',Time(TimeInd),'ThresholdPercentage',20,'DoubleSided',1);
        ClusterCorrectionTJB(CorrObsv(1,TimeInd)',CorrShuff(:,TimeInd)',{[],[2,4,8]},'Time',Time(TimeInd),'ThresholdPercentage',20,'DoubleSided',0);

        sgtitle({RepSpace(ThisFileNameSpecs);CorrMethod})
        SourceFig=gcf;
        %  fp.FormatAllAxesFig(gcf)

        %% create a summery plot for each of the stat specs for each metric so we can look at all of them together
        SourceSpSet={[2,4,6],[2,4,4],[2,4,1],[2,4,7],[2,4,8]}; % referint to subplots inside the figure
        SourceSpSetName={'STD','uncorrpval','Cluster','ClustTJBDblsided','ClusterTJBSingsided'};
        if nThisCorrMethod==1
            h=fp.RenderFigure(length(SourceSpSet),[]);
            [Htot,TargetSps]=cellfun(@(x) fp.RenderSubplots([],[],x,nMethods),h,'UniformOutput',0);
        end

        % loop on figures and methods
        for tsp=1:length(SourceSpSet)
            figure(SourceFig);
            SourceSp=subplot(SourceSpSet{tsp}(1),SourceSpSet{tsp}(2),SourceSpSet{tsp}(3)); % choose the relevant Sp
            fp.CopySubplot(SourceSp,TargetSps{tsp}(nThisCorrMethod)); % copy into target
            fp.AppendFigTitles(TargetSps{tsp}(nThisCorrMethod),{CorrMethod},'AppendTitles',1)
            sgtitle(RepSpace(ThisFileNameSpecs))
        end
     
    end
end

if AvgDim2;Ext='AvgDim2';else;Ext='';end
    % save figures
    fprintf('\nSaving current figures ...')
    fp.SaveCurrentFigs(['AllStats_' NStgTxt Ext ThisFileNameSpecs],[AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep]);
   
    % save specific figures in one file 
    for kk=[1 2 4]
        fp.SaveFigSeries(['AllStats_' NStgTxt Ext MetricName '_' SourceSpSetName{kk}],[AnalysisOpts.ResultsSavePath filesep 'Stattests' filesep],Htot{kk});       
    end
end

