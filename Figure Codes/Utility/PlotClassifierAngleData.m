function h=PlotClassifierAngleData(FigName,DataPath)

FigParams=fig_params;
ManData=ManipulateData;

load([DataPath 'AnalysisOpts.mat']);
load([DataPath 'PFC_ClassifierAngleData.mat']);

h=figure;
Title=FigName;
t=1;


SubPlot={[1 2],[3],[4]};
nSp=3;

for d=1:2
    dd=ones(1,nSp);
    for s=1:nSp
        subplot(3,nSp,s);
        hold on;axis square;
        Ang3d=cell2mat(EncodeAng_deg(d).Ang);
        ndd=find(SubPlot{s}==d);
        if ndd
            [CircMean]=arrayfun(@(x) (ManData.MovingCircularMeanSTD(deg2rad(Ang3d(:,x)'),2,15)),1:size(Ang3d,2),'uniformoutput',0);
            CircMean=cellfun(@rad2deg,CircMean,'UniformOutput',0);


            arrayfun(@(x) FigParams.PlotMeanStd(ClassifierOpts.AnalysisOpts.Time,CircMean{x},[],'Time','Angle(Deg)',ndd,...
                1,Title,'AppendTitles',1,'NormalizebyMax',0,'WidthSmoothing',1,'SmoothingMethod','movmean','include_n',0),1:size(Ang3d,2));

            FigParams.plot_significance_level(StatTest(d).clusters{t},...
                StatTest(d).statsummery{t},ClassifierOpts.AnalysisOpts.Time,'auto','k',[],[],'WidthSmoothing',15,'SmoothingMethod','movmean');

            % add shuffle data
            [CircMeanShSmoothed]=arrayfun(@(x) (ManData.MovingCircularMeanSTD(deg2rad(EncodeAng_deg(d).Ang_Sh{x}'),2,15)),1:size(EncodeAng_deg(d).Ang_Sh,2),'uniformoutput',0);
            CircMeanShSmoothed=cellfun(@rad2deg,CircMeanShSmoothed,'UniformOutput',0);

            MeanCircMeanShSmoothed=ManData.circ_mean_deg(CircMeanShSmoothed{1},1);
            StdCircMeanShSmoothed=ManData.circ_std_deg(CircMeanShSmoothed{1},1);

            FigParams.PlotMeanStd(ClassifierOpts.AnalysisOpts.Time,MeanCircMeanShSmoothed,...
                StdCircMeanShSmoothed,[],[],[0.5 0.5 0.5],1,[],...
                'p_line_style','none','p_marker_size',15,'p_marker','.','STD_method','bootstrap','WidthSmoothing',1,'SmoothingMethod','movmean');

            va=axis;
            plot([ClassifierOpts.SpkCountPeriod(1,1) ClassifierOpts.SpkCountPeriod(1,1) ],[va(3) va(4)])
            plot([ClassifierOpts.SpkCountPeriod(1,2) ClassifierOpts.SpkCountPeriod(1,2) ],[va(3) va(4)])
            xlim([-0.2 0.55])

            % plot average angle
            Sp=subplot(3,nSp,nSp+s);
            FigParams.BarPlot(ndd,EncodeAng_deg(d).AngSpkCnt,ndd,'Trial','Angle(Deg)','',...
                'LegendTxt','AvgAng','IsthisAxisTime',0,'include_n',0);

            FigParams.AddDetailedSignificanceStar(ndd,StatTest(d).pval_AngSpkCnt_Sh,ndd,Sp,'SigStar_fontsize',20);

            % plot bar for the shuffle angels
            MeanSpkSh=ManData.circ_mean_deg(EncodeAng_deg(d).AngSpkCnt_Sh{t},2);
            StdSpkSh=ManData.circ_std_deg(EncodeAng_deg(d).AngSpkCnt_Sh{t},2);
            FigParams.PlotMeanStd(ndd,MeanSpkSh,...
                StdSpkSh,[],[],[0.5 0.5 0.5],0,[],...
                'p_line_style','none','p_marker_size',15,'p_marker','.','STD_method','bootstrap');
            ylim([50 100])


            dd(s)=dd(s)+1;
        end
    end
end