GetFig
SM=@(x,w) smoothdata(x,2,'movmean',10);
Col={'r','g','g','r'};
width=[3 1 3 1];
subplot(221)
hold on;arrayfun(@(x) plot(AnalysisOpts.Time,SM(mean(ColorAxisTime{1}{x},1)),'Color',Col{x},'LineWidth',width(x)),1:4)
title('Color Encoding Axis Trn2/Tst3 - Learning step 1')
subplot(222)
hold on;arrayfun(@(x) plot(AnalysisOpts.Time,SM(mean(ShapeAxisTime{1}{x},1)),'Color',Col{x},'LineWidth',width(x)),1:4)
title('Shape Encoding Axis Trn1/Tst3 - Learning step 1')

subplot(223)
hold on;arrayfun(@(x) plot(AnalysisOpts.Time,SM(mean(ColorAxisTime{16}{x},1)),'Color',Col{x},'LineWidth',width(x)),1:4)
title('Color Encoding Axis Trn2/Tst3 - Learning step last')

subplot(224)
hold on;arrayfun(@(x) plot(AnalysisOpts.Time,SM(mean(ShapeAxisTime{16}{x},1)),'Color',Col{x},'LineWidth',width(x)),1:4)
hold on;arrayfun(@(x) plot(AnalysisOpts.Time,SM(mean(ShapeAxisTime{1}{x},1)),'.-','Color',Col{x},'LineWidth',width(x)),1:4)

title('Shape Encoding Axis Trn1/Tst3 - Learning step last')

for i=1:4;subplot(2,2,i);legend({'RB','GT','GB','RT'});end

figure
TimeInd=AnalysisOpts.Time>=-0.1 & AnalysisOpts.Time<=0.6;
for i=1:10
    Sp=subplot(2,5,i);
    hold on;
    arrayfun(@(x) plot(AnalysisOpts.Time(TimeInd),SM(EncodingDist.([EncodingVars{i} 'Avg'])(x,TimeInd)),'color',Cop(x,:),'linewidth',2),1:16);
    title(EncodingVars{i})   
    axis square
    axis tight
    obj.FigParams.PlotVerticalLine(0,Sp,'k')
   % obj.FigParams.PlotHorizontalLine(0.5,Sp,'k')
end
figure
colormap(copper(16))
for i=1:3
    Sp=subplot(1,3,i);
    hold on;
    arrayfun(@(x) plot(AnalysisOpts.Time(TimeInd),SM(log(CompressionEncoding.TrlAvg.(CompIndexVars{i})(x,TimeInd))),'color',Cop(x,:),'linewidth',2),1:16);
    colorbar
    title(sprintf('Compression Index %s',CompIndexVars{i}))  
    axis square
    axis tight
    obj.FigParams.PlotVerticalLine(0,Sp,'k')
end

Sp=GetFig;
ObjTxt={'RB','GT','GB','RT'};
Col=copper(nTrialRange);
for TrlRng=[1:nTrialRange]%floor(nTrialRange/2)
    sp=subplot(4,4,TrlRng);
    hold on
    arrayfun(@(ThisObj) text(ColorAxis{TrlRng}(ThisObj),ShapeAxis{TrlRng}(ThisObj),ObjTxt{ThisObj},...
        'Color',Col(TrlRng,:),'FontSize',14),1:4)
    xlim([0.2 0.85]);ylim([0.2 0.85])
    axis square
    xlabel('Color Enc Axis');ylabel('Shape Enc Axis');
    obj.FigParams.PlotVerticalLine(0.5,sp,[0.5 0.5 0.5],'p_line_style','--')
    obj.FigParams.PlotHorizontalLine(0.5,sp,[0.5 0.5 0.5],'p_line_style','--')
    title(sprintf('TrialRange:%i',TrialRange(TrlRng)));
end
xlabel('Color Encoding Axis');
ylabel('Shape Encoding Axis');
