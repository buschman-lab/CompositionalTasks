function PlotRMI4BetaResp(obj,BetaAxis1,BetaAxis2,Ts1,Ts2)
global AnalysisOpts


BetaAxis1Abs=cellfun(@abs,BetaAxis1,'UniformOutput',0);
BetaAxis2Abs=cellfun(@abs,BetaAxis2,'UniformOutput',0);

BetaAxisAbsAvg1=mean(obj.ManData.ReshapeCell2Mat(BetaAxis1Abs,3),3);
BetaAxisAbsAvg2=mean(obj.ManData.ReshapeCell2Mat(BetaAxis2Abs,3),3);
Time=AnalysisOpts.Time;
nTim=length(Time);
NNeu=size(BetaAxisAbsAvg1,2);
TimeInd1=Time>=0.3 & Time<=0.5;

%% plot
figure
% plot orthogonal representation for 
subplot(131)
hold on
AngBetEncodingAxes=obj.CalculateAngleBetEncodingAxes(BetaAxis1,BetaAxis2);

obj.FigParams.PlotMeanStd(Time,AngBetEncodingAxes',[],'Time','Angle(Deg)',d,...
    obj.MeanStdPlotType,'','AppendTitles',1,'NormalizebyMax',0,'WidthSmoothing',Navg,'SmoothingMethod','movmean','include_n',0,'LegendTxt',[Ts1 ' ' Ts2]);
 
% now plot selectivity of axis 
RMI=cell2mat(arrayfun(@(Neu) arrayfun(@(t) obj.ManData.CalResponseModulationIndex([BetaAxisAbsAvg1(t,Neu) BetaAxisAbsAvg2(t,Neu)]),...
    1:nTim)',1:NNeu,'UniformOutput',0));

% sort the RMI 
AvgRMI=mean(RMI(TimeInd1,:),1);
[~,SortRMIInd]=sort(AvgRMI);
Sp=subplot(132);

% obj.FigParams.Image(Time,1:NNeu,RMI(:,SortRMIInd)',...
%     {AnalysisOpts.Xlabel},{AnalysisOpts.Xlabel},'RMI','Seelctivityindex',Sp,'image_colormap',AnalysisOpts.ClassifierAccuracyColormap,...
%     'WidthSmoothing',obj.WidthSmoothing,'SmoothingMethod','movmean');

imagesc(Time,1:NNeu,RMI(:,SortRMIInd)')
xlabel(AnalysisOpts.Xlabel);ylabel('Neuron #');
colormap(AnalysisOpts.ClassifierAccuracyColormap);
colorbar

Sp=subplot(133);
BetaAxisAbsAvg1Tim=mean(BetaAxisAbsAvg1(TimeInd1,:),1);
BetaAxisAbsAvg2Tim=mean(BetaAxisAbsAvg2(TimeInd1,:),1);

hold on
obj.FigParams.PlotMeanStd(1:NNeu,BetaAxisAbsAvg1Tim(SortRMIInd),[],'Neuron #','RMI',4,...
    obj.MeanStdPlotType,'','AppendTitles',1,'NormalizebyMax',0,'WidthSmoothing',1,'SmoothingMethod','movmean','include_n',0,'LegendTxt','Beta 1','p_line_width',1);

obj.FigParams.PlotMeanStd(1:NNeu,BetaAxisAbsAvg2Tim(SortRMIInd),[],'Neuron #','RMI',5,...
    obj.MeanStdPlotType,'','AppendTitles',1,'NormalizebyMax',0,'WidthSmoothing',1,'SmoothingMethod','movmean','include_n',0,'LegendTxt','Beta 2','p_line_width',1);
 
obj.FigParams.PlotMeanStd(1:NNeu,AvgRMI(SortRMIInd),[],'Neuron #','RMI',d,...
    obj.MeanStdPlotType,'','AppendTitles',1,'NormalizebyMax',0,'WidthSmoothing',1,'SmoothingMethod','movmean','include_n',0,'LegendTxt','Avg RMI');
