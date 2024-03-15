% common utility codes go here so we don't look for them 

%% plotting 
CoreMotifsFig=cell(1,3);
[CoreMotifsFig(1:3)]=obj.FigParams.RenderFigure(3,[0 0 0.8 0.8]);
           
varargout(1:4)=obj.FigParams.RenderFigure(4,[0 0 0.8 0.8]); % create figures
[hp]=FigParams.Plot(obj,X,Y,Col,Xlbl,Ylbl,Title,varargin);
[ClustExmplFig{:}]=obj.PlotMotifClusts(Motifs_Aligned,TemplateMotif,CoreMotifs,clust_ind_Aligned,TSNEMotifsDTW,PCAMotifs);
DendroFig=cell(1);
[DendroFig{1}]=obj.PlotDendrogramwithMotifs(Z,cutoff,clust_ind,CoreMotifs);
[hp]=BarPlot(obj,X,Y,Col,Xlbl,Ylbl,Title,varargin) % plots whatever
hp=HistogramPlot(obj,X,edges,Col,Xlbl,Ylbl,Title,varargin) % plots histrogram

%legend 
 % Plot AUC and Accuracy 
yyaxis left
hp{1}=obj.PlotSVMPerfMetric(ClassifierResults,ClassifierOpts,PerfMetrics{1},Sp,Cond,'ThisColor',Cols(1,:),'MeanStdPlotType',3);
hp{2}=obj.PlotSVMPerfMetric(ClassifierResults,ClassifierOpts,PerfMetrics{2},Sp,Cond,'ThisColor',Cols(2,:),'MeanStdPlotType',3);
yyaxis right
hp{3}=obj.PlotSVMPerfMetric(ClassifierResults,ClassifierOpts,PerfMetrics{3},Sp,Cond,'ThisColor',Cols(3,:),'MeanStdPlotType',3);            

legend(PerfMetrics(~cellfun(@isempty,hp)),'Location','Best');
legend('boxoff');

% properties 'p_line_width', 'p_line_style','p_marker_size'
 [hp,Y,YSTD]=PlotMeanStd(obj,X,Y,YSTD,Xlbl,Ylbl,Col,Shaded,Title,varargin)
 obj.FigParams.PlotMeanStd(opts.TimSeqAvg,AvgHChs.AvgTim,[],'Time to Stim',AnalysisFields{f},FieldCOL(f,:),1,'')                            
 helperCWTTimeFreqPlot(Mean_NewMotif_padd,1:size(Mean_NewMotif_padd,2),1:size(Mean_NewMotif_padd,1),'justplot1',[''],'Time','Frequency(Hz)',0)

% axis label
xticklabels(arrayfun(@(x) num2str(XtickLbl(x),3),1:3:length(XtickLbl),'UniformOutput',0));

%Saving plots 
 AnalysisOpts.CurrentCh='';
[~,~,AnalysisData.SelectivityFig]=obj.ManData.GetFileName([],'SelectivityTargOn','SaveInResults',1,'WantedDate','ALL');
obj.FigParams.SaveFigSeries([],AnalysisData.SelectivityFig,[Figs])
          
%Varibales
obj.ManData.IsVarExistinFile(XCorrFileName,'DistMat') & ~obj.OverWrite
%Manipulate data
out=ReshapeCell2Mat(obj,data,Dim,varargin) 
out=BinData(~,binsiz,nbin,str,stp)
SaveVar(obj,AnalysisPathName,VarVal,VarName,ExtraTxt,varargin) % saves a specific variable into file and folders

% cluster funcs
ProcessingStreamLine(0, 0, [], 6, 1, 'MotifAnalysis.Phenograph_K',5);

% clustring 
Dist=1-AnalysisData.MotifSimilarityMat;
yOut = squareform(Dist,'tovector');

X=TSNEMotifs;
idx=dbscan(yOut,0.001,30,'Distance','precomputed');
%idx=kmeans(X,20);
ClustInd=idx;
Nclust=length(unique(ClustInd));
UniqueClusts=unique(ClustInd);
figure;hold on;Col=distinguishable_colors(Nclust);arrayfun(@(x) plot(X(ClustInd==x,1),X(ClustInd==x,2),'.','color',Col(UniqueClusts==x,:)),unique(ClustInd))
arrayfun(@(x) text(2*x,50,num2str(x),'color',Col(x,:)),1:Nclust)




%How do I vary color along a 2D line?
%https://www.mathworks.com/matlabcentral/answers/5042-how-do-i-vary-color-along-a-2d-line
x = 0:.05:2*pi;
y = sin(x);
z = zeros(size(x));
col = x;  % This is the color, vary with x in this case.
surface([x;x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);