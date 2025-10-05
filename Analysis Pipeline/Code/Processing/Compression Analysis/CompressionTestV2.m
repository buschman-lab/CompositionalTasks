%% plot Xcorr values
for i=1:16
    Color=obj.ManData.ReshapeCell2Mat(ColorAxisTime{i},1)';
    Shape=obj.ManData.ReshapeCell2Mat(ShapeAxisTime{i},1)';
    Color=obj.ManData.SmoothData(Color,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
    Shape=obj.ManData.SmoothData(Shape,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);

    % calculate covariance at each time point 
    Cov=arrayfun(@(x) cov(Color(:,x),Shape(:,x)),1:141,'UniformOutput',0);
    CovCol(:,i)=cellfun(@(x) x(1,1),Cov);
    CovSh(:,i)=cellfun(@(x) x(2,2),Cov);
    CovColSh1(:,i)=cellfun(@(x) x(1,2),Cov);  
    CovColSh2(:,i)=cellfun(@(x) x(2,1),Cov);     
end
h1=obj.FigParams.RenderFigure(1,[]);
Col=copper(16);MM=10; 
subplot(221);hold on
arrayfun(@(x) plot(AnalysisOpts.Time,CovCol(:,x),'color',Col(x,:)),1:16);title('color');xlabel('Time from Sample On')
subplot(222);hold on
arrayfun(@(x) plot(AnalysisOpts.Time,CovSh(:,x),'color',Col(x,:)),1:16);title('shape');xlabel('Time from Sample On')
subplot(223);hold on
arrayfun(@(x) plot(AnalysisOpts.Time,CovColSh1(:,x),'color',Col(x,:)),1:16);title('colorshape');xlabel('Time from Sample On')
subplot(224);hold on
arrayfun(@(x) plot(AnalysisOpts.Time,CovColSh2(:,x),'color',Col(x,:)),1:16);title('colorshape');xlabel('Time from Sample On')



%% plot superimpose rectangles during learning 
h2=obj.FigParams.RenderFigure(1,[]);

Col=copper(16);
k=1;
for T=[0:0.05:0.4]
    Ta=find(AnalysisOpts.Time>=T,1,'first');
        subplot(3,3,k)
    for i=[1:16]
        
        data_x=obj.ManData.ReshapeCell2Mat(ColorAxisTime{i},1)';
        data_y=obj.ManData.ReshapeCell2Mat(ShapeAxisTime{i},1)';
        
        data_x=obj.ManData.SmoothData(data_x,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
        data_y=obj.ManData.SmoothData(data_y,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);

        Cov=arrayfun(@(x) cov(data_x(:,x),data_y(:,x)),1:141,'UniformOutput',0);
        
        data_x=data_x(:,Ta);
        data_y=data_y(:,Ta);
        cov_matrix = Cov{Ta};
        CovRectangle(data_x,data_y,cov_matrix,Ta,Col,i,0)
        % CovEllipse(data_x,data_y,cov_matrix,T,Col,i)
       % CovParralelogram(data_x,data_y,cov_matrix,Ta,Col,i,0)
    end
    k=k+1;
end

%% plot encoding for each object 
h3=obj.FigParams.RenderFigure(16,[]);
Col=copper(16);
for i=[1:16 ]
    k=1;figure(h3{i})
    for T=[0:0.05:0.4]
        Ta=find(AnalysisOpts.Time>=T,1,'first');
                
        subplot(3,3,k)
        data_x=obj.ManData.ReshapeCell2Mat(ColorAxisTime{i},1)';
        data_y=obj.ManData.ReshapeCell2Mat(ShapeAxisTime{i},1)';
        
        data_x=obj.ManData.SmoothData(data_x,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
        data_y=obj.ManData.SmoothData(data_y,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);

        Cov=arrayfun(@(x) cov(data_x(:,x),data_y(:,x)),1:141,'UniformOutput',0);
        
        data_x=data_x(:,Ta);
        data_y=data_y(:,Ta);
        cov_matrix = Cov{Ta};
        CovRectangle(data_x,data_y,cov_matrix,Ta,Col,i,2)
        %CovParralelogram(data_x,data_y,cov_matrix,Ta,Col,i,1)
        k=k+1;
    end
end

%% plot encoding for each object with distribution of encodings
Col=copper(16);
h4=obj.FigParams.RenderFigure(16,[]);

for i=[1:16 ]
    k=1;figure(h4{i})
    for T=[0:0.05:0.4]
        Ta=find(AnalysisOpts.Time>=T,1,'first');
                
        subplot(3,3,k)
        data_x=obj.ManData.ReshapeCell2Mat(ColorAxisTime{i},1)';
        data_y=obj.ManData.ReshapeCell2Mat(ShapeAxisTime{i},1)';
        
        data_x=obj.ManData.SmoothData(data_x,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);
        data_y=obj.ManData.SmoothData(data_y,obj.WidthSmoothing,'SmoothingMethod',obj.SmoothingMethod,'DimSmoothing',2);

        Cov=arrayfun(@(x) cov(data_x(:,x),data_y(:,x)),1:141,'UniformOutput',0);
        
        data_x=data_x(:,Ta);
        data_y=data_y(:,Ta);
        cov_matrix = Cov{Ta};
        CovRectangle(data_x,data_y,cov_matrix,Ta,Col,i,2)
        %CovParralelogram(data_x,data_y,cov_matrix,Ta,Col,i,1)
        k=k+1;
    end
end

h=[h1 h2 h3 h4];