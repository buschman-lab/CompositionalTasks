%Warping test
MeanPoints=[-1 -1;1 1;-1 1;1 -1];
X=cell2mat(arrayfun(@(x) 0.2*rand(1,100)+x,MeanPoints(:,1)','UniformOutput',0));
Y=cell2mat(arrayfun(@(y) 0.2*rand(1,100)+y,MeanPoints(:,2)','UniformOutput',0));
scatter(X,Y);axis square

VarX=var(X);
VarY=var(Y);
CovXY=cov(X,Y);

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
figure
Col=copper(16);MM=10;
subplot(221);hold on
arrayfun(@(x) plot(AnalysisOpts.Time,CovCol(:,x),'color',Col(x,:)),1:16);title('color')
subplot(222);hold on
arrayfun(@(x) plot(AnalysisOpts.Time,CovSh(:,x),'color',Col(x,:)),1:16);title('shape')
subplot(223);hold on
arrayfun(@(x) plot(AnalysisOpts.Time,CovColSh1(:,x),'color',Col(x,:)),1:16);title('colorshape')
subplot(224);hold on
arrayfun(@(x) plot(AnalysisOpts.Time,CovColSh2(:,x),'color',Col(x,:)),1:16);title('colorshape')




% Define your covariance matrix
i=16;
figure
for T=1:5:141
    cla
    data_x=obj.ManData.ReshapeCell2Mat(ColorAxisTime{i},1)';
    data_y=obj.ManData.ReshapeCell2Mat(ShapeAxisTime{i},1)';
    Cov=arrayfun(@(x) cov(data_x(:,x),data_y(:,x)),1:141,'UniformOutput',0);
    
    data_x=data_x(:,T);
    data_y=data_y(:,T);
    cov_matrix = Cov{T};
    
    % Calculate eigenvalues and eigenvectors
    [eigenvec, eigenval] = eig(cov_matrix);
    
    % Determine major and minor axes
    major_axis = sqrt(max(diag(eigenval)));
    minor_axis = sqrt(min(diag(eigenval)));
    
    % Calculate the angle of rotation for the ellipse
    angle = atan2(eigenvec(2, 1), eigenvec(1, 1));
    
    % Create an array of angles for the ellipse
    theta = linspace(0, 2 * pi, 100);
    
    % Calculate the points on the ellipse
    x = major_axis * cos(theta);
    y = minor_axis * sin(theta);
    
    % Rotate the points to match the orientation of the major axis
    x_rotated = x * cos(angle) - y * sin(angle);
    y_rotated = x * sin(angle) + y * cos(angle);
    
    % Calculate the center of the ellipse (mean values of your data)
    center_x = mean(data_x);
    center_y = mean(data_y);
    
    % Translate the points to the center
    x_final = center_x + x_rotated;
    y_final = center_y + y_rotated;
    
    % Plot the ellipse and data points
    plot(x_final, y_final, 'b');  % Plot the ellipse in blue
    hold on;
    scatter(data_x, data_y, 'ro');  % Scatter plot of data points in red
    xlabel('Variable X');
    ylabel('Variable Y');
    title('Confidence Ellipse');
    axis equal;
    grid on;
    title(sprintf('T:%0.2f',AnalysisOpts.Time(T)))
    pause
end


% Define your covariance matrix
i=16;
figure
for T=1:5:141
    cla
    data_x=obj.ManData.ReshapeCell2Mat(ColorAxisTime{i},1)';
    data_y=obj.ManData.ReshapeCell2Mat(ShapeAxisTime{i},1)';
    Cov=arrayfun(@(x) cov(data_x(:,x),data_y(:,x)),1:141,'UniformOutput',0);
    
    data_x=data_x(:,T);
    data_y=data_y(:,T);
    cov_matrix = Cov{T};
    % Calculate the eigenvalues and eigenvectors of the covariance matrix
    [eigenvectors, eigenvalues] = eig(cov_matrix);
    
    % Extract the major and minor eigenvalues
    major_eigenvalue = max(diag(eigenvalues));
    minor_eigenvalue = min(diag(eigenvalues));
    
    % Determine the angle of rotation for the rectangle
    angle = atan2(eigenvectors(2, 1), eigenvectors(1, 1));
    
    % Calculate the standard deviations along each axis
    std_major = sqrt(major_eigenvalue);
    std_minor = sqrt(minor_eigenvalue);
    
    % Calculate the half-lengths of the sides of the rectangle
    half_width = 2 * std_minor;
    half_height = 2 * std_major;
    
    % Calculate the center of the rectangle (mean values of your data)
    center_x = mean(data_x);
    center_y = mean(data_y);
    
    % Calculate the coordinates of the four corners of the rotated rectangle
    x1 = center_x + half_width * cos(angle) - half_height * sin(angle);
    y1 = center_y + half_width * sin(angle) + half_height * cos(angle);
    x2 = center_x - half_width * cos(angle) - half_height * sin(angle);
    y2 = center_y - half_width * sin(angle) + half_height * cos(angle);
    x3 = center_x - half_width * cos(angle) + half_height * sin(angle);
    y3 = center_y - half_width * sin(angle) - half_height * cos(angle);
    x4 = center_x + half_width * cos(angle) + half_height * sin(angle);
    y4 = center_y + half_width * sin(angle) - half_height * cos(angle);
    
    % Plot the rotated rectangle and data points
    plot([x1, x2, x3, x4, x1], [y1, y2, y3, y4, y1], 'b');  % Plot the rotated rectangle in blue
    hold on;
    scatter(data_x, data_y, 'ro');  % Scatter plot of data points in red
    xlabel('Variable X');
    ylabel('Variable Y');
    title('Rotated Rectangle Representation');
    axis equal;
    grid on;
    title(sprintf('T:%0.2f',AnalysisOpts.Time(T)))
    pause
end



%% plot superimpose rectangles during learning 
figure
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

Col=copper(16);
for i=[1:16 ]
    k=1;figure
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
%%

