function CovRectangle(data_x,data_y,cov_matrix,T,Col,i,PlotScatter)
global AnalysisOpts
nTrls=length(data_x)/4;
objID=cell2mat(arrayfun(@(x) x*ones(1,nTrls),1:4,'UniformOutput',0));
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
plot([x1, x2, x3, x4, x1], [y1, y2, y3, y4, y1], 'color',Col(i,:));  % Plot the rotated rectangle in blue
hold on;
Marker={'B','T','B','T'};%{'o','+','o','+'};%[red bunny, green tee, green bunny, red tee]
%AnalysisOpts.QuadrantColors=[AnalysisOpts.QuadrantColors(1:2,:);AnalysisOpts.QuadrantColors(1:2,:)];
if PlotScatter==1
    %arrayfun(@(x) scatter(data_x(objID==x), data_y(objID==x), 36,AnalysisOpts.QuadrantColors(x,:),Marker{x}),1:4);  % Scatter plot of data points in red
    arrayfun(@(x) text(data_x(objID==x), data_y(objID==x), Marker{x},'color',AnalysisOpts.QuadrantColors(x,:),...
        'FontSize',9),4:-1:1);  % Scatter plot of data points in red   
elseif PlotScatter==2 % plot center as the mean
    arrayfun(@(x) text(mean(data_x(objID==x)), mean(data_y(objID==x)), Marker{x},'color',AnalysisOpts.QuadrantColors(x,:),...
        'FontSize',20),1:4);  % Scatter plot of data points in red   
end
xlabel('Color encoding axis');
ylabel('Shape encdoing axis');
axis equal;
xlim([-0.2 1.2])
ylim([-0.2 1.2])
grid on;
title(sprintf('T:%0.2fsecs-Train:%i',AnalysisOpts.Time(T),i))
end
