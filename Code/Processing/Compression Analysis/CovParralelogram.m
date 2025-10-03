function CovParralelogram(data_x,data_y,cov_matrix,T,Col,i,PlotScatter)
global AnalysisOpts
% Calculate the eigenvalues and eigenvectors of the covariance matrix

% Calculate the eigenvalues and eigenvectors of the covariance matrix
[eigenvectors, eigenvalues] = eig(cov_matrix);

% Determine the angle of rotation for the parallelogram (in degrees)
angle = atan2(eigenvectors(2, 1), eigenvectors(1, 1)) * 180 / pi; % Convert to degrees

% Calculate the standard deviations along each axis
std_major = sqrt(max(diag(eigenvalues)));
std_minor = sqrt(min(diag(eigenvalues)));

% Define the scaling factors for the parallelogram sides
scaling_factor_major = 2 * std_minor; % Adjust the scaling as needed
scaling_factor_minor = 2 * std_major; % Adjust the scaling as needed

% Calculate the coordinates of the four corners of the parallelogram
center_x = mean(data_x);
center_y = mean(data_y);

% Define the vertices of the parallelogram
vertex1_x = center_x + scaling_factor_major * cosd(angle);
vertex1_y = center_y + scaling_factor_major * sind(angle);

vertex2_x = center_x - scaling_factor_minor * sind(angle);
vertex2_y = center_y + scaling_factor_minor * cosd(angle);

vertex3_x = center_x - scaling_factor_major * cosd(angle);
vertex3_y = center_y - scaling_factor_major * sind(angle);

vertex4_x = center_x + scaling_factor_minor * sind(angle);
vertex4_y = center_y - scaling_factor_minor * cosd(angle);

% Plot the parallelogram and data points
 plot([vertex1_x, vertex2_x, vertex3_x, vertex4_x, vertex1_x], ...
    [vertex1_y, vertex2_y, vertex3_y, vertex4_y, vertex1_y], 'color',Col(i,:)); % Plot the parallelogram in blue
hold on;
if PlotScatter
    scatter(data_x, data_y, 'ro');  % Scatter plot of data points in red
end
ylabel('Color encoding axis');
xlabel('Shape encdoing axis');
axis equal;
grid on;
title(sprintf('T:%0.2f',AnalysisOpts.Time(T)))
end
