function CovParralelogramMaxFit(data_x,data_y,cov_matrix,T,Col,i,PlotScatter)
global AnalysisOpts
% Calculate the eigenvalues and eigenvectors of the covariance matrix

% Calculate the eigenvalues and eigenvectors of the covariance matrix
[eigenvectors, eigenvalues] = eig(cov_matrix);

% Find the major and minor eigenvalues and their corresponding eigenvectors
major_eigenvalue = max(diag(eigenvalues));
minor_eigenvalue = min(diag(eigenvalues));

major_eigenvector = eigenvectors(:, find(diag(eigenvalues) == major_eigenvalue));
minor_eigenvector = eigenvectors(:, find(diag(eigenvalues) == minor_eigenvalue));

% Calculate the angle of rotation for the parallelogram (in degrees)
angle = atan2(major_eigenvector(2), major_eigenvector(1)) * 180 / pi;

% Calculate the half-lengths along major and minor axes
half_length_major = sqrt(major_eigenvalue);
half_length_minor = sqrt(minor_eigenvalue);

% Calculate the coordinates of the four corners of the parallelogram
center_x = mean(data_x);
center_y = mean(data_y);

vertex1_x = center_x + half_length_major * cosd(angle) + half_length_minor * sind(angle);
vertex1_y = center_y + half_length_major * sind(angle) - half_length_minor * cosd(angle);

vertex2_x = center_x - half_length_major * cosd(angle) + half_length_minor * sind(angle);
vertex2_y = center_y - half_length_major * sind(angle) - half_length_minor * cosd(angle);

vertex3_x = center_x - half_length_major * cosd(angle) - half_length_minor * sind(angle);
vertex3_y = center_y - half_length_major * sind(angle) + half_length_minor * cosd(angle);

vertex4_x = center_x + half_length_major * cosd(angle) - half_length_minor * sind(angle);
vertex4_y = center_y + half_length_major * sind(angle) + half_length_minor * cosd(angle);

% Plot the minimum bounding parallelogram and data points
 plot([vertex1_x, vertex2_x, vertex3_x, vertex4_x, vertex1_x], ...
      [vertex1_y, vertex2_y, vertex3_y, vertex4_y, vertex1_y], 'color',Col(i,:)); % Plot the parallelogram in blue
% data = [data_x, data_y];
% k=convhull(data);
% plot(data(k, 1), data(k, 2), 'color',Col(i,:)); % Plot the convex hull in blue
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
