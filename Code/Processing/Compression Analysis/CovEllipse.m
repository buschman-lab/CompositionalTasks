function CovEllipse(data_x,data_y,cov_matrix,T,Col,i)
global AnalysisOpts 

% Calculate the eigenvalues and eigenvectors of the covariance matrix
[eigenvectors, eigenvalues] = eig(cov_matrix);

% Extract the major and minor eigenvalues
major_eigenvalue = max(diag(eigenvalues));
minor_eigenvalue = min(diag(eigenvalues));

% Determine the angle of rotation for the ellipse (in degrees)
angle = atan2(eigenvectors(2, 1), eigenvectors(1, 1)) * 180 / pi; % Convert to degrees

% Calculate the standard deviations along each axis
std_major = sqrt(major_eigenvalue);
std_minor = sqrt(minor_eigenvalue);

% Calculate the semi-axes lengths for the rotated ellipse
semi_major = 2 * std_major; % Double the standard deviation for semi-major axis
semi_minor = 2 * std_minor; % Double the standard deviation for semi-minor axis

% Calculate the center of the ellipse (mean values of your data)
center_x = mean(data_x);
center_y = mean(data_y);

% Generate points on the rotated ellipse
t = linspace(0, 2 * pi, 100);  % 100 points around the ellipse
x = center_x + semi_major * cosd(angle) * cos(t) - semi_minor * sind(angle) * sin(t);
y = center_y + semi_major * sind(angle) * cos(t) + semi_minor * cosd(angle) * sin(t);

% Plot the rotated ellipse and data points
plot(x, y, 'color',Col(i,:));  % Plot the rotated ellipse in blue
hold on;
scatter(data_x, data_y, 'ro');  % Scatter plot of data points in red
xlabel('Color encoding axis');
ylabel('Shape encdoing axis');
axis square;
grid on;
end
