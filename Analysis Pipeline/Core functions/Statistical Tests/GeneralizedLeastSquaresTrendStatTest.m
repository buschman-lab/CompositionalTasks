function [beta_gls,p_values]=GeneralizedLeastSquaresTrendStatTest(t,y)
%t=(1:16)';    
% % Time points
% y = [2.1, 2.4, 2.7, 3.0, 3.1, ...  % Replace with your data
%      3.3, 3.7, 4.0, 4.1, 4.4, ...
%      4.6, 4.8, 5.0, 5.3, 5.5, 5.7]';

% Step 1: Fit OLS for Comparison
X = [ones(size(t)), t];            % Design matrix (intercept and time)
beta_ols = (X' * X) \ (X' * y);    % OLS coefficients
residuals = y - X * beta_ols;      % OLS residuals

% Step 2: Estimate Autocorrelation
rho = corr(residuals(1:end-1), residuals(2:end)); % Lag-1 autocorrelation

% Step 3: Construct Covariance Matrix (AR(1) Example)
n = length(t);
Sigma = rho .^ abs((1:n)' - (1:n)); % Covariance matrix for AR(1)

% Step 4: Cholesky Decomposition
L = chol(Sigma, 'lower');          % Lower triangular matrix
invL = inv(L);                     % Inverse of L

% Step 5: Transform Data
y_trans = invL * y;                % Transform response
X_trans = invL * X;                % Transform predictors

% Step 6: GLS Regression
beta_gls = (X_trans' * X_trans) \ (X_trans' * y_trans); % GLS coefficients
y_pred_gls = X * beta_gls;         % Predicted values

% Step 7: Statistical Significance
sigma_gls = inv(X_trans' * X_trans);
std_errors = sqrt(diag(sigma_gls));
t_values = beta_gls ./ std_errors;
p_values = 2 * (1 - tcdf(abs(t_values), n - size(X, 2))); % Two-tailed test

% Results
% disp('GLS Coefficients:');
% disp(beta_gls);
% disp('P-Values:');
% disp(p_values);
% 
% % Plot
% figure;
% plot(t, y, 'o', t, y_pred_gls, '-r', 'LineWidth', 1.5);
% xlabel('Time'); ylabel('Observed Values');
% title('GLS Trend Detection');
% legend('Observed Data', 'GLS Trend Line');