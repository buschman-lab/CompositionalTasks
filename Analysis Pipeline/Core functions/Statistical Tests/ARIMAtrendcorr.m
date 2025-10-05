
function [rho,p,SensSlop,residuals]=ARIMAtrendcorr(y,t,model,CorrType)
% Step 1: Fit an ARIMA Model to Prewhiten
if isempty(model)
    model = arima('ARLags', 1, 'D', 0, 'Constant', 0); % AR(1) model
end
UseAutoCorr=0;

if UseAutoCorr==1   
    %Step 1: Compute Autocorrelations for Lags 1 to 5
    [acf, lags] = autocorr(y, 'NumLags', 5); % Compute ACF up to lag 5
    rho = acf(2:6);                         % Extract lag-1 to lag-5 autocorrelations

    % Step 2: Prewhiten the Data Using Lags 1 to 5
    n = length(y);
    y_prewhitened = nan(n, 1);               % Initialize prewhitened data
    for i = 6:n                              % Start from the 6th point to account for all lags
        y_prewhitened(i) = y(i) - sum(rho .* flip(y(i-5:i-1))'); % Adjust based on all 5 lags
    end
    residuals = y_prewhitened(6:end);   % Remove the first 5 points (incomplete adjustment)
    t_prewhitened = t(6:end);               % Adjust time vector
sens_slope

elseif UseAutoCorr==0
    fit = estimate(model, y,'display','off');                          % Fit model to data
    residuals = infer(fit, y);                         % Extract prewhitened residuals
    t_prewhitened=t;

elseif UseAutoCorr==2 % use a simple differencing model 
    residuals=diff(y);
    t_prewhitened=t(1:end-1);
end


% Step 2: Perform Spearman/Kendall Correlation on Prewhitened Data
[rho,p] = corr(residuals, t_prewhitened, 'Type', CorrType);      % Spearman correlation
SensSlop=sens_slope(t_prewhitened,residuals);
 
