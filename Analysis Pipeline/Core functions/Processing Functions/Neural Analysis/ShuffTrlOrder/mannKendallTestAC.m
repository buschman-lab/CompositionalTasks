function [tau, pValue, z,h] = mannKendallTestAC(x, alpha)
    % MANNKENDALLTESTAC - Mann-Kendall test with autocorrelation adjustment
    %
    % Inputs:
    %   x     - Time series data (vector)
    %   alpha - Significance level (e.g., 0.05 for a 95% confidence level)
    %
    % Outputs:
    %   h      - Test decision (1 if trend is significant, 0 otherwise)
    %   tau    - Kendall's tau (correlation coefficient)
    %   pValue - P-value for the test
    %   z      - Standardized test statistic (z-score)
    
    % Check inputs
    if nargin < 2
        alpha = 0.05; % Default significance level
    end
    x = x(:); % Ensure input is a column vector
    n = length(x);
    if n < 2
        error('Time series must have at least two observations.');
    end

    % Compute the Mann-Kendall statistic (S)
    S = 0;
    for i = 1:n-1
        for j = i+1:n
            S = S + sign(x(j) - x(i));
        end
    end

    % Compute Kendall's tau
    tau = S / (n * (n - 1) / 2);

    % Adjust for autocorrelation (Effective Sample Size, ESS)
    lagMax = n - 1;
    acf = autocorr(x, 'NumLags', lagMax); % Compute autocorrelation function
    % Use only non-zero lags for the autocorrelation adjustment
    rhoSum = sum(acf(2:end) .* (1 - (1:lagMax)' / n)); % Weighted autocorrelations
    ESS = n / (1 + (2 * rhoSum)); % Effective Sample Size (single scalar)

    % Compute the variance of S adjusted for ESS
    uniqueVals = unique(x);
    ties = histcounts(x, uniqueVals); % Count ties
    tieCorrection = sum(ties .* (ties - 1) .* (2*ties + 5));
    varS = (ESS * (ESS - 1) * (2*ESS + 5) - tieCorrection) / 18;

    % Compute the z-score
    if S > 0
        z = (S - 1) / sqrt(varS);
    elseif S < 0
        z = (S + 1) / sqrt(varS);
    else
        z = 0;
    end

    % Compute the two-tailed p-value
    pValue = 2 * (1 - normcdf(abs(z)));

    % Decision rule
    h = pValue < alpha;

%     % Display results
%     fprintf('Mann-Kendall Test with Autocorrelation Adjustment:\n');
%     fprintf('--------------------------------------------------\n');
%     fprintf('Effective Sample Size (ESS): %.2f\n', ESS);
%     fprintf('Test Statistic (S): %d\n', S);
%     fprintf('Variance (Var(S)): %.4f\n', varS);
%     fprintf('Z-Score: %.4f\n', z);
%     fprintf('P-Value: %.4f\n', pValue);
%     fprintf('Kendall''s Tau: %.4f\n', tau);
%     fprintf('Trend is significant: %d (1 = Yes, 0 = No)\n', h);
end
