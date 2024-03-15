function [params, fit_curve, slope_sign] = fitSigmoidAutoSlope(x, y)
    % Fit a sigmoid function to the data (x, y) and automatically determine slope sign

    % Sigmoid function model
    sigmoid_model = @(params, x) params(1) + (params(2) - params(1)) ./ (1 + exp(params(4) * (x - params(3))));

    % Initial guess for parameters
    initial_params = [min(y), max(y), median(x), 1];

    % Perform nonlinear least squares fit
    params_pos = lsqcurvefit(sigmoid_model, initial_params, x, y);

    % Change the sign of the slope
    params_neg = params_pos;
    params_neg(4) = -params_neg(4);

    % Compute error for both positive and negative slopes
    error_pos = sum((sigmoid_model(params_pos, x) - y).^2);
    error_neg = sum((sigmoid_model(params_neg, x) - y).^2);

    % Choose the parameters with lower error
    if error_pos < error_neg
        params = params_pos;
        slope_sign = 1; % Positive slope
    else
        params = params_neg;
        slope_sign = -1; % Negative slope
    end

    % Generate fitted curve using the obtained parameters
    fit_curve = sigmoid_model(params, x);

    % Plot the data and the fitted curve (optional)
    figure;
    plot(x, y, 'o', 'DisplayName', 'Data');
    hold on;
    plot(x, fit_curve, 'r-', 'DisplayName', 'Fit');
    xlabel('X-axis');
    ylabel('Y-axis');
    legend('show');
    title(['Sigmoid Fit (Slope Sign: ' num2str(slope_sign) ')']);

end
