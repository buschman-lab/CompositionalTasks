function [params, fit_curve] = fitSigmoid(x, y)
    % Fit a sigmoid function to the data (x, y)

    % Sigmoid function model
    sigmoid_model = @(params, x) params(1) + (params(2) - params(1)) ./ (1 + exp(-(x - params(3)) / params(4)));

    % Initial guess for parameters
    initial_params = [min(y), max(y), median(x), 1];

    % Perform nonlinear least squares fit
    params = lsqcurvefit(sigmoid_model, initial_params, x, y);

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
    title('Sigmoid Fit');

end
