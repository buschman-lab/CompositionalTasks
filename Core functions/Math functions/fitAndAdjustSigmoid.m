% The sigmoid function typically has four parameters. The parameters define the shape and position of the sigmoid curve. The general form of the sigmoid function is given by:
% 
% \[ f(x) = a + \frac{b - a}{1 + e^{-(x - c)/d}} \]
% 
% Here are the meanings of the parameters:
% 
% - \(a\): The lower asymptote, which is the minimum value that the function approaches as \(x\) goes to negative infinity.
% 
% - \(b\): The upper asymptote, which is the maximum value that the function approaches as \(x\) goes to positive infinity.
% 
% - \(c\): The x-value (horizontal shift) of the sigmoid's midpoint. It represents the inflection point of the curve.
% 
% - \(d\): The slope or steepness of the sigmoid curve. It influences how quickly the sigmoid transitions from the lower asymptote to the upper asymptote.
% 
% When fitting a sigmoid function to data, these parameters are estimated to best match the shape of the sigmoid to the observed data. The process of estimating these parameters is often done using optimization techniques, such as least squares fitting.
% 
% In a MATLAB implementation, these parameters are often represented as a vector:
% 
% \[ \text{params} = [a, b, c, d] \]
% 
% When using tools like `lsqcurvefit` in MATLAB, you provide an initial guess for these parameters, and the function iteratively adjusts them to minimize the difference between the predicted sigmoid curve and the observed data.
% 
% It's worth noting that in some contexts, variations of the sigmoid function might have additional parameters or constraints based on the specific application. Adjustments may be made to handle cases where the sigmoid needs to be shifted, scaled, or constrained differently.
function [slop,asymptoteDiff,params] = fitAndAdjustSigmoid(x, y,PlotFlag)
    % Fit a sigmoid function to the data (x, y) and adjust for ascending or descending trend
    

    % Sigmoid function model
    sigmoid_model = @(params, x) params(1) + (params(2) - params(1)) ./ (1 + exp(-(x - params(3)) / params(4)));

    % Initial guess for parameters
    initial_params = [min(y), max(y), median(x), 1];

   % opts = optimset('Display', 'off');
    opts = statset('nlinfit');
 %   opts.RobustWgtFun = 'bisquare';
    opts.Display='Off';
   % modelstr='y~ params1 + (params2 - params1) ./ (1 + exp(-(x - params3) / params4))';
    % Perform nonlinear least squares fit
    %params = lsqcurvefit(sigmoid_model, initial_params, x, y,[],[],options);
    params = fitnlm(x, y, sigmoid_model, initial_params, 'Options',opts);
    params=cell2mat(table2cell(params.Coefficients));

    % Generate fitted curve using the obtained parameters
 %   fit_curve = sigmoid_model(params, x);

    % Check the trend (ascending or descending)
   slope = params(2,1) - params(1,1);  % Slope of the sigmoid curve
    if slope < 0
        % If the slope is negative, adjust the sigmoid for a descending trend
        descending_sigmoid_model = @(params, x) params(1) + (params(2) - params(1)) ./ (1 + exp((x - params(3)) / params(4)));
     %   params = lsqcurvefit(descending_sigmoid_model, params, x, y,[],[],options);
        params = fitnlm(x, y, descending_sigmoid_model, initial_params, 'Options',opts);
  %      fit_curve = descending_sigmoid_model(params, x);
        params=cell2mat(table2cell(params.Coefficients));

        params(4,1)=-1*params(4,1);
    end

    if PlotFlag
        % Plot the data and the fitted curve (optional)
        figure;
        plot(x, y, 'o', 'DisplayName', 'Data');
        hold on;
        plot(x, fit_curve, 'r-', 'DisplayName', 'Adjusted Fit');
        xlabel('X-axis');
        ylabel('Y-axis');
        legend('show');
        title('Adjusted Sigmoid Fit');
    end
    slop=params(4,1);
    asymptoteDiff=params(1,1)-params(2,1);
end

%% use the code below to test this function 
% % Example data
% x = linspace(1, 10, 100);
% y = 2 - 3 ./ (1 + exp(-(x - 5) / 1.5)) + 0.2 * randn(size(x));
% 
% % Call the fitAndAdjustSigmoid function
% [params, fit_curve] = fitAndAdjustSigmoid(x, y);
% 
% % Display the fitted parameters
% disp('Fitted Parameters:');
% disp(params);

