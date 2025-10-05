% Function to calculate Sen's slope
function sen_slope = sens_slope(x, y)
    % Input: x - independent variable (e.g., time)
    %        y - dependent variable (e.g., measurements)
    
    n = length(x); % Number of data points
    slopes = [];    % Initialize an empty array to store slopes
    
    % Loop over all pairs of points (i, j) where i < j
    for i = 1:n-1
        for j = i+1:n
            % Calculate the slope for the pair (i, j)
            slope = (y(j) - y(i)) / (x(j) - x(i));
            slopes = [slopes; slope]; % Append the slope to the array
        end
    end
    
    % Calculate the median of the slopes
    sen_slope = median(slopes);
end
