function [slope,intercept]=LinearRegression(time, measure)

coefficients = polyfit(time, measure, 1);
slope = coefficients(1);
intercept = coefficients(2);

end
