function results = analysis_func(x, y)
% This function takes two vectors of data, and performs some simple
% analysis of them.
% It displays a plot of the data with some features highlighted, and
% returns a struct with the following fields for the main results:
% - x_mean, x_var, y_mean, y_var: mean and variance of each vector
% - correlation: correlation coefficient between the two data variables
% - regression: coefficients of the linear regression line, suitable for
%   passing to polyval
% - distances: distance of each (x,y) point to the regression line

results = struct;

% Compute summary statistics
results.x_mean = sum(x) / length(x);
results.y_mean = mean(y);

results.x_var = var(x);
results.y_var = var(y);

R = corrcoef(x, y);
results.correlation = R(1,2);

% Do a linear regression
p = polyfit(x, y, 1);
results.regression = p;
line_x = 1:19;
line_y = polyval(p, line_x);

% Plot the data
plot(x, y, '+b', line_x, line_y, '-k');
xlim([1 19]);
ylim([min([2 min(y)]) max([14 max(y)])]);

% Look at distances of points from the regression line
results.distances = abs(-p(1)*x + y - p(2)) / sqrt(p(1)^2 + 1^2);
[~, min_i] = min(results.distances);
[~, max_i] = max(results.distances);
hold on
plot(x(min_i), y(min_i), 'og');
plot(x(max_i), y(max_i), 'or', 'MarkerSize', 15);
hold off

end
