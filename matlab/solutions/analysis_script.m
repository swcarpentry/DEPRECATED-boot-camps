%% Exploring data example

% First we load some sample data, containing x and y variables
load('../data/sample1');

%% Compute summary statistics
x_mean = sum(x) / length(x);
y_mean = mean(y);

x_var = var(x);
y_var = var(y);

R = corrcoef(x, y);

disp(['Means: x=' num2str(x_mean) ', y=' num2str(y_mean)]);
disp(['Variances: x=' num2str(x_var) ', y=' num2str(y_var)]);
disp(['Correlation: ' num2str(R(1,2))]);

%% Do a linear regression
p = polyfit(x, y, 1);
disp(['Linear regression: y = ' num2str(p(1)) 'x + ' num2str(p(2))]);
line_x = 1:19;
line_y = polyval(p, line_x);

%% Plot the data
plot(x, y, '+b', line_x, line_y, '-k');
xlim([1 19]);
ylim([min([2 min(y)]) max([14 max(y)])]);

%% Look at distances of points from the regression line
% The formula for a line $ ax + by + c = 0 $ and point $(x_0, y_0)$ is
% $$ d = \frac{|ax_0 + by_0 + c|}{\sqrt{a^2 + b^2}} $$
distances = abs(-p(1)*x + y - p(2)) / sqrt(p(1)^2 + 1^2);
[min_d, min_i] = min(distances);
[max_d, max_i] = max(distances);
hold on
plot(x(min_i), y(min_i), 'og');
plot(x(max_i), y(max_i), 'or', 'MarkerSize', 15);
hold off
disp(['Maximum distance to regression line = ' num2str(max_d)]);
disp(['Average ........................... = ' num2str(mean(distances))]);
disp(['Variance .......................... = ' num2str(var(distances))]);
