% Generate sample data files for the analysis script

% Preallocate array of structures
data(6) = struct('x', 0, 'y', 0);

% The first four datasets come from Anscombe, F. J. (1973).
% "Graphs in Statistical Analysis". American Statistician 27 (1): 17–21.
anscombe = [
   10.0000    8.0400   10.0000    9.1400   10.0000    7.4600    8.0000    6.5800;
    8.0000    6.9500    8.0000    8.1400    8.0000    6.7700    8.0000    5.7600;
   13.0000    7.5800   13.0000    8.7400   13.0000   12.7400    8.0000    7.7100;
    9.0000    8.8100    9.0000    8.7700    9.0000    7.1100    8.0000    8.8400;
   11.0000    8.3300   11.0000    9.2600   11.0000    7.8100    8.0000    8.4700;
   14.0000    9.9600   14.0000    8.1000   14.0000    8.8400    8.0000    7.0400;
    6.0000    7.2400    6.0000    6.1300    6.0000    6.0800    8.0000    5.2500;
    4.0000    4.2600    4.0000    3.1000    4.0000    5.3900   19.0000   12.5000;
   12.0000   10.8400   12.0000    9.1300   12.0000    8.1500    8.0000    5.5600;
    7.0000    4.8200    7.0000    7.2600    7.0000    6.4200    8.0000    7.9100;
    5.0000    5.6800    5.0000    4.7400    5.0000    5.7300    8.0000    6.8900];

for i=1:4
    data(i).x = anscombe(:, i*2-1);
    data(i).y = anscombe(:, i*2);
end

% We also include a couple of randomly generated datasets
n_points = length(data(1).x);
xmin = 1; xmax = 19;

% Linearly distributed
data(5).x = xmin + (xmax-xmin)*rand(n_points, 1);
data(5).y = 0.66 * data(5).x + 3.3 + 1.5*randn(n_points, 1);

% Following a quadratic curve
data(6).x = xmin + (xmax-xmin)*rand(n_points, 1);
data(6).y = 1.1 + ((data(6).x - 10)/3).^2 + 1.5*randn(n_points, 1);


% Write data to file
save('all_samples', 'data');
for i=1:6
    fname = ['sample' num2str(i)];
    x = data(i).x;
    y = data(i).y;
    save(fname, 'x', 'y');
end
