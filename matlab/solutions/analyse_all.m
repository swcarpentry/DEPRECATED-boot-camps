%% Analyse all data sets

load('../data/all_samples');

for i=1:length(data)
    figure;
    results = analysis_func(data(i).x, data(i).y);
    disp(results);
end
