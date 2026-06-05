files = dir(fullfile('../../app/data/', '*.csv'));
filenames = string({files.name});
filenames = erase(filenames, ".csv");

for f = filenames
    graph_results(f);
end