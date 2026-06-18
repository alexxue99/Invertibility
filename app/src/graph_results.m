function [] = graph_results(filename)
A = readmatrix("../../app/data/" + filename + ".csv");
V = A(1,1);      % exact value
X = A(:,2:end);

fig=figure;
fig.Units = 'inches';
fig.Position = [1 1 6 4];

if ~startsWith(filename, 'diff')
    yline(V,'black--','LineWidth',.5);
else
    means = mean(X, 2);
    plot(1:length(means), means, 'ko', 'MarkerFaceColor', 'k');
end

hold on;

if endsWith(filename, "swap")
    shortM = endsWith(filename, "shortM_swap");
else
    shortM = endsWith(filename, "shortM");
end

if shortM
    Ns = [3, 4, 5, 6];
    xaxis = '$\log N$';
else
    Ns = [25, 100, 400];
    xaxis = '$N$';
end

boxplot(X.', ...
    'Labels', string(Ns), ...
    'Widths', 0.7, ...
    'Whisker', 1.5, ...
    'Symbol', '');

d = dictionary();
d('a') = '$a$'; d('beta') = '$\beta$'; d('gamma') = '$\gamma$'; d('diff') = 'difference'; d('m') = '$M$'; d('nu') = '$\nu t$'; 
d('pwd') = '$\mu t_{uv}$'; d('wd') = '$\mu t_w$';

ytitle = d(extractBefore(filename, "_")); 
xlabel(xaxis,'Interpreter','latex');
ylabel(ytitle,'Interpreter','latex');

set(gca,'FontSize',14);
set(gcf,'Renderer','painters');

box on;
hold off;

h = findobj(gca, 'Tag', 'Upper Whisker');
upper = max(arrayfun(@(x) max(x.YData), h));
h = findobj(gca, 'Tag', 'Lower Whisker');
lower = min(arrayfun(@(x) min(x.YData), h));
pad = 0.05*(upper - lower);

if pad == 0
    pad = 0.5;
end

ylim([lower-pad upper+pad]);

ax = gca;
ax.Units = 'normalized';
ax.Position = [0.13 0.15 0.78 0.72];

fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
fig.PaperSize = [6 4];
fig.PaperPositionMode = 'manual';

print(fig, "../../app/draws/" + filename + ".eps", ...
    '-depsc', '-vector', '-loose')
%close all;
end