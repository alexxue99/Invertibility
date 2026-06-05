function [] = graph_results(filename)
A = readmatrix("../../app/data/" + filename + ".csv");
V = A(1,1);      % exact value
X = A(:,2:end);

fig=figure;

if ~startsWith(filename, 'diff')
    yline(V,'black--','LineWidth',.5);
else
    means = mean(X, 2);
    plot(1:length(means), means, 'ko', 'MarkerFaceColor', 'k');
end

hold on;

shortM = endsWith(filename, "shortM");
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

lower = max(X(:));
upper = min(X(:));

for i = 1:size(X, 1)
Q = quantile(X(i, :), [0.25 0.75]);
IQR = Q(2) - Q(1);
lower = min(lower, Q(1) - 1.5*IQR);
upper = max(upper, Q(2) + 1.5*IQR);
end

lower = min(X(X >= lower));
upper = max(X(X <= upper));
pad = 0.05*(upper - lower);

if pad == 0
    pad = 0.5;
end

ylim([lower-pad upper+pad]);

saveas(gcf, "../../app/draws/" + filename + ".eps", 'epsc')
close all;
end