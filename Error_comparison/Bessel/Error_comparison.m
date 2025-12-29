clear all
clc;
mu=input('value of mu=');


    % Get data from each method
    [x1, y1] = error_bes_fom(mu);
    [x2, y2] = error_bes_som(mu);         
    [x3 ,y3] = error_bes_mnm(mu);
    [x4, y4] = error_bes_tom(mu);          


%%

figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM-B');  % Red
hold on;
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'SOM-B');  % Black

plot(x3, y3, '-', 'Color', 'g', 'LineWidth', 1, 'DisplayName', 'MNM-B');  % Green
plot(x4, y4, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TOM');  % Blue

% Labels and grid
grid on;
xlabel('zeros of J_{10000}(x)');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of cylinder (mu = %d)', mu));
legend({'FOM-B','SOM-B','MNM-B','TOM-B'},'Location', 'best');
set(gca, 'FontSize', 12);

