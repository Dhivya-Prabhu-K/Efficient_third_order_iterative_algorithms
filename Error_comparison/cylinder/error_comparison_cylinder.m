clear all
clc;
mu=input('value of mu=');
alpha=input('value of alpha=');

    % Get data from each method
    [x1, y1] = error_cyl_fom(mu, alpha); 
    [x2, y2] = error_cyl_som(mu, alpha);          % Second: Black
    [x3 ,y3] = error_cyl_mnm(mu, alpha); % Third: Green (returns only y)
    [x4, y4] = error_cyl_tom(mu, alpha);          % Fourth: Blue



figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM');  % Red
hold on;
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'SOM-C');  % Black

plot(x3, y3, '-', 'Color', 'g', 'LineWidth', 1, 'DisplayName', 'MNM-C');  % Green
plot(x4, y4, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TOM');  % Blue

% Labels and grid
grid on;
xlabel('zeros of C_{10000}(x)');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of cylinder (mu = %d)', mu));
legend('Location', 'best');
set(gca, 'FontSize', 12);
