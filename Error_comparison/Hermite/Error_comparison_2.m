clear all
clc;
n=10000;
    % Get data from each method
    [x1, y1] = r_error_fom_2(n);          % First: Red
    [x2, y2] = r_error_som_2(n);          % Second: Black
    [x3 ,y3] = hermite_errors_asy(n); % Third: Green (returns only y)
    [x4, y4] = r_error_tom_2(n);          % Fourth: Blue

    % Define index for y3 if x not available
   

    % Plot each curve with specific color and marker
    figure;
    semilogy(x1, y1, 'o-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'FOM');
    hold on;
    semilogy(x2, y2, 's-', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', 'SOM');
    %semilogy(x3,  y3, '^-', 'Color', 'g', 'LineWidth', 1.5, 'DisplayName', 'ASY');
    semilogy(x4, y4, 'd-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'TOM');

    % Labels and grid
    grid on;
    xlabel('zeros of H_{10000}(x)');
    ylabel('Relative error');
    % title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
    legend('Location', 'best');
    set(gca, 'FontSize', 12);
%%

figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM-H');  % Red
hold on;
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'SOM-H');  % Black
%plot(x3, y3, '-', 'Color', 'g', 'LineWidth', 1, 'DisplayName', 'ASY');  % Green
plot(x4, y4, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TOM');  % Blue

% Labels and grid
grid on;
xlabel('zeros of H_{10000}(x)');
ylabel('Relative error');
% title(sprintf('Relative Error Comparison of Hermite Zeros (n = %d)', n));
legend('Location', 'best');
set(gca, 'FontSize', 12);
