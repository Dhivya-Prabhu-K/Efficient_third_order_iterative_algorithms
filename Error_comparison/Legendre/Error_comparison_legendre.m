clear all
clc;
n=10000;

% Storing values
[x1, y1] = compare_leg_fourth_1(n);         
[x2, y2] = compare_leg_second_1(n);         
[x3, y3] = compare_leg_third_1(n);           


figure;
plot(x1, y1, '-', 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'FOM'); hold on;  % Red
plot(x2, y2, '-', 'Color', 'k', 'LineWidth', 1, 'DisplayName', 'SOM');  % Black
plot(x3, y3, '-', 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TOM');  % Blue


grid on;
xlabel('zeros of P_{10000}(x)');
ylabel('Relative error');
legend('Location', 'best');
set(gca, 'FontSize', 12);




