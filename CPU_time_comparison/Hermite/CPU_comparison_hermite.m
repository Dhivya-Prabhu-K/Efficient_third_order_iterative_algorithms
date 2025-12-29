clear all
clc;
n_values = linspace(1000000, 1300000, 15);
num_runs = 10;

time_method1 = zeros(size(n_values));
time_method2 = zeros(size(n_values));
time_method3 = zeros(size(n_values));
time_method4 = zeros(size(n_values));

for j = 1:length(n_values)
    n = round(n_values(j));
    time_method1(j) = cpu_time_fom_2(n, num_runs);
    time_method2(j) = cpu_time_som_2(n, num_runs);
    time_method3(j) = cpu_time_asy(n, num_runs);
    time_method4(j) = cpu_time_tom_2(n, num_runs);
   
end

%% plotting the results


figure;
plot(n_values, time_method1, '-sr', 'LineWidth', 2); hold on;
plot(n_values, time_method2, '-dk', 'LineWidth', 2);
plot(n_values, time_method3, '-^g', 'LineWidth', 2);
plot(n_values, time_method4, '-ob', 'LineWidth', 2);

xlabel('n', 'FontSize', 14);
ylabel('Average CPU Time (seconds)', 'FontSize', 14);

% Create the legend
h_legend = legend('FOM-H', 'SOM-H', 'ASY', 'TOM'); 
set(h_legend, 'Location', 'NorthWest', 'FontSize', 12);

title('Hermite Polynomial', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);


% Inset axes
inset_pos = [0.55, 0.27, 0.3, 0.3]; % [left, bottom, width, height]
inset_ax = axes('Position', inset_pos);
plot(inset_ax, n_values, time_method1, '-sr', 'LineWidth', 1.5); hold(inset_ax, 'on');
plot(inset_ax, n_values, time_method2, '-dk', 'LineWidth', 1.5);
plot(inset_ax, n_values, time_method4, '-ob', 'LineWidth', 1.5);
title(inset_ax, ' FOM-H vs SOM-H vs TOM', 'FontSize', 10);
% xlabel(inset_ax, 'n', 'FontSize', 10);
% ylabel(inset_ax, '', 'FontSize', 10);
grid(inset_ax, 'on');
set(inset_ax, 'FontSize', 10);
box(inset_ax, 'on');


