clear all
clc;
n_values = linspace(1000000, 1300000, 15);
num_runs = 10;

time_method1 = zeros(size(n_values));
time_method2 = zeros(size(n_values));
time_method3 = zeros(size(n_values));


for j = 1:length(n_values)
    n = round(n_values(j));
    time_method1(j) = cpu_time_legen_fourth_1(n, num_runs);
    time_method2(j) = cpu_time_legen_second_1(n, num_runs); 
    time_method3(j) = cpu_time_legen_third_1(n, num_runs);
end

% Plotting the results
figure;
plot(n_values, time_method1, '-sr', 'LineWidth', 2); hold on;
plot(n_values, time_method2, '-dk', 'LineWidth', 2);
plot(n_values, time_method3, '-ob', 'LineWidth', 2);

xlabel('n', 'FontSize', 14);
ylabel('Average CPU Time (seconds)', 'FontSize', 14);

% Create the legend
h_legend = legend('FOM-L','SOM-L','TOM'); 
set(h_legend, 'Location', 'NorthWest', 'FontSize', 12);

title('Legender Polynomial', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);

