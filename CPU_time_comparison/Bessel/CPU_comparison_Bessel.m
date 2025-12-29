clear; clc;
n_values = linspace(10000, 110000, 15);
num_runs = 10;

time_method1 = zeros(size(n_values));
time_method2 = zeros(size(n_values));
time_method3 = zeros(size(n_values)); 
time_method4 = zeros(size(n_values));

for j = 1:length(n_values)
    n = round(n_values(j));
    time_method1(j) = cpu_time_bes_fom(n, num_runs);
    time_method2(j) = cpu_time_bes_som(n, num_runs);
    time_method3(j) = cpu_time_bes_mnm(n, num_runs);
    time_method4(j) = cpu_time_bes_tom(n, num_runs);   
end

%% Plotting the results


figure;
plot(n_values, time_method3, '-^g', 'LineWidth', 2); hold on;
plot(n_values, time_method1, '-sr', 'LineWidth', 2);
plot(n_values, time_method2, '-dk', 'LineWidth', 2);
plot(n_values, time_method4, '-ob', 'LineWidth', 2);

xlabel('\mu', 'FontSize', 14);
ylabel('Average CPU Time (seconds)', 'FontSize', 14);

% Create the legend
h_legend = legend('MNM-B','FOM-B','SOM-B','TOM'); 
set(h_legend, 'Location', 'NorthWest', 'FontSize', 12);

% title('Bessel zeros', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);
