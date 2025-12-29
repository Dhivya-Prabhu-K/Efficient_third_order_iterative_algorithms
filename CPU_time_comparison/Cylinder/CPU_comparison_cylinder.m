clear; clc;
n_values = linspace(10000, 110000, 15);
alp=0.75;
num_runs = 10;
time_method1 = zeros(size(n_values));
time_method2 = zeros(size(n_values));
time_method3 = zeros(size(n_values));
time_method4 = zeros(size(n_values));

for j = 1:length(n_values)
    n = round(n_values(j));
    time_method1(j) = cpu_time_cyl_fom(n, alp, num_runs);
    time_method4(j) = cpu_time_cyl_tom(n, alp, num_runs);
    time_method2(j) = cpu_time_cyl_som(n, alp, num_runs);
    time_method3(j) = cpu_time_cyl_mnm(n, alp, num_runs);
   
   
end


%%
figure;

plot(n_values, time_method3, '-^g', 'LineWidth', 2);hold on;
plot(n_values, time_method1, '-sr', 'LineWidth', 2);
plot(n_values, time_method2, '-dk', 'LineWidth', 2);
plot(n_values, time_method4, '-ob', 'LineWidth', 2);

xlabel('\mu', 'FontSize', 14);
ylabel('Average CPU Time (seconds)', 'FontSize', 14);

% Create the legend
h_legend = legend('MNM-C','FOM-C','SOM-C','TOM-C'); 
set(h_legend, 'Location', 'NorthWest', 'FontSize', 12);

% title('Cylinder function \alpha=0.5', 'FontSize', 16);
grid on;
set(gca, 'FontSize', 12);
