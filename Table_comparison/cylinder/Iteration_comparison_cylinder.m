clear; clc;
n_values = [1000, 3000, 6000, 9000, 11000];
alp=0.75;
num_runs = 10;
time_method1 = zeros(size(n_values));
time_method2 = zeros(size(n_values));
time_method3 = zeros(size(n_values));
time_method4 = zeros(size(n_values));

for j = 1:length(n_values)
    n = round(n_values(j));
    time_method1(j) = cpu_time_cyl_fom(n, alp, num_runs);
    time_method2(j) = cpu_time_cyl_mnm(n, alp, num_runs);
    time_method3(j) = cpu_time_cyl_som(n, alp, num_runs);
    time_method4(j) = cpu_time_cyl_tom(n, alp, num_runs);
   
   
end

s1=time_method1;
s2=time_method2;
s3=time_method3;
s4=time_method4;
s1 = s1'; 
s2 = s2';
s3 = s3';
s4 = s4';


for j = 1:length(n_values)
    n = round(n_values(j));    
    [~, I1(j)] = iter_cyl_fom(n, alp);
    k1(j) = I1(j);
    [~, I2(j)] = iter_cyl_mnm(n, alp);
    k2(j) = I2(j);
    [~, I3(j)] = iter_cyl_som(n, alp);
    k3(j) = I3(j);
    [~, I4(j)] = iter_cyl_tom(n, alp);
    k4(j) = I4(j);
end  
% n1 = sum(I1);
% n2 = sum(I2);
% n3 = sum(I3);

%% plotting the results

t1=k1;
t2=k2;
t3=k3;
t4=k4;
t1 = t1';
t2 = t2';
t3 = t3';
t4 = t4';


%%
format long g
n_values = n_values';

data = [n_values, t1, s1, t2, s2, ...
           t3, s3, t4, s4];

% ---- Print header ----
fprintf('%8s | %15s %15s | %15s %15s | %15s %15s | %15s %15s\n', ...
        'mu', 'FOM-C','', 'MNM-C','', 'SOM-C','', 'TOM','');

fprintf('%8s | %15s %15s | %15s %15s | %15s %15s | %15s %15s\n', ...
        '', 'I. Count','A. Time', 'I. Count','A. Time', 'I. Count','A. Time', 'I. Count','A. Time');

fprintf('%s\n', repmat('-',1,145));

% ---- Print data ----
for i = 1:size(data,1)
    fprintf('%8d | %15.10g %15.10g | %15.10g %15.10g | %15.10g %15.10g | %15.10g %15.10g\n', ...
        data(i,1), ...
        data(i,2), data(i,3), ...
        data(i,4), data(i,5), ...
        data(i,6), data(i,7), ...
        data(i,8), data(i,9));
end
