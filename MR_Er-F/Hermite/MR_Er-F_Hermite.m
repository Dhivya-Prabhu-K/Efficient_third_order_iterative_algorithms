clear all
clc;
n_values = [1000000, 1050000, 1100000, 1200000, 1300000];



for j = 1:length(n_values)
    n = n_values(j);    
    x4 = her_fom(n);
 
    x2 = her_som(n);
    k2(j) = max (abs((x2-x4)./x4));

    x3 = her_tom(n);
    k3(j) = max (abs((x3-x4)./x4));


end 


%% Create numeric table (no strings)

T = table(n_values', ...
          k2', ...
          k3', ...
    'VariableNames', {'n', ...
    'SOM_MR_error', ...
    'TOM_MR_error'});

%% Display format
format short e
disp(T)