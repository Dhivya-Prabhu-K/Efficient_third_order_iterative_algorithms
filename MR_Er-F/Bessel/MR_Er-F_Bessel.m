clear; clc;
n_values = [1000, 3000, 6000, 9000, 11000];

for j = 1:length(n_values)
    n = n_values(j);    
    x4 = bes_fom(n);
    
    x2 = bes_som(n);
    if length(x4)>length(x2)
        x4 = x4(2:1:length(x4));
    end
    m4(j) = max(x4);
    k2(j) = max (abs((x2-x4)./x4));

    x3 = bes_tom(n);
    k3(j) = max (abs((x3-x4)./x4));

    x1 = bes_mnm(n);
    k1(j) = max (abs((x1-x4)./x4));

end 


%% Create numeric table (no strings)

T = table(n_values', ...
          k1', ...
          k2', ...
          k3', ...
    'VariableNames', {'n', ...
    'MNM_MR_error', ...
    'SOM_MR_error', ...
    'TOM_MR_error'});

%% Display format
format short e
disp(T)