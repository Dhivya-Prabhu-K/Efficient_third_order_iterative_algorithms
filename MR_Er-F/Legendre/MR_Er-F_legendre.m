clear all
clc;
n_values = [1000000, 1050000, 1100000, 1200000, 1300000];



for j = 1:length(n_values)
    n = n_values(j);    
    x4 = fsweep_fourth_1(n);
 
    x2 = fsweep_second_1(n);
    r2 = abs((x2-x4)./x4);
    n2(j) = no_terms(r2);
    k2(j) = max (r2);

    x3 = fsweep_third_1(n);
    r3 = abs((x3-x4)./x4);
    k3(j) = max (r3);
    n3(j) = no_terms(r3);


end 


%% Create table with mixed formatting

format short e   % for scientific notation columns

T = table( ...
    n_values', ...
    k2',  n2', ...
    k3',  n3', ...
    'VariableNames', {'n', ...
    'SOM_MR_error','SOM_NA', ...
    'TOM_MR_error','TOM_NA'});

%% Convert integer columns to full integer display
T.n       = int64(T.n);
T.SOM_NA  = int64(T.SOM_NA);
T.TOM_NA  = int64(T.TOM_NA);

disp(T)