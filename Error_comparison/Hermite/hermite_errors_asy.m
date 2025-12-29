function [x_value, y_value] = hermite_errors_asy(n)
    % Ensure the Symbolic Math Toolbox is available
    assert(license('test','Symbolic_Toolbox'), 'Symbolic Toolbox is required');

    % Double precision ASY
    x_double = hermite_asy_2(n);  

    % Extended precision ASY
    x_vpa = hermite_asy_2_vpa(n);    


    x_vpa_dbl = (x_vpa);  

    % Ensure dimensions match
    assert(length(x_double) == length(x_vpa_dbl), 'Output lengths differ');

    % Compute relative errors
    rel_error = abs(1 - x_double ./ x_vpa_dbl);
    x_value = x_double;  % your x values
for i=1:length(x_value)
    if rel_error(i)>0
     y_value(i)=log(rel_error(i));
    else
      y_value(i)=-50;
    end
end

end
