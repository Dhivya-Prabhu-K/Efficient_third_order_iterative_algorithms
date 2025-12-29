function [x_value, y_value]=compare_leg_second_1(n)
    % Error calculation of Legendre polynomial using TOM
    x_double=fsweep_second_1(n);
    x_vpa=fsweep_second_1_vpa(n);

    % Relative calculation
    rel_error = abs(1 - x_double ./ x_vpa);
    % x value
    x_value = x_double; 


   for i=1:length(x_value)
        if rel_error(i)>0
            y_value(i)=log(rel_error(i));
        else
            y_value(i)=-50;
        end
   end
   % y_value = rel_error;
end