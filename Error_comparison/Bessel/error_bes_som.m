function [x_value, y_value] = error_bes_som(mu_input)


   mu = mu_input;
   s=mu;
   f1=@(x) besselj(mu,x);
   f2=@(x) besselj(mu-1,x);
   h=@(x) f1(x)/f2(x);
   tol = 1e-10;
   x0 = mu + 10000;
   k = h(x0);
   x1 = x0 - atan(k);
   while abs(x1-x0)>=tol
         x0=x1;
         k = h(x0);
         x1 = x0 - atan(k);
   end
   v1 = [];
   v1(1) = x1;
   i = 1;
   while x1 > s && i < 3
         i = i + 1;
         x0 = x1 - pi / 2;
         k = h(x0);
         x1 = x0 - (atan(k));
         while abs(x1 - x0) >= tol && x1 > s
               x0 = x1;
               k = h(x0);
               x1 = x0 - (atan(k));
        end
        v1(i) = x1;
   end

    while i >= 3 && x1 >s
        i = i + 1;
        x0 = x1 - (v1(i - 2) - v1(i - 1));
        k = h(x0);
        x1 = x0 - (atan(k));
        while abs(x1 - x0) >= tol && x1 > s
              x0 = x1;
              k = h(x0);
              x1 = x0 - (atan(k));
        end
        if x1 > s
            v1(i) = x1;
        end
    end


    %% -------- Third-Order Method (VPA 32-digit) --------
    digits(32);
    mu = vpa(mu_input);
    s_vpa=mu;
    f1=@(x) besselj(mu,x);
    f2=@(x) besselj(mu-1,x);
    h_vpa=@(x) f1(x)/f2(x);
    % h_vpa = @(x) besselj(mu, x) / besselj(mu - 1, x);
    eps = vpa(1e-11);
    x0 = mu + vpa(10000);
    k = h_vpa(x0);
    x1 = x0 - (atan(k));
    while abs(x1-x0)>=eps
          x0=x1;
          k = h_vpa(x0);
          x1 = x0 - atan(k);
    end
    vq = sym([]);
    vq(1) = x1;
    i = 1;

    while x1 > s_vpa && i < 3
          i = i + 1;
          x0 = x1 - vpa(pi) / 2;
          k = h_vpa(x0);
          x1 = x0 - (atan(k));
        while abs(x1 - x0) >= eps && x1 > s_vpa
              x0 = x1;
              k = h_vpa(x0);
              x1 = x0 - (atan(k));
        end
        vq(i) = x1;
    end

    while i >= 3 && x1 > s_vpa
          i = i + 1;
          x0 = x1 - (vq(i - 2) - vq(i - 1));
          k = h_vpa(x0);
          x1 = x0 - (atan(k));
          while abs(x1 - x0) >= eps && x1 > s_vpa
                x0 = x1;
                k = h_vpa(x0);
                x1 = x0 - (atan(k));
         end
         if x1 > s_vpa
            vq(i) = x1;
         end
    end

    x_value=v1;
    %% -------- Error Calculation --------
    len = min(length(v1), length(vq));  % truncate to smallest
    vq_dbl = vq(1:len);
    v1_dbl = v1(1:len);                 % already in double

    error = abs(1 - v1_dbl ./ vq_dbl);
    for i=1:len
        y(i)=log(error(i));
    end
    y_value=y;

end
