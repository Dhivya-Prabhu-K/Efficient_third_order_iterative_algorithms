function  v1 = bes_som(mu)

    s=mu;

        %% -------- Atan(h(x)) Method --------   
        v1 = [];
        h = @(x) (besselj(mu, x)) / (besselj(mu - 1, x));
        tic;
        x0 = mu + 100000;
        x0=x0-(pi/2)*((1-sign(h(x0)))/2);
        x1 = x0 - atan(h(x0));
        while abs(x1 - x0) >= 1e-10
                x0 = x1;
                x1 = x0 - atan(h(x0));
        end
        v1(1) = x1;
        i = 1;
        while x1 > s && i<3
                i = i + 1;
                x0 = x1 - pi / 2;
                x1 = x0 - atan(h(x0));
                while abs(x1 - x0) >= 1e-10 && x1 > s
                        x0 = x1;
                        x1 = x0 - atan(h(x0));
                end
                v1(i) = x1;
        end
        while i>=3 && x1>s
                i=i+1;
                x0 = x1 - (v1(i - 2) - v1(i - 1));
                 x1 = x0 - atan(h(x0));
                while abs(x1 - x0) >= 1e-10 && x1 > s
                        x0 = x1;
                        x1 = x0 - atan(h(x0));
                end
                v1(i) = x1;
        end

        zn = length(v1);
    v1 = v1(1:1:zn-1);
    v1 = v1';
end