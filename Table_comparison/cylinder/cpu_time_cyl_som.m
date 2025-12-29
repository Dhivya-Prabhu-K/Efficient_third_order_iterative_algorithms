function avg_time = cpu_time_cyl_som(mu, alp, num_runs)

    s=mu;
    times = zeros(1, num_runs);

    for run = 1:num_runs
        %% -------- Atan(h(x)) Method --------   
        v1 = [];
        f1=@(x) besselj(mu,x)*cos(alp)-bessely(mu,x)*sin(alp);
        f2=@(x) besselj(mu-1,x)*cos(alp)-bessely(mu-1,x)*sin(alp);
        h=@(x) f1(x)/f2(x);
        tic;
        x0 = mu + 100000;
        k = h(x0);
        x0=x0-(pi/2)*((1-sign(k))/2);
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
        times(run) = toc;
    end
    avg_time = mean(times);
        zn = length(v1);
v1 = v1(1:1:zn-1);
end