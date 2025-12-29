function  [avg_time, no] = cpu_time_bes_som(mu, num_runs)

    s=mu;
    times = zeros(1, num_runs);

    for run = 1:num_runs
        %% -------- Atan(h(x)) Method --------   
        v1 = [];
        no=[];
        h = @(x) (besselj(mu, x)) / (besselj(mu - 1, x));
        tic;
        x0 = mu + 100000;
        x0=x0-(pi/2)*((1-sign(h(x0)))/2);
        x1 = x0 - atan(h(x0));
        iter=1;
        while abs(x1 - x0) >= 1e-10
                x0 = x1;
                x1 = x0 - atan(h(x0));
                iter=iter+1;
        end
        v1(1) = x1;
        no(1) = iter;
        i = 1;
        while x1 > s && i<3
                i = i + 1;
                x0 = x1 - pi / 2;
                x1 = x0 - atan(h(x0));
                 iter=1;
                while abs(x1 - x0) >= 1e-10 && x1 > s
                        x0 = x1;
                        x1 = x0 - atan(h(x0));
                         iter=iter+1;
                end
                v1(i) = x1;
                no(i) = iter;
        end
        while i>=3 && x1>s
                i=i+1;
                x0 = x1 - (v1(i - 2) - v1(i - 1));
                 x1 = x0 - atan(h(x0));
                 iter=1;
                while abs(x1 - x0) >= 1e-10 && x1 > s
                        x0 = x1;
                        x1 = x0 - atan(h(x0));
                        iter=iter+1;
                end
                v1(i) = x1;
                 no(i) = iter;
        end
        times(run) = toc;
    end
    avg_time = mean(times);
end