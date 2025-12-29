function [avg_time] = cpu_time_bes_mnm(mu, num_runs)

    s=mu;
    times = zeros(1, num_runs);

    for run = 1:num_runs
      %% -------- New Iterative Scheme --------    
        v1 = [];
        h = @(x) (besselj(mu, x)) / (besselj(mu - 1, x));
        % b = @(x) (h(x) ./ abs(h(x))) .* min(1, abs(h(x)));
        tic;
        x0 = mu + 100000;
        k = h(x0);
        x0=x0-(pi/2)*((1-sign(k))/2);
        k= h(x0);
        x1 = x0 - 1/ (1/k + k);
        iter=1;
        while abs(x1 - x0) >= 1e-10
              x0 = x1;
              k = h(x0);
              x1 = x0 - 1/ (1/k + k);
              iter=iter+1;
        end
        v1(1) = x1;
        i = 1;
        while x1 > s && i<3
                i = i + 1;
                x0 = x1 - pi / 2;
                k = h(x0);
                x1 = x0 - 1/ (1/k + k);
                while abs(x1 - x0) >= 1e-10 && x1 > s
                        x0 = x1;
                        k = h(x0);
                        x1 = x0 - 1/ (1/k + k);
                end
            v1(i) = x1;
        end
        while i>=3 && x1>s
                i=i+1;
                x0 = x1 - (v1(i - 2) - v1(i - 1));
                k = h(x0);
                x1 = x0 - 1/ (1/k + k);
                while abs(x1 - x0) >= 1e-10 && x1 > s
                        x0 = x1;
                        k = h(x0);
                        x1 = x0 - 1/ (1/k + k);
                end
            v1(i) = x1;
        end
        
        times(run) = toc;
    end
    avg_time = mean(times);
    zn = length(v1);
    v1 = v1(1:1:zn-1);
end