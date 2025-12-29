function avg_time = cpu_time_asy(n, num_runs)
times = zeros(1, num_runs);

for run = 1:num_runs
    tic;
    hermite_asy_2(n);
    times(run) = toc;
end

avg_time = mean(times);
end
