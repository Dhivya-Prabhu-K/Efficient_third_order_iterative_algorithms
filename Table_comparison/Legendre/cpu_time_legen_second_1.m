function avg_time = cpu_time_legen_second_1(n, num_runs)

   times = zeros(1, num_runs);

    for run = 1:num_runs
        tic;
        fsweep_second_1(n);
        times(run) = toc;
    end
    avg_time = mean(times);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%function fixedzser
