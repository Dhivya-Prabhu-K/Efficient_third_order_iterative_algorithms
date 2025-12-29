function avg_time = cpu_time_cyl_fom(mu, alpha, num_runs)
% epsil = 1e-10;
% Nn = floor(n/2);
% ind = n - 2*Nn;
 s=mu;
times = zeros(1, num_runs);
for run = 1:num_runs
     %% -------- Fourth-Order Method --------
    v1 = [];
    % f = @(x) (1/x)*(0.5 - mu)*besselj(mu, x)+besselj(mu - 1, x);
    f1=@(x) besselj(mu,x)*cos(alpha)-bessely(mu,x)*sin(alpha); % Normal transform
    f2 = @(x) (1/x)*(0.5-mu)*f1(x)+(cos(alpha)*besselj(mu-1,x)-sin(alpha)*bessely(mu-1,x));
    h = @(x) f1(x)/f2(x);
    A = @(x) 1 +(0.25 - mu^2)/(x^2);
    W = @(x) sqrt(A(x));

    %% Start timing fourth-order
    tic;
    x0 = mu + 100000;
    k=W(x0);
    x1 = x0 -(1/ k)*atan(k*h(x0));
    while abs(x1 - x0) >= 1e-10
        x0 = x1;
    k=W(x0);
     x1 = x0 - (1/ k)* atan(k* h(x0));
    end
    v1(1) = x1;
    i = 1;
    while x1 > s
        i = i + 1;
        x0 = x1 - pi/sqrt(A(x1));
        k=W(x0);
        x1 = x0 - (1/ k)* atan(k* h(x0));
        while abs(x1 - x0) >= 1e-10 && x1 > s
            x0 = x1;
            k=W(x0);
            x1 = x0 -(1/ k)*atan(k*h(x0));
        end
            v1(i) = x1; 
     end
    times(run) = toc;
end
avg_time = mean(times);
zn = length(v1);
v1 = v1(1:1:zn-1);
end