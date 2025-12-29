function [avg_time, su] = iter_bes_fom(mu)
% epsil = 1e-10;
% Nn = floor(n/2);
s=mu;
num_runs = 1;
    times = zeros(1, num_runs);

for run = 1:num_runs
     %% -------- Fourth-Order Method --------
    v1 = [];
    no=[];
    f = @(x) (1/x)*(0.5 - mu)*besselj(mu, x)+besselj(mu - 1, x);
    h = @(x) besselj(mu, x)/f(x);
    A = @(x) 1 +(0.25 - mu^2)/(x^2);
    W = @(x) sqrt(A(x));

    %% Start timing fourth-order
    tic;
    x0 = mu + 100000;
    k=W(x0);
    x1 = x0 -(1/ k)*atan(k*h(x0));
    iter=1;
    while abs(x1 - x0) >= 1e-10
        x0 = x1;
    k=W(x0);
     x1 = x0 - (1/ k)* atan(k* h(x0));
     iter=iter+1;
    end
    v1(1) = x1;
    no(1)=iter;
    i = 1;
    while x1 > mu
        i = i + 1;
        x0 = x1 - pi/sqrt(A(x1));
        k=W(x0);
        x1 = x0 - (1/ k)* atan(k* h(x0));
        iter=1;
        while abs(x1 - x0) >= 1e-10 && x1 > mu
            x0 = x1;
            k=W(x0);
            x1 = x0 -(1/ k)*atan(k*h(x0));
            iter=iter+1;
        end
            v1(i) = x1;
            no(i) = iter;
     end
    times(run) = toc;
end
avg_time = mean(times);
su = sum(no);
zn = length(v1);
v1 = v1(1:1:zn-1);
end