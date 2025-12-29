function [x_value, y_value] = error_cyl_fom(mu_input, alpha)
mu = mu_input;
s = 2 / sqrt(3) * sqrt(mu^2 - 0.25);

%% -------- Fourth-Order Method (double precision) --------
f1 = @(x) besselj(mu, x) * cos(alpha) - bessely(mu, x) * sin(alpha);
f2 = @(x) (1 ./ x) .* (0.5 - mu) .* f1(x) + ...
          (cos(alpha) * besselj(mu - 1, x) - sin(alpha) * bessely(mu - 1, x));
h = @(x) f1(x) ./ f2(x);
A = @(x) 1 + (0.25 - mu^2) ./ (x.^2);
W = @(x) sqrt(A(x));

x0 = mu + 10000;
k = W(x0);
m = h(x0);
x1 = x0 - (1 / k) * atan(k * m);
while abs(x1-x0)>=1e-10 && x1>mu
  x0=x1;
  k = W(x0);
m = h(x0);
x1 = x0 - (1 / k) * atan(k * m);
end
v1 = [];
v1(1) = x1;
i = 1;

while x1 > mu
    i = i + 1;
    x0 = x1 - pi / sqrt(A(x1));
    k = W(x0);
    m = h(x0);
    x1 = x0 - (1 / k) * atan(k * m);
    while abs(x1 - x0) >= 1e-10 && x1 > mu
        x0 = x1;
        k = W(x0);
        m = h(x0);
        x1 = x0 - (1 / k) * atan(k * m);
    end
  if x1 > mu
      v1(i) = x1;
 end
end

%% -------- Fourth-Order Method (VPA 32-digit) --------
digits(32);
mu_sym = vpa(mu_input);
alpha_sym = vpa(alpha);

f1_vpa = @(x) besselj(mu_sym, x) * cos(alpha_sym) - ...
              bessely(mu_sym, x) * sin(alpha_sym);
f2_vpa = @(x) (1 ./ x) .* (vpa(0.5) - mu_sym) .* f1_vpa(x) + ...
              (cos(alpha_sym) * besselj(mu_sym - 1, x) - ...
               sin(alpha_sym) * bessely(mu_sym - 1, x));
h_vpa = @(x) f1_vpa(x) ./ f2_vpa(x);
A_vpa = @(x) 1 + (vpa(0.25) - mu_sym^2) ./ (x.^2);
W_vpa = @(x) sqrt(A_vpa(x));

x0 = mu_sym + vpa(10000);
k = W_vpa(x0);
m = h_vpa(x0);
x1 = x0 - (1 / k) * atan(k * m);
while abs(x1-x0)>=vpa(1e-10) 
   x0=x1;
   k = W(x0);
   m = h(x0);
   x1 = x0 - (1 / k) * atan(k * m);
end
vq = sym([]);
vq(1) = x1;
i = 1;

while x1 > mu_sym
    i = i + 1;
    x0 = x1 - vpa(pi) / sqrt(A_vpa(x1));
    k = W_vpa(x0);
    m = h_vpa(x0);
    x1 = x0 - (1 / k) * atan(k * m);
    while abs(x1 - x0) >= vpa(1e-10) && x1 > mu_sym
        x0 = x1;
        k = W_vpa(x0);
        m = h_vpa(x0);
        x1 = x0 - (1 / k) * atan(k * m);
    end
    if x1 > mu_sym
       vq(i) = x1;
    end
end

%% -------- Error Calculation --------
len = min(length(v1), length(vq));
vq_dbl = vq(1:len);
v1_dbl = v1(1:len);

error = abs(1 - v1_dbl ./ vq_dbl);
y = log(error);

x_value = v1;
y_value = y;
end
