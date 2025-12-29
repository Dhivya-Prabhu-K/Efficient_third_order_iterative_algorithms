function [x_value, y_value] = r_error_fom_2(n)
% Double precision FOM
epsil = (10)^(-10);
Nn = floor(n/2);
ind = n - 2*Nn;

xp = 0;
delta = (1 + ind)/2 * pi / sqrt(2*n + 1);
xn = xp + delta;

if ind > 0
    u0 = 0; up1 = 1; i = 1;
    xc = zeros(1, Nn + ind);
else
    u0 = 1; up1 = 0; i = 0;
    xc = zeros(1, Nn);
end

while i < Nn + ind
    i = i + 1;
    err = 1 + epsil;
    while err > epsil
        t = xn - xp;
        delta = t;
        um1 = 0; um2 = 0;
        suma = u0 + up1 * t;
        sumad = up1;
        k = -1;
        errod = 1;
        tol = 1e-25;
        max_k = 50;
        while (errod > tol && k < max_k) || k < 5
            k = k + 1;
            up2 = (xp^2 - 2*n - 1)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
            um2 = um1; um1 = u0; u0 = up1; up1 = up2;
            ntd = up2 * t; sumad = sumad + ntd;
            if i > 1 + ind
                errod = abs(ntd / sumad);
            end
            t = t * delta / (k + 2);
            nt = up2 * t;
            suma = suma + nt;
        end
        u0 = suma; up1 = sumad;
        hs = u0 / up1;
        ws = sqrt(2*n + 1 - xn^2);
        delta = -atan(ws * hs) / ws;
        xp = xn;
        xn = xp + delta;
        err = abs(delta / xn);
    end
    xc(i) = xn;
    delta = pi / sqrt(2*n + 1 - xn^2);
    
    xn = xn + delta;
end

x16 = xc;  % Store 16-digit result

% Extended precision FOM
digits(32)
n = vpa(n);
epsil = vpa((10)^(-11));
Nn = floor(n/2);
ind = n - 2*Nn;

xp = vpa(0);
delta = ((1 + ind)/2) * vpa(pi) / sqrt(2*n + 1);
xn = xp + delta;

if ind > 0
    u0 = vpa(0); up1 = vpa(1); i = 1;
    xc = sym(zeros(1, double(Nn + ind)));
else
    u0 = vpa(1); up1 = vpa(0); i = 0;
    xc = sym(zeros(1, double(Nn)));

end

while i < Nn + ind
    i = i + 1;
    err = vpa(1) + epsil;
    while err > epsil
        t = xn - xp;
        delta = t;
        um1 = vpa(0); um2 = vpa(0);
        suma = u0 + up1 * t;
        sumad = up1;
        k = -1;
        errod = vpa(1);
        tol = vpa(1e-80);
        max_k = 100;
        while (errod > tol && k < max_k) || k < 5
            k = k + 1;
            up2 = (xp^2 - 2*n - 1)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
            um2 = um1; um1 = u0; u0 = up1; up1 = up2;
            ntd = up2 * t; sumad = sumad + ntd;
            if i > 1 + ind
                errod = abs(ntd / sumad);
            end
            t = t * delta / (k + 2);
            nt = up2 * t;
            suma = suma + nt;
        end
        u0 = suma; up1 = sumad;
        hs = u0 / up1;
        ws = sqrt(2*n + 1 - xn^2);
        delta = -vpa(atan(ws * hs)) / ws;
        xp = xn;
        xn = xp + delta;
        err = abs(delta / xn);
    end
    xc(i) = xn;
    delta = vpa(pi) / sqrt(2*n + 1 - xn^2);
    
    xn = xn + delta;
end

x32 = xc;  % 32-digit result

% --- Compute relative error ---
error = zeros(1, double(Nn + ind));
for j = 1:double(Nn + ind)
    error(j) = abs(1 - x16(j) / (x32(j)));
end
x_value = x16;  % your x values
for i=1:length(xc)
    if error(i)>0
     y_value(i)=log(error(i));
    else
      y_value(i)=-50;
    end
end 
end
