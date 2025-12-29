function [x_value, y_value]=r_error_som_2(n)
% Double precision SOM
epsil = 10^(-10);
Nn = floor(n/2);
ind = n - 2*Nn;
xp = 0;  
 delta = pi / (2*sqrt(2*(n+1)));
xn = xp + delta;
iter=zeros(1, Nn+ind);
ycd =(zeros(1, double(Nn+ind)));
if ind > 0
    u0 = 0;
    up1 = 1;
    i = 1;
    xc(i) = 0;
    deri(i) = 1;
else
    u0 = 1;
    up1 = 0;
    i = 0;
end

while i < Nn + ind
    i = i + 1;
    err = 1 + epsil;
    while err > epsil
        iter(i)=iter(i)+1;

        t = xn - xp;
        delta = t;
        um1 = 0;
        um2 = 0;
        suma = u0 + up1*t;
        sumad = up1;
        k = -1;
        errod = 1;
        tol = 1e-25;
        max_k = 50;

        while (errod > tol && k < max_k) || k < 5
            k = k + 1;
            up2 = (xp*xp-2*n-1)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
            um2 = um1;
            um1 = u0;
            u0 = up1;
            up1 = up2;
            ntd = up2*t;
            sumad = sumad + ntd;
            if i > 1+ind
                errod = abs(ntd/sumad);
            end
            t = t*delta/(k+2);
            nt = up2*t;
            suma = suma + nt;
        end

        u0 = suma;
        up1 = sumad;
        hs = -sqrt(2*(n+1))*(u0/(2*xn*u0-up1));
      
        delta =  -atan(hs)/(sqrt(2*(n + 1)));
        xp = xn;
        xn = xp + delta;
        err = abs(delta/xn);
    end

    ycd(i) = xn;
    xc(i) = xn;

    if i < 3
        xn = xc(i) + pi/ (2*sqrt(2*(n+1)));
    else
        xn = xc(i) + (xc(i-1) - xc(i-2));
    end
end
ycd';

%% Extended precision SOM
digits(32)

n = vpa(n);
epsil = vpa((10)^(-11));
Nn = floor(n/2);
ind = n - 2*Nn;
xp = vpa(0);
delta = vpa(pi) / (2*sqrt(2*(n+1)));
xn = xp + delta;

%%% Fix array allocation!
xc = sym(zeros(1, double(Nn+ind)));
ycq = sym(zeros(1, double(Nn+ind)));

if ind > 0
    u0 = vpa(0);
    up1 = vpa(1);
    i = 1;
    xc(i) = vpa(0);
    deri(i) = vpa(1);
else
    u0 = vpa(1);
    up1 = vpa(0);
    i = 0;
end

while i < Nn + ind
    i = i + 1;
    err = vpa(1) + epsil;
    while err > epsil
       
        t = xn - xp;
        delta = t;
        um1 = vpa(0);
        um2 = vpa(0);
        suma = u0 + up1*t;
        sumad = up1;
        k = -1;
        errod = vpa(1);
        tol = vpa(1e-80);
        max_k = 100;

        while (errod > tol && k < max_k) || k < 5
            k = k + 1;
            up2 = (xp*xp-2*n-1)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
            um2 = um1;
            um1 = u0;
            u0 = up1;
            up1 = up2;
            ntd = up2*t;
            sumad = sumad + ntd;
            if i > 1+ind
                errod = abs(ntd/sumad);
            end
            t = t*delta/(k+2);
            nt = up2*t;
            suma = suma + nt;
        end

        u0 = suma;
        up1 = sumad;
        hs = -sqrt(2*(n+1))*(u0/(2*xn*u0-up1));
 
        delta =  -vpa(atan(hs))/(sqrt(2*(n + 1)));
        xp = xn;
        xn = xp + delta;
        err = abs(delta/xn);
    end

    ycq(i) = xn;
    xc(i) = xn;

    if i < 3
        xn = xc(i) + vpa(pi) / (2*sqrt(2*(n+1)));
    else
        xn = xc(i) + (xc(i-1) - xc(i-2));
    end
end


% Output vector from extended precision
ycq';

%% Compute relative error
error = zeros(1, double(Nn+ind));

for i = 1:Nn+ind
    error(i) = abs(1 - ycd(i) / ycq(i)); % <--- convert ycd(i) back to double for comparison
end
x_value = ycd;  % your x values
for i=1:length(ycd)
    if error(i)>0
     y_value(i)=log(error(i));
    else
      y_value(i)=-50;
    end
end % your y values
end

