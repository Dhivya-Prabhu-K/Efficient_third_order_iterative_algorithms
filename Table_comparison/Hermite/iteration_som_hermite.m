function q = iteration_som_hermite(n)
epsil = 1e-10;
Nn = floor(n/2);
ind = n - 2*Nn;

    xp = 0;
    delta = pi / (2*sqrt(2*(n+1)));
    xn = xp + delta;
    if ind > 0
        u0 = 0; up1 = 1; i = 1;
        xc = zeros(1, Nn + ind);
         I = zeros(1, Nn + ind);
    else
        u0 = 1; up1 = 0; i = 0;
        xc = zeros(1, Nn);
        I = zeros(1, Nn);
    end
    while i < 3
        i = i + 1;
        I(i) = 0;
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
            hs = -sqrt(2*(n + 1))*(u0 / (xn*u0 - up1));
            delta = -atan(hs)/(sqrt(2*(n + 1)));
            xp = xn;
            xn = xp + delta;
            err = abs(delta / xn);
            I(i) = I(i)+1;
        end
        xc(i) = xn;
        xn = xn + pi / (2*sqrt(2*(n+1)));
    end
    while i >= 3 && i < Nn + ind
        i = i + 1;
        I(i) = 0;
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
            hs = -sqrt(2*(n + 1))*(u0 / (xn*u0 - up1));
            delta = -atan(hs)/(sqrt(2*(n + 1)));
            xp = xn;
            xn = xp + delta;
            err = abs(delta / xn);
            I(i) = I(i)+1;
        end
        xc(i) = xn;
        xn = xn + (xc(i-1) - xc(i-2));
    end


q = sum(I);
end
