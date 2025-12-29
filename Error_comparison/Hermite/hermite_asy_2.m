function [x] = hermite_asy_2(n)
% For computing nodes for the hermite polynomial using NM 
% initial guess and functions evaluation are done with Asymptotic formula 

    x0 = initialguess(n);
    t0 = x0./sqrt(2*n+1);
    the= acos(t0);
   
    for it =1:20
        [f, df] = asy_airy_her(n, the);
        dthe = -f./(sqrt(2)*sqrt(2*n+1).*df.*sin(the));
        the = the - dthe; %Newton step
        if norm(dthe, inf) < sqrt(eps)/10
            break;
        end
    end
    t0 = cos(the);
    x = sqrt(2*n+1).*t0; %change of variable
end

% Evaluating the function and its derivative

function [f, df] = asy_airy_her(n, the)
    mu2  = 2*n+1;
    cosv = cos(the);
    sinv = sin(the);
    sin2v = 2.*cosv.*sinv;
    et = 0.5.*the-0.25.*sin2v;
    ch = -(3.*et./2).^(2/3);
    ph = (-ch./ sinv.^2).^(1/4);
    con = 2.* sqrt(pi).* mu2.^(1/6).*ph;
    z = mu2.^(2/3).*ch;
    Airy = real(airy(0, z));
    DAiry = real(airy(1, z));

    % some constants
    a0 = 1;
    b0 = 1;
    a1 = 15/144;
    b1 = -7/5 * a1;
    a2 = 5*7*9*11/2/144^2;
    b2 = -13/11 * a2;
    a3 = 7*9*11*13*15*17/6/144^3;
    b3 = -19/17 * a3;

    % Polynomials need to evalute the coefficients of Expansions
    u0 = 1;
    u1 = (cosv.^3 - 6 .* cosv) ./ 24;
    u2 = (-9 .* cosv.^4 + 249 .* cosv.^2 + 145) ./ 1152;
    u3 = (-4042 .* cosv.^9 + 18189 .* cosv.^7 - 28287 .* cosv.^5 - 151995 .* cosv.^3 - 259290 .* cosv) ./ 414720;
    
    % Evaluation of the series for f
    A0 = 1;
    f = A0 .* Airy;
    B0 = -(a0 .* u1 .* ph.^6 + a1 .* u0) ./ ch.^2;
    f = f + B0 .* DAiry ./ mu2.^(4/3);
    A1 = (b0 .* u2 .* ph.^12 + b1 .* u1 .* ph.^6 + b2 .* u0) ./ ch.^3;
    f = f + A1 .* Airy ./ mu2.^2;
    B1 = -(u3 .* ph.^18 + a1 .* u2 .* ph.^12 + a2 .* u1 .* ph.^6 + a3 .* u0) ./ ch.^5;
    f = f + B1 .* DAiry ./ mu2.^(10/3);

    f = con .* f;

    % Evaluation of derivative values
    con1= sqrt(2*pi).*mu2.^(1/3)./ ph;
   
    % Polynomials need to evaluate the coefficients for Df
    v0 = 1;
    v1 = (cosv.^3 + 6 .* cosv) ./ 24;
    v2 = (15 .* cosv.^4 - 327 .* cosv.^2 - 143) ./ 1152;
    v3 = (259290 .* cosv + 238425 .* cosv.^3 - 36387 .* cosv.^5 + 18189 .* cosv.^7 - 4042 .* cosv.^9) ./ 414720;

    % Evaluation of the series for df

    C0 = -(b0 .* ph.^6 .* v1 + b1 .* v0) ./ ch;
    df = C0 .* Airy ./ mu2.^(2/3);
    D0 = a0 .* v0;
    df = df + D0 .* DAiry;
    C1 = -(v3 .* ph.^18 + b1 .* v2 .* ph.^12 + b2 .* v1 .* ph.^6 + b3 .* v0) ./ ch.^4;
    df = df + C1 .* Airy ./ mu2.^(8/3);
    D1 = (a0 .* v2 .* ph.^12 + a1 .* v1 .* ph.^6 + a2 .* v0) ./ ch.^3;
    df = df + D1 .* DAiry ./ mu2.^2;

    df = con1 .* df;
end

% Function file for computing initial guess

function x_in = initialguess(n)
    if n < 20
        error('n must be >= 20.');
    end
    if mod(n, 2) == 1  % odd
        nh = (n - 1)/2;
        a = 0.5;
    else  % even
        nh = n/2;
        a = -0.5;
    end
     mu = 4 * nh + 2 * a + 2;

     T = @(t) t.^(2/3) .* (1 + (5/48) * t.^(-2) - (5/36) * t.^(-4) + (77125/82944) * t.^(-6) ...
        - (108056875/6967296) * t.^(-8) + (162375596875/334430208) * t.^(-10));

    airyrts = -T(3 * pi / 8 * (4 * (1:nh) - 1));

    airy_roots = [-2.338107410459767, -4.087949444130970, -5.520559828095551, ...
                       -6.786708090071759, -7.944133587120853, -9.022650853340980, ...
                       -10.04017434155809, -11.00852430373326, -11.93601556323626, ...
                       -12.82877675286576];
    airyrts(1:min(10, nh)) = airy_roots(1:min(10, nh));

    x_airy = sqrt(abs( ...
        mu + (2.^(2/3)) .* airyrts .* mu.^(1/3) ...
        + (1/5 * 2.^(4/3)) .* airyrts.^2 .* mu.^(-1/3) ...
        + (11/35 - a^2 - 12/175) .* airyrts.^3 ./ mu ...
        + ((16/1575) .* airyrts + (92/7875) .* airyrts.^4) .* 2.^(2/3) .* mu.^(-5/3) ...
        - ((15152/3031875) .* airyrts.^5 + (1088/121275) .* airyrts.^2) .* 2.^(1/3) .* mu.^(-7/3) ...
    ));
    x_airy = flip(x_airy);

    % initial guess from non-linear problem

    init = (pi/2) * ones(1, nh);
    cons2 = ((4 * nh + 3) - 4 * (1:nh)) / mu * pi;
    for k = 1:7
        N = init - sin(init) - cons2;
        dN = 1 - cos(init);
        dinit = N ./ dN;
        init = init - dinit;
    end

    tnk = cos(init/2).^2;
    x_sin = sqrt(mu .* tnk - ((tnk + 1/4) ./ (tnk - 1).^2 + (3*a^2 - 1)) ./ (3 * mu));

    p = 0.4985 + eps;
    x_in = [x_sin(1:floor(p * n)), x_airy(ceil(p * n):end)];

    if mod(n, 2) == 1
        x_in = [0, x_in];
        x_in = x_in(1:(nh + 1));
    else
        x_in = x_in(1:nh);
    end
end