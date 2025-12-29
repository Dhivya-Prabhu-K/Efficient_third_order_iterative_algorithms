function [x]=fsweep_fourth_1_vpa(n)
% To find all zeros of Legendre polynomial using FOM
% This code is a modification of the code available in 
% `https://personales.unican.es/segurajj/gaussian.html'
    digits(32);
    ind=n-2*floor(n/2);
    x = sym(zeros(floor(n/2)+ind, 1));
    L=vpa(2*n+1);
    nor=vpa(1.e-130);
    ind=n-2*floor(n/2);
        if ind>0
            y0=0;y1=nor; i=1; x(1)=0;
        else
            y0=nor;y1=0; i=0;
        end
    xi= vpa(tanh(((ind+1)*pi)/sqrt(L^2-1)));
    xt=0;
    t=xi-xt;
    [y0, y1] = seriestay_1_vpa(n, xt, t, y0, y1);
    [xt, x0, y0, y1] = fixedzser_fourth_1_vpa(n, t, y0, y1);
    i = i + 1;
    x(i) = xt;
    sqde=sqrt(0.25*((L*L-1)*(1-xt*xt)));
    de=vpa(-tanh(pi/sqde));
    xi=(xt-de)/(1-xt*de);
    t=xi-x0;
    [y0,y1]=seriestay_1_vpa(n,x0,t,y0,y1);
        while i<floor(n/2)+ind
            [xt,x0,y0,y1]=fixedzser_fourth_1_vpa(n,xi,y0,y1);
            i=i+1;
            x(i)=xt;   
            sqde=sqrt(0.25*((L*L-1)*(1-xt*xt)));
            de=vpa(-tanh(pi/sqde));
            xi=(xt-de)/(1-xt*de);
            t=xi-x0;
            [y0,y1]=seriestay_1_vpa(n,x0,t,y0,y1);
        end
end