function [x,x0,y0,y1]=fixedzser_second_1(n,x0,y0,y1)
% fixed point on z and using Taylor for (1-x)^(a+1)/2*(1+x)^(b+1)/2*P
    epsil=10^(-15);
    h=(n+1)*(y0/(-n*x0*y0+(1-x0*x0)*y1));
    argu=(1/(n+1))*atan(h);
    de=tanh(argu);
    x=(x0-de)/(1-x0*de);   
    erro=abs(1-x/x0);
        while erro>epsil
            t=x-x0;
            [y0,y1]=seriestay_1(n,x0,t,y0,y1);
            x0=x;
            h=(n+1)*(y0/(-n*x0*y0+(1-x0*x0)*y1));
            argu=(1/(n+1))*atan(h);
            de=tanh(argu);
            x=(x0-de)/(1-x0*de);   
            erro=abs(1-x/x0);
        end
end