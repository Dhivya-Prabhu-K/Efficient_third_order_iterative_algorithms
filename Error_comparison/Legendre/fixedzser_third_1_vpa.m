function [x,x0,y0,y1]=fixedzser_third_1_vpa(n,x0,y0,y1)
% TOM on z and using Taylor for (1-x)^(a+1)/2*(1+x)^(b+1)/2*P
    digits(32);
    epsil=vpa(10^(-16));
    h=(n+1)*(y0/(-n*x0*y0+(1-x0*x0)*y1));
    argu=vpa(2*h/((n+1)*(2+h^2+2*x0*h)));
    de=vpa(tanh(argu));
    x=(x0-de)/(1-x0*de);   
    erro=abs(1-x/x0);
        while erro>epsil
            t=x-x0;
            [y0,y1]=seriestay_1(n,x0,t,y0,y1);
            x0=x;
            h=(n+1)*(y0/(-n*x0*y0+(1-x0*x0)*y1));
            argu=vpa(2*h/((n+1)*(2+h^2+2*x0*h)));
            de=vpa(tanh(argu));
            x=(x0-de)/(1-x0*de);   
            erro=abs(1-x/x0);
        end
end