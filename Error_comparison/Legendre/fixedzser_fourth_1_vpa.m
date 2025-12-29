function [x,x0,y0,y1]=fixedzser_fourth_1_vpa(n,x0,y0,y1)
% FOM on z and using Taylor for (1-x)^(1)/2*(1+x)^(1)/2*P 
    digits(32);
    L=vpa(2*n+1);
    epsil=vpa(10^(-16));
    sq=0.25*((L*L-1)*(1-x0*x0));
    sqde=sqrt(abs(sq));
    h=y0/(y1*(1-x0*x0)+x0*y0);
   
    argu=vpa(atan(sqde*h)/sqde);       
   de=vpa(tanh(argu));
   x=(x0-de)/(1-x0*de);   
   erro=abs(1-x/x0);

    while erro>epsil
        t=x-x0;
        [y0,y1]=seriestay_1(n,x0,t,y0,y1);
        x0=x;
        sq=0.25*((L*L-1)*(1-x0*x0));
        sqde=sqrt(abs(sq));
        h=y0/(y1*(1-x0*x0)+x0*y0);
        
        argu=vpa(atan(sqde*h)/sqde);        
        de=vpa(tanh(argu));
        x=(x0-de)/(1-x0*de);   
        erro=abs(1-x/x0);
    end
end


