function [x,x0,y0,y1]=fixedzser_fourth_1(n,x0,y0,y1)
% fixed point on z and using Taylor for (1-x)^(a+1)/2*(1+x)^(b+1)/2*P

    L=2*n+1;
    epsil=10^(-15);
    sq=0.25*((L*L-1)*(1-x0*x0));
    sqde=sqrt(abs(sq));
    h=y0/(y1*(1-x0*x0)+x0*y0);
   
    argu=atan(sqde*h)/sqde;    
    if ~h<0
       argu=argu-pi/sqde; 
   end    
   de=tanh(argu);
   x=(x0-de)/(1-x0*de);   
   erro=abs(1-x/x0);

    while erro>epsil
        t=x-x0;
        [y0,y1]=seriestay_1(n,x0,t,y0,y1);
        x0=x;
        sq=0.25*((L*L-1)*(1-x0*x0));
        sqde=sqrt(abs(sq));
        h=y0/(y1*(1-x0*x0)+x0*y0);
        
        argu=atan(sqde*h)/sqde;       
        de=tanh(argu);
        x=(x0-de)/(1-x0*de);   
        erro=abs(1-x/x0);
    end
end

