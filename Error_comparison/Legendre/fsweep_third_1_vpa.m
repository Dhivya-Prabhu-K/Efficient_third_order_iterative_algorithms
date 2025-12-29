function [x]=fsweep_third_1_vpa(n)
% This code is to find all zeros of Legendre polynomial using TOM
    digits(32);
    ind=n-2*floor(n/2);
    x = sym(zeros(floor(n/2)+ind, 1));
    xi=vpa(tanh(pi/(2*(n+1))));
    nor=1.e-130;
        if ind>0
            y0=0;y1=nor;x(1)=0; i=1;
        else
            y0=nor;y1=0;i=0;
        end
    xt=vpa(0);
    t=xi;
    [y0,y1]=seriestay_1_vpa(n,xt,t,y0,y1);
      while i< 3
            [xt,x0,y0,y1]=fixedzser_third_1_vpa(n,xi,y0,y1);
            i=i+1;
            x(i)=xt;
            de=vpa(-tanh(pi/(2*(n+1))));
            xi=(xt-de)/(1-xt*de);
            t=xi-x0;
            [y0,y1]=seriestay_1_vpa(n,x0,t,y0,y1);
      end
     while i>=3 && i<floor(n/2)+ind
            [xt,x0,y0,y1]=fixedzser_third_1_vpa(n,xi,y0,y1);
            i=i+1;
            x(i)=xt;
            de=vpa(-(x(i-1)-x(i-2)));
            xi=(x(i)*(1-x(i-1)*x(i-2))-de)/(1-x(i-1)*x(i-2)-x(i)*de);
            t=xi-x0;
            [y0,y1]=seriestay_1_vpa(n,x0,t,y0,y1);
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%function fixedzser
   