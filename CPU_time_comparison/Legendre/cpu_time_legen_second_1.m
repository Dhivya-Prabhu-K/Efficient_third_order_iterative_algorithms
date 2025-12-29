function avg_time = cpu_time_legen_second_1(n, num_runs)

   times = zeros(1, num_runs);

    for run = 1:num_runs
        tic;
        fsweep_second_1(n);
        times(run) = toc;
    end
    avg_time = mean(times);

    
function x=fsweep_second_1(n)
% This code is to find all zeros of Legendre polynomial using SOM
     ind=n-2*floor(n/2);
     x = zeros(floor(n/2)+ind, 1);
     xi=tanh(pi/(2*(n+1)));
     nor=1.e-130;
    
    if ind>0
       y0=0;y1=nor;x(1)=0; i=1;
    else
       y0=nor;y1=0;i=0;
    end
    xt=0;
    t=xi;
    [y0,y1]=seriestay_1(n,xt,t,y0,y1);
  while i< 3
       [xt,x0,y0,y1]=fixedzser_second_1(n,xi,y0,y1);
       i=i+1;
       x(i)=xt;
       de=-tanh(pi/(2*(n+1)));
       xi=(xt-de)/(1-xt*de);
       t=xi-x0;
       [y0,y1]=seriestay_1(n,x0,t,y0,y1);
  end
  while i>=3 && i<floor(n/2)+ind
         [xt,x0,y0,y1]=fixedzser_second_1(n,xi,y0,y1);
         i=i+1;
         x(i)=xt;
          de=-(x(i-1)-x(i-2));
         xi=(x(i)*(1-x(i-1)*x(i-2))-de)/(1-x(i-1)*x(i-2)-x(i)*de);
         t=xi-x0;
        [y0,y1]=seriestay_1(n,x0,t,y0,y1);
  end


function [x,x0,y0,y1]=fixedzser_second_1(n,x0,y0,y1)
% SOM on z and using Taylor for (1-x)^(a+1)/2*(1+x)^(b+1)/2*P
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

function [y0,yp1]=seriestay_1(n,xp,t,y0,yp1)
% TAYLOR EN X para (1-x)^((a+1)/2)*(1+x)^((b+1)/2)*P(n,a,b,x)
% for the case of legendre a=b=0
 epso=1.0e-19;
 ym1=0;ym2=0; 
 delta=t;
 y0=y0*1.e-30;
 yp1=yp1*1.e-30;
 suma=y0+yp1*t;
 sumad=yp1;
 j=-1;
 errod=1; 
 xp2=xp*xp;
 
 L=2*n+1;L2=L*L;
 while (errod>epso)&&(j<100)
     j=j+1;
     c1=4*(j+2)*(j+1)*(1-xp2)^2;
     c2=-16*(j+1)*j*(1-xp2)*xp;
     c3=0.5*j*(j-1)*(48*xp2-16)+(L2-1)*(1-xp2)+(2)*(1-xp)+(2)*(1+xp);
     c4=16*(j-1)*(j-2)*xp-2*(L2-1)*xp;
     c5=4*(j-2)*(j-3)+(1-L2);
     yp2= -1/c1*(c2*yp1+c3*y0+c4*ym1+c5*ym2);
     ym2=ym1;   ym1=y0;     y0=yp1;     yp1=yp2;  
     ntd=(j+2)*yp2*t;
     sumad=sumad+ntd;
     t=t*delta;
     nt=yp2*t;
     suma=suma+nt;
     if (j>20) && (y0*yp1~=0)
        errod=max(abs((ntd/sumad)),abs((nt/suma)));
     end
 end
 y0=suma;yp1=sumad;
 y0=y0*1.e30;
 yp1=yp1*1.e30;


