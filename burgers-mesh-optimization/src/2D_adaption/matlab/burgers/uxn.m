%burgers1d functions

function [uout] = uxn(x,y,RE,Lref,flag)
%   !up to 5th derivative of burgers' exact solution
%   !tested using finite difference, derivative 1-5 are verified accurate.

%     ! f = c1*tanh(c2x)
    c1 = -2;
    c2 = RE/(2*Lref);

    f = c1*tanh(c2*(x+y));
    if (flag==0)
      uout=f;
      return
    end

    f1 = -c2/c1*(f.^2-c1.^2);
    if (flag==1)
      uout=f1;
      return
    end

    f2 = -c2/c1*(2.*f.*f1);
    if (flag==2)
      uout=f2;
      return
    end

    f3 = -2*c2/c1*(f1.*f1+f.*f2);
    if (flag==3)
      uout=f3;
      return
    end

    f4 = -2*c2/c1*(3*f1.*f2+f.*f3);
    if (flag==4)
      uout=f4;
      return
    end
    
    f5 = -2*c2/c1*(3*(f2.^2+f1.*f3)+f1.*f3+f.*f4);
    if (flag==5) 
      uout=f5;
      return
    end


    uout=0;



end