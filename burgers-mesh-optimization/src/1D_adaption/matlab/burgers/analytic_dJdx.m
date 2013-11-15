% !  
% !      N-1  
% !      --
% !  J = \  f(i)
% !      /
% !      --
% !      i=1
% !
% !          1 /                     \/               \
% !  f(i) = ---| TE(i+1)^2 + TE(i)^2 || x(i+1) - x(i) |
% !          2 \                     /\               /
% !
% !  dJ     -1  /                       \ 
% ! ----- = --- | TE(i+1)^2 - TE(i-1)^2 | 
% ! dx(i)    2  \                       /
% !
% !          /         dTE(i+1)        dTE(i)  \/             \
% !        + | TE(i+1) --------- + TE(i)------ || x(i+1)-x(i) |
% !          \           dx              dx    /\             /
% !
% !          /        dTE(i)           dTE(i-1) \/             \
% !        + | TE(i) -------- + TE(i-1)-------- || x(i)-x(i-1) |
% !          \          dx                dx    /\             /
% !
% !          /          dTE(i+1)  \/               \
% !        + | TE(i+1) ---------- || x(i+2)-x(i+1) |
% !          \            dx      /\               /
% !
% !          /          dTE(i-1)  \/               \
% !        + | TE(i-1) ---------- || x(i-1)-x(i-2) |
% !          \            dx      /\               /


function djdx = analytic_dJdx(x,TE,dTEi_dx,dTEip1_dx,dTEim1_dx,imax)


pdf_pdx = zeros(imax-2);


  for i = 2:imax-1

    pdf_pdx(i) =-(TE(i+1).^2-TE(i-1).^2)/2 ...
                +(TE(i+1)*dTEip1_dx(i)+TE(i)*dTEi_dx(i))*(x(i+1)-x(i))...
                +(TE(i)*dTEi_dx(i)+TE(i-1)*dTEim1_dx(i))*(x(i)-x(i-1)) ; 

    if (i>2)
      pdf_pdx(i) = pdf_pdx(i)+TE(i-1)*dTEim1_dx(i)*(x(i-1)-x(i-2));
    end
    
    if (i<imax-1)
      pdf_pdx(i) = pdf_pdx(i)+TE(i+1)*dTEip1_dx(i)*(x(i+2)-x(i+1));
    end

  end

  djdx = zeros(imax,1);
  djdx(2:imax-1) = pdf_pdx(2:imax-1);


end 