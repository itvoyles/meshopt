function [intF Area] = trap_sum2d(x,y,f)

  %setup indexing
  [imax,jmax] = size(x);
  i = 1:imax-1;
  j = 1:jmax-1;
  
  %vectorize cell boundaries
  v1(i,j,1) = x(i+1,j)-x(i,j);
  v1(i,j,2) = y(i+1,j)-y(i,j);

  v2(i,j,1) = x(i,j+1)-x(i,j);
  v2(i,j,2) = y(i,j+1)-y(i,j);

  v3(i,j,1) = x(i,j+1)-x(i+1,j+1);
  v3(i,j,2) = y(i,j+1)-y(i+1,j+1);

  v4(i,j,1) = x(i+1,j)-x(i+1,j+1);
  v4(i,j,2) = y(i+1,j)-y(i+1,j+1);

  %calculates area by calculating the area of two half-quads 
  %A = (v1 x v2)/2+(v3 x v4)/2  [v1 x v2 => crossproduct]

  Area = 0.5.*abs(v1(i,j,1).*v2(i,j,2)-v1(i,j,2).*v2(i,j,1)) ...
       + 0.5.*abs(v3(i,j,1).*v4(i,j,2)-v3(i,j,2).*v4(i,j,1));

  fave = (f(i,j)+f(i+1,j+1)+f(i,j+1)+f(i+1,j))./4.;
  intF= sum(sum(Area.*fave));


end
