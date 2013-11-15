function uxx = d2udx2(f,dx)

uxx = zeros(size(f));

imax = length(f);
i = 2:imax-1;
uxx(i) = (f(i+1)-2*f(i)+f(i-1))/(dx^2);

i=1;
uxx(i) = (2*f(i)-5*f(i+1)+4*f(i+2)-f(i+3))/(dx^2);

i=imax;
uxx(i) = ( -f(i-3)+4*f(i-2)-5*f(i-1)+2*f(i) )/(2*dx);




end