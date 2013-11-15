function ux = dudx(f,dx)

ux = zeros(size(f));

imax = length(f);
i = 2:imax-1;
ux(i) = (f(i+1)-f(i-1))/(2*dx);

i=1;
ux(i) = (-3*f(i)+4*f(i+1)-f(i+2))/(2*dx);

i=imax;
ux(i) = (f(i-2)-4*f(i-1)+3*f(i))/(2*dx);




end