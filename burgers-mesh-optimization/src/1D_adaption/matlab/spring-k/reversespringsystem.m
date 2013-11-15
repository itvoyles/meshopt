function k = reversespringsystem(x,E)


imax = length(x)-1;

A = zeros(imax,imax);
b = zeros(imax,1);

dx = x(2:end)-x(1:end-1);
for ii = 1:imax-1
   A(ii,ii) = dx(ii);
   A(ii,ii+1) = -dx(ii+1);
   b(ii) = 0;
end


if (E==0) || (nargin == 1);
    A(imax,imax) = 1;
    b(imax) = 1;
else

    %specify total energy in spring system
    A(imax,:) = 1/2*dx.^2;
    b(imax) = E;

    %specify sum of spring constants
    % A(imax,:) = 1;
    % b(imax) = E;

end





k = A\b;

% k = abs(k)./min(abs(k));


