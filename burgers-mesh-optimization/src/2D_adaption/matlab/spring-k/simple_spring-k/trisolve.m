function x = trisolve(LL,D,UU,B)


m = length(D);

for i= 1:m-1
   k = LL(i)/D(i);
%    LL(i) = LL(i)-k*D(i);
   D(i+1) = D(i+1)-k*UU(i);
   B(i+1) = B(i+1)-k*B(i);    
end

x = zeros(m,1);
x(m) = B(m)/D(m);
for col = m-1:-1:1
    x(col) = ( B(col) - UU(col).*x(col+1) )./D(col);  
end