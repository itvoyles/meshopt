function x=gauss_seidel(A,b,x0,imax)


x = x0;
n = size(A,1);

nmax = 10000;
tol = 1e-10;

for k = 1:nmax
    y = x;
    
    %first iteration
    i = 1;
    index = [i+1:min(i+imax,n)];
    summation = b(i) - A(i,index)*x(index);
    x(i) = summation/A(i,i);
    
    %middle iterations
    for i = 2:n-1
%        for j = 1:i-1
%            sum = sum-A(i,j)*x(j);
%        end
        index = [max(1,i-imax):i-1];
        summation = b(i)-A(i,index)*x(index);
%        for j = i+1:n
%            sum = sum-A(i,j)*x(j);
%        end
       index = [i+1:min(i+imax,n)];
       summation = summation - A(i,index)*x(index);
       
       x(i) = summation/A(i,i);
    
    end

    i = n;
    index = [max(1,i-imax):i-1];
    summation = b(i)-A(i,index)*x(index);
    x(i) = summation/A(i,i);
    
    
    resid = A*x-b;
    if max(abs(resid)) < tol
        return
    end
end