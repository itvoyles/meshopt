function x = pentasolve(A,b,imax)

    %diagonalize A matrix and record reduction history
    m = size(A,1);
    Aold = A;
    Anew = A;
    B = b;
    
    
for col = 1:m-1

    rowmat = [col:min(imax+col,m)];
    
%     rowmat = max(1,rowmat);
    for r = 2:length(rowmat)

        row = rowmat(r);
            k = Aold(row,col)/Aold(col,col);
            Aold(row,rowmat) = Aold(row,rowmat)-k*Aold(col,rowmat);
            B(row) = B(row)-k*B(col);
    end
%     Aold = Anew;

end

% opts.UT = true;
% xx = linsolve(Anew,B,opts);
% xx1 = A\b;


x(m) = B(m)/Anew(m,m);
for col = m-1:-1:1
    rowmat = [col+1:min(m,imax+col)];
    
    x(col) = ( B(col) - sum(Aold(col,rowmat).*x(rowmat)) )/Aold(col,col);  
%     x(rowmat(1)) = x(rowmat(1))/Anew(col,col);
    
    
end



