function [ dadx_cc ] = dadx( x )
%Function "dadx"
%   Takes cell-center locations, returns the area variations at
%    those positions.

n = length(x);
dadx_cc = zeros(n,1);

for i = 1:n,
    
    if abs(x(i)) <= 1,
        
        dadx_cc(i) = 0.4*pi*cos( pi*( x(i) - 0.5 ) );        
        
    elseif abs(x(i))>1 && abs(x(i))<=2,
        
        dadx_cc(i) = 0.0;
        
    end
    
end

end


