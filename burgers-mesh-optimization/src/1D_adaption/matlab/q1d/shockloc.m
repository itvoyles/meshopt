function [ xs, As2 ] = shockloc( pb,p01)
%Function "shockloc"
%   Iterates to find the shock location. This function modifies the
%     value of A*_2 [ As2 ]
global astar toler gmma gm1 gp1 xthroat rightbndry


Ms = 1.5; % Initial guess for pre-shock Mach number
syms g real % Symbolic variable for symbolic differentiation
dfcn = matlabFunction(diff(funcloc(g,pb,p01)));
% keyboard
count = 0;
fcn = funcloc(Ms,pb,p01);

while fcn > toler && count < 500,
    
    delmach = -fcn/dfcn(Ms);
    
    Ms = Ms + delmach;
    
    fcn = funcloc(Ms,pb,p01);
    count = count + 1;
end

% for j=1:500,
%     
%     fcn = funcloc(Ms,pb,p01);
%     
%     if abs(fcn) < toler,
%         break
%     end
%     
%     delmach = -fcn/dfcn(Ms);
%     
%     Ms = Ms + delmach;
%     
% end
Ms = abs(Ms);

% Use Ms to determine the shock area, As:

aratshk = (1/Ms) * ( 2*(1 + Ms^2 * gm1/2)/gp1)^(-1/2 + gmma/gm1); % Pre-shock area ratio

As = aratshk * astar(1);

% Use the shock area to determine the shock location:

% xL = x(k);
xL = xthroat;
xR = rightbndry;


% Bisection method:

% for j=1:1000,
%     
%     xs = (xL - xR)/2;
%     Ag = geom(xs); % Guessed area
%     
%     Adiff = Ag - As;
%     
%     if abs(Adiff) < toler,
%         break
%     elseif Adiff > 0,
%         xR = xs;
%     elseif Adiff < 0,
%         xL = xs;
%     end
%     
%     xs = (xR - xL)/2;   
%     
% end
count = 0;
xs = (xL - xR)/2;
Ag = geom(xs); % Guessed area

Adiff = Ag - As;

while abs(Adiff) > toler && count < 1000,
    
    if Adiff > 0,
        xR = xs;
    elseif Adiff < 0,
        xL = xs;
    end
    
    xs = (xL - xR)/2;
    Ag = geom(xs); % Guessed area
    
    Adiff = Ag - As;
    count = count + 1;
end

xs = abs(xs);
% Calculate downstream A*: A2*
As2 = astar(1) * (1 + 2*gmma*(Ms^2 - 1)/gp1)^(1/gm1) ...
    * ( (2+gm1*Ms^2)/(gp1*Ms^2) )^(gmma/gm1);

end
