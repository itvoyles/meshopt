function [y]=gridmake(imax,c)

% imax = 65;
% c = 1.5;

x = linspace(-1,1,imax);

y = 2*tan(c*x)/max(tan(c*x));

% dy = y(2:end) - y(1:end-1);

% figure(1)
% subplot(2,1,1), plot(x,y,'-*');
% subplot(2,1,2), plot(x(1:end-1),dy,'-*');
end