global Lref iter

Lref = 8;
iter = 0;
hold off
x = linspace(-4,4,10)';
xcold = (x(2:end)+x(1:end-1))/2;
for i = 1:10
k = -xcold.^2+17;
x = springsystem(k);
xc = (x(2:end)+x(1:end-1))/2;
max(xc-xcold)
xcold = xc;
plot(xc,k,'o-')
hold on
end
% E = sum(1/2*(x(2:end)-x(1:end-1).^2))
% ktest = reversespringsystem(x,E);
% max(k - ktest')


figure(1)
plot(x,k,'o-');
dx = (x(2:end)-x(1:end-1));
xc = (x(2:end)+x(1:end-1))/2;
figure(2)
plot(xc,dx,'o-');