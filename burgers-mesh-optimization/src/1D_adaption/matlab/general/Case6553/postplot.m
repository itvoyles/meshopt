function postplot
global fin
% Function "postplot.m"
%   Reads from output files to create plots for examination.

[solu] = load('exact.dat'); % solu = [rhoc uc pc]

[xdxf] = load('faces.dat'); % xdxf = [x dxf]

x = xdxf(:,1); dxf = xdxf(:,2);

[TEuni] = load('TEuniform.dat'); % TEuni = [xcu dxu TE1u TE2u TE3u]

xcu = TEuni(:,1); dxu = TEuni(:,2); TE1u = TEuni(:,3);
TE2u = TEuni(:,4); TE3u = TEuni(:,5);

[TEfin] = load('TE.dat'); % TEfin = [xc dx TE1 TE2 TE3];

[~, hold1] = system('head -n 1 opt-hist-1d.out');
[~, hold2] = system('tail -n 1 opt-hist-1d.out');
fid5 = fopen('temp.dat','w');
fprintf(fid5,hold1);
fprintf(fid5,hold2);
[fin] = load('temp.dat');

xc = TEfin(:,1); dx = TEfin(:,2); TE1 = TEfin(:,3);
TE2 = TEfin(:,4); TE3 = TEfin(:,5);

figure(1)
subplot(2,2,1), plot(x,dxf,'--*r',xc,dx,'-ok')
xlabel('x-location')
ylabel('Cell-size (dx)')
h1 = legend('Faces','Cell-Centers');
set(h1,'FontSize',5)

subplot(2,2,2), plot(xcu,TE1u,'--*',xc,TE1,'-o')
xlabel('x-location')
ylabel('Cont TE')
h2 = legend('Uniform','Optimized','Location','SouthWest');
set(h2,'FontSize',5)

subplot(2,2,3), plot(xcu,TE2u,'--*',xc,TE2,'-o')
xlabel('x-location')
ylabel('X-Mtm TE')
% h3 = legend('Uniform','Optimized');
% set(h3,'FontSize',5)

subplot(2,2,4), plot(xcu,TE3u,'--*',xc,TE3,'-o')
xlabel('x-location')
ylabel('Energy TE')
% h4 = legend('Uniform','Optimized');
% set(h4,'FontSize',5)

figure(2)
subplot(2,2,1), plot(xc,solu(:,1))
subplot(2,2,2), plot(xc,solu(:,2))
subplot(2,2,3), plot(xc,solu(:,3))

end