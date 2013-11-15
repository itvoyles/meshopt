function [J,djdx] = burgers_functional_calc(x)
global scaling

%setup executable path
% exe_path = '~/Software/burgers-mesh-adaption/library/bin';
% path(path,exe_path);

%writes grid to file
y = zeros(size(x));
z = zeros(size(x));
plot3d_write('grid.grd',x,y,z)

%calculates sensitivies
system('burgers1d-te.exe <burgers.inp');
% !burgers1d-te.exe <burgers.inp


%reads in sensitivies
[x,J,djdx,x_xi,dxxi_dx]=import_burgers1d_te_data('griddata.dat','functional.dat');

% load scalingdata.mat
% scaling = 0.01/J;
J=J*scaling;
djdx = djdx*scaling;

% fprintf('%12.4e %12.4e\n',J,max(abs(djdx)));
end