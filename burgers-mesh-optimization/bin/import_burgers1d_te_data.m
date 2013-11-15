function [x,J,djdx,x_xi,dxxi_dx]=import_burgers1d_te_data(filename1,filename2)

filename = filename1;
functional=filename2;


% fid = fopen(functional,'r');
J = dlmread(functional);
% J = importdata(functional);



% data = importdata(filename);
data = dlmread(filename);
imax = data(1);
datamat = data(2:end,:);
% datavec = data(1:end);

% datamat = reshape(datavec,6,imax)';

x = datamat(:,1);
djdx = datamat(:,2);
x_xi = datamat(:,3);

dxxi_dx1 = diag(datamat(2:end,4),-1);
dxxi_dx2 = diag(datamat(1:end,5),0);
dxxi_dx3 = diag(datamat(1:end-1,6),1);

dxxi_dx = dxxi_dx1+dxxi_dx2+dxxi_dx3;
dxxi_dx(1,:) = 0;
dxxi_dx(end,:) = 0;

dxxi_dx = sparse(dxxi_dx);

end