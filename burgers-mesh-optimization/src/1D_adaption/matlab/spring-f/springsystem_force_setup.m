function springsystem_force_setup(param)
global SSAinv SSA Lref cells_hold

% param=imax-cells_hold+1;
%create matrix
Ad = diag(2*ones(param,1),0);
Au = -diag(ones(param-1,1),1);
Al = -diag(ones(param-1,1),-1);
A = Ad+Au+Al;

%setup boundary conditions
A(1,:) = 0;
A(1,1) = 1;
A(param,:) = 0;
A(param,param) = 1;

SSA = A;

%invert
SSAinv = A^-1;
