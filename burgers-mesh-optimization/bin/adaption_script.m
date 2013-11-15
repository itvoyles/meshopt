function adaption_script
global scaling reynoldsnumber nu Lref alpha TEmax fidhistory
%setup executable path
exe_path = '~/Software/burgers-mesh-adaption/library/bin';
path(path,exe_path);
exe_path = '~/Software/burgers-mesh-adaption/library/src/burgers-1d/matlab';
path(path,exe_path);

tic

%% burgers problem setup
reynoldsnumber = 128;
alpha = reynoldsnumber/2;
Lref = 8;
uref = 2;
nu = uref*Lref/reynoldsnumber;

%options
interp = 0; %interpolate a coarse grid to fine grid for quicker convergence
nodesgoal = 129; %if interp = 1, number of nodes in interpolated grid (must be factor of 2)
run_initial_grid = 0; %runs solver on initial grid before adaption
scaling = 1000; %scaling factor for J and dJ/dx incase sensitivies are too low.

%small number limiting nearness of nodes
eps = 1e-5;


%% read grid
fprintf('Reading grid...');
% x0 = plot3d_read('grid.grd');
% to generate evenly spaced mesh --------
imax = 129;
x0 = linspace(-Lref/2,Lref/2,imax);
% x0(2:imax-1) = linspace(-1,1,imax-2);
% x0(1)=-Lref/2;
% x0(imax) = Lref/2;
%---------------------------------------
fprintf('Done\n');
nodes = length(x0);

fidhistory = fopen('history.dat','w');
fprintf(fidhistory,'Variables="x""x_xi"\n')

%% calculate truncation error on evenly spaced mesh for comparison
xeq = x0;
xeq = linspace(-Lref/2,Lref/2,nodes);
dxi = 1;
x_xi = dudx(xeq,dxi);
x_xixi = d2udx2(xeq,dxi);

TEmax = max(max(abs(Burgers_TE(xeq,x_xi,x_xixi,dxi,nu,alpha,Lref))));



%% interpolation section coarse-->fine. 
%USED IF STARTING FROM A COARSER MESH. MUST BE A FACTOR OF 2 DIFFERENCE( I.E.
%REFINEMENT FACTOR = 2)
%example: adapt mesh on 33 mesh points, read in 33 node mesh and add points to
%match nodesgoal value.

if interp == 1
if nodes ~= nodesgoal
    x0old = x0;
    n = nodes;
    while n < nodesgoal
        r = (nodesgoal-1)/(n-1);
        if mod(r,2) == 0
           x0new(1:2:(n-1)*2+1) = x0old;
           x0new(2:2:end-1) = (x0new(1:2:end-2)+x0new(3:2:end))/2;
           x0old = x0new;  
        else
          break
        end
        n = length(x0old);
    end
    x0 = x0new';
    nodes = length(x0);
else
    fprintf('Incorrect grid size!\n')
    return
end
end




%% Runs solver on input grid
if run_initial_grid == 1
    fprintf('Running solver...');
    burgers_solver(nodes,reynoldsnumber);
    fprintf('Done\n');
    system('cp output.dat orig_grid_output.dat')
end



%% Run optimization
fprintf('Adapting Mesh...\n');
% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------
dx = 1/(nodes-1);1e-8; % Typical Delta X
tolerance = 1e-30; % default 1e-6;
 
options = optimset('Algorithm','sqp',...
                   'GradObj','on',...
                   'MaxFunEvals',10^10,...
                   'MaxIter',100000,...
                   'FinDiffType','central',...
                   'RelLineSrchBnd',0.05,...
                   'TolFun',tolerance,...
                   'TolCon',tolerance,...
                   'TolX',tolerance,...
                   'PlotFcns',@dataoutput,...
                   'TypicalX',ones(size(x0))*dx);
%                    'UseParallel','always',...
options1 = optimset('Algorithm','sqp',...
                   'GradObj','off',...
                   'MaxFunEvals',10^10,...
                   'MaxIter',100000,...
                   'FinDiffType','central',...
                   'RelLineSrchBnd',0.05,...
                   'TolFun',tolerance,...
                   'TolCon',tolerance,...
                   'TolX',tolerance,...
                   'PlotFcns',@dataoutput,...
                   'UseParallel','always',...
                   'TypicalX',ones(size(x0))*dx);  

               
% sets up linear constraints
A1 = diag(ones(nodes,1),0);
A2 = diag(-ones(nodes-1,1),1);
A = A1+A2;
A(end,:) = 0;

b = -eps*ones(nodes,1);
b(end) = 0;

Aeq = zeros(nodes,nodes);
Aeq(1) = 1;
Aeq(end) = 1;

beq = zeros(nodes,1);
beq(1) = -4;
beq(end) = 4;



% -------------------------------------------------------------------------
% Optimization Call
% -------------------------------------------------------------------------               

exitflag = 0;

while exitflag == 0

[x,fval,exitflag,output,lambda,grad,hessian] = ...
 fmincon(@burgers1d_functional_calc,x0,A,b,Aeq,beq,[],[],[],options1);

x0 = x;

end

fprintf('Done');
toc



dxi = 1;
x_xi = dudx(x,dxi);
x_xixi = d2udx2(x,dxi);

TE = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref);
fid = fopen('TE.dat','w');
for i = 1:length(x)
  fprintf(fid,'%23.15e  %23.15e\n',x(i),TE(i)); 
end
fclose(fid);


%% runs solver on adapted mesh 

fprintf('Running solver...');
 burgers_solver(nodes,reynoldsnumber);
fprintf('Done\n');

%% Write new grid to file
y = zeros(size(x));
z = zeros(size(x));
fprintf('Writing Mesh...');
plot3d_write('grid.grd',x,y,z)
fprintf('Done\n');



function stop = dataoutput(x, optimValues, state)
global TEmax nu alpha Lref fidhistory

% subplot(2,1,1)
% semilogy(optimValues.iteration,optimValues.fval,'ok')

%print to screen ---------------------------------------------------------------
dxi = 1;
x_xi = dudx(x,dxi);
x_xixi = d2udx2(x,dxi);

TE = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref);
TEratio = max(max(abs(TE)))/TEmax;

if optimValues.iteration==1
    fprintf('\n Iter F-count            J(x)              TE\n')
    printquest = 1;    
end

if optimValues.iteration>=1
fprintf('%5.0f  %5.0f  %15.6e %15.6e\n',...
        optimValues.iteration,optimValues.funccount,optimValues.fval, 1/TEratio)
end
%-------------------------------------------------------------------------------

%write tecplot history
if mod(optimValues.iteration,5) == 0
    fprintf(fidhistory,'Zone t="%10.0f"\n',optimValues.iteration);
    fprintf(fidhistory,'i=%4.0f\n',length(x));
    for i = 1:length(x)
    fprintf(fidhistory,'%23.15e %23.15e \n',x(i),x_xi(i));
    end
end

%plot grid spacing
if mod(optimValues.iteration,1) == 0
plot((x(1:end-1)+x(2:end))/2,x(2:end)-x(1:end-1),'-ok')
end



% write to history file --------------------------------------------------------
fid = fopen('opt-hist-1d.out','a');
fprintf(fid,'%5.0f  %5.0f  %15.6e %15.6e\n',...
        optimValues.iteration,optimValues.funccount,optimValues.fval, 1/TEratio);
stop = false;
fclose(fid);


