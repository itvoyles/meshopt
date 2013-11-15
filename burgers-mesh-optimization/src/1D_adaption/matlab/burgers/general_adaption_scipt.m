function general_adaption_scipt
global reynoldsnumber nu Lref alpha TEmax curiter curfeval fidhistory ...
       leftbndry rightbndry method freq_history forward_mapping TEtolerance

% clear
% clc
   
%% Setup executable path
working_dir = '/home/roycfd/Tyrone/burgers-mesh-optimization/';
exe_path = [working_dir,'library/bin'];
path(path,exe_path);
exe_path = [working_dir,'library/src/burgers-1d/matlab'];
path(path,exe_path);
exe_path = [working_dir,'/library/src/springsystem_force'];
path(path,exe_path);
exe_path = [working_dir,'/library/src/springsystem'];
path(path,exe_path);
%% burgers problem setup 

% input and general configuration options
restart = 0;            % restart optimization from saved config.mat file

% following controls are saved to config.mat file
imax = 129;                      % Number of grid nodes
reynoldsnumber = 32;            % Reynolds number of Burgers' equation
alpha = reynoldsnumber/2;       % Scaling factor for Burgers' equation
Lref = 8;                       % Reference domain size
uref = 2;                       % Reference velocity
nu = uref*Lref/reynoldsnumber;  % Viscosity calculation




method = 'spring-f';        % Choose adaption method:
                        % 'pure'     adaption with analytic TE calculations
                        % 'spring-k' spring adaption using spring constants
                        % 'spring-f' spring adaption using net force

% optimization options
maxiter = 1000;         % maximum number of adaption iterations
opt_restart = 1;        % restat optimization after every maxiter
gradobj = 'off';         % user supplied gradient (only for method='pure')
algorithm = 'sqp';      % optimization algorithm
findifftype = 'forward';% finite difference type if gradobj='off'

tolerance = 1e-40; % default 1e-6;
tolfun_tol = 1e-10;

func_tolerance = 1e-8; % change in max function min between outer iterations

% output
freq_history = 1;       % Frequency to write grid to history file

% small number limiting nearness of nodes
eps_small = 1e-6;




% following controls are not saved to config.mat file

readgrid = 0;           % read grid file 'grid.grd' = 1, create grid = 0
interp = 0;             % Interpolate a coarse grid to fine grid
nodesgoal = 129;        % Target fine grid

run_initial_grid = 0;   % Runs initial grid before adaption
halfdomain = 0;         % -1) Run negative half of domain
                        %  1) Run positive half of domain
                        %  0) Run full domain










if restart == 0
    save config.mat imax reynoldsnumber alpha Lref uref nu method ...
         freq_history eps_small maxiter opt_restart gradobj algorithm findifftype ...
         tolerance tolfun_tol
     
elseif restart == 1
    load config.mat
else
    error('Variable: restart must have value 0 or 1!')
end

% Grid input --------------------------------------------------------------
if readgrid == 1
    fprintf('Reading grid...');
    x0 = plot3d_read('grid.grd');
    if length(x0)~=imax && interp==0
        warning('Grid size does not match specified imax variable. Adjusting imax accordingly.')
    else
        fprintf('Done\n');
    end
    imax = length(x0);
else
    fprintf('Created evenly spaced grid.\n');
    x0 = linspace(-Lref/2,Lref/2,imax);
end
    x0 = reshape(x0,length(x0),1);
    


    
% grid interpolation ------------------------------------------------------

% USED IF STARTING FROM A COARSER MESH. MUST BE A FACTOR OF 2 DIFFERENCE 
% ( I.E. REFINEMENT FACTOR = 2)
% example: adapt mesh on 33 mesh points, read in 33 node mesh and add points to
% match nodesgoal value.
if interp == 1 && imax~=nodesgoal
    if mod(nodesgoal-1,imax-1)==0
        x0old = x0;
        n = imax;
        while n < nodesgoal
            r = (nodesgoal-1)/(n-1);
            if mod(r,2) == 0
               x0new(1:2:(n-1)*2+1) = x0old;
               x0new(2:2:end-1) = (x0new(1:2:end-2)+x0new(3:2:end))/2;
               x0old = x0new;  
            end
            n = length(x0old);
        end
        x0 = x0new';
        imax = length(x0);
    else
        error('Variable: nodesgoal must use a refinement factor divisible by 2!\n')

    end
end





% boundary conditions -----------------------------------------------------
if halfdomain == 0;
    leftbndry = -Lref/2;
    rightbndry = Lref/2;

elseif halfdomain == -1;
    leftbndry = -Lref/2;
    rightbndry = 0;
    imax = (imax-1)/2+1;
    x0 = x0(1:imax);
   
elseif halfdomain == 1;
    leftbndry= 0;
    rightbndry = Lref/2;
    imax = (imax-1)/2+1;
    x0 = x0(imax:end);

else
    error('Variable: Halfdomain must equal -1, 1, or 0!')
end

%% Initial calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start timer
tic 

% Calculate truncation error on equally spaced mesh
xeq = linspace(leftbndry,rightbndry,imax);
dxi = 1;
x_xi = dudx(xeq,dxi);
x_xixi = d2udx2(xeq,dxi);

TEmax = max(max(abs(Burgers_TE(xeq,x_xi,x_xixi,dxi,nu,alpha,Lref))));


% Setup history file to record adaption history
fidhistory = fopen('history.dat','w');
fprintf(fidhistory,'Variables="x""x_xi"\n');
% fprintf(fidhistory,'DATAPACKING="block"\n')



%% Runs solver on input grid
if run_initial_grid == 1
    fprintf('Running solver...');
    burgers_solver(imax,reynoldsnumber);
    fprintf('Done\n');
    system('cp output.dat orig_grid_output.dat')
end







%% Setup and define options for optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Setup method
% -------------------------------------------------------------------------
if strcmp(method,'pure')
    
    [A,b,Aeq,beq] = pure_linear_constraints(imax,eps_small,leftbndry,rightbndry);
    reverse_mapping = @(x)x;
    forward_mapping = @(x)x;
    
    
elseif strcmp(method,'spring-k')
    [A,b,Aeq,beq] = spring_k_linear_constraints(imax,eps_small);
    reverse_mapping = @spring_k_reverse_mapping;
    forward_mapping = @springsystem;
    
    
    
elseif strcmp(method,'spring-f')
    springsystem_force_setup(imax);
    [A,b,Aeq,beq] = spring_f_linear_constraints(imax,eps_small,leftbndry,rightbndry);
    reverse_mapping = @reversespringsystem_force;
    forward_mapping = @springsystem_force;    
    
else
    
    error('Variable: method is incorrectly defined!')
    
end



% -------------------------------------------------------------------------
% Initialize 
% -------------------------------------------------------------------------
DV = reverse_mapping(x0);

% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------

curiter = 0;
curfeval = 0;

options = optimset('Algorithm',algorithm,...
                   'GradObj',gradobj,...
                   'MaxFunEvals',10^10,...
                   'MaxIter',maxiter,...
                   'FinDiffType',findifftype,...
                   'RelLineSrchBnd',0.4,...
                   'TolFun',tolfun_tol,...
                   'TolCon',1e-10,...
                   'TolX',tolerance,...
                   'UseParallel','always',...
                   'PlotFcns',@dataoutput,...
                   'TypicalX',ones(size(DV,1),1)*1/(imax-1));






%% Call Optimization routine
fprintf('Adapting Mesh...\n');

% -------------------------------------------------------------------------
% Optimization Call
% -------------------------------------------------------------------------               

fvalold = 0;
exitflag = 0;
while exitflag == 0

    [DV,fval,exitflag,output,lambda,grad,hessian] = ...
     fmincon(@burgers1d_functional_calc,DV,A,b,Aeq,beq,[],[],[],options);

    
    
    
    if abs(fval-fvalold)>func_tolerance && exitflag ~= -1
        exitflag = 0;
    end
    fvalold = fval;
    
    if opt_restart == 0, exitflag = 1; end

end

% map design variable to grid nodes
x = forward_mapping(DV);

fprintf('Done\n');

%% Post processing

% if halfdomain double the double the domain
if halfdomain==-1 || halfdomain==1
    x = reshape(x,1,length(x));
    if halfdomain==-1
        xnew = [x,-fliplr(x(1:end-1))];
    elseif halfdomain==1
        xnew = [-fliplr(x(2:end)),x];
    end
    x = xnew;
    imax = (imax-1)*2+1;
end




% Calculate truncation error
dxi = 1;
x_xi = dudx(x,dxi);
x_xixi = d2udx2(x,dxi);
TE = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref);

% Write final truncation error to file
fid = fopen('TE.dat','w');
for i = 1:imax
  fprintf(fid,'%23.15e  %23.15e %23.15e\n',x(i),x_xi(i),TE(i)); 
end
fclose(fid);



%% Save info for final adapted grid using external fortran program and solves 
%Burgers' equation on the new mesh
% [J,djdx] = burgers_functional_calc(x);
% 
% fprintf('Running solver...');
%  burgers_solver(imax,reynoldsnumber);
% fprintf('Done\n');

%% Write new grid to file
% -------------------------------------------------------------------------
% Write new grid to file
% -------------------------------------------------------------------------



    
    


y = zeros(size(x));
z = zeros(size(x));
fprintf('Writing Mesh...');
plot3d_write('grid.grd',x,y,z)
fprintf('Done\n');

toc



%% Define output function
% -------------------------------------------------------------------------
% OPTIMIZATION ITERATION OUTPUT HISTORY FUNCTION
% -------------------------------------------------------------------------
function stop = dataoutput(DV, optimValues, state)
global TEmax nu alpha Lref curiter curfeval iter fidhistory forward_mapping ...
       freq_history TEtolerance maxTEold leftbndry rightbndry


% map design variable to grid nodes
x = forward_mapping(DV);
% x = springsystem_force(DV);


iter = optimValues.iteration;

% print to screen ---------------------------------------------------------

% Calculate truncation error
dxi = 1;
x_xi = dudx(x,dxi);
x_xixi = d2udx2(x,dxi);
TE = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref);

% Calculate ratio of max truncation error to max truncation on evenly
% spaced grid
TEratio = max(max(abs(TE)))/TEmax;

% If first iteration print header
if optimValues.iteration==1
    fprintf('\n Iter F-count            J(x)              TE\n')
    printquest = 1;    
end

% print values to screen
if optimValues.iteration>=1
fprintf('%5.0f  %5.0f  %15.6e %15.6e %15.6e\n',...
    optimValues.iteration,optimValues.funccount,optimValues.fval, 1/TEratio,abs(maxTEold-max(abs(TE))))
end



% write tecplot history ---------------------------------------------------

if mod(optimValues.iteration,freq_history) == 0
    fprintf(fidhistory,'Zone t="%10.0f"\n',optimValues.iteration);
    fprintf(fidhistory,'i=%4.0f\n',length(x));
    for i = 1:length(x)
    fprintf(fidhistory,'%23.15e %23.15e \n',x(i),x_xi(i));
    end
end




% plot mesh distribution --------------------------------------------------
plot(x,x_xi,'-ok')
axis([leftbndry,rightbndry,0,1.2*max(x_xi)])


% Record adaption history to file -----------------------------------------
fid = fopen('opt-hist-1d.out','a');
fprintf(fid,'%5.0f  %5.0f  %15.6e %15.6e\n',...
        curiter,curfeval,optimValues.fval, 1/TEratio);
curiter = curiter + 1;
curfeval = curfeval + 1;
stop = false;
fclose(fid);


% Check stopping tolerance based on dTE between iterations
if optimValues.iteration<1 || isempty(maxTEold)
    maxTEold = max(abs(TE));
   
else
    if abs(maxTEold-max(abs(TE)))<eps
%         stop = true;
    end
    maxTEold = max(abs(TE));
end



%% Define adaption method dependent functions
% Must define:
% Method to return linear constraints [A,b,Aeq,Beq] 
% Forward mapping function if neccisary to transform design variable to grid nodes
% Reverse mapping function if neccisary to transform grid nodes to design
%   variables

% If no mapping functions are defined then define a blank function

%% pure functions
    function [A,b,Aeq,beq] = pure_linear_constraints(imax,eps_small,leftbndry,rightbndry)
        
        % sets up linear inequality constraints
        % specifies that x_(i+1) > x_(i)+eps_small using the linear system Ax<=b,
        % specifies that 
        A1 = diag(ones(imax,1),0);
        A2 = diag(-ones(imax-1,1),1);
        A = A1+A2;
        A(end,:) = 0;

        b = -eps_small*ones(imax,1);
        b(end) = 0;

        % setups up linear equality constraints
        Aeq = zeros(2,imax);
        Aeq(1) = 1;
        Aeq(2,end) = 1;

        beq = zeros(2,1);
        beq(1) = leftbndry;
        beq(2) = rightbndry;


%% spring-k functions

    function [A,b,Aeq,beq] = spring_k_linear_constraints(imax,eps_small)

        % Linear inequality constraints
        % Ak <= b
        % defines all spring constants must be greater than small number (eps_small)
        A = -diag(ones(imax-1,1),0);
        b = -eps_small*ones(imax-1,1);

        % Linear equality constraints
        % Ak = b
        Aeq = [];
        beq = [];

    function k = spring_k_reverse_mapping(x)
    % This functions defines the energy in the spring system equal to the
    % energy in an equally distributed spring system in order to define a
    % unique value for spring constants.
        
        % define equally spaced system and calculate spring constants
        xeq = linspace(x(1),x(end),length(x));
        k0 = reversespringsystem(xeq,0);
        
        % check that the size of x is a column vector
        x = reshape(x,length(x),1);

        % calculate energy of the equally spaced spring system
        dx0 = x(2:end)-x(1:end-1);
        E0 = sum(1/2*k0.*dx0.^2);
        
        % calculate spring constants for input grid 
        k = reversespringsystem(x,E0);
        
    % k = spring_k_forward_mapping is defined in separate file
    
    
%% spring-f functions
    function [A,b,Aeq,beq] = spring_f_linear_constraints(imax,eps_small,leftbndry,rightbndry)
    global SSAinv 

        % linear inequality constraints
        % specifies that x_(i+1) > x_(i)+eps_small using the linear system KF=x,
        % where SSAinv = k^(-1): x_(i)-x_(i+1) = (K_(i) - K_(i+1))F <= -eps_small
        A = -(SSAinv(2:end,:)-SSAinv(1:end-1,:));
        b = -eps_small*ones(imax-1,1);

        % linear equality constraints
        % specifies that x_(L) = right boundary and x_(1) = left boundary
        Aeq = zeros(2,imax);
        Aeq(1,1) = 1;
        Aeq(2,end)=1;

        beq = zeros(2,1);
        beq(1) = leftbndry;
        beq(2) = rightbndry;

    