function general_adaption_scipt
global reynoldsnumber nu Lref alpha TEmax curiter curfeval fidhistory ...
       leftbndry rightbndry method freq_history forward_mapping TEtolerance
global k xthroat toler astar tzero clls imax eqn
global pb pzero sol order flux_type lim_type kappa shloc As2 flag functional
global gmma gm1 gp1 xg xgm1 xgp1 xg2m1 gxgm1 gxg2m1 gm1xgp1 gp1xgm1 cv cp

% clear
% clc
   
%% Setup executable path
%working_dir = 'C:\Users\cjroy\Dropbox\AFRL-Summer-Visits\2013\Repository\burgers-mesh-optimization\';
working_dir = 'C:\Users\Limited\Documents\Mesh Optimization\Repository\burgers-mesh-optimization\'; 
% working_dir = 'Y:\2013-Voyles\Q1D Matlab 072513\Repository\burgers-mesh-optimization\';
% working_dir = 'C:\Users\alyanaej\Documents\MeshOptimization\Repository\burgers-mesh-optimization\';
exe_path = [working_dir,'bin'];
path(path,exe_path);
exe_path = [working_dir,'src/1D_adaption/matlab/burgers'];
path(path,exe_path);
exe_path = [working_dir,'src/1D_adaption/matlab/q1d'];
path(path,exe_path);
exe_path = [working_dir,'src/1D_adaption/matlab/spring-f'];
path(path,exe_path);
exe_path = [working_dir,'src/1D_adaption/matlab/spring-k'];
path(path,exe_path);
%% Problem setup 

% input and general configuration options
restart = 0;            % restart optimization from saved config.mat file

% following controls are saved to config.mat file
if numel(imax) ~= 0,
    disp('Using predefined imax')
    disp(imax)
else
    imax = 65;
end

% Number of grid nodes/ FVM faces (choose odd)
eqn = 'q1d';                       % Case used:
                                   %    'q1d'     = Quasi-1D Euler equations
                                   %    'burgers' = Burgers' equation

if strcmp(eqn,'burgers'),
    % Burgers' Equation:
    reynoldsnumber = 32; %16;            % Reynolds number of Burgers' equation
    alpha = reynoldsnumber/2;       % Scaling factor for Burgers' equation
    Lref = 8;                       % Reference domain size
    uref = 2;                       % Reference velocity
    nu = uref*Lref/reynoldsnumber;  % Viscosity calculation
    
elseif strcmp(eqn,'q1d'),
    % Q1D Euler Equations:
    Lbnd = -2;  % Domain boundary: 'Lbnd' is the inflow
    Rbnd = 2;   % Domain boundary: 'Rbnd' is the outflow
    
    if numel(sol) ~= 0,
        disp('Using predefined sol var')
        disp(sol)
    else
        sol = 'sub'; % Specify exact solution type: 'sub' -- Isentropic Subsonic-subsonic
                       %                              'super' -- Isentropic Subsonic-supersonic
                       %                              'shock' -- Internal shock
    end 
    
    
    if numel(flux_type) ~= 0,
        disp('Using predefined flux_type var')
        disp(flux_type)
    else
        flux_type = 'vanleer'; % Available:  'roe'     = Roe FDS
    end              %             'vanleer' = Van Leer FVS
    
    if numel(lim_type) ~= 0,
        disp('Using predefined lim_type var')
        disp(lim_type)
    else                  
    lim_type = 'vanleer'; % Available: 'ospre','minmod','superbee','sweby','vanalbada','vanleer'              
    end
    
    kappa = -1; % MUSCL constant
    
    if numel(functional) ~= 0,
        disp('Using predefined functional var')
        disp(functional)
    else 
        functional = 'cont'; % Type of functional used: 'cont'  : continuity
                             %                          'xmtm'  : x-momentum
                             %                          'energy': energy
                             %                          'rss'   : root sum square
    end
    
    pb = 297000; % Pa, Back pressure (sub: 297330, shock : <sub, 297000)
    
    switch sol
        case 'sub'
            pb = 297330;
        case 'super'
            pb = 293000; % not used
        case 'shock'
            pb = 297000;
    end
    
    tzero = 300.0; % K, Stagnation temperature
    pzero = 300000.0; % Pa, Chamber pressure
    order = 3; % Order of the quadrature
    toler = 1E-12; % Iterative tolerance
    
    % Geometry related:
    k = ceil( (imax+1)/2); % Index of throat (symmetric nozzle, uniform grid)
    x0 = linspace(Lbnd,Rbnd,imax); % Vector of uniformly-spaced face locations
    xthroat = x0(k); % location of throat
    clear x0 % Clear the above x0 definition
    clls = imax - 1; % Number of cells is one fewer than the number of faces
    As2 = -99.9; % Unused unless sol = 'shock' 
    shloc = -99.9; % Unused unless sol = 'shock'
    flag = 0; % Activates the use of two A* values (sol = 'shock' only)
    
    % Gamma-related constants:
    gmma = 1.4; % Ratio of specific heats
    R = 286.9; % Metric gas constant for air
    gm1  = gmma - 1;
    gp1  = gmma + 1;
    xg   = 1/gmma;
    xgm1 = 1/(gmma - 1);
    xgp1 = 1/(gmma + 1);
    xg2m1 = 1/(gmma*gmma - 1);
    gxgm1 = gmma / (gmma - 1);
    gxg2m1 = gmma / (gmma*gmma - 1);
    gm1xgp1 = gm1/gp1;
    gp1xgm1 = gp1/gm1;
    
    cv = R*xgm1;
    cp = R*gxgm1;
    
    % Other
    astar = zeros(2,1);
    
else
    error('Enter a valid case for "eqn"')
end


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

tolerance = 1e-6; % default 1e-6;
tolfun_tol = 1e-10;

func_tolerance = 1e-12; % change in max function min between outer iterations

% output
freq_history = 1;       % Frequency to write grid to history file

% small number limiting nearness of nodes
eps_small = 1e-6;



% following controls are not saved to config.mat file

readgrid = 1;           % read grid file 'grid.grd' = 1, create grid = 0
interp = 0;             % Interpolate a coarse grid to fine grid
nodesgoal = 129;        % Target fine grid

run_initial_grid = 0;   % Runs initial grid before adaption
halfdomain = 0;         % -1) Run negative half of domain
                        %  1) Run positive half of domain
                        %  0) Run full domain
                        %  - Always run full domain for Q1D equations










if restart == 0
    
    if strcmp(eqn,'burgers'),
        save config.mat imax reynoldsnumber alpha Lref uref nu method ...
            freq_history eps_small maxiter opt_restart gradobj algorithm findifftype ...
            tolerance tolfun_tol
    elseif strcmp(eqn,'q1d'),
        save config.mat imax Lbnd Rbnd sol pb pzero tzero order toler method lim_type flux_type...
            freq_history eps_small maxiter opt_restart gradobj algorithm findifftype ...
            tolerance tolfun_tol
    end
    
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
    
    if strcmp(eqn,'burgers')
        x0 = linspace(-Lref/2,Lref/2,imax);
    elseif strcmp(eqn,'q1d')
        x0 = linspace(Lbnd,Rbnd,imax);
    end
    
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

    if (exist('Lbnd','var') && exist('Rbnd','var')),
        leftbndry = Lbnd;   % For the Q1D case
        rightbndry = Rbnd;
    else
        leftbndry = -Lref/2;
        rightbndry = Lref/2;
    end
        

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

%% Q1D Preliminary Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(eqn,'q1d') && strcmp(sol,'shock')
    
    flag = 1; % Activate the use of two A* values
    astar(:) = geom(xthroat); % A* for upstream, A* = A_t
    
    [ shloc, As2 ] = shockloc(pb,pzero); % Calculate shock location
                                            % and downstream A*
                                         
end
    
fprintf('%7.5f \n',shloc)

%% Initial calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start timer
tic 

% Calculate truncation error on equally spaced mesh
xeq = linspace(leftbndry,rightbndry,imax);

if strcmp(eqn,'burgers'),
    dxi = 1;
    x_xi = dudx(xeq,dxi);
    x_xixi = d2udx2(xeq,dxi);
    
    TEmax = max(max(abs(Burgers_TE(xeq,x_xi,x_xixi,dxi,nu,alpha,Lref))));
    
    
    % Setup history file to record adaption history
    fidhistory = fopen('history.dat','w');
    fprintf(fidhistory,'Variables="x""x_xi"\n');
    % fprintf(fidhistory,'DATAPACKING="block"\n')
    
elseif strcmp(eqn,'q1d'),
    xc = 0.5 .* ( xeq(2:end) + xeq(1:end-1) );
    dx = xeq(2:end) - xeq(1:end-1);
    [TE] = q1d_TE( imax,xeq,pb,pzero,sol,order,kappa,flag,flux_type,lim_type,shloc,As2 );
    fiduniform = fopen('TEuniform.dat','w');
    for i = 1:imax-1,
        fprintf(fiduniform,'%23.15e  %23.15e %23.15e %23.15e %23.15e\n',xc(i),dx(i),TE(i+1,1),TE(i+1,2),TE(i+1,3) );
    end
    fclose(fiduniform);
    
    TEmax = max( abs( q1d_TE( imax,xeq,pb,pzero,sol,order,kappa,flag,flux_type,lim_type,shloc,As2 ) ) ); 
                         % Currently returns max TE in each
                         % equation, [cont x-mtm energy]
                         
    % Set up history file to record adaptation history
    fidhistory = fopen('history.dat','w');
    fprintf(fidhistory,'Variables="x""dx"\n');
    % FIXME: 
    %        -> appropriate history file
    
end


%% Runs solver on input grid
if run_initial_grid == 1
    fprintf('Running solver...');
    burgers_solver(imax,reynoldsnumber);
    fprintf('Done\n');
    system('cp output.dat orig_grid_output.dat')
end


% FIXME: need to add this capability for Q1D equations




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

if strcmp(eqn,'burgers'),
    
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
    
elseif strcmp(eqn,'q1d'),
    
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
        'PlotFcns',@dataoutq1d,...
        'TypicalX',ones(size(DV,1),1)*1/(imax-1));
    
end


%% Call Optimization routine
fprintf('Adapting Mesh...\n');

% -------------------------------------------------------------------------
% Optimization Call
% -------------------------------------------------------------------------               

fvalold = 0;
exitflag = 0;

if strcmp(eqn,'burgers'),
    
    while exitflag == 0
        
        [DV,fval,exitflag,output,lambda,grad,hessian] = ...
            fmincon(@burgers1d_functional_calc,DV,A,b,Aeq,beq,[],[],[],options);
           
        if abs(fval-fvalold)>func_tolerance && exitflag ~= -1
            exitflag = 0;
        end
        fvalold = fval;
        
        if opt_restart == 0, exitflag = 1; end
        
    end
    
elseif strcmp(eqn,'q1d'),
    
    while exitflag == 0
        
        [DV,fval,exitflag,output,lambda,grad,hessian] = ...
            fmincon(@q1d_functional_calc,DV,A,b,Aeq,beq,[],[],[],options);
        
      
        
        if abs(fval-fvalold)>func_tolerance && exitflag ~= -1
            exitflag = 0;
        end
        fvalold = fval;
        
        if opt_restart == 0, exitflag = 1; end
        
    end
    
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

if strcmp(eqn,'burgers'),
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
    
elseif strcmp(eqn,'q1d'),
    
    xc = 0.5 .* ( x(2:end) + x(1:end-1) );
    dx = x(2:end) - x(1:end-1);
    
    dxf = zeros(imax,1);
    dxf(2:end-1) = 0.5*( dx(1:end-1) + dx(2:end) );
    dxf(1) = dx(1) - ( xc(1)-x(1) )/( xc(2)-xc(1) ) * ( dx(2)-dx(1) );
    dxf(end) = dx(end) + ( x(end)-xc(end) )/( xc(end)-xc(end-1) ) * ( dx(end)-dx(end-1) );
    
    TE = q1d_TE( imax,x,pb,pzero,sol,order,kappa,flag,flux_type,lim_type,shloc,As2 );
    
    % Write final truncation error to file
    fid = fopen('TE.dat','w');
    for i = 1:imax-1,
        fprintf(fid,'%23.15e  %23.15e %23.15e %23.15e %23.15e\n',xc(i),dx(i),TE(i+1,1),TE(i+1,2),TE(i+1,3) );
    end
    fclose(fid);   
    
    fid2 = fopen('faces.dat','w');
    for i = 1:imax,
        fprintf(fid2,'%23.15e %23.15e\n',x(i),dxf(i) );
    end
    
end

[ solu ] = exactsol(x,imax,sol,pb,pzero,order,shloc,flag,As2 );

save('exact.dat','solu','-ascii');
% 
% fid3 = fopen('exact.dat','w');
% for i = 1:i

% FIXME: 1) Arguments for q1d_TE
%        2) TE printing for Q1D equations



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
% global pb pzero sol order flux_type lim_type kappa shloc As2 flag
% global eqn imax
% FIXME: How do I deal with the TE ratio in the multi-equation system?

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

    