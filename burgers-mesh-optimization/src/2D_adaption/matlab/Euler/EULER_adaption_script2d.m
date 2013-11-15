function EULER_adaption_script2d
global RE Lref imax jmax eps_small indxi indxip1 indxim1 indxjp1 indxjm1 indyi indyip1 indyim1 indyjp1 indyjm1 ...
    indxip1jp1 indyip1jp1 djdxmat scale TEmax minangle nextiter iter TEL2 fidhistory solver grid_save_freq func_method ...
    forward_mapping reverse_mapping freq_history xgs ygs method graditer TEfunc argin x y x_xsi y_xsi x_eta y_eta J exact1 exact2 exact3 exact4 parameters F1 G1 TEMatrix
% clear
clc
   close all
%% Setup executable path
%working_dir = '/home/grad3/bpickeri/Mesh_Opt_Working_Copy/burgers-mesh-optimization';
%setup executable path
%exe_path = [working_dir,'/library/bin'];
%path(path,exe_path);
%exe_path = [working_dir,'/library/src/2D_adaption/matlab/general'];
%path(path,exe_path);
%exe_path = [working_dir,'/library/src/2D_adaption/matlab/burgers'];
%path(path,exe_path);
%exe_path = [working_dir,'/library/src/2D_adaption/matlab/spring-f'];
%path(path,exe_path);
%exe_path = [working_dir,'/library/src/2D_adaption/matlab/spring-k'];
%path(path,exe_path);
%exe_path = [working_dir,'/library/src/2D_adaption/matlab/spring-k/coupled_spring-k'];
%path(path,exe_path);
%exe_path = [working_dir,'/library/src/2D_adaption/matlab/spring-k/simple_spring-k'];
%path(path,exe_path);
%% problem setup 

% input and general configuration options
readgrid = 0;           % read grid file 'grid.grd' = 1, create grid = 0
restart = 0;            % restart optimization from saved config.mat file
load_DV = 0;            % load the saved design variable for initialization
DV_file = ['DV-200.mat'];   % Select design variable file to load if load_DV = 1

% following controls are saved to config.mat file
m = 17;                          %grid nodes square domain
imax = m;                       % Number of grid nodes in i direction
jmax = m;                       % Number of grid nodes in j direction
RE = 8;            % Reynolds number of Burgers' equation
Lref = 4;                       % Reference domain size (meters)

%**************************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***************************************************************************************
%Setup for expansion fan
%************************EDIT THESE PARAMETERS**************************************
m=m; %dimension of square domain
xdim=Lref; %length of domain, meters
ydim=Lref; %width of domain, meters

M1 = 1.5;    %upstream mach number
angle = 12.5;   %turn angle in degrees
P1 = 101.325;   %upstream pressure in kPa
T1 = 273.16;    %upstream temperature in K

%MUSCL extrapolation variables
epsilon = 0.0; % 0 -> 1st order, 1 -> 2nd order 
kappa = 1.0; % -1=full upwind; 0=upwind bias; 1=central diff; 1/2=leonards quick; 1/3=3rd order 
%***********************************************************************************

P1 = P1*1000; %convert from kPa to Pa
[ M2, T2, P2] = PM_exp_fan(M1, T1, P1, angle );  %compute downstream conditions

%calculate the angles of the mach lines (radians)
u1 = asin(1/M1);
u2 = asin(1/M2);
%convert to angle counterclockwise from the bottom of the domain
u1 = u1 + pi*angle/180;

%calculate primitive variables for upstream flow
totalvel = M1*sqrt(1.4*287.04*T1);
uvel1 = totalvel*cos(pi*angle/180);
vvel1 = totalvel*sin(pi*angle/180);
rho1 = P1/(287.04*T1);

%calculate primitive variables for downstream flow
uvel2 = M2*sqrt(1.4*287.04*T2);
vvel2 = 0;
rho2 = P2/(287.04*T2);

%calculate prand meyer function for upstream conditions
v1 = sqrt( (1.4+1)/(1.4-1) )*atan( sqrt( (1.4-1)*(M1^2-1)/(1.4+1) ) ) - atan( sqrt( M1^2-1 ) );

%put exact solution parameters into a vector to be passed to the mex function 
parameters = [rho1,uvel1,vvel1,P1,rho2,uvel2,vvel2,P2,u1,u2,epsilon,kappa,v1];

%Generate a uniform grid, x and y vary from 0 to 4 meters
x=repmat(linspace(0,xdim,m)',1,m);
y=repmat(linspace(0,ydim,m),m,1);

%Generate flux matrix
F1=repmat(linspace(0,xdim,m)',1,(m-1));
G1=repmat(linspace(0,xdim,(m-1))',1,m);

%Generate the neccessary Matlab matrices and data structures for the exact
%solutions to four equations
 exact1=repmat(linspace(-2,2,(m+1))',1,(m+1));
 exact2=repmat(linspace(-2,2,(m+1))',1,(m+1));
 exact3=repmat(linspace(-2,2,(m+1))',1,(m+1));
 exact4=repmat(linspace(-2,2,(m+1))',1,(m+1));
 
 exact1=exact1*0;
 exact2=exact2*0;
 exact3=exact3*0;
 exact4=exact4*0;

[x_xsi, y_xsi, x_eta, y_eta, J] = gridMetricsCM(x,y);   %create matrices of correct size to store grid info
    
x_xsi=x_xsi*0;
y_xsi=y_xsi*0;
x_eta=x_eta*0;
y_eta=y_eta*0;
J=J*0;

TEMatrix = J*0;
%End setup for expansion fan
%*************************************************************************************


method = 'coupled-spring-f';        % Choose adaption method:
                        % 'pure' adaption with original system
                        % 'coupled-spring-k' spring adaption using spring constants
                        % 'simple-spring-k' spring adaption using spring constants
                        % 'coupled-spring-f' spring adaption using net force
                        
solver = 1;     % Choose solver type for coupled linear system
                        % 1 - MATLAB left divide
                        % 2 - If coupled-spring-k - penta-diagonal solver (more efficient for larger systems ~> 50x50
                        %     If simple-spring-k  - tri-diagonal solver
                        % 3 - Gauss-Seidel solver (not very fast)

% optimization options
maxiter = 100;       % maximum number of adaption iterations
opt_restart = 0;        % restat optimization after every maxiter
algorithm = 'interior-point';      % optimization algorithm
                                            % 'interior-point'
                                            % 'trust-region-reflective' - generates a warning but runs much faster than 'interior-point'
hessianapprox = 'lbfgs';  % hessian approximation
                                  % use 'lbfgs' for 'interior-point'
                                  % use 'fin-diff-grads' for 'trust-region-reflective'
gradon='off';
findifftype = 'forward';% finite difference type if gradobj='off'

tolerance = 0.0; % default 1e-6;
tolfun_tol = 1e-40;

func_tolerance = 1e-8; % change in max function min between outer iterations

% output
freq_history = 50;       % Frequency to write grid to history file

% small number limiting nearness of nodes
eps_small = 1e-6;

grid_save_freq = 200;   % Frequency to save intermediate addaption grids and restart file

func_method= 1;         % Select which functional to use
                        % 1 - Integral of TE^2
                        % 2 - Integral of TE^2 + TE_(rh)^2



if restart == 0
    save config.mat imax jmax RE Lref method ...
         freq_history eps_small maxiter opt_restart algorithm findifftype ...
         tolerance tolfun_tol grid_save_freq func_method
     
elseif restart == 1
    load config.mat
else
    error('Variable: restart must have value 0 or 1!')
end

% Grid input --------------------------------------------------------------
if readgrid == 1
    fprintf('Reading grid...');
    [x0,y0] = plot3d_read('grid.grd');
    if length(x0)~=imax || length(y0)~=jmax
        warning('Grid size does not match specified imax variable. Adjusting imax accordingly.')
    else
        fprintf('Done\n');
    end
    imax = length(x0);
else
    fprintf('Created evenly spaced grid.\n');
    x0 = repmat(linspace(0,Lref,imax)',1,jmax);     %edit to change domain size and offset
    y0 = repmat(linspace(0,Lref,jmax),imax,1);
end
    
% ----------------------
% Define TE function
% ----------------------
TEfunc = @Euler2D_TE;
argin = [RE,Lref];

    

%% Initial calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start timer
tic 

% Calculate truncation error on equally spaced mesh
xeq = repmat(linspace(0,Lref,imax)',1,jmax);        %edit to change domain size and offset
yeq = repmat(linspace(0,Lref,jmax),imax,1);

TE = TEfunc(x,y,x_xsi,y_xsi,x_eta,y_eta,J,exact1,exact2,exact3,exact4,parameters,F1,G1,TEMatrix);
TEorig = TE;
TEmax = max(max(abs(TEorig)));
TEL2 = sqrt( sum(sum(TEorig).^2)/numel(TEorig) );

% Setup history file to record adaption history
fidhistory = fopen('adaption-history.dat','w');
fprintf(fidhistory,'Variables="x""y""TE"\n');


graditer = -1; %variable to allow gradient calculation only at begining of iteration when GradObj='on'



%% Setup and define options for optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Setup method
% -------------------------------------------------------------------------
if strcmp(method,'pure')
    
    [A,b,Aeq,beq] = pure_linear_constraints(imax,jmax,eps_small,-Lref/2,Lref/2);   %edit to change domain size and offset
    reverse_mapping = @pure_reversemapping;
    forward_mapping = @pure_forwardmapping;
    
    
elseif strcmp(method,'coupled-spring-k')
    springsystem2d_initialize;
    [A,b,Aeq,beq] = spring_k_linear_constraints(imax,jmax,eps_small);
    reverse_mapping = @reversespringsystem2d_coupled;
    forward_mapping = @springsystem2d_coupled;
    
    if solver == 3
      xgs = reshape(x0,imax,jmax); %just to initialize starting matrix if the gauss-seidel solver is chosen
      ygs = reshape(y0,imax,jmax);
    end
    
elseif strcmp(method, 'simple-spring-k')
    springsystem2d_initialize;
    [A,b,Aeq,beq] = spring_k_linear_constraints(imax,jmax,eps_small);
    reverse_mapping = @reversespringsystem2d_simple;
    forward_mapping = @springsystem2d_simple;
    
elseif strcmp(method,'coupled-spring-f')
    springsystem2d_force_setup(imax,jmax);
    [A,b,Aeq,beq] = spring_f_linear_constraints(imax,jmax,eps_small,0,Lref);        %edit to change domain size and offset
    reverse_mapping = @reversespringsystem2d_force;
    forward_mapping = @springsystem2d_force;    
    
    
else
    
    error('Variable: method is incorrectly defined!')
    
end






% -------------------------------------------------------------------------
% Initialize 
% -------------------------------------------------------------------------
if load_DV == 1
    load(DV_file)
else
    DV = reverse_mapping(x0,y0);
    iter = 0;
end

% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------

options = optimset('Algorithm',algorithm,...
                   'Hessian',hessianapprox,...
                   'GradObj',gradon,...
                   'MaxFunEvals',10^10,...
                   'MaxIter',maxiter,...
                   'FinDiffType',findifftype,...
                   'RelLineSrchBnd',0.4,...
                   'TolFun',tolfun_tol,...
                   'TolCon',1e-10,...
                   'TolX',tolerance,...
                   'UseParallel','always',...
                   'PlotFcns',@dataoutput,...
                   'SubproblemAlgorithm','cg',...
                   'TypicalX',ones(size(DV,1),1)*1/(imax-1));

%                    'DerivativeCheck','on',...
%                    'Display','iter-detailed',...

%% Call Optimization routine
fprintf('Adapting Mesh...\n');

% -------------------------------------------------------------------------
% Optimization Call
% -------------------------------------------------------------------------               

fvalold = 0;
exitflag = 0;
while exitflag == 0

    [DV,fval,exitflag,output,lambda,grad,hessian] = ...
     fmincon(@EULER_functional,DV,[],[],Aeq,beq,[],[],[],options);

    
    if abs(fval-fvalold)>func_tolerance && exitflag ~= -1
        exitflag = 0;
    end
    fvalold = fval;
    
    if opt_restart == 0, exitflag = 1; end

end

% map design variable to grid nodes
[x,y] = forward_mapping(DV);

fprintf('Done\n');

%% Post processing




% Calculate truncation error
TE = TEfunc(x,y,x_xsi,y_xsi,x_eta,y_eta,J,exact1,exact2,exact3,exact4,parameters,F1,G1,TEMatrix);
%Extend TE matrix
TE(imax+1,:)=0;
TE(:,jmax+1)=0;
% Write final truncation error to file
fid = fopen('TE.dat','w');
for i = 1:numel(x)
  fprintf(fid,'%23.15e  %23.15e %23.15e\n',x(i),y(i),TE(i)); 
end
fclose(fid);


%% Write new grid to file
% -------------------------------------------------------------------------
% Write new grid to file
% -------------------------------------------------------------------------

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
global TEmax Lref iter fidhistory forward_mapping ...
       freq_history TEL2 grid_save_freq RE imax jmax TEfunc argin x y x_xsi y_xsi x_eta y_eta J exact1 exact2 exact3 exact4 parameters F1 G1 TEMatrix 

   
% update iteration counter
if optimValues.iteration > 0
    iter = iter + 1;
end   
   
% map design variable to grid nodes
[x,y] = forward_mapping(DV);


% print to screen ---------------------------------------------------------

% Calculate truncation error
TE = TEfunc(x,y,x_xsi,y_xsi,x_eta,y_eta,J,exact1,exact2,exact3,exact4,parameters,F1,G1,TEMatrix);
% TE1storder_terms = Burgers2d_TE1storder_terms(x,y,RE,Lref);
%Extend TE matrix

TE2 = sqrt( sum(sum( TE.^2))/(imax*jmax) );

% Calculate ratio of max truncation error to max truncation on evenly
% spaced grid
TEratio = TEmax/max(max(abs(TE)));
TEratio2 = TEL2/TE2;

% If first iteration print header
if optimValues.iteration==1
    fprintf('\n Iter F-count            J(x)              TE\n')
end

if optimValues.iteration>=1
fprintf('%5.0f  %5.0f  %15.6e %15.6e %15.6e\n',...
        iter,optimValues.funccount,optimValues.fval, TEratio, TEratio2)
end



% write tecplot history ---------------------------------------------------

if mod(iter,freq_history) == 0 || iter==0
    fprintf(fidhistory,'Zone t="%10.0f"\n',iter);
    fprintf(fidhistory,'i=%4.0f\n',size(x,1));
    fprintf(fidhistory,'j=%4.0f\n',size(x,2));
    fprintf(fidhistory,'DATAPACKING="block"\n');
    
    fprintf(fidhistory,'%23.15e %23.15e %23.15e %23.15e\n',x);
    if mod(imax*jmax,4)~=0; fprintf(fidhistory,'\n');end
    
    fprintf(fidhistory,'%23.15e %23.15e %23.15e %23.15e\n',y);
    if mod(imax*jmax,4)~=0; fprintf(fidhistory,'\n');end
    
    fprintf(fidhistory,'%23.15e %23.15e %23.15e %23.15e\n',TE);
    if mod(imax*jmax,4)~=0; fprintf(fidhistory,'\n');end
end

% fclose(fidhistory)


% plot mesh distribution --------------------------------------------------

% d = (TE(1:end-1,1:end-1)+TE(2:end,1:end-1)...
%     +TE(1:end-1,2:end)+TE(2:end,2:end))/4;
%xc=(x(2:end,2:end)+x(2:end,1:end-1)+x(1:end-1,2:end)+x(1:end-1,1:end-1))/4;
%yc=(y(2:end,2:end)+y(2:end,1:end-1)+y(1:end-1,2:end)+y(1:end-1,1:end-1))/4;
d=0*x;
surf(x,y,d,TE)
axis([0,Lref,0,Lref])      %edit to change domain size and offset 
view([0,0,90])

% Record adaption history to file -----------------------------------------

if iter==0
    fid = fopen('opt-hist-2d.dat','w');
    fprintf(fid,'Variables="iteration""Function count""Function value""max TE/max TEi""L2 TE/L2 TEi"\n');
    fclose(fid);
else
    fid = fopen('opt-hist-2d.dat','a');
    fprintf(fid,'%5.0f  %5.0f  %15.6e %15.6e %15.6e\n',...
            iter,optimValues.funccount,optimValues.fval, TEratio, TEratio2);
    fclose(fid);
end


% Save intermediate grid levels
if mod(iter,grid_save_freq)==0
    z = zeros(size(x));
    plot3d_write(['grid-',num2str(iter),'.grd'],x,y,z)
    
    dvout = ['DV-',num2str(iter),'.mat'];
    save(dvout, 'DV','iter')
end



stop = false;


%% Define adaption method dependent functions
% Must define:
% Method to return linear constraints [A,b,Aeq,Beq] 
% Forward mapping function if neccisary to transform design variable to grid nodes
% Reverse mapping function if neccisary to transform grid nodes to design
%   variables

% If no mapping functions are defined then define a blank function

%% pure functions
    function [A,b,Aeq,beq] = pure_linear_constraints(imax,jmax,eps_small,leftbndry,rightbndry)
        
        % sets up linear inequality constraints
        % specifies that x_(i+1) > x_(i)+eps_small using the linear system Ax<=b,
        % specifies that
        
        %setup indicies
        xivec = 1:imax*jmax;
        yivec = xivec+imax*jmax;

        imat = reshape(xivec,imax,jmax);
        jmat = reshape(yivec,imax,jmax);
        
        n = imax*jmax;
        A = zeros(2*n,2*n);
        b = zeros(2*n,1);
        %set xi<xi+1
        indxi = reshape(imat(1:end-1,1:end),(imax-1)*(jmax),1);
        indxip1 = reshape(imat(2:end,1:end),(imax-1)*(jmax),1);
        for j = 1:length(indxi)
            A(indxi(j),indxi(j)) = 1;
            A(indxi(j),indxip1(j)) = -1;
            b(indxi(j)) = -eps_small;
        end

        %set xj<xj+1
        indyj = reshape(imat(1:end,1:end-1),(imax)*(jmax-1),1)+n;
        indyjp1 = reshape(imat(1:end,2:end),(imax)*(jmax-1),1)+n;
        for j = 1:length(indyj)
            A(indyj(j),indyj(j)) = 1;
            A(indyj(j),indyjp1(j)) = -1;
            b(indyj(j)) = -eps_small;
        end

        A = sparse(A);
        b = sparse(b);
        % spy(A);

        %set equality contraints %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Aeq = zeros(2*n,2*n);
        beq = zeros(2*n,1);

        %set x1 and ximax are stationary
        indxl = reshape(imat(1,:),jmax,1);
        indxr = reshape(imat(imax,:),jmax,1);

        indyb = reshape(imat(:,1),imax,1)+n;
        indyt = reshape(imat(:,jmax),imax,1)+n;

        indborder = [indxl;indxr;indyb;indyt];
        for i = 1:imax
            Aeq(indxl(i),indxl(i)) = 1;
            beq(indxl(i)) = leftbndry;

            Aeq(indxr(i),indxr(i)) = 1;
            beq(indxr(i)) = rightbndry;
        end
        
        for i = 1:jmax
            Aeq(indyb(i),indyb(i)) = 1;
            beq(indyb(i)) = leftbndry;
            
            Aeq(indyt(i),indyt(i)) = 1;
            beq(indyt(i)) = rightbndry;
        end
        

        Aeq = sparse(Aeq);
        beq = sparse(beq);

        
    function [x,y] = pure_forwardmapping(X)
            % convertes design vector to two 2D arrays of nodes
            global imax jmax
            
            n = imax*jmax;
            
            x = reshape(X(1:n),imax,jmax);
            y = reshape(X(n+1:2*n),imax,jmax);
            
    function [X] = pure_reversemapping(x,y)
            % convertes two 2D arrays of nodes to vector
            [imax,jmax]=size(x);
            n = imax*jmax;
            
            X = [reshape(x,n,1);reshape(y,n,1)];

%% coupled spring-k functions

    function [A,b,Aeq,beq] = spring_k_linear_constraints(imax,jmax,eps_small)
        % Linear inequality constraints
        % Ak <= b
        % defines all spring constants must be greater than small number (eps_small)
        A = -diag(ones( (imax-1)*jmax+(jmax-1)*imax,1));
        b = -eps_small*ones( (imax-1)*jmax+(jmax-1)*imax,1);
        A = sparse(A);
        b = sparse(b);

        % Linear equality constraints
        % Ak = b
        Aeq = [];
        beq = [];

    function springsystem2d_initialize
    % initilizes indexing variables for efficient matrix setup
    global east west north south ...
           ixrowsleft ixrowsright ixrowsbottom ixrowstop...
           ixelemleft ixelemright ixelembottom ixelemtop imax jmax
       
    [east west north south ixrowsleft ixrowsright ixrowsbottom ixrowstop ...
       ixelemleft ixelemright ixelembottom ixelemtop] = setup_springsystem(imax,jmax);
    
%% spring-f functions
    function [A,b,Aeq,beq] = spring_f_linear_constraints(imax,jmax,eps_small,leftbndry,rightbndry)
    global SSAxinv SSAyinv

        % Set up inequality constraints
        ind = reshape([1:imax*jmax],imax,jmax);
        indx = reshape(ind(2:end,:),(imax-1)*jmax,1);
        SSAxjj = SSAxinv(indx,:);
        SSAxjjm1 = SSAxinv(indx-1,:);
        Ax = -(SSAxjj-SSAxjjm1);

        indyi = reshape(ind(:,2:end),imax*(jmax-1),1);
        indyim1=reshape(ind(:,1:end-1),imax*(jmax-1),1);
        SSAyjj = SSAyinv(indyi,:);
        SSAyjjm1 = SSAyinv(indyim1,:);
        Ay = -(SSAyjj-SSAyjjm1);
        
        Zx = zeros((imax-1)*jmax,imax*jmax);
        Zy = zeros((jmax-1)*imax,imax*jmax);
        A = [Ax,Zx;Zy,Ay];
        b = -eps_small*ones(size(A,1),1);

%         A = sparse(A);
%         b = sparse(b);


        % Set up equality constraints
        Aeq = zeros(2*imax+2*jmax,2*imax*jmax);
        beq = zeros(2*imax+2*jmax,1);

        cnt = 1;
        for i = 1:imax:imax*jmax
            Aeq(cnt,i) = 1;
            beq(cnt) = leftbndry;
            cnt = cnt + 1;
        end

        for i = imax:imax:imax*jmax
            Aeq(cnt,i) = 1;
            beq(cnt) = rightbndry;
            cnt = cnt + 1;
        end

        for i = imax*jmax+1:imax*jmax+imax
            Aeq(cnt,i) = 1;
            beq(cnt) = leftbndry;
            cnt = cnt + 1;
        end

        for i = imax*jmax+imax*(jmax-1)+1:imax*jmax+imax*(jmax)
            Aeq(cnt,i) = 1;
            beq(cnt) = rightbndry;
            cnt = cnt + 1;
        end

        % Set up equality constraints
%         Aeq = zeros(2*imax+2*jmax,2*imax*jmax);
%         beq = zeros(2*imax+2*jmax,1);
% 
%         cnt = 1;
%         for i = 1:imax:imax*jmax
%             Aeq(cnt,1:imax*jmax) = SSAxinv(i,:);
%             beq(cnt) = leftbndry;
%             cnt = cnt + 1;
%         end
% 
%         for i = imax:imax:imax*jmax
%             Aeq(cnt,1:imax*jmax) = SSAxinv(i,:);
%             beq(cnt) = rightbndry;
%             cnt = cnt + 1;
%         end
% 
%         for i = 1:imax
%             Aeq(cnt,imax*jmax+1:2*imax*jmax) = SSAyinv(i,:);
%             beq(cnt) = leftbndry;
%             cnt = cnt + 1;
%         end
% 
%         for i = imax*(jmax-1)+1:imax*jmax
%             Aeq(cnt,imax*jmax+1:2*imax*jmax) = SSAyinv(i,:);
%             beq(cnt) = rightbndry;
%             cnt = cnt + 1;
%         end
        
    