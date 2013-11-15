function [x,f,Converged,output,lambda]=slp_trust(fun,X0,Options,vlb,vub,grd,varargin)
% usage: [x,f,exitflag,output,lambda]=slp_trust(fun,x0,options,vlb,vub,grd,P1,P2,...)
%
% SLP       Finds the constrained minimum of a function of several
%           variables using Sequential Linear Programming.
%
% inputs:   fun     - function handle (or string) which returns the
%                     value of the objective function and a vector of
%                     constraints (i.e. [f,g]=fun(x)).  f is minimized
%                     such that g<zeros(g).
%           X0      - initial vector of design variables
%           Options - (optional) a structure according to optimset
%           vlb     - (optional) vector of lower bounds on design variables
%           vub     - (optional) vector of upper bounds on design variables
%           grd     - (optional) function handle that returns a vector of 
%                     function gradients and a matrix of constraint gradients
%                     (i.e. [fp,gp]=grd(x)).
%           Pn      - (optional) variables directly passed to fun and grd
%
%                     Note: optional inputs can be skipped by inputing []
%
% outputs:  x       - vector of design variables at the optimal solution
%           f       - final value of objective function
%           exitflag- Convegence flag
%           output  - Structure of output results (iteration history)
%           lambda  - Structure of Lagrange multipliers
            
%
% Written by:    Robert A. Canfield
%                Air Force Institute of Technology
%                AFIT/ENY, Bldg. 640
%                2950 Hobson Way
%                WPAFB, OH  45433-7665
%
% Created:        4/13/06
% Last modified: 10/30/11
%
% The function format is based on the MATLAB function constr.m written
% by Andy Grace of MathWorks, 7/90, supplemented with optimization toolbox 2
% structures for options and output. 
%---------------------------------------------------------------------
%
% Check inputs
%
if nargin<6, grd=[];
if nargin<5, vub=[];
if nargin<4, vlb=[];
if nargin<3, Options=[];
   if nargin<2
      disp('usage: [x,f,exitflag,output,lambda]=slp(fun,x0,options,vlb,vub,grd,P1,P2,...)')
      return
   end; end; end; end; 
end;
fd_gradients = isempty(grd);
numberofvariables=length(X0(:));
%
% Process options
%
if isfield(Options,'MoveLimit'),     MoveLimit     = Options.MoveLimit;     else MoveLimit = 0.2; end
if isfield(Options,'MoveReduction'), MoveReduction = Options.MoveReduction; else MoveReduction = 0.8; end
if isfield(Options,'TrustRegion'),   TrustRegion   = Options.TrustRegion;   else TrustRegion = 'on'; end
if isfield(Options,'TypicalX'),      TypicalX      = Options.TypicalX;      else TypicalX = X0(:); end
if ischar(TypicalX), TypicalX=eval(TypicalX); end
trust = strcmp(TrustRegion,'on'); if trust, contracted=0; end
options = optimset(optimset('fmincon'),Options);
TolFun  = optimget(options,'TolFun');
TolCon  = optimget(options,'TolCon');
TolX    = optimget(options,'TolX');
Display = optimget(options,'Display');
MaxIter = optimget(options,'MaxIter');
if isempty(MaxIter)
    MaxIter = 50*length(X0);
end
%
% Check lower and upper bounds: vlb and vub
%
[x,vlb,vub] = checkbounds(X0,vlb,vub,numberofvariables); % Optimization toolbox routine
MoveLimit = min(MoveLimit, max([(x-vlb);(vub-x)]));
if isempty(TypicalX)
   delx = max(TolX, MoveLimit*abs(x));
else
   delx = max(TolX, MoveLimit*abs(TypicalX(:)));
end
%
% Set up function and gradient calls
%
fcnstr = 'feval(fun,xshape(x),varargin{:})';
if all(size(x)==size(X0))
   xshape = inline('x');
else
   xshape = inline(['reshape(x(:),' int2str(size(X0,1)) ',' int2str(size(X0,2)) ')']);
end;
if fd_gradients
   DiffMaxChange = optimget(options,'DiffMaxChange');
   DiffMinChange = optimget(options,'DiffMinChange');
   fdstr = 'fdgrd(fun,x,f,g,DiffMinChange,DiffMaxChange,xshape,varargin{:})';
end
if fd_gradients
   grdstr = fdstr;
else
   grdstr = 'feval(grd,xshape(x),varargin{:})';
end
%
% Initial function evaluations.
%
[f, g]         = eval(fcnstr); g=g(:);
[gradf, gradg] = eval(grdstr); gradf=gradf(:);
[mg,mj]=max(g);
ndv  = length(x);
if ndv ~= length(gradf), error('Objective gradient length does not match # design variables'), end
if trust
   penalty = norm(gradf) ./ max( eps, sqrt( diag(gradg'*gradg) ) );
%  penalty = guessLM( g, 0, -10*TolCon, gradf, gradg ); %MCA-specific penalty
   Merit = f + penalty'*max(0,g);
end
if nargout>3
   output.f(1) = f;
   output.g(:,1) = g;
end
Iter=0;
if strcmpi(Display,'iter')
   disp(' ')
   disp('         Sequential Linear Programming Iteration History')
   disp('Iteration   Objective   Constraint_max  Index   Step-size  Merit       MoveLimit  TrustRatio')
      disp(sprintf('%9.0f  %11.4g  %14.6g  %5.0f  %10.4g  %9.4g',[Iter f mg mj 0 Merit]))
end
%
% ------------------------------------------------------------------------
% SLP loop
% ------------------------------------------------------------------------
%
optLP = optimset(optimset('linprog'),'Display','off');
%optLP.TolFun = min([TolFun; abs(f/TolFun); TolCon; TolX]); %MCA-specific TolFun
stat=0;
Converged=0;
while ~Converged && (stat>=0 || stat<=-2) && Iter<MaxIter %(allow infeasible to continue)
   Iter = Iter + 1;
   x0 = x;
   f0 = f;
   g0 = g;
   mg0 = mg;
%  delx = MoveLimit*TypicalX;     %MCA-specific delx
   dxlb = max( vlb-x0, -delx );
   dxub = min( vub-x0,  delx );
   MoveLimit0 = MoveLimit;
   [dx,df,stat,out,lambda]=linprog(gradf,gradg',-g,[],[],dxlb,dxub,[],optLP);
   x = max(vlb, min(vub, x0 + dx));
   [f, g] = eval(fcnstr); g=g(:);         % Evaluate functions
   [mg,mj] = max(g);

   % Update move limit.
   bound = lambda.lower>TolFun & x>vlb+TolX ...
         | lambda.upper>TolFun & x<vub-TolX;   
   Bound = any( bound );   
   if trust
      Merit0 = Merit;
      fapx = f0 + df;
      gapx = g0 + gradg'*dx;
      penalty = max( lambda.ineqlin, (penalty+lambda.ineqlin)/2 );
      Merit = f + penalty'*max(0,g);
      MeritApx = fapx + penalty'*max(0,gapx);
      MoreViolated = mg > max(TolCon,mg0);
      if Merit<Merit0 && MoreViolated % Increase penalty if violation increased
         penalty(mj) = penalty(mj) + (Merit0-Merit)/mg;
         Merit = Merit0;               % Reset merit to force contraction
      end
      FCD = MeritApx < Merit0; % Fractional Cauchy Decrease
      if FCD
         TrustRatio = (Merit0-Merit) / (Merit0-MeritApx);
      else
         TrustRatio = NaN;
      end
      rejected = (Merit >= Merit0 && MoreViolated) || (f >= f0 && mg >= max(0,mg0));
      if TrustRatio < 0.10 || ~FCD || rejected
         MoveLimit = 0.5*MoveLimit0;
         delx      = 0.5*delx;
      elseif TrustRatio > 0.75 && ~MoreViolated && ~contracted && Bound
         MoveLimit   = 2.0*MoveLimit0;
         delx(bound) = 2.0*delx(bound);
      elseif ~Bound %MCA-specific (finely tuned delx)
         steplength = max(0.5,max(abs(x-x0)./delx));
         if steplength<1
            MoveLimit = steplength*MoveLimit;
            delx      = steplength*delx;
         end
      end
      contracted = MoveLimit < MoveLimit0;
   else
      if Bound || mg>=mg0
         MoveLimit = MoveReduction*MoveLimit;
      end
      rejected=0;
   end
   
   % Check convergence.
   if ~rejected
      [gradf, gradg] = eval(grdstr); gradf=gradf(:); % Evaluate gradients
      Lagrangian = gradf + gradg*lambda.ineqlin ...
                         + (lambda.upper - lambda.lower).*(~bound);
      Slowed = abs(f-f0)    < TolFun && max(abs(x-x0)) < TolX && ...
               abs(f-Merit) < TolFun && stat>=0; %&& ~Bound;
      Converged = (Slowed || norm(Lagrangian,inf) < TolFun) ...
                && max(g) < TolCon;
   else
      Lagrangian = NaN;
   end

   % Output
   if strcmpi(Display,'iter') % Print intermediate results
      BoundStr = [];
      if Converged
         if Bound, BoundStr='Bound'; else BoundStr='Unbound'; end
      elseif trust
         if rejected, BoundStr='Rejected'; elseif ~FCD, BoundStr='~FCD'; end
      end
      disp([sprintf('%9.0f  %11.4g  %14.6g  %5.0f  %10.4g  %10.4g  %9.4g  %10.4g  ',...
           [Iter f mg mj max(abs(x-x0)) Merit MoveLimit0 TrustRatio]) BoundStr])
   end
   if nargout>3
      output.f(end+1) = f;
      output.g(:,end+1) = g;
   end
   
   if rejected % Reset to previous point
      x=x0; f=f0; g=g0; mg=mg0; Merit=Merit0;
   end
end
if nargout>3, output.iterations = Iter; output.message=stat; end
if nargout>4, lambda.ineq = lambda.ineqlin; end
%
% ------------------------------------------------------------------------
%
% Print final results
%
if strcmpi(Display,'iter')
   disp('           ----------   --------------         ----------')
   disp(sprintf('   Criteria %9.4g  %14.6g         %10.4g  %10.4g',[TolFun TolCon TolX]))
end
if ~strcmpi(Display,'off')
   if Converged
      disp(['SLP converged. Final objective function value = ',num2str(f)])
      disp(['               Lagrangian gradient   2-norm = ',num2str(norm(Lagrangian))])
      disp(['               Lagrangian gradient inf-norm = ',num2str(norm(Lagrangian,inf))])
      if strcmp(Display,'Iter')
         disp( '               Lagrange Multipliers   (j)')
         TolLM = sqrt(eps)*max([1;lambda.ineqlin]);
         disp(sprintf('%34.4g  %4.0f\n',[lambda.ineqlin(lambda.ineqlin>TolLM),...
                                         find(lambda.ineqlin>TolLM)]'))
         bound = find( lambda.lower>TolFun | lambda.upper>TolFun );
         if any(bound)
             disp( '               Lower    Upper         (i)')
             disp(sprintf('%22.4g  %10.4g  %4.0f\n',...
                 [lambda.lower(bound),...
                  lambda.upper(bound), bound]'))
         end
      end
   else
      disp(['SLP did NOT converge in ',num2str(Iter),' iterations.'])
   end
end
if stat<0 && ~Converged
   error('slp_trust:linprog',out.message)
end
%--------------------------------------------------------------------------






%--------------------------------------------------------------------------
function [df,dg] = fdgrd( fcn, x0, f0, g0, xmin, xmax, sx, varargin )
% PFDGRD      Calculates first forward finite difference gradients
%             of objective and constraint functions.
%
% usage:        [df,dg]=sqpfdgrd(fcn, x0, f0, g0, xmin, xmax, x1,P1,...,P15)
%
% inputs:       fcn      - function evaluation call
%               x0       - current design variable vector
%               f0       - objective value at x0
%               g0       - constraint values at x0
%               xmin     - minimum finite difference step
%               xmax     - maximum finite difference step
%               sx       - inline function to re-shape x
%               Pn       - optional variables directly passed to fcnstr
%
% outputs:      df       - finite difference objective gradient vector
%               dg       - finite difference constraint gradients matrix
%
% Written by:   Robert A. Canfield
%               AFIT/ENY, Bldg. 640
%               2950 Hobson Way
%               WPAFB, OH  45433-7665
%
% Created:      4/14/06
% Modified:      5/5/08

% Local variables
%
% dx....... Finite difference step
% f........ Perturbed objective value
% g........ Perturbed constraint values
% i........ Loop variable for current design variable perturbation

%--BEGIN
%
if nargin<5 || isempty(xmin), xmin=1e-8; end
if nargin<6 || isempty(xmax), xmax=1e-1; end
if nargin<7 || isempty(sx),   sx=inline('x'); end
% Less stringent relative change in x than 1.e-8 may be 
% needed for implicit functions that require numeric evaluation.
dx = min( max(1.e-8*abs(x0(:)),xmin), xmax );
dg = zeros(length(dx),length(g0(:)));
df = zeros(size(dg,1),1);

% Forward Finite Difference loop.
for i=1:length(x0(:)),
   x = x0;
   x(i) = x(i) + dx(i); 
   [f,g] = feval(fcn,sx(x),varargin{:});
   df(i,1) = (f - f0) / dx(i);
   if ~isempty(g)
      dg(i,:) = (g(:) - g0(:))'/ dx(i);
   end;
end
%--------------------------------------------------------------------------