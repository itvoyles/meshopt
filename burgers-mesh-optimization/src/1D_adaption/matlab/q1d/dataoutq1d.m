function stop = dataoutq1d(DV, optimValues, state)
global TEmax nu alpha Lref curiter curfeval iter fidhistory forward_mapping
global freq_history TEtolerance maxTEold leftbndry rightbndry
global pb pzero sol order flux_type lim_type kappa shloc As2 flag
global eqn param functional

% map design variable to grid nodes
x = forward_mapping(DV);
% x = springsystem_force(DV);


iter = optimValues.iteration;

% print to screen ---------------------------------------------------------

% Calculate truncation error

xc = 0.5 .* ( x(2:end) + x(1:end-1) );
dx = x(2:end) - x(1:end-1);

dxf = zeros(param,1);
dxf(2:end-1) = 0.5*( dx(1:end-1) + dx(2:end) );
dxf(1) = dx(1) - ( xc(1)-x(1) )/( xc(2)-xc(1) ) * ( dx(2)-dx(1) );
dxf(end) = dx(end) + ( x(end)-xc(end) )/( xc(end)-xc(end-1) ) * ( dx(end)-dx(end-1) );

TE = q1d_TE( param,x,pb,pzero,sol,order,kappa,flag,flux_type,lim_type,shloc,As2 );

% Select index for functional type

switch functional
    case 'cont'
        g = 1;
    case 'xmtm'
        g = 2;
    case 'energy'
        g = 3;
    case 'rss'
        g = 4;
end
    
% Calculate ratio of max truncation error to max truncation for each
% equation on evenly spaced grid

TEratio1 = max(abs(TE(:,1)))/TEmax(1); % continutity
TEratio2 = max(abs(TE(:,2)))/TEmax(2); % x-momentum
TEratio3 = max(abs(TE(:,3)))/TEmax(3); % energy

% If first iteration print header
if optimValues.iteration==1
    fprintf('\n Iter F-count            J(x)              TE1             TE2             TE3\n')
    printquest = 1;    
end

if g ~= 4,
    maxTE = max( abs((TE(:,g))) );
elseif g == 4,
    maxTE = sqrt(  ( max( abs(TE(:,1)) )/TEmax(1) ).^2 ...
        + ( max( abs(TE(:,2)) )/TEmax(2) ).^2 ...
        + ( max( abs(TE(:,3)) )/TEmax(3) ).^2 );
end

% print values to screen
if optimValues.iteration>=1
fprintf('%5.0f  %5.0f  %15.6e %15.6e %15.6e %15.6e %15.6e\n',...
    optimValues.iteration,optimValues.funccount,optimValues.fval, 1/TEratio1, 1/TEratio2, 1/TEratio3, abs(maxTEold-maxTE) )
end
if optimValues.fval < 0,
    keyboard
end

% write tecplot history ---------------------------------------------------

dx2 = [dx;0]; % For writing history, since length(dx)~=length(x) 

if mod(optimValues.iteration,freq_history) == 0
    fprintf(fidhistory,'Zone t="%10.0f"\n',optimValues.iteration);
    fprintf(fidhistory,'i=%4.0f\n',length(x));
     for i = 1:length(x)
     fprintf(fidhistory,'%23.15e %23.15e\n',x(i),dx2(i));
     end
end

% FIXME: Face locations need to be noted (they are the ones being moved),
%   but what is printed/plotted with them? Perhaps plot both face and
%   cell-center locations against dx (averages for the faces)?


% plot mesh distribution --------------------------------------------------
% figure(1)
plot(x,dxf,'--*r',xc,dx,'-ok')
axis([leftbndry,rightbndry,0,1.2*max( max(dxf),max(dx) )])
xlabel('Location')
ylabel('Cell Size (dx)')
legend('Face','Cell Center')


% Record adaption history to file -----------------------------------------
fid = fopen('opt-hist-1d.out','a');
fprintf(fid,'%5.0f  %5.0f  %15.6e %15.6e %15.6e %15.6e\n',...
        curiter,curfeval,optimValues.fval, 1/TEratio1,1/TEratio2,1/TEratio3);
curiter = curiter + 1;
curfeval = curfeval + 1;
stop = false;
fclose(fid);


% Check stopping tolerance based on dTE between iterations

% if optimValues.iteration<1 || isempty(maxTEold)
%     maxTEold = max(abs(TE(:,g)));
%    
% else
%     if abs(maxTEold-max(abs(TE(:,g))))<eps
% %         stop = true;
%     end
%     maxTEold = max(abs(TE(:,g)));
% end

if optimValues.iteration<1 || isempty(maxTEold)
    maxTEold = maxTE;
   
else
    if abs(maxTEold-maxTE)<eps
%         stop = true;
    end
    maxTEold = maxTE;
end
