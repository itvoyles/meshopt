function djdxmat = calc_djdx(x,y)
% global djdxmat
global TEfunc argin
  [imax,jmax]=size(x);

  err_lim = 0.0001; %max error [fraction] in djdx
  djdxmat = zeros(imax,jmax,2);
  
%   djdxlim = sum(sum(sum(abs(djdxmat))))/(imax*jmax);
  
  TE = TEfunc(x,y,argin);

  for j = 1:jmax
  for i = 1:imax
  
%   if any(abs(djdxmat(i,j,:))>=djdxlim)
      
%   dx = (max(x(:,j))-min(x(:,j)))./1000000.;

% it was found that this value of dx results in the lowest error for 2D
% Burgers' Equation. Values tested that were larger and smaller result in
% a greater error in the approximation due to increasing truncation error
% for larger values and increasing round-off error for smaller values.
  dx = 8e-7; 


  %find index for submatrix to calculate djdx (submatrix = 9x9)
   len = 5;
    if (i<1+len) 
      il = max(i-len,1);
      ih = il+2*len;
      iloc = i;
      ex(1) = 0;
      ex(2) = 2;
    elseif (i>imax-len)
      ih = min(i+len,imax);
      il = ih-2*len;
      iloc = i-il+1;
      
      ex(1) = 2;
      ex(2) = 0;
   else
      il = i-len;
      ih = i+len;
      iloc = 1+len;
      ex(1) = 1;
      ex(2) = 1;
   end

   if (j<1+len) 
      jl = max(j-len,1);
      jh = jl+2*len;
      jloc = j;
      ex(3) = 0;
      ex(4) = 2;
   elseif (j>jmax-len)
      jh = min(j+len,jmax);
      jl = jh-2*len;
      jloc = j-jl+1;
      ex(3) = 2;
      ex(4) = 0;
   else
      jl = j-len;
      jh = j+len;
      jloc = 1+len;
      ex(3) = 1;
      ex(4) = 1;
   end

  if (il < 1)
      il = 1;
      ex(1) = 0;
      iloc = i;
  end
  if (ih > imax)
      ih = imax;
      ex(2) = 0;
      iloc=i;
  end
  if(jh >jmax)
      jh = jmax;
      ex(4) = 0;
      jloc=j;
  end
  if (jl<1)
      jl=1;
      ex(3)=0;
      jloc=j;
  end
  
  %override and use entire matrix
%   il=1; ih=imax;
%   jl=1; jh=jmax;
%   iloc=i; jloc=j;
%   ex=zeros(1,4);
  
  
  xloc = x(il:ih,jl:jh);
  yloc = y(il:ih,jl:jh);
  TEsmall = TE(il:ih,jl:jh);
 
  
  Jsmall = trap_sum2d(xloc(1+ex(1):end-ex(2),1+ex(3):end-ex(4)),...
                      yloc(1+ex(1):end-ex(2),1+ex(3):end-ex(4)),...
                      TEsmall(1+ex(1):end-ex(2),1+ex(3):end-ex(4) ).^2 );


  djdxmat(i,j,1) = djdx_loc_full(xloc,yloc,iloc,jloc,dx,ex,Jsmall);
  djdxmat(i,j,2) = djdy_loc_full(xloc,yloc,iloc,jloc,dx,ex,Jsmall);
  

  
  % Check for NaN
  if isnan(djdxmat(i,j,1))
      djdxmat(i,j,1) = 0;
  end
  if isnan(djdxmat(i,j,2))
      djdxmat(i,j,2) = 0;
  end
  
  
  
  end
  end


  
%   djdxmat(:,:,1) = djdxloc;
%   djdxmat(:,:,2) = djdyloc;
  
  djdxmat(1,1:jmax,1) = 0.;
  djdxmat(imax,1:jmax,1) = 0.;

  djdxmat(1:imax,1,2) = 0.;
  djdxmat(1:imax,jmax,2) = 0.;
  
%   djdxmat(:,:,1) = djdxmat(:,:,1).*TE;
%   djdxmat(:,:,2) = djdxmat(:,:,2).*TE;



  %write(*,'(9e12.4)') djdx(1:9,1:9)

end