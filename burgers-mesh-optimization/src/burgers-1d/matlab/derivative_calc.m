function der = dnfdxn(f,dx,dorder,order)
%calculate an arbitrary derivative of arbitrary order
% x : nodes
% f : value
% dorder : derivative order
% order : order of accuracy



% nterms = max(2*floor( ((dorder-1) + (order-1)+1)/2 ),1); %number of TS
% expansions
nterms = (dorder-1) + (order-1)+1;
sten_width = nterms/2; % numer of nodes in stencil minus 1

N = length(f);
if sten_width+1 > N
    fprintf('Error, not enough points to calculate derivative!\n')
end

%create stencils
i = 1;
stencil(i,:) = (i-1):(i-1)+nterms;


for i = 1:nterms
    stencil(i+1,:) = stencil(i,:)-1;
end


%solve Taylor Series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(stencil(:,1));
   clear A
    
   I  = find(stencil(i,:) ~= 0);
   I0 = find(stencil(i,:) == 0);
   cur_stencil = stencil(i,I)';
   
   pwr = [1:nterms+2];
   
   lcursten = length(cur_stencil);
   trunc = zeros(lcursten,length(pwr)-lcursten);
   for j = 1:lcursten
      mat = cur_stencil(j).^pwr./factorial(pwr);
      A(j,:) = mat(1:lcursten);
      trunc(j,:) = mat(lcursten+1:end);
   end
    A = A';
    trunc = trunc';
   
   b = zeros(size(A,1),1);
   b(dorder) = 1;
   
   multipliers = A\b;
   coef(i,I) = multipliers;
   coef(i,I0) = -sum(coef(i,I));
   
   trunc_terms = sum(trunc.*repmat(multipliers,1,size(trunc,1))',2);
   Ipwr = find(trunc_terms ~= 0);
   final_order(i) = pwr(lcursten+Ipwr(1))-dorder;
   final_trunc(i) = trunc_terms(Ipwr(1));
       
end

coef = coef/dx^dorder;

%Calculate derivative %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Icenter = ceil((nterms+1)/2);
for i = 1:length(f)
   phy_stencil = stencil + i;
   
   if min(phy_stencil(Icenter,:)) < 1
       nshift = 1-min(phy_stencil(Icenter,:));
       Isten = Icenter-nshift;
   elseif max(phy_stencil(Icenter,:)) > length(f)
       nshift = max(phy_stencil(Icenter,:) - length(f));
       Isten = Icenter + nshift;
   else
       Isten = Icenter;
   end
   
   
   cur_coef = coef(Isten,:);

 
   der(i) = cur_coef*f(phy_stencil(Isten,:));
   
    
end

    der = reshape(der,size(f));

end