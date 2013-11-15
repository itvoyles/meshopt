function [varargout] = EULER_functional(DV)
global x y x_xsi y_xsi x_eta y_eta J exact1 exact2 exact3 exact4 parameters F1 G1 TEMatrix forward_mapping grad


[x,y] = forward_mapping(DV);

exact1=exact1*0;
exact2=exact2*0;
exact3=exact3*0;
exact4=exact4*0;
J_func = truncErrorU1(x,y,x_xsi,y_xsi,x_eta,y_eta,J,exact1,exact2,exact3,exact4,parameters,F1,G1,TEMatrix);
        

varargout{1} = J_func;
varargout{2} = grad;

end