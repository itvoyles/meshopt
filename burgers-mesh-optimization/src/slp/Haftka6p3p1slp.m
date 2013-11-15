% Script to solve Haftka Exercise 6.3.1 using SLP

% Initialize guess for design variables and move limit bounds
x0  = [ 1; 1];
xlb = [ 0; 0];
xub = [10; 10];

% Initialize termination criteria tolerances
options=optimset('TolX',0.01,'TolCon',1e-3,'Display','iter');
options.MoveLimit=0.5;

[xopt,fval] = slp_trust(@fHaftka6p3p1,x0,options,xlb,xub,@gHaftka6p3p1)