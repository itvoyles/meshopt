%% AOE5204_H7_IanVoyles
% Problem 1 Script
% Part (a)

format long

A = [1 2 0 0; 3 4 5 0; 0 1 2 3; 0 4 -5 -6];

[E,L] = eig(A,'nobalance')

% Check that A = ELE^-1:
format short
Acheck = E*L/E

% Form the block-diagonal eigenvalue matrix:
format long
Eb = [E(:,1) real(E(:,2)) abs(imag(E(:,2))) E(:,4)]

% Form the block-diagonal matrix:

Lb = Eb\A*Eb

% Check that A = Eb Lb Eb^-1:
format short
Acheck2 = Eb*Lb/Eb

% Part (b)

% Part (c)



