clear;close;clc

%Optimization of Area of skirt and plenum
m = 0.7;
g = 9.81;
Qfmax = 0.453069552/60;
% Pfmax = 1*g/Area_Plenum;

syms H l w lambda1 lambda2 lambda3
f = l*w + H*w*0.5; %Objective Function -> Area
%Constraint 1 -> Minimum pressure to hover
p = (1*g/(l*w + H*w*0.5))*(1-(3.525e-4/Qfmax))-(m*g)/(l*w) == 0;
%Dimensional Constraints from rulebook
g = (l+H)^2 - 0.60^2 == 0; 
h = w^2 - 0.27^2 == 0;

L = f + lambda1*lhs(g) + lambda2*lhs(h)+ lambda3*lhs(p);

dLdH = diff(L,H) == 0;
dLdw = diff(L,w) == 0;
dLdl = diff(L,l) == 0;
dLdlambda1 = diff(L,lambda1) == 0;
dLdlambda2 = diff(L,lambda2) == 0;
dLdlambda3 = diff(L,lambda3) == 0;

Sol = solve([dLdH;dLdw;dLdl;dLdlambda1;dLdlambda2;dLdlambda3], [H w l lambda1 lambda2 lambda3], 'Real', true);

w = abs(Sol.w(1));
l = abs(Sol.l(1));
H = abs(Sol.H(1));

length_of_rectangular_area= l
Total_Width = w
Length_of_triangular_area = H

% lambda1 = abs(Sol.lambda1(1))
% lambda2 = abs(Sol.lambda2(1))