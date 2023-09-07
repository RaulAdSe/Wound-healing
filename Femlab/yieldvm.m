function f = yieldvm(S,G,dL,Sy)
%***********************************************  
% yield:
%   Calculates the value of the von Mises
%   yield function for plane stress.
%   Reference: Krenk (1993): Nonlinear analysis
%   with finite elements, p.93.
% Syntax:
%   f = yieldvm(S,G,dL,Sy)
% Input:
%   S   : stress vector = [  S11 
%                            S22 
%                            S12 ]  
%   G   : material properties = [ E nu S0 H ]
%   dL  : plastic multiplier increment.
%   Sy  : current (intermediate) yield stress
% Output:
%   f   : value of yield function.
% Date:
%   Version 1.0    04.05.95
%*********************************************** 

% translate material properties
E  = G(1);     
nu = G(2);
S0 = G(3);
H  = G(4);

% define pseudo-variables
E1 = 2*H + E/(1-nu);
E2 = 2*H + 3*E/(1+nu);

S1 = S(1)+S(2);
S2 = S(1)-S(2);
S3 = S(3);

% calculate 3 contributions to f
f1 = S1^2/(2*Sy+dL*E1)^2; 
f2 = (3*S2^2)/(2*Sy+dL*E2)^2;
f3 = (12*S3^2)/(2*Sy+dL*E2)^2;

f = f1 + f2 + f3 - 1;
