function df = dyieldvm(S,G,dL,Sy) 
%***********************************************
% dyieldvm:
%   Calculates the value of the gradient of the 
%   von Mises yield function for plane stress.
%   Reference: Krenk (1993): Nonlinear analysis
%   with finite elements, p.93.
% Syntax:
%   df = dyieldvm(S,G,dL,Sy)
% Input:
%   S   : stress vector = [  S11 
%                            S22 
%                            S12 ]  
%   G   : material properties = [ E nu S0 H ]
%   dL  : plastic multiplier increment.
%   Sy  : current (intermediate) yield stress
% Output:
%   df  : gradient of yield function.
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

xi1 = 2*Sy + E1*dL;  
xi2 = 2*Sy + E2*dL; 

S1 = S(1)+S(2);
S2 = S(1)-S(2);
S3 = S(3);

% calculate the two contributions to the gradient
df1 = -2*E1*S1^2/(xi1^3);  
df2 = -2*E2*(3*S2^2+12*S3^2)/(xi2^3); 

df = df1 + df2;

