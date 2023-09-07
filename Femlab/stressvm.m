function [S,dL] = stressvm(S,G,Sy)
%***********************************************  
% stressvm:
%   Calculates the updated stress and plastic
%   multiplier of used in a radial return 
%   algorithm for von Mises materials in plane 
%   stress.
% Syntax:
%   [S,dL] = stressvm(S,G,Sy)
% Input:
%   S   : elastic predictor stress = [ S11 
%                                      S22 
%                                      S12 ]  
%   G   : material properties = [ E nu S0 H ]
%   Sy  : current (intermediate) yield stress
% Output:
%   S   : updated stress vector
%   dL  : plastic multiplier increment.
% Date:
%   Version 1.0    04.05.95
%*********************************************** 

E  = G(1);
nu = G(2);
H  = G(4);

% initial value of yield function
dL = 0;
f = yieldvm(S,G,dL,Sy);  

% Newton-Raphson iterations 
while abs(f) > 1e-6
  df  = dyieldvm(S,G,dL,Sy);
  ddL = -f/df;
  dL  = dL + ddL;
  f   = yieldvm(S,G,dL,Sy);
end 

% update yield stress  
Sy = Sy + H*dL;

% diameter and center of 'circle'
E1 = E/(1-nu);        
E2 = 3*E/(1+nu);
S1 = (S(1)+S(2))/(1+0.5*dL*E1/Sy);   
S2 = (S(1)-S(2))/(1+0.5*dL*E2/Sy);

% total stress
S(1) = 0.5*(S1+S2);  
S(2) = 0.5*(S1-S2); 
S(3) = S(3)/(1+0.5*dL*E2/Sy);

