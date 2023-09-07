function [S,dL] = stressdp(S,G,Sy0,dE,dS)
%***********************************************  
% stressdp:
%   Calculates the updated stress and plastic
%   multiplier of used in a backward Euler
%   algorithm for associated Drucker-Prager
%   materials in plane stress.
% Syntax:
%   [S,dL] = stressdp(S,G,Sy)
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

% set material parameters
E  = G(1);  
nu = G(2);
S0 = G(3);
H  = G(4);
phi = G(5);

% elastic flexibility matrix
C  = 1/E ...   
      *[   1  -nu    0
	  -nu   1    0 
	   0    0   2*(1+nu) ];

% deviator, mean and equivalent stress 
[Sd,Sm] = devstress(S);
Seq = eqstress(S);

% D-P yield condition
f = Seq + phi*Sm - Sy0; 

% gradient of yield function 
sd = [Sd(1),Sd(2),2*Sd(3)]'; 
mp = [1 1 0]';
df = 3/(2*Seq)*sd + phi/3*mp; 

% initialization...
R  = zeros(3,1);
deltaS = zeros(3,1);
dL=0;

FTOL = 1e-6; RTOL=1e-3*norm(dE); 
while (norm(R) > RTOL | abs(f) > FTOL)
     
  % calculate 2nd derivative of yield function
  d2f1 = 3/(2*Seq)*diag([1 1 2]);
  d2f2 = 9/(4*Seq^3)*sd*sd';
  d2f =  d2f1 - d2f2; 
      
  % tangent stiffness  
  K  = [ C + dL*d2f  df
         df'         -H ];
    
  % evaluate subincrement
  delta = K\[R;-f];

  deltaS = deltaS + delta(1:3); 
  dL = dL + delta(4);

  % new stress values
  [Sd,Sm] = devstress(S+deltaS);
  Seq = eqstress(S+deltaS);

  % update yield stress
  Sy = Sy0 + dL*H;

  % D-P yield condition
  f = Seq + phi*Sm - Sy;

  % strain residual 
  sd = [Sd(1),Sd(2),2*Sd(3)]'; 
  df = 3/(2*Seq)*sd + phi/3*mp;
  R  = dE - C*(dS+deltaS) - dL*df;

end

% update stresses with plastic correction
S = S + deltaS;
