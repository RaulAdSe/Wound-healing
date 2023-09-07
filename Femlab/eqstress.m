function eqs = eqstress(S)
%***************************************************
% eqstress:
%   Evaluates equivalent (von Mises) stress of a
%   stress vector S. The type of stress assumption
%   is indicated by the number of stress components.
%       Plane stress = 3 components   
%       Plane strain = 4 components   
%       3D stress    = 6 components   
% Syntax:
%   eqs = eqstress(S)
% Input:
%   S   :  stress vector.
% Output: 
%   eqs :  equivalent stress. 
% Date:
%   Version 1.0    04.05.95
%***************************************************

if rows(S) == 3      % plane stress
  eqs = S(1)^2 + S(2)^2 - S(1)*S(2) + 3*S(3)^2;
else                 % plane strain or 3D
  % normal stress contribution
  eqs = 0.5*( (S(1)-S(2))^2 + (S(2)-S(3))^2 + (S(3)-S(1))^2 );   
  % shear stress contribution
  eqs = eqs + sum(3*S(4:rows(S)).^2); 
end
eqs = sqrt(eqs);    
