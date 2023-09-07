function Ke = ket3p(Xe,Ge)
%***************************************************
% KeT3P: 
%   Creates the combined element conductivity and
%   dissipation matrix of potential triangular 
%   3-node element.
% Syntax:
%   Ke = ket3p(Xe,Ge)
% Input:
%   Xe   :  coordinates Xe = [x1 y1; x2 y2; x3 y3].
%   Ge   :  material parameters Ge = [k b] or Ge = k.
% Output:
%   Ke   :  element conductivity/dissipation matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Form side vectors a(i,:), i = 1,2,3.
a  = [ Xe(3,:)-Xe(2,:)
       Xe(1,:)-Xe(3,:)
       Xe(2,:)-Xe(1,:) ];

% Find triangle area A.
A  = 0.5*abs(det(a(1:2,1:2)));

% Set conductivity coefficient.
k  = Ge(1);   

% isotropic conductivity matrix
D = eye(2)*k;

% gradient matrix
B = 1/(2*A)*[-a(:,2) a(:,1)]';
 
% Element conductivity matrix.
Ke = A*B'*D*B; 

% Add conditional contribution from dissipation b = Ge(2).
if cols(Ge) == 2
  b = Ge(2);

% Dissipation matrix.
  Kb = (b*A/12)*[ 2  1  1
                  1  2  1
                  1  1  2 ]; 

% Add dissipation contribution to Ke.
  Ke = Ke + Kb;  

end
