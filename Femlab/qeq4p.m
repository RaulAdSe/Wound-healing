function [qe,Se,Ee] = qeq4p(Xe,Ge,Ue)
%***************************************************
% QeQ4P: 
%   Evaluates gradient and flux for the current
%   potentials. Creates the element internal force 
%   vector of a potential quadrilateral 4-node element.
% Syntax:
%   [qe,Se,Ee] = qeq4p(Xe,Ge,Ue)
% Input:
%   Xe   :  coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4].
%   Ge   :  material parameters Ge = [k b] or Ge = k.
%   Ue   :  element nodal potentials vector [u1; u2; u3; u4].
% Output: 
%   qe   :  element internal force vector.
%   Se   :  element flux (stress) matrix.
%   Ee   :  element gradient (strain) matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Gauss abscissae and weights.
r = [-1 1]/sqrt(3);
w = [ 1 1]; 

% Set isotropic (conductivity) stiffness.
k  = Ge(1);
D  = [k 0
      0 k];

qe = zeros(4,1);
% Gauss integration of internal force vector.
for i = 1:2 
for j = 1:2
  % organize the Gauss points as the element nodes
  gp = i + 3*(j-1) - 2*(i-1)*(j-1);

  % set up gradient matrix B
  dN = [ -(1-r(j))  (1-r(j))  (1+r(j)) -(1+r(j)) 
         -(1-r(i)) -(1+r(i))  (1+r(i))  (1-r(i)) ]/4;
  Jt = dN*Xe;
  B  = Jt\dN;

  % evaluate gradient and flux
  Ee(gp,:) = (B*Ue)';
  Se(gp,:) = Ee(gp,:)*D;

  % evaluate internal forces
  qe = qe + w(i)*w(j)*( B'*Se(gp,:)')*det(Jt);
end
end



