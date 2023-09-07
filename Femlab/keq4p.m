function [Ke,A] = keq4p(Xe,Ge)
%***************************************************
% KeQ4P: 
%   Creates the combined element conductivity and
%   dissipation matrix of potential quadrilateral
%   4-node element.
% Syntax:
%   Ke = keq4p(Xe,Ge)
% Input:
%   Xe   :  coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4].
%   Ge   :  material parameters Ge = [k b] or Ge = k.
% Output:
%   Ke   :  element conductivity/dissipation matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Gauss abscissae and weights.
r = [-1 1]/sqrt(3);
w = [ 1 1]; 

% Set isotropic conductivity matrix
k  = Ge(1);
D  = [k 0
      0 k];

% Set isotropic dissipation parameter.
if cols(Ge) == 2
  b = Ge(2);
else
  b = 0.0;
end

% Initialize system matrix.
Ke = zeros(4);

% Gauss integration of system matrix.
for i = 1:2
for j = 1:2
  N  = [ (1-r(i))*(1-r(j)) (1+r(i))*(1-r(j)) ...
         (1+r(i))*(1+r(j)) (1-r(i))*(1+r(j)) ]/4;
  dN = [ -(1-r(j))  (1-r(j))  (1+r(j)) -(1+r(j)) 
         -(1-r(i)) -(1+r(i))  (1+r(i))  (1-r(i)) ]/4;
  Jt = dN*Xe;
  B  = Jt\dN;
  Ke = Ke + w(i)*w(j)*( B'*D*B + N'*b*N )*det(Jt);
end
end



