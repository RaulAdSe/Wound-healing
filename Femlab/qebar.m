function [qe,Se,Ee] = qebar(Xe0,Xe1,Ge)
%***************************************************
% QeBar: 
%   Evaluates the strain and stress for the current
%   displacements. Creates the element internal
%   force vector for elastic bar element.
% Syntax:
%   qe = qebar(Xe0,Xe1,Ge)
% Input:
%   Xe0 : initial coordinates = [x1    y1    (z1)
%                                x2    y2    (z3)   ].
%   Xe1 : current coordinates = [x1+u1 y1+v1 (z1+w1)
%                                x2+u2 y2+v2 (z3+w2)].
%   Ge  : element properties Ge = [A E] or Ge = AE.
% Output:
%   qe  : element nodal force vector.
%   Se  : element axial stress.
%   Ee  : element axial strain.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Form initial element (column) vector a0(:).
a0 = (Xe0(2,:)-Xe0(1,:))';
l0 = sqrt(a0'*a0);

% Form current element (column) vector a1(:).
a1 = (Xe1(2,:)-Xe1(1,:))';
l1 = sqrt(a1'*a1);

% Get element properties.
A  = Ge(1);
if  cols(Ge) == 1
  E = 1;
else
  E = Ge(2);
end

% Find axial strain and stress 
Ee  = 0.5*(l1^2 - l0^2)/l0^2; 
Se = E*Ee;

% Find internal element force in current position.
N  = A*Se;
qe =  (N/l0)*[-a1
               a1];
