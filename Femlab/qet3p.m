function [qe,Se,Ee] = qet3p(Xe,Ge,Ue)
%***************************************************
% QET3P: 
%   Evaluates gradient and flux for the current
%   potentials. Creates the element internal force 
%   vector of a potential triangular 3-node element.
% Syntax:
%   [q,Se,Ee] = qet3p(Xe,Ge,Ue)
% Input:
%   Xe   :  coordinates Xe = [x1 y1; x2 y2; x3 y3].
%   Ge   :  material parameters Ge = [k b] or Ge = k.
%   Ue   :  nodal potentials
% Output:
%   qe   :  element internal force vector.
%   Se   :  element flux (stress) vector.
%   Ee   :  element gradient (strain) vector.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Form side vectors L (3 rows). 
a  = [ Xe(3,:) - Xe(2,:)
       Xe(1,:) - Xe(3,:)
       Xe(2,:) - Xe(1,:) ];

% Find triangle area A. 
A  = 0.5*abs(det(a(1:2,1:2)));

% Form strain matrix B.
B  = ((1/(2*A))*[-a(:,2)  a(:,1)])';

% Set isotropic (conductivity) stiffness. 
k  = Ge(1);
D  = eye(2)*k;

% Form gradient and flux vectors
Ee = (B*Ue)';
Se = Ee*D;

% Integrate internal force
qe = (B'*Se')*A;

