function Ke = ket3e(Xe,Ge)
%***************************************************
% KeT3E: 
%   Creates the element stiffness matrix of elastic
%   triangular 3-node element in plane stress or 
%   plane strain.
% Syntax:
%   Ke = ket3e(Xe,Ge)
% Input:
%   Xe  :  corner coordinates (3 rows).
%   Ge  :  element material data: [E , nu , (type)]
%          optional: type = 1 : plane stress (default)
%                    type = 2 : plane strain
% Output:
%   Ke :  element stiffness matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Form side vectors a (3 rows).
a  = [ Xe(3,:) - Xe(2,:)
       Xe(1,:) - Xe(3,:)
       Xe(2,:) - Xe(1,:) ];

% Find triangle area A.
A  = 0.5*abs(det(a(1:2,1:2)));

% Form gradient vectors dN (3 columns).
dN  = (1/(2*A))*[-a(:,2)  a(:,1)]';

% Form strain matrix B by expanding the gradient matrix dN.
B  = [dN(1,1)      0     dN(1,2)      0     dN(1,3)      0  
	  0    dN(2,1)       0    dN(2,2)       0    dN(2,3) 
      dN(2,1)  dN(1,1)   dN(2,2)  dN(1,2)   dN(2,3)  dN(1,3)];

% Elasticity matrix 
E  = Ge(1);  
nu = Ge(2);
if cols(Ge) == 2                  % default: plane stress 
  type = 1; 
else 
  type = Ge(3); 
end

if type == 1                      % plane stress
  D  = E/(1-nu^2) ... 
	*[  1  nu     0
	   nu   1     0
	    0   0  (1-nu)/2 ];  
else                              % plane strain
  D  = E/((1+nu)*(1-2*nu)) ...   
	*[ 1-nu   nu       0
	    nu   1-nu      0
	     0     0  (1-2*nu)/2 ];  
end

% Integrate element stiffness matrix.
Ke = (B'*D*B)*A;
           
