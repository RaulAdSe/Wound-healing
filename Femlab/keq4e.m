function Ke = keq4e(Xe,Ge)
%***************************************************
% KeQ4E: 
%   Creates the element stiffness matrix of elastic
%   quadrilateral 4-node element in plane stress or 
%   plane strain.
%   If Xe contains 8 node coordinate set the stiff-
%   ness matrix of 8-node element with quadratic 
%   shape functions is evaluated using reduced 
%   integration.
% Syntax:
%   Ke = keq4e(Xe,Ge)
% Input:
%   Xe   : coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4]
%   Ge   : element material data: [E , nu , (type)]
%          optional: type = 1 : plane stress (default)
%                    type = 2 : plane strain
% Output:
%   Ke   : element stiffness matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************  

% Gauss abscissae and weights.
r = [-1 1]/sqrt(3);
w = [ 1 1]; 

% Set isotropic elasticity matrix
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

% determine number of nodes per element
nnodes = rows(Xe);

% Initialize stiffness matrix.
Ke = zeros(2*nnodes);

% Gauss integration of stiffness matrix.
for i = 1:2
for j = 1:2

  % Parametric derivatives:
  dN = [ -(1-r(j))  (1-r(j))  (1+r(j)) -(1+r(j)) 
	 -(1-r(i)) -(1+r(i))  (1+r(i))  (1-r(i)) ]/4;

  if nnodes==8
    % evaluate the quadratic terms for the midside nodes
    dN8 = ... 
     [ -r(i)*(1-r(j))   0.5*(1-r(j)^2) -r(i)*(1+r(j))  -0.5*(1-r(j)^2) 
       -0.5*(1-r(i)^2) -r(j)*(1+r(i))   0.5*(1-r(i)^2) -r(j)*(1-r(i)) ];  

    % modify corner nodes 
    dN(:,1) = dN(:,1) - 0.5*dN8(:,1) - 0.5*dN8(:,4);  
    dN(:,2) = dN(:,2) - 0.5*dN8(:,2) - 0.5*dN8(:,1);  
    dN(:,3) = dN(:,3) - 0.5*dN8(:,3) - 0.5*dN8(:,2);  
    dN(:,4) = dN(:,4) - 0.5*dN8(:,4) - 0.5*dN8(:,3);  
    
    % expand gradient matrix  
    dN = [dN, dN8];
  end 

  % transform to global coordinates
  Jt = dN*Xe;
  dN = Jt\dN;

  % set up 4 node part of the gradient matrix 
  B  = [  dN(1,1)    0    dN(1,2)   0     dN(1,3)    0    dN(1,4)   0
	     0    dN(2,1)    0    dN(2,2)    0    dN(2,3)    0    dN(2,4)
	  dN(2,1) dN(1,1) dN(2,2) dN(1,2) dN(2,3) dN(1,3) dN(2,4) dN(1,4) ];  

  if nnodes == 8
    % set up gradient matrix for midside nodes
    B8 = [  dN(1,5)    0    dN(1,6)   0     dN(1,7)    0    dN(1,8)   0  
	       0    dN(2,5)    0    dN(2,6)    0    dN(2,7)    0    dN(2,8)
	    dN(2,5) dN(1,5) dN(2,6) dN(1,6) dN(2,7) dN(1,7) dN(2,8) dN(1,8) ];   
    
    % expand gradient matrix 
    B = [B,B8]; 
  end 

  Ke = Ke + w(i)*w(j)*( B'*D*B )*det(Jt);

end
end



