function K = assmk(K,Ke,Te,dof)
%***************************************************
% AssmK: 
%   Assembles system matrix by adding element 
%   matrix to existing global matrix. The system
%   matrix may be a stiffness matrix, mass matrix,
%   conductivity matrix, etc. 
% Syntax:
%   K = assmk(K,Ke,Te)
%   K = assmk(K,Ke,Te,dof)
% Input:
%   K  :  global matrix.
%   Ke :  element matrix.
%   Te :  element topology vector.
%   dof:  degrees of freedom per node.
% Output:
%   K  :  updated global matrix. 
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Default no=umber of dof = 1
if nargin == 3 
  dof = 1;
end 

% Set number of element nodes.
ne = cols(Te) - 1;

% Define global address vector ig( ) for element dofs.
ig  = zeros(1,ne*dof);
for i = 1:ne
  for j = 1:dof
    ig((i-1)*dof+j) = (Te(i)-1)*dof + j;
  end
end

% Add element contributions.
for i = 1:ne*dof
  for j = 1:ne*dof
    K(ig(i),ig(j)) = K(ig(i),ig(j)) + Ke(i,j);
  end
end
     
