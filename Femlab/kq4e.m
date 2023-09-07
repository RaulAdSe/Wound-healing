function K = kq4e(K,T,X,G)
%***************************************************
% KQ4E: 
%   Creates and assembles stiffness matrix for
%   a group of elastic quadrilateral 4-node 
%   elements in plane stress or plane strain.
% Syntax:
%   K = kq4e(K,T,X,G)
% Input:
%   K  :  initial global stiffness matrix.
%   T  :  element topology matrix.
%   X  :  node coordinate matrix. 
%   G  :  material property matrix. 
% Output:
%   K  :  new global stiffness matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% determine number of nodes per element  
nnodes = cols(T)-1; 

for j = 1:rows(T)

  % define element arrays
  Xe = X(T(j,1:nnodes),:);
  Ge = G(T(j,nnodes+1),:); 

  % evaluate element stiffness
  Ke = keq4e(Xe,Ge);

  % assemble element stiffness into global stiffness
  K  = assmk(K,Ke,T(j,:),2);
end

           
