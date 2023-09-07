function K = kt3e(K,T,X,G)
%***************************************************
% KT3E: 
%   Creates and assembles stiffness matrix for
%   a group of elastic triangular 3-node
%   elements in plane stress or plane strain.
% Syntax:
%   K = kt3e(K,T,X,G)
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

% Generate and assemble T3E elements.
for j = 1:rows(T)

  % define element arrays
  Xe = X(T(j,1:3),:);
  Ge = G(T(j,4),:);

  % evaluate element stiffness
  Ke = ket3e(Xe,Ge);

  % assemble element stiffness into global stiffness
  K  = assmk(K,Ke,T(j,:),2);
end

           
