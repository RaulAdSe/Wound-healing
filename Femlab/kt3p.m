function K = kt3p(K,T,X,G)
%***************************************************
% KT3P: 
%   Creates and assembles combined conductivity 
%   and dissipation matrix for a group of 
%   potential triangular 3-node elements.
% Syntax:
%   K = kt3p(K,T,X,G)
% Input:
%   K  :  initial global system matrix.
%   T  :  element topology matrix.
%   X  :  node coordinate matrix. 
%   G  :  material property matrix. 
% Output:
%   K  :  new global system matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

for j = 1:rows(T)

  % define element arrays
  Xe = X(T(j,1:3),:);
  Ge = G(T(j,4),:); 

  % evaluate element matrix
  Ke = ket3p(Xe,Ge);
  
  % assemble into global matrix
  K  = assmk(K,Ke,T(j,:));
end

           
