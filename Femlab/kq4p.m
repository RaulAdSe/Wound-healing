function K = kq4p(K,T,X,G)
%***************************************************
% KQ4P: 
%   Creates and assembles combined conductivity 
%   and dissipation matrix for a group of 
%   potential quadrilateral 4-node elements.
% Syntax:
%   K = kq4p(K,T,X,G)
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
  Xe = X(T(j,1:4),:);
  Ge = G(T(j,5),:); 

  % evaluate element matrix
  Ke = keq4p(Xe,Ge);

  % assemble global matrix
  K  = assmk(K,Ke,T(j,:));
end

           
