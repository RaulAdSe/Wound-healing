function K = kq4eps(K,T,X,G,S,E,type)
%***************************************************
% KQ4EPS: 
%   Creates and assembles tangent stiffness matrix
%   for a group of elasto-plastic quadrilateral
%   4-node elements in plane stress.
% Syntax:
%   K = kq4eps(K,T,X,G,S,E)
%   K = kq4eps(K,T,X,G,S,E,type)
% Input:
%   K    : initial global stiffness matrix.
%   T    : element topology matrix.
%   X    : node coordinate matrix. 
%   G    : material property matrix. 
%   S    : current stress matrix.
%   E    : current strain matrix.
%   type : parameter setting material model
%          type = 1 -> Von Mises
%          type = 2 -> Drucker-Prager
% Output:
%   K    : new global tangent stiffness matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% set default material model to Von Mises
if nargin==6
  type = 1;
end

% if not defined - expand S 
if cols(S) ~= 16
  S(1,16) = 0;
end
if cols(E) ~= 16 
  E(1,16) = 0;
end

for j = 1:rows(T)
 
  % extract element information from global arrays
  Xe = X(T(j,1:4),:);
  Ge = G(T(j,5),:); 
  % select row j and reshape into element format
  Se = reshape(S(j,:),4,4)';
  Ee = reshape(E(j,:),4,4)';

  % evaluate element stiffness 
  Ke = keq4eps(Xe,Ge,Se,Ee,type);

  % assemble element stiffness into global stiffness
  K  = assmk(K,Ke,T(j,:),2);
end

           
