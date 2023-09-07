function K = kbar(K,T,X,G,u)
%***************************************************
% KBar: 
%   Creates and assembles tangent stiffness matrix
%   of elastic bar elements.
% Syntax:
%   K = kbar(K,T,X,G,u)
%   K = kbar(K,T,X,G)
% Input:
%   K    :  existing global stiffness matrix.
%   T    :  topology matrix for bar elements.
%   X    :  initial node coordinate matrix. 
%   G    :  bar element property matrix. 
%   u    :  displacement vector.
% Output:
%   K    :  new global tangent stiffness matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Define current node cordinate array.
if  nargin == 4
  X1 = X;                                
else
  X1 = X + reshape(u,cols(X),rows(X))'; 
end

for j = 1:rows(T)

  % define element arrays
  Xe0 = X(T(j,1:2),:);              
  Xe1 = X1(T(j,1:2),:);              
  Ge  = G(T(j,3),:);               
 
  % evaluate element stiffness
  Ke  = kebar(Xe0,Xe1,Ge);

  % assemble element stiffness into global stiffness
  K   = assmk(K,Ke,T(j,:),cols(X)); 
end

           
