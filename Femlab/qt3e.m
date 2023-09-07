function [q,S,E] = qt3e(q,T,X,G,u)
%***************************************************
% QT3E: 
%   Creates and assembles internal force vector,
%   stress matrix and strain matrix for a group of
%   elastic triangular 3-node elements in plane
%   stress or plane strain.
% Syntax:
%   [q,S,E] = qt3e(q,T,X,G,u)
% Input:
%   q  :  initial global force vector.
%   T  :  element topology matrix.
%   X  :  node coordinate matrix. 
%   G  :  material property matrix. 
%   u  :  global displacement vector.
% Output:
%   q  :  new global internal force vector.
%   S  :  global stress matrix.
%   E  :  global strain matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% reshape global displacement vector
U = reshape(u,cols(X),rows(X))';

for j = 1:rows(T)

  % define element arrays 
  Xe =  X(T(j,1:3),:);              
  Ue =  U(T(j,1:3),:);              
  Ue = reshape(Ue',6,1);
  Ge  = G(T(j,4),:);
 
  % evaluate internal forces for element
  [qe,Se,Ee]  = qet3e(Xe,Ge,Ue);   

  % assemble into global arrays
  q = assmq(q,qe,T(j,:),cols(X));  
  S(j,:) = Se;
  E(j,:) = Ee; 
 
end

           

