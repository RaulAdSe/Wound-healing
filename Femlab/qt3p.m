function [q,S,E] = qt3p(q,T,X,G,u)
%***************************************************
% QT3P: 
%   Creates and assembles internal force vector,
%   gradient matrix and flux matrix for a group
%   of potential triangular 3-node elements.
% Syntax:
%   [q,S,E] = qt3p(T,X,G,u)
% Input:
%   T  :  element topology matrix.
%   X  :  node coordinate matrix. 
%   G  :  material property matrix.  
%   u  :  global potential vector.
% Output:
%   q  :  internal force vector.
%   S  :  global flux matrix.
%   E  :  global gradient matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

for j = 1:rows(T)  

  % define element arrays
  Xe = X(T(j,1:3),:);
  Ge = G(T(j,4),:);
  Ue = u(T(j,1:3)); 

  % evaluate element internal force
  [qe,Se,Ee] = qet3p(Xe,Ge,Ue); 

  % assemble into global arrays
  q = assmq(q,qe,T(j,:));
  S(j,:) = Se;  
  E(j,:) = Ee;   
end

              
