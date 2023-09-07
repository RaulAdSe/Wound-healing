function [q,S,E] = qbar(q,T,X,G,u)
%***************************************************
% QBar: 
%   Creates and assembles current internal force 
%   vector, stress matrix and strain matrix for 
%   a group of elastic bar elements.
% Syntax:
%   q = qbar(q,T,X,G,u)
% Input:
%   q    :  existing global internal force vector.
%   T    :  topology matrix for bar elements.
%   X    :  initial node coordinate matrix. 
%   G    :  bar element property matrix. 
%   u    :  displacement vector.
% Output:
%   q    :  internal forces from bar element group.
%   S    :  axial stress 
%   E    :  axial strain
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Define current node cordinate array.
if  nargin == 4
  X1 = X;                                
else
  X1 = X + reshape(u,cols(X),rows(X))'; 
end

% Generate and assemble elastic bar elements.
for j = 1:rows(T)

  % define element arrays
  Xe0 = X(T(j,1:2),:);              
  Xe1 = X1(T(j,1:2),:);              
  Ge  = G(T(j,3),:);               

  % evaluate element internal force
  [qe,Se,Ee] = qebar(Xe0,Xe1,Ge);   

  % assemble into global arrays
  q = assmq(q,qe,T(j,:),cols(X));  
  S(j,1) = Se;
  E(j,1) = Ee;
end

           
