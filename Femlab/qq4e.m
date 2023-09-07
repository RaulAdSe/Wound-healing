function [q,S,E] = qq4e(q,T,X,G,u)
%***************************************************
% QQ4E: 
%   Creates and assembles internal force vector,
%   stress matrix and strain matrix for a group of
%   elastic quadrilateral 4-node elements in plane
%   stress or plane strain.
% Syntax:
%   [q,S,E] = qq4e(q,T,X,G,u)
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

% number of dof per node
dof = cols(X);

% determine number of nodes per element
nnodes = cols(T)-1;

% reshape global displacement vector
U = reshape(u,cols(X),rows(X))';

for j = 1:rows(T)

  % define element arrays 
  Xe =  X(T(j,1:nnodes),:);              
  Ue =  U(T(j,1:nnodes),:);              
  Ue = reshape(Ue',nnodes*dof,1);
  Ge = G(T(j,nnodes+1),:);
 
  % evaluate internal forces for element
  [qe,Se,Ee]  = qeq4e(Xe,Ge,Ue);   

  % assemble into global arrays
  q = assmq(q,qe,T(j,:),dof);  
  S(j,:) = reshape(Se',1,12); 
  E(j,:) = reshape(Ee',1,12); 
 
end

           
