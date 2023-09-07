function [q,S,E] = qq4p(q,T,X,G,u)
%***************************************************
% QQ4P: 
%   Creates and assembles internal force vector,
%   gradient matrix and flux matrix for a group
%   of potential quadrilateral 4-node elements.
% Syntax:
%   [q,S,E] = qq4p(q,T,X,G,u)
% Input:
%   q  :  global internal force vector
%   T  :  element topology matrix.
%   X  :  node coordinate matrix. 
%   G  :  material property matrix.  
%   u  :  global potential vector.
% Output:
%   q  :  updated internal force vector
%   S  :  global flux matrix in the following format:
%          [ Sx^1 Sy^1 ... Sx^4 Sy^4 ]
%         with one row for each element.
%   E  :  global gradient matrix in the following format:
%          [ ex^1 ey^1 ... ex^4 ey^4 ]
%         with row for each element.
% Date:
%   Version 1.0    04.05.95
%***************************************************

for j = 1:rows(T)  

  % define element arrays
  Xe = X(T(j,1:4),:);
  Ge = G(T(j,5),:);
  Ue = u(T(j,1:4)); 

  % evaluate element internal force
  [qe,Se,Ee] = qeq4p(Xe,Ge,Ue);

  % assemble global arrays
  q = assmq(q,qe,T(j,:),1);
  S(j,:) = reshape(Se',1,8);
  E(j,:) = reshape(Ee',1,8);
end

              
