function [q,Sn,En] = qq4epe(q,T,X,G,u,S,E,type)
%***************************************************
% QQ4EPE: 
%   Creates and assembles current internal force 
%   vector, stress matrix and strain matrix for 
%   a group of elasto-plastic quadrilateral 
%   4-node elements in plane strain.
% Syntax:
%   [q,Sn,En]  = qq4epe(q,T,X,G,u,S,E)
%   [q,Sn,En]  = qq4epe(q,T,X,G,u,S,E,type)
% Input:
%   q    : current global internal force vector.
%   T    : topology matrix for elements.
%   X    : initial node coordinate matrix. 
%   G    : element property matrix
%   u    : global displacement vector.
%   S    : current stress matrix
%   E    : current strain matrix
%   type : parameter setting material model
%          type = 1 -> Von Mises
%          type = 2 -> Drucker-Prager
% Output:
%   q    : internal forces from 4-node element group.
%   Sn   : updated stress matrix.
%   En   : updated strain matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************
  
% set default material model to Von Mises
if nargin==7
  type = 1;
end

% if not defined - expand S and E
if cols(S) ~= 20
  S(1,20) = 0;
end
if cols(E) ~= 20
  E(1,20) = 0;
end

% reshape global displacement vector
U = reshape(u,cols(X),rows(X))';

for j = 1:rows(T)   % We are looping the elements:

  % define element arrays 
  Xe =  X(T(j,1:4),:);              
  Ue =  U(T(j,1:4),:);              
  Ue = reshape(Ue',8,1);
  Ge  = G(T(j,5),:);

  % select row j and reshape into element format
  Se = reshape(S(j,:),5,4)';
  Ee = reshape(E(j,:),5,4)';
 
  % evaluate internal forces for element
  [qe,Sen,Een]  = qeq4epe(Xe,Ge,Ue,Se,Ee,type);     % Each element has a contribution, whicht then is assembled into global

  % assemble into global arrays
  q = assmq(q,qe,T(j,:),cols(X));  
  Sn(j,:) = reshape(Sen',1,20); 
  En(j,:) = reshape(Een',1,20); 
 
end

           
