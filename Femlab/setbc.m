function [K,p,ks] = setbc(K,p,C,dof)
%***************************************************
% SetBC: 
%   Sets boundary conditions by imposing diagonal 
%   springs with 'large' stiffness ks.
% Syntax:
%   [K,p,ks] = setbc(K,p,C) 
%   [K,p,ks] = setbc(K,p,C,dof) 
%   [K,p] = setbc(K,p,C) 
%   [K,p] = setbc(K,p,C,dof) 
% Input:
%   K   :  original global stiffness matrix.
%   p   :  global load vector.
%   C   :  constraint matrix, C = [ node1 dof1 u1 
%                                   node2 dof2 ...].
%          or for dof=1,      C = [ node1 u1
%                                   node2 ...].
%   dof :  number of dof pr. node.
% Output:
%   K   :  global stiffness matrix with springs. 
%   p   :  load vector including spring loads.
%   ks  :  stiffness of constraining springs. 
% Date:
%   Version 1.0    04.05.95
%***************************************************

% default number of dof = 1
if nargin < 4 
  dof = 1;
end 

% Set spring stiffness.
ks = 10^6*max(abs(diag(K)));

% Introduce constraining spring stiffness and loads.
for i = 1:rows(C)
  if dof == 1
    j = C(i,1);
  else
    j = (C(i,1)-1)*dof + C(i,2);
  end
  K(j,j) = K(j,j) + ks;
  p(j)   = p(j) + ks*C(i,cols(C));
end

     
