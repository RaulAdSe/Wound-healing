function [K,p,q] = init(nn,dof)
%***************************************************
% init:
%   Initializes global stiffness matrix, global
%   load vector and global internal force vector.
% Syntax:
%   K = init(nn,dof)
%   [K,p] = init(nn,dof)
%   [K,p,q] = init(nn,dof)
% Input:
%   nn  :  number of global nodes.
%   dof :  number of degrees of freedom per node.
% Output:
%   K   :  initialized global stiffness matrix.
%   p   :  initialized global load vector.
%   q   :  initialized global internal force vector.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% square matrix K(nn*dof,nn*dof)
K = zeros(nn*dof);

% column vector p(nn*dof)
if nargout>1 
  p = zeros(nn*dof,1);
end

% column vector q(nn*dof)
if nargout == 3
  q = zeros(nn*dof,1);
end
