function q = assmq(q,qe,Te,dof)
%***************************************************
% AssmQ: 
%   Assembles global force vector by adding element 
%   nodal forces to existing global force vector.
% Syntax:
%   q = assmq(q,qe,Te)
%   q = assmq(q,qe,Te,dof)
% Input:
%   q  :  global vector.
%   qe :  element vector.
%   Te :  element topology vector.
%   dof:  degrees of freedom per node.
% Output:
%   q  :  updated global vector. 
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Default number of dof = 1
if nargin == 3 
  dof = 1;
end 

% Set number of element nodes.
ne = cols(Te) - 1;

% Define global address vector ig( ) for element dofs.
ig  = zeros(1,ne*dof);
for i = 1:ne
  for j = 1:dof
    ig((i-1)*dof+j) = (Te(i)-1)*dof + j;
  end
end

% Add element contributions.
for i = 1:ne*dof
  q(ig(i)) = q(ig(i)) + qe(i);
end
     
