function p = setload(p,P)
%***************************************************
% SetLoad: 
%   Sets the nodal loads in global vector form. 
% Input:
%   p   :  initial global load vector.
%   P   :  loads given as [NodeNo Load(1) .. Load(dof)].
% Output:
%   p    :  loads in global column vector. 
% Date:
%   Version 1.0    04.05.95
%***************************************************

% determine number of dof per node
dof = cols(P) - 1;

% Put loads into global load vector.
for i = 1:rows(P)
  j = (P(i,1)-1)*dof; 
  p(j+1:j+dof) = P(i,1+1:1+dof)';
end

     
