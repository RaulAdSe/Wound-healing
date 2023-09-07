function n = rnorm(f,C,dof)
%...................................................
% rnorm:
%   Evaluates the reduced Euclidian norm of the
%   vector, f. Constrained terms in f are identified
%   from the matrix C specifies the prescribed dof. 
% Syntax:
%   n = rnorm(f,C,dof)
% Input:
%   f   : force vector
%   C   : constraint matrix
%   dof : number of dof pr. element node
% Output:
%   n   : Reduced Euclidian norm of f
% Date:
%   Version 1.0    04.05.95
%...................................................

% set up array 'fix' marking constrained terms in f
fix = zeros(rows(f),1);
for i=1:rows(C)
  dof_no = dof*(C(i,1)-1) + C(i,2);
  fix(dof_no) = 1;
end

% evaluate reduced norm
n = 0;
for i=1:rows(f)
  if fix(i) ~= 1
    n = n + f(i)^2;
  end
end
n = sqrt(n);
