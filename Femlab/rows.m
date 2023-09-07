function m = rows(X)
%***************************************************
% rows:
%   Determines number of rows in matrix X.
% Syntax:
%   m = rows(X)
% Input:
%   X :  matrix.
% Output:
%   m :  number of rows in X.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Call MATLAB function size.
[m,n] = size(X);
