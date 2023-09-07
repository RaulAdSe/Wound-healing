function n = cols(X)
%***************************************************
% cols:
%   Determines number of columns in matrix X.
% Syntax:
%   n = cols(X)
% Input:
%   X :  matrix.
% Output:
%   n :  number of columns in X.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Call MATLAB function size.
[m,n] = size(X);
