function plotu(T,X,u)
%***************************************************
% PlotU:
%   Plots a 2D or 3D color contour plot for the 
%   nodal value u in all elements, held in the topology 
%   matrix T, using the coordinate matrix X. Use
%   COLORMENU to add a colormap menu to the figur.
% Syntax:
%   plotu(T,X,u)
% Input:
%   T         :  element topology matrix.
%   X         :  node coordinate matrix. (2D or 3D)
%   u         :  nodal value vector.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% number of elements in T
nel=rows(T);

% number of nodes per element
nno=cols(T)-1; 
order = [1:nno,1]; 
if nno==8
 order = [1 5 2 6 3 7 4 8 1];
end

% set spatial dimension
nx=cols(X);

% initial graphics commands
colormap(jet);
% define view point
if nx==2
  view(2);
else
  view(-30,max(X(:,3))*1.5);
end

% plot contours for all elements
for j=1:nel
  % define patch vertices and vertex intensities
  x(1:(nno+1),:) = X(T(j,order),:);  
  z(1:(nno+1))   = u(T(j,order));
  if nx==2 
    patch(x(:,1),x(:,2),z)
  else
    patch(x(:,1),x(:,2),x(:,3),z)
  end
end

axis('equal');
axis('off');
