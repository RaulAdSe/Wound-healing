function plotq4(T,X,S,comp)
%***************************************************
% PlotQ4:
%   Plots a 2D color contour plot for the stress or 
%   strain S of 4 node quadrilateral elements, held 
%   in the topology matrix T, using the coordinate
%   matrix X.
% Syntax:
%   plotq4(T,X,S,comp)
% Input:
%   T         :  element topology matrix.
%   X         :  node coordinate matrix. (2D or 3D)
%   S         :  global stress/strain matrix.
%   comp      :  stress/strain component to plot
% Date:
%   Version 1.0    04.05.95
%***************************************************

% number of elements in T
nel = rows(T);
% number of available stress/strain components
ncomp = cols(S)/4;

if comp > ncomp
  fprintf('\n\nRequested output component No. %d does not exist\n',comp);
  fprintf('Press any key to break plotting routine.\n\n');
  pause;
  return;
end;

% Gauss abscissae and weights.
r = [-1 1]*sqrt(3);
w = [ 1 1]; 
% Shape functions
for i=1:2
for j=1:2
  % organize the Gauss points as the element nodes
  gp = i + 3*(j-1) - 2*(i-1)*(j-1);
  N(gp,:)  = [ (1-r(i))*(1-r(j)) (1+r(i))*(1-r(j)) ...
               (1+r(i))*(1+r(j)) (1-r(i))*(1+r(j)) ]/4;
end
end

% Extrapolation from Gauss values to nodes
for i=1:nel
  s(1:4) = S(i,((1:4)-1)*ncomp+comp);
  for j=1:4
    Snodes(i,j) = N(j,:)*s';
  end 
end

% initial graphics commands
colormap(jet); 
view(2);

% plot contours for all elements
for i=1:nel
  hold on
  x(1:4,:) = X(T(i,1:4),:); 
  x(5,:)   = X(T(i,1),:);
  s(1:4)   = Snodes(i,1:4);
  s(5)     = Snodes(i,1);
  patch(x(:,1),x(:,2),s)
end 
axis('equal');
axis('off');

hold off;


