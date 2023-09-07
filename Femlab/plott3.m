function plott3(T,X,S,comp)
%***************************************************
% PlotT3:
%   Plots a 2D color contour plot for the stress or
%   strain S of 3 node triangle elements, held in the 
%   topology matrix T, using the coordinate matrix X.
% Syntax:
%   plott3(T,X,S,comp)
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
ncomp = cols(S);

if comp > ncomp
  fprintf('\n\nRequested output component No. %d does not exist\n',comp);
  fprintf('Press any key to break plotting routine.\n\n');
  pause;
  break;
end;

% initial graphics commands
colormap(jet);
view(2);

% plot contours for all elements
for i=1:nel
  x(1:3,:) = X(T(i,1:3),:);
  x(4,:)   = X(T(i,1),:);
  s(1:3)   = zeros(3,1)+S(i,comp);  % artificial way to avoid loop
  s(4)     = S(i,comp);
  patch(x(:,1),x(:,2),s)
end
axis('equal');
axis('off');

hold off;


