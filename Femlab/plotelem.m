function plotelem(T,X,str,nonum)
%***************************************************
% PlotElem:
%   Plots elements in topology matrix T with
%   coordinate matrix X. Uses linear line segment
%   between all nodes.
% Syntax:
%   plotelem(T,X)
%   plotelem(T,X,'linetype')
%   plotelem(T,X,'linetype',nonum)
% Input:
%   T         :  element topology matrix.
%   X         :  node coordinate matrix.
%   'linetype':  string defining linetype ('y-','c:',...).
%   nonum     :  if nonum = 1 -> node numbers will be plotted 
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Close element contours for 3 or more nodes
if  cols(T) == 3
  T(:,3) = T(:,2);
else
  T(:,cols(T)) = T(:,1);
end

% Convert 1D problem to 2D.
if cols(X) == 1
  X = [X zeros(rows(X),1)]
end

% define line style and color
if nargin == 2
  str1 = 'yo';
  str2 = 'y-';
else
  if str(1) == ':' | str(1) == '-'  % check if line color is defined
    str1 = 'yo'; 
    str2 = ['y' str];
  else
    str1 = [str(1) 'o'];
    str2 = str;
  end
end
 
nnodes = cols(T)-1;

% Plot 2D elements by calling function 'plot'
if  cols(X) == 2
  order = [1:nnodes,1];
  if nnodes == 8
    % define node numbering for 8 node elements
    order = [1 5 2 6 3 7 4 8 1];
  end;
% Plot nodes to scale geometry  
  plot(X(:,1),X(:,2),str1)
  hold on

% Plot node numbers
  if nargin==4 
    if nonum==1
      for I=1:rows(X)
        text(X(I,1),X(I,2),int2str(I))
      end
    end
  end

% Plot 2D elements
  for j = 1:rows(T)
    plot(X(T(j,order),1),X(T(j,order),2),str2)
  end
end

% Plot 3D elements by calling function 'plot3'
if  cols(X) == 3
% Plot nodes to scale geometry 
  plot3(X(:,1),X(:,2),X(:,3),str1)
  hold on
% Plot 3D elements
  for j = 1:rows(T)
    plot3(X(T(j,:),1),X(T(j,:),2),X(T(j,:),3),str2)
  end
end

axis('equal')    % for use in ver. 4.0
axis('off') 

% enable the plot to be overwritten
hold off 
