%......................................... 
% File:flowq4.m 
%
% Driver for linear analysis of potential
% problem with Q4P elements. 
%
% Required input:
%  T1:  Topology for 4-node elements
%  X:   Node coordinates
%  G:   Material properties
%  C:   Prescribed displacements
%  P:   Prescribed nodal loads
%  dof: No. of dof per node
%.........................................   

% Create new figure window 
figure

% Plot elements in topology matrix T1.
subplot(2,2,1); 
plotelem(T1,X); 
title('Element mesh');

% Initialization:  
%   K: Global conductivity
%   p: Global load vector
%   q: Global internal force vector 
[K,p,q] = init(rows(X),dof);

% Assembly of element conductivity matrices.
K = kq4p(K,T1,X,G);
% use sparse storage 
K = sparse(K); 

% Set nodal loads.
if exist('P')
  p = setload(p,P); 
end

% Set boundary conditions.
[K,p] = setbc(K,p,C);

% Solve equations: 
%  u: potential vector
u = K\p;

% Write potential. 
u

% plot potential contours 
subplot(2,2,2)
plotu(T1,X,u);
title('Potentials');

% Postprocessing: 
%  q: internal force vector
%  S: flux matrix
%  E: gradient matrix
[q,S,E] = qq4p(q,T1,X,G,u);

% contour plot of gradient components
subplot(2,2,3);
plotq4(T1,X,E,1)
title('Ex');
subplot(2,2,4);
plotq4(T1,X,E,2)
title('Ey');

