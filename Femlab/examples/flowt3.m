%......................................... 
% File:flowt3.m 
%
% Driver for linear analysis of potential
% problem with T3P elements. 
%
% Required input:
%  T2:  Topology for 3-node elements
%  X:   Node coordinates
%  G:   Material properties
%  C:   Prescribed displacements
%  P:   Prescribed nodal loads
%  dof: No. of dof per node
%.........................................   

% Create new figure window 
figure

% Plot elements in topology matrix T2.
subplot(2,2,1); 
plotelem(T2,X); 
title('Element mesh');

% Initialization:  
%   K: Global conductivity
%   p: Global load vector
%   q: Global internal force vector 
[K,p,q] = init(rows(X),dof);

% Assembly of element conductivity matrices.
K = kt3p(K,T2,X,G);
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
plotu(T2,X,u);
title('Potentials');

% Postprocessing: 
%  q: internal force vector
%  S: flux matrix
%  E: gradient matrix
[q,S,E] = qt3p(q,T2,X,G,u);

% contour plot of gradient components
subplot(2,2,3);
plott3(T2,X,E,1)
title('Ex');
subplot(2,2,4);
plott3(T2,X,E,2)
title('Ey');
