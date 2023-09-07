%.........................................
% File:elastic.m 
%
% Driver for linear analysis of elastic 
% problem with Q4E elements. 
%
% Required input:
%  T:   Topology
%  X:   Node coordinates
%  G:   Material properties
%  C:   Prescribed displacements
%  P:   Prescribed nodal loads
%  dof: No. of dof per node
%.........................................   

% Initialization:  
%   K: Global stiffness matrix 
%   p: Global load vector
%   q: Global internal force vector 
[K,p,q] = init(rows(X),dof);

% Assembly of element stiffness matrices.
K = kq4e(K,T,X,G);

% Set nodal loads:
p = setload(p,P); 

% Set boundary conditions. 
[K,p] = setbc(K,p,C,dof);

% Solve equations:
%  u: displacement vector
u = K\p;

% Postprocessing: 
%  q: internal force vector
%  S: stress matrix
%  E: strain matrix
[q,S,E] = qq4e(q,T,X,G,u);

% Open figure window No. 1
figure(1)
% Plot elements and nodes of topology T. 
plotelem(T,X)
hold on
% Plot displaced configuration X+U (scaled)
U = reshape(u,cols(X),rows(X))';
plotelem(T,X+5*U,'c:')  
hold off

% Open figure window No. 2
figure(2)
% Contour plot of stress component No. 1: Sx
plotq4(T,X,S,1)

% extract reaction forces
H = reaction(q,C,dof,1);
V = reaction(q,C,dof,2);
