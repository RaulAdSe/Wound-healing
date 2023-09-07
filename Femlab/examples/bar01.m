%.......................................
% File bar01.m
%
% Example: Nonlinear 2D truss structure  
% Element type: BAR
% No. of Nodes: 3  
% No. of Elements : 2
%.......................................

% Coordinates of 3 nodes,
X = [ -1.00  0.00 
       1.00  0.00 
       0.00  0.10 ];

% Topology matrix T(node1,node2,propno),
T = [ 1  3   1
      2  3   1 ];
      
% Element property matrix G = [ A E ],
G = [ 1.0  1.0 ];

% Prescribed loads
P = [ 3   0.000 -0.0001 ];

% Boundary conditions   
C = [ 1  1  0.0
      1  2  0.0
      2  1  0.0
      2  2  0.0 ];

% control parameters for non-linear analysis
no_loadsteps   = 20;
% iteration limits
i_max = 8;
i_d   = 3;
% convergence threshold
TOL = 1E-3;

% plot parameters used in nlbar.m
% limits for element plot
elaxis = [-1 1 -.2 .2];
% axis load-displacement curve
plotaxis  = [0 0.25 -1e-3 1e-3];
% dof to be plotted in load-displacement curve
plotdof = 6;
