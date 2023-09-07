%.......................................
% File bar02.m
%
% Example: Nonlinear truss structure  
% Element type: BAR
% No. of Nodes: 4  
% No. of Elements : 3
%.......................................

% Coordinates of 4 nodes,
X = [  0.00  0.00  0.00
       5.00  0.00  0.00
       2.50  4.33  0.00
       2.50  2.165 0.1 ];

% Topology matrix T(node1,node2,propno),
T = [ 1  4   1
      2  4   1 
      3  4   1 ];
      
% Element property matrix G = [ A E ],
G = [ 1.0  1.0 ];

% Prescribed loads
P = [ 4   0.00  0.00 -0.00001 ];

% Boundary conditions   
C = [ 1  1  0.0
      1  2  0.0
      1  3  0.0
      2  1  0.0
      2  2  0.0
      2  3  0.0 
      3  1  0.0
      3  2  0.0 
      3  3  0.0 ];


% control parameters for non-linear analysis
no_loadsteps   = 15;
% iteration limits
i_d   = 3;
i_max = 10; 
% convergence threshold
TOL = 1E-3;

% plot parameters used in nlbar.m
% limits for element plot
elaxis = [0 5 0 5 -.2 .2];
% axis load-displacement curve
plotaxis  = [0 0.25 -4e-5 4e-5];
% dof to be plotted in load-displacement curve
plotdof = 12;
