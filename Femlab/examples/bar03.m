%.......................................
% File bar03.m
%
% Example: 12-bar truss structure  
% Element type: BAR
% No. of Nodes: 9  
% No. of Elements : 12
%.......................................

% Coordinates of nodes,
X = [  -1.697 -1.0  0.0 
       -1.697  1.0  0.0
        0.0   -1.0  0.0
        0.0    1.0  0.0
        1.697 -1.0  0.0
        1.697  1.0  0.0
       -1.414  0.0  1.0
        0.0    0.0  1.0
        1.414  0.0  1.0  ];

% Topology matrix T(node1,node2,propno),
T = [ 1  7   1
      1  8   1 
      2  7   1 
      2  8   1
      3  8   1
      4  8   1
      5  8   1
      5  9   1
      6  8   1
      6  9   1
      7  8   1
      8  9   1  ];

% Element property matrix G = [ A E ],
G = [ 1.0  1.0 ];

% Prescribed loads
P = [ 7   0.00  0.00 -0.0225
      8   0.00  0.00 -0.015
      9   0.00  0.00 -0.0225 ];

% Boundary conditions   
C = [ 1  1  0.0
      1  2  0.0
      1  3  0.0
      2  1  0.0
      2  2  0.0
      2  3  0.0 
      3  1  0.0
      3  2  0.0
      3  3  0.0 
      4  1  0.0
      4  2  0.0
      4  3  0.0 
      5  1  0.0
      5  2  0.0
      5  3  0.0 
      6  1  0.0
      6  2  0.0 
      6  3  0.0 ];

%control parameters for non-linear analysis
no_loadsteps = 100;
% iteration limits
i_max = 10;
i_d   =  4;
% convergence threshold
TOL = 1E-3;

% plot parameters used by nlbar.m
% limits for element plot
elaxis = [-2 2 -1 1 -1.2 1.2];
% dof to be plotted in load-displacement curve
plotdof=21;
% axis load-displacement curve
plotaxis = [0 2 -.2 .2];
