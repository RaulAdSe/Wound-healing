%.........................................
% File:square.m 
%
% Example : Elastic or elasto-plastic 
%           tension test
% Element Type: Q4E, Q4EPE, Q4EPS 
% No. of Nodes: 9
% No. of Elements: 4
%......................................... 

% number of degrees of freedom per element node
dof = 2;

% Node coordinates
X = [   
      -1  -1
       0  -1
       1  -1
      -1   0
       0   0
       1   0
      -1   1
       0   1
       1   1 
    ];

% Element topology
T = [
       1 2 5 4   1
       2 3 6 5   1
       4 5 8 7   1
       5 6 9 8   1
    ];

% Prescribed displacements
C = [
       1  1  0
       1  2  0
       4  1  0
       4  2  0
       7  1  0
       7  2  0
    ];

% Material properties
%       E  nu  S0   H    phi
G = [ 100 0.3 1.00 0.00 0.50];

% Reference loading 
P = [ 
       3 0.025   0.00
       6 0.05    0.00
       9 0.025   0.00
    ];
elastic;
% control parameters for nonlinear analysis
% total No. of load steps
no_loadsteps = 30; 
% max No. of iterations before restart with smaller increment
i_max = 20;
% desired No. of iterations
i_d = 8;
% convergence threshold
TOL = 1e-2;

% limits for plots
elaxis = [-2 2 -2 2];
% dof to monitor as (u,f)
plotdof = 11;
% interval for stress/strain component e.g. equivalent stress
stressaxis = [-0.00001,1.1];
strainaxis = [-0.00001,0.11];
