%...........................................
% File: hole.m:
%
% Example : Plate with hold (1/4 modelled).
% Element Types: Q4E, Q4EPE, Q4EPS 
% No. of Nodes: 65
% No. of Elements: 48
%........................................... 

% reference nodal load 
P0 = 0.02;

% number of degrees of freedom per element node
dof = 2;

% material data: 
%      E      nu     S0   H=0.1*E   phi(DP)
G = [100.0   0.3    1.0    10.0      0.5  ];

% symmetry constraints
C = [ 
       1  2   0.0
      10  2   0.0
      19  2   0.0
      28  2   0.0
      37  2   0.0
      46  2   0.0
      51  2   0.0
      56  2   0.0 
      61  2   0.0
       9  1   0.0
      18  1   0.0
      27  1   0.0
      36  1   0.0
      45  1   0.0
    ];

% reference load
P = [  61   0.5*P0   0.0
       62   P0       0.0
       63   P0       0.0
       64   P0       0.0
       65   0.5*P0   0.0 ];


% parametric definition of node coordinates

a  = 1;
ba = 1.25*a;
b  = 1.5*a;
bb = 1.75*a;
c  = 2*a;

X = [ 
      a                      0
      cos(pi/16)*a           sin(pi/16)*a
      cos(pi/8)*a            sin(pi/8)*a
      cos(3*pi/16)*a         sin(3*pi/16)*a
      cos(pi/4)*a            sin(pi/4)*a
      cos(5*pi/16)*a         sin(5*pi/16)*a
      cos(3*pi/8)*a          sin(3*pi/8)*a
      cos(7*pi/16)*a         sin(7*pi/16)*a
      0                      a
      ba                     0
      cos(pi/16)*ba          sin(pi/16)*ba
      cos(pi/8)*ba*1.025     sin(pi/8)*ba*1.025
      cos(3*pi/16)*ba*1.05   sin(3*pi/16)*ba*1.05
      cos(pi/4)*ba*1.075     sin(pi/4)*ba*1.075
      cos(5*pi/16)*ba*1.05   sin(5*pi/16)*ba*1.05
      cos(3*pi/8)*ba*1.025   sin(3*pi/8)*ba*1.025
      cos(7*pi/16)*ba        sin(7*pi/16)*ba
      0                      ba
      b                      0
      cos(pi/16)*b           sin(pi/16)*b
      cos(pi/8)*b*1.05       sin(pi/8)*b*1.05
      cos(3*pi/16)*b*1.1     sin(3*pi/16)*b*1.1
      cos(pi/4)*b*1.175      sin(pi/4)*b*1.175
      cos(5*pi/16)*b*1.1     sin(5*pi/16)*b*1.1
      cos(3*pi/8)*b*1.05     sin(3*pi/8)*b*1.05
      cos(7*pi/16)*b         sin(7*pi/16)*b
      0                      b
      bb                     0
      cos(pi/16)*bb          sin(pi/16)*bb
      cos(pi/8)*bb*1.05      sin(pi/8)*bb*1.05
      cos(3*pi/16)*bb*1.15   sin(3*pi/16)*bb*1.15
      cos(pi/4)*bb*1.275     sin(pi/4)*bb*1.275
      cos(5*pi/16)*bb*1.15   sin(5*pi/16)*bb*1.15
      cos(3*pi/8)*bb*1.05    sin(3*pi/8)*bb*1.05
      cos(7*pi/16)*bb        sin(7*pi/16)*bb
      0                      bb
      c                      0
      c                      a/2*0.8
      c                      a*0.85
      c                      3*a/2*0.9
      c                      c
      3*a/2*0.9              c
      a*0.85                 c
      a/2*0.8                c
      0                      c 
      3*a                    0
      3*a                    a/2         
      3*a                    a  
      3*a                    3*a/2       
      3*a                    2*a         
      4*a                    0
      4*a                    a/2         
      4*a                    a  
      4*a                    3*a/2       
      4*a                    2*a         
      5*a                    0
      5*a                    a/2         
      5*a                    a  
      5*a                    3*a/2       
      5*a                    2*a         
      6*a                    0
      6*a                    a/2         
      6*a                    a  
      6*a                    3*a/2       
      6*a                    2*a         
    ];

% topology 
T = [1 10 11  2    1 
     2 11 12  3    1
     3 12 13  4    1
     4 13 14  5    1
     5 14 15  6    1
     6 15 16  7    1
     7 16 17  8    1
     8 17 18  9    1

    10 19 20 11    1
    11 20 21 12    1
    12 21 22 13    1
    13 22 23 14    1
    14 23 24 15    1
    15 24 25 16    1
    16 25 26 17    1
    17 26 27 18    1

    19 28 29 20    1
    20 29 30 21    1
    21 30 31 22    1
    22 31 32 23    1
    23 32 33 24    1
    24 33 34 25    1
    25 34 35 26    1
    26 35 36 27    1

    28 37 38 29    1
    29 38 39 30    1
    30 39 40 31    1
    31 40 41 32    1 
    32 41 42 33    1
    33 42 43 34    1
    34 43 44 35    1
    35 44 45 36    1

    37 46 47 38    1
    38 47 48 39    1
    39 48 49 40    1
    40 49 50 41    1
    
    46 51 52 47    1
    47 52 53 48    1
    48 53 54 49    1
    49 54 55 50    1

    51 56 57 52    1
    52 57 58 53    1
    53 58 59 54    1
    54 59 60 55    1
 
    56 61 62 57    1
    57 62 63 58    1
    58 63 64 59    1
    59 64 65 60    1
 ];

% control parameters for nonlinear analysis
% total No. of load steps
no_loadsteps = 20; 
% max No. of iterations before restart with smaller increment
i_max = 20;
% desired No. of iterations
i_d = 8;
% convergence threshold
TOL = 1e-2;

% limits for plots
elaxis = [0 7 0 2];
% dof to monitor as (u,f)
plotdof = 123;
% interval for stress/strain component e.g. equivalent stress
stressaxis = [-0.00001,1.1];
strainaxis = [-0.00001,0.11];


