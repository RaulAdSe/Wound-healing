%.........................................
% File:canti.m 
%
% Example : Bending of cantilever beam. 
% Element Type: Q4E 
% No. of Nodes: 27
% No. of Elements: 16
%......................................... 

% Topology for 4-node element    
%     Node1 Node2 Node3 Node4    PropNo
T = [  1     4     5     2         1 
       2     5     6     3         1
       4     7     8     5         1
       5     8     9     6         1
       7    10    11     8         1
       8    11    12     9         1 
      10    13    14    11         1
      11    14    15    12         1 
      13    16    17    14         1
      14    17    18    15         1
      16    19    20    17         1 
      17    20    21    18         1
      19    22    23    20         1
      20    23    24    21         1
      22    25    26    23         1
      23    26    27    24         1  ];

% Node coordinates  
%      X      Y
X = [  0      0 
       0      0.5
       0      1
       0.5    0 
       0.5    0.5 
       0.5    1
       1      0 
       1      0.5
       1      1 
       1.5    0 
       1.5    0.5
       1.5    1 
       2      0 
       2      0.5
       2      1 
       2.5    0 
       2.5    0.5
       2.5    1 
       3      0
       3      0.5 
       3      1 
       3.5    0
       3.5    0.5 
       3.5    1 
       4      0 
       4      0.5
       4      1     ];

% No. of dof per node (Q4E element)  
dof = 2;

% Material properties (Elastic element)
%       E     nu   type = plane stress
G = [  100    0.3    1  ];

% Constrained dof 
%    NodeNo   DofNo   U
C = [  1        1     0.0
       2        1     0.0
       2        2     0.0 
       3        1     0.0 ];

% Nodal loads 
%    NodeNo    Px    Py
P = [  25      0.0   -0.05   
       26      0.0   -0.10
       27      0.0   -0.05 ];
