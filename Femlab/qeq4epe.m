function [qe,Se,Ee] = qeq4epe(Xe,Ge,Ue,Se,Ee,type)
%***************************************************
% QeQ4EPE: 
%   Evaluates the stress and strain for the current
%   displacements. Creates the element internal force 
%   vector for elasto-plastic quadrilateral 4-node
%   element in plane strain. 
%   The material model is von Mises or associated
%   Drucker-Prager with linear isotropic hardening. 
% Syntax: 
%   [qe,Se,Ee] = qeq4epe(Xe,Ge,Ue,Se,Ee)
%   [qe,Se,Ee] = qeq4epe(Xe,Ge,Ue,Se,Ee,type)
% Input:
%   Xe   : Coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4]
%   Ge   : Element material data: [E , nu, Sy, H, phi]
%   Ue   : Element nodal displacements [u1; v1; u2; v2...]
%   Se   : Element stress matrix: [ Sx Sy Sz Sxy Seq 
%                                  .....            ]
%   Ee   : Element strain matrix: [ Ex Ey Ez Exy Ep
%                                  .....            ]
%          with one row for each Gauss point. Seq is the
%          consistent equivalent stress. Ep is the total
%          equivalent plastic strain. 
%   type : parameter setting material model
%          type = 1 -> Von Mises (default)
%          type = 2 -> Drucker-Prager
% Output:  
%   qe   : Element internal force vector.
%   Se   : Updated element stresses.
%   Ee   : Updated element strains.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% Define material models
VonMises = 1;
DruckerPrager = 2;

% default material model: von Mises
if nargin==5
  type = VonMises;
end

% Gauss abscissae and weights.
r = [-1 1]/sqrt(3);
w = [ 1 1]; 

% "reduced" integration of dilatation part of B (r=(0,0))
dN = [ -1  1  1 -1  
       -1 -1  1  1 ]/4;
Jt = dN*Xe; 
dN = Jt\dN;
W = [dN(1,1) dN(2,1) dN(1,2) dN(2,2) ... 
     dN(1,3) dN(2,3) dN(1,4) dN(2,4) ]; 
m = [1 1 0 0]';
I = eye(4);

% Set material parameters for isotropic plane strain
E  = Ge(1);  
nu = Ge(2);
S0 = Ge(3);
H  = Ge(4);
phi = 0.0;
if type==DruckerPrager & cols(Ge) == 5
  phi = Ge(5);
end

% elastic stiffness
D  = E/((1+nu)*(1-2*nu)) ...   
      *[ 1-nu   nu    nu      0
	  nu   1-nu   nu      0 
          nu    nu   1-nu     0
	  0     0     0   (1-2*nu)/2 ];  
% shear and bulk modulus
mu = E/(2*(1+nu));
k = E/(3*(1-2*nu));

% Initialize force vector
qe = zeros(8,1);

% Gauss integration of internal force vector.
for i = 1:2
for j = 1:2
  
  gp = i + 3*(j-1) - 2*(i-1)*(j-1);

  % evaluate strain matrix B
  dN = [ -(1-r(j))  (1-r(j))  (1+r(j)) -(1+r(j))  
	 -(1-r(i)) -(1+r(i))  (1+r(i))  (1-r(i)) ]/4;
  Jt = dN*Xe;
  dN = Jt\dN;
  B  = [  dN(1,1)    0    dN(1,2)   0     dN(1,3)    0    dN(1,4)   0 
	     0    dN(2,1)    0    dN(2,2)    0    dN(2,3)    0    dN(2,4)
             0       0       0      0        0       0       0      0
	  dN(2,1) dN(1,1) dN(2,2) dN(1,2) dN(2,3) dN(1,3) dN(2,4) dN(1,4) ];
  B = (I - 0.5*m*m')*B + 0.5*m*W;

  % evaluate total strain, En, and strain increment, dE.
  e = B*Ue;
  dE = e - Ee(gp,1:4)';
  Ee(gp,1:4) = e';

  % evaluate elastic trial stress 
  dS = D*dE;    %dS: stress increment, D:tensor of elastic moduli, dE: increment of strains. We assume elasticity: thus the update of stresses i proportional to increment of strains
  S = Se(gp,1:4)' + dS(1:4);

  % updated yield stress
  Ep = Ee(gp,5);
  Sy = S0 + Ep*H;

  % deviatoric, mean and equivalent stress 
  [Sd,Sm] = devstress(S);
  Seq = eqstress(S);

  % evaluate yield condition (phi=0 -> von Mises)
  f = Seq + phi*Sm - Sy;

  % determine new stress and plastic multiplier
  dL = 0;
  if f >= 0
    dL = f/(H + 3*mu + k*phi^2);
    mp = [1 1 1 0]';
    S = S - dL*(3*mu*Sd/Seq + k*phi*mp);    
  end
  
  % update total stresses
  Se(gp,1:4) = S(1:4)';
  % update equivalent stress
  Se(gp,5) = eqstress(S); 

  % update total equivalent plastic strain
  Ee(gp,5) = Ee(gp,5) + dL;

  % update internal force
  qe = qe + w(i)*w(j)*( B'*S )*det(Jt);   
  
end
end



