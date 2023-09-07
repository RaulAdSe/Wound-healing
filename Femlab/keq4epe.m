function Ke = keq4epe(Xe,Ge,Se,Ee,type)
%***************************************************
% KeQ4EPE: 
%   Creates the element tangent stiffness matrix 
%   of elasto-plastic quadrilateral 4-node element
%   in plane strain. The strain gradient B is re-
%   placed by a modified B using reduced integration
%   of the dilatational components. 
%   The material model is von Mises or associated 
%   Drucker-Prager with linear isotropic hardening. 
% Syntax: 
%   Ke = keq4epe(Xe,Ge,Se,Ee)
%   Ke = keq4epe(Xe,Ge,Se,Ee,type)
% Input:
%   Xe   : Coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4].
%   Ge   : Element material data Ge = [E nu Sy H phi]
%   Se   : Element stress matrix Se = [Sx Sy Sz Sxy Seq 
%                                      ......          ].
%   Ee   : Element strain matrix Ee = [Ex Ey Ez 2Exy Ep
%                                      ......             ].
%          with one row for each Gauss point. Seq is the 
%          equivalent stress, Ep is the total equivalent 
%          plastic strain.
%   type : parameter setting material model
%          type = 1 -> Von Mises (default)
%          type = 2 -> Drucker-Prager
% Output:
%   Ke   : Element stiffness matrix.
% Date:
%   Version 1.0    04.05.95
%***************************************************

% material model identifiers
VonMises = 1;
DruckerPrager = 2;

% default material model: von Mises
if nargin==4
  type = VonMises;
end

% Gauss abscissae and weights.
r = [-1 1]/sqrt(3);
w = [ 1 1]; 

% "reduced" integration of dilatation part of B
% r = (0,0)
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
k  = E/(3*(1-nu));

% Initialize stiffness matrix.
Ke = zeros(8);

% Gauss integration of stiffness matrix.
for i = 1:2
for j = 1:2
  
  % organize the Gauss points as the element nodes
  gp = i + 3*(j-1) - 2*(i-1)*(j-1);

  % set up gradient matrix
  dN = [ -(1-r(j))  (1-r(j))  (1+r(j)) -(1+r(j))  
	 -(1-r(i)) -(1+r(i))  (1+r(i))  (1-r(i)) ]/4;
  Jt = dN*Xe; 
  dN = Jt\dN;
  B  = [  dN(1,1)    0    dN(1,2)   0     dN(1,3)    0    dN(1,4)   0 
	     0    dN(2,1)    0    dN(2,2)    0    dN(2,3)    0    dN(2,4)
             0       0       0      0        0       0       0      0
	  dN(2,1) dN(1,1) dN(2,2) dN(1,2) dN(2,3) dN(1,3) dN(2,4) dN(1,4) ];
  B = (I - 0.5*m*m')*B + 0.5*m*W;

  % updated yield stress
  Ep = Ee(gp,5);
  Sy = S0 + Ep*H;

  % deviatoric, mean and equivalent stress
  [Sd,Sm] = devstress(Se(gp,1:4)');
  Seq = Se(gp,5);

  % evaluate yield condition (phi=0 -> von Mises)
  f = Seq + phi*Sm - Sy; 

  % evaluate elasto-plastic stiffness...
  if  f < 0

    Dep = D;

  else

    % evaluate deviatoric part 
    Dp = 9*mu^2*(Sd*Sd')/(Seq^2);

    % evaluate mean stress part (D-P materials)
    if phi ~= 0.0
      mp = [1 1 1 0]';
      Dp = Dp + (k*phi)^2*(mp*mp');
    end 

    % scaling factor: (df/ds)'*D*(df/ds)
    factor = 3*mu + k*phi^2;

    % elasto-plastic stiffness
    Dep = D - Dp/(H+factor);

  end 

  Ke = Ke + w(i)*w(j)*( B'*Dep*B )*det(Jt);

end
end



