function Ke = keq4eps(Xe,Ge,Se,Ee,type)
%***************************************************
% KeQ4EPS: 
%   Creates the element tangent stiffness matrix 
%   of elasto-plastic quadrilateral 4-node element
%   in plane stress. 
%   The material  model is von Mises or associated 
%   Drucker-Prager with linear isotropic hardening. 
% Syntax: 
%   Ke = keq4eps(Xe,Ge,Se,Ee)
%   Ke = keq4eps(Xe,Ge,Se,Ee,type)
% Input:
%   Xe   : Coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4].
%   Ge   : Element material data Ge = [E nu S0 H phi]
%   Se   : Element stress matrix Se = [Sx Sy Sxy Seq 
%                                      ......        ].
%          with one row for each Gauss point.
%   Ee   : Element stress matrix Se = [Ex Ey 2Exy Ep 
%                                      ......        ].
%          with one row for each Gauss point. Seq is the
%          equivalent stress.  Ep is the total equivalent
%          plastic strain. 
%   type : parameter setting material model
%          type = 1 -> Von Mises (default)
%          type = 2 -> Drucker-Prager
% Output:
%   Ke   : Element tangent stiffness matrix.
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

% Set material parameters for isotropic plane stress
E  = Ge(1);  
nu = Ge(2);
S0 = Ge(3);
H  = Ge(4);
phi = 0.0;
if type==DruckerPrager & cols(Ge) == 5
  phi = Ge(5);
end 

% elastic stiffness
D  = E/(1-nu^2) ...   
      *[  1   nu   0
	  nu  1    0 
	  0   0  (1-nu)/2 ];  
% shear and bulk modulus
mu = E/(2*(1+nu));
k  = E/(3*(1-nu));

% Initialize stiffness matrix.
Ke = zeros(8);

% Gauss integration of tangent stiffness matrix.
for i = 1:2
for j = 1:2
  
  % organize the Gauss points as the element nodes
  gp = i + 3*(j-1) - 2*(i-1)*(j-1);

  % set up gradient matrix B
  dN = [ -(1-r(j))  (1-r(j))  (1+r(j)) -(1+r(j))  
	 -(1-r(i)) -(1+r(i))  (1+r(i))  (1-r(i)) ]/4;
  Jt = dN*Xe; 
  dN = Jt\dN;
  B  = [  dN(1,1)    0    dN(1,2)   0     dN(1,3)    0    dN(1,4)   0 
	     0    dN(2,1)    0    dN(2,2)    0    dN(2,3)    0    dN(2,4)
	  dN(2,1) dN(1,1) dN(2,2) dN(1,2) dN(2,3) dN(1,3) dN(2,4) dN(1,4) ];

  % updated yield stress
  Ep = Ee(gp,4);
  Sy = S0 + Ep*H;

  % deviatoric, mean and equivalent stress
  [Sd,Sm] = devstress(Se(gp,1:3)');
  Seq = Se(gp,4);

  % evaluate yield condition (phi=0 -> von Mises)
  f = Seq + phi*Sm - Sy; 

  % evaluate elasto-plastic stiffness...
  if  f < 0 

    Dep = D;

  else

    % evaluate vector a = D*d(Seq)/ds 
    a = [Sd(1)+nu*Sd(2),nu*Sd(1)+Sd(2),(1-nu)*Sd(3)]';
    a = (3*E/(2*Seq*(1-nu^2)))*a;

    % use simple version for von Mises or D-P with phi=0
    if phi == 0.0

      % scaling factor: (df/ds)'*D*(df/ds)
      factor = 3*mu*(1 - (1-2*nu)*(3*Sm)^2/(6*(1-nu)*Seq^2));

    else

      % a = D*df/ds = D*(d(Seq)/ds + phi*d(sm)/ds)
      a = a + k*phi*[1 1 0]';

      % gradient b = df/ds (strain format)
      b = (3/(2*Seq))*[Sd(1) Sd(2) 2*Sd(3)]' + phi/3*[1 1 0]';

      % scaling factor: (df/ds)'*D*(df/ds)
      factor = b'*a;

    end

    % elasto-plastic stiffness
    Dep = D - (a*a')/(H+factor);

  end 

  Ke = Ke + w(i)*w(j)*( B'*Dep*B )*det(Jt);
  
end
end



