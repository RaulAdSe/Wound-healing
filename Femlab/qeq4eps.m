function [qe,Se,Ee] = qeq4eps(Xe,Ge,Ue,Se,Ee,type)
%***************************************************
% QeQ4EPS: 
%   Evaluates the stress and strain for the current
%   displacements. Creates the element internal force 
%   vector for elasto-plastic quadrilateral 4-node 
%   element in plane stress. 
%   The material model is von Mises or associated
%   Drucker-Prager with linear isotropic hardening. 
% Syntax: 
%   [qe,Se,Ee] = qeq4eps(Xe,Ge,Ue,Se,Ee)
%   [qe,Se,Ee] = qeq4eps(Xe,Ge,Ue,Se,Ee,type)
% Input:
%   Xe   : Coordinates Xe = [x1 y1; x2 y2; x3 y3; x4 y4]
%   Ge   : Element material data: [E , nu, S0, H, phi]
%   Ue   : Element nodal displacements [u1; v1; u2; v2...]
%   Se   : Element stress matrix: [ Sx Sy Sxy Seq
%                                   .....          ]
%   Ee   : Element strain matrix: [ Ex Ey 2Exy Ep
%                                   .....          ]
%          with one row for each Gauss point. Seq is the
%          equivalent stress. Ep is the total equivalent
%          plastic strain.
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

% Set material parameters for isotropic plane stress
E  = Ge(1);  
nu = Ge(2);
S0 = Ge(3);
H  = Ge(4);
phi = 0.0;
if type==DruckerPrager & cols(Ge) == 5
  phi = Ge(5);
end

% elastic stiffness matrix
D  = E/(1-nu^2) ...   
      *[   1  nu    0
	  nu   1    0 
	   0   0   (1-nu)/2 ];  

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
	  dN(2,1) dN(1,1) dN(2,2) dN(1,2) dN(2,3) dN(1,3) dN(2,4) dN(1,4) ];

  % evaluate total strain, En, and strain increment, dE.
  e = B*Ue;
  dE = e - Ee(gp,1:3)';
  Ee(gp,1:3) = e';

  % evaluate elastic trial stress 
  dS = D*dE; 
  S = Se(gp,1:3)' + dS(1:3);

  % updated yield stress
  Ep = Ee(gp,4);
  Sy = S0 + Ep*H;

  % deviatoric, mean and equivalent stress
  [Sd,Sm] = devstress(S);
  Seq = eqstress(S);

  % evaluate yield condition (phi=0 -> von Mises)
  f = Seq + phi*Sm - Sy; 

  % determine new stress and plastic multiplier
  dL = 0;
  if f >= 0
    if phi == 0.0       % von Mises 
      [S,dL] = stressvm(S,Ge,Sy);
    else                % Drucker-Prager
      [S,dL] = stressdp(S,Ge,Sy,dE,dS);
    end 
  end
  
  % update total stresses
  Se(gp,1:3) = S(1:3)';
  % update equivalent stress
  Se(gp,4) = eqstress(S); 

  % update total equivalent plastic strain 
  Ee(gp,4) = Ee(gp,4) + dL; 

  % update internal force
  qe = qe + w(i)*w(j)*( B'*S )*det(Jt);   

end
end
