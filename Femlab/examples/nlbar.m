%.........................................
% File: nlbar.m 
%
% Driver for nonlinear analysis of elastic 
% problem with BAR elements. The orthogonal 
% residual method is used for equilibrium
% iterations, see e.g. Krenk (1993): Non-
% linear analysis with finite elements.  
%
% Required input:
%  T:   Topology
%  X:   Node coordinates
%  G:   Material properties
%  C:   Prescribed displacements
%  P:   Prescribed nodal loads
%  dof: No. of dof per node
%
%  no_loadsteps: Total No. of loadsteps 
%  i_max:        Max. No. of iterations before restart
%  i_d:          Desired No. of iterations
% 
%  plotaxis:  Limits for load-displacement curve 
%  plotdof:   No. of dof to be monitored 
%  elaxis:    Limits for element plot
%.........................................    

% Macro definitions for messages displayed during analysis   
format1 = 'Convergence NOT obtained in loadstep %d\n';    
format2 = 'Convergence obtained in loadstep %d in %d iterations\n';

% plots...

% remove old versions of plot variables
clear X1 U F

% Define title on figures
txt = sprintf('Loadstep No. 0'); 

% Plot load-displacement curve 
figure(1)
axis(plotaxis);
title(txt);
hold off;

% Plot geometry
figure(2);  
plotelem(T,X);
axis(elaxis);
title(txt);
hold off;

% Allow the user to position windows
fprintf('\n\n\nPosition figures on screen\n');   
fprintf('and press any key to continue\n\n'); 
pause;

% Initialize load and displacement vectors
dof = cols(X); 
ndof = cols(X)*rows(X);
u  = zeros(ndof,1);
du = zeros(ndof,1);
f  = zeros(ndof,1);
df = zeros(ndof,1);
df = setload(df,P);   

% reset counters
n = 1;
i = i_d;

% start loadsteps - Orthogonal residual method
while (n<=no_loadsteps)
  
  if i<i_max  

    % new loadstep

    K = zeros(ndof);
    K = kbar(K,T,X,G,u);
    [Kt,df] = setbc(K,df,C,dof);
    du0 = Kt\df;

    % check for new loading direction
    if du'*du0 < 0
      df = -df;
      du0 = -du0;
    end
 
    % define step size measures
    if n==1
      l0=norm(du0);
      l = l0;
      l_max = 2*l0;
    else
      l  = norm(du);
      l0 = norm(du0);
    end;
  end;  

  % step size adjustment
  if i_d<=i & i<i_max                % normal convergence
    du = min(l/l0,l_max/l0)*du0;
  elseif i<i_d                       % fast convergence
    du = min(2*l/l0,l_max/l0)*du0;
  else                               % restart
    du0 = 0.5*du0;
    du  = du0;
  end; 

  % equilibrium iterations...
  for i=1:i_max
    q = zeros(ndof,1);
    q = qbar(q,T,X,G,u+du);
    dq = q - f;
    xi = (dq'*du)/(df'*du);
    r = - dq + xi*df;
    if rnorm(r,C,dof) < TOL*rnorm(df,C,dof)
      break
    else
      [Kt,r] = setbc(K,r,C,dof);
      delta_u = Kt\r;
      du = du + delta_u;
    end
  end

  if i>=i_max
    % convergence is not reached - restart
    fprintf(format1,n);
  else
    % convergence is reached - update u,f,n
    fprintf(format2,n,i);
    f = f + xi*df;
    u = u + du;
    
    % plots...    

    X1 = X + reshape(u,cols(X),rows(X))';
    F(n+1) = f(plotdof);
    U(n+1) = u(plotdof);

    txt = sprintf('Loadstep No. %d',n);
    figure(1); 
    plot(-U,-F,-U(n+1),-F(n+1),'cx'); 
    axis(plotaxis);
    title(txt);
    hold off; 
    
    figure(2); 
    plotelem(T,X1);
    axis(elaxis);
    title(txt);
    hold off;

    n = n+1;  
  end;
end; 

