%...........................................................
% plastpe.m:
%   Driver for elasto-plastic plane strain analysis with
%   the orthogonal residual method.
%   
%   Input format: see hole.m
%   
%   NB.: The plotting facilities requires the full MATLAB 
%   installation and 256-color graphics driver.
%...........................................................

% Macro definitions for messages displayed during analysis   
format1 = 'Convergence NOT obtained in loadstep %d\n';    
format2 = 'Convergence obtained in loadstep %d in %d iterations\n';

% clear old version of plot variables
clear U F

% Initialize load and displacement vectors
dof = cols(X);
ndof = cols(X)*rows(X);
f = zeros(ndof,1);
df = zeros(ndof,1);
df = setload(df,P);
u = zeros(ndof,1);
du = zeros(ndof,1);

% Initialize stress and strain matrices
nelem  = rows(T);
S = zeros(nelem,1);
E = zeros(nelem,1);

% Position plot windows
figure(1)
figure(2)
fprintf('\n\n\nPosition figure windows on screen\n');   
fprintf('and press any key to continue\n\n'); 
pause;

% reset and set material model - default:von Mises 
clear mattype
fprintf('\n\n')
mattype = input('Choose material model (von Mises=1, Drucker-Prager=2): ');
fprintf('\n\n');
if mattype==[] | mattype > 2  
  mattype=1;
end

% reset counters
n = 1;
i = i_d;
i_tot=0;

% start loadsteps - Orthogonal residual method
while (n<=no_loadsteps)
  
  if i<i_max

    % new loadstep

    K = zeros(ndof);
    K = kq4epe(K,T,X,G,S,E,mattype);
    K = sparse(K);
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
      l = norm(du);
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
    [q,Sn,En] = qq4epe(q,T,X,G,u+du,S,E,mattype);
    dq = q - f;
    xi = (dq'*du)/(df'*du);
    r = - dq + xi*df;

    if rnorm(r,C,dof) < TOL*rnorm(df,C,dof)
      break
    else
      fprintf('Iteration No. %d \nResidual norm = %f\n',i,rnorm(r,C,dof));
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
    S = Sn;
    E = En;

    % plots...
    
    % load-displacement curve
    U(n+1) = u(plotdof);  
    F(n+1) = f(plotdof);
    figure(1)
    plot(U,F,U,F,'yx',U(n+1),F(n+1),'gx')
    
    % plot of equivalent plastic strain in deformaed geometry
    figure(2), clg;  
    U1 = reshape(u,cols(X),rows(X))';
    plotq4(T,X+U1,E,5);
    caxis(strainaxis);
    axis(elaxis);
    hold off

    i_tot = i_tot+i;
    n = n+1;
  end;
end; 

fprintf('Total No. of Iterations = %d\n',i_tot);
