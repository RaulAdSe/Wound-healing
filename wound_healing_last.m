 %..................................................
% Linear Elasticity example with a single element
%.................................................. 
clear all; close all; clc;
%%
addpath '..\Femlab'; % Path to functions
% Number of dof/node
dof = 2;
E = logspace(0,2,20);
%
% GEOMETRY DEFINITION
%
% Nodal coordinates
t = 0.03;
L = 1;
X = [];
elx = 50; 
ely = 50;
for ii = 0:elx
    for jj = 0:ely
        X = [X; L/2*ii/elx L/2*jj/ely];
    end
end

% Topology, element definition
% (nodes, material)
T = [];
for ii = 1:elx
    for jj = 1:ely
        T = [T; ((ii-1)*(ely+1))+jj ii*(ely+1)+jj ii*(ely+1)+jj+1 ((ii-1)*(ely+1))+jj+1 1];
    end
end

%% Variable Young modulus
for kk = 1:20
%
%
% MATERIAL PROPERTIES
%       E     nu   hypothesis (plane stress=1, plane strain=2)
G = [  E(kk)    0.3    1  ];


%
% BOUNDARY CONDITIONS
%
% Prescribed displacements
%(node, dof , prescribed value)
C = [];
for ii = 1:length(X)
     if X(ii,1) > 0.4999
        if X(ii,2) > L/4
            C = [C; ii 1 0.0];
        end
     end
     if X(ii,2) < 0.00001
        C = [C; ii 2 0.0];
        if X(ii,1) > 0.4999
            nslit = ii;
        end
     end
end
        

% Nodal loads
%    node    load
P = [];
for ii = 1:length(X)
    if X(ii,1) < 0.0000001
        if X(ii,2) > 0.499 || X(ii,2) < 0.0001 
            P = [P; ii -t/2/ely 0.0];
        else
            P = [P; ii -t/ely 0.0];
        end
    end
end
        
% SOLVE: equilbirum
%
% Initialization
%   K: stiffness matrix
%   p: external load vector
%   q: internal load vector 
[K,p,q] = init(rows(X),dof);

% Build global stiffness matrix
K = kq4e(K,T,X,G);

% Applied loads
p = setload(p,P); 

% Impose prescribed values
Penalty=false;

if Penalty
    [K,p] = setbc(K,p,C,dof);
else
    cdof=(C(:,1)-1)*dof+C(:,2);
    p=p-K(:,cdof)*C(:,3);
    K(cdof,:)=[];
    p(cdof)=[];
    K(:,cdof)=[];
end

% Solve system
%  u: displacement vector
u = K\p;
if ~Penalty
    ndof=1:length(q);
    ndof(cdof)=[];
    uf=u;
    u=0*q;
    u(ndof)=uf;
    u(cdof)=C(:,3);
end
U = reshape(u,cols(X),rows(X))';

d(kk) = U(nslit,1);
end
figure()
plot(E,-2*d, 'LineWidth', 1.5)
xlabel('$\textbf{\textit{E}}$', 'Interpreter', 'latex')
ylabel('$\textbf{\textit{d}}$', 'Interpreter', 'latex')
% xlabel("Young Modulus")
% ylabel("Wound opening")
%title("Wound opening as a function of E")

%% Variable boundary traction
t = linspace(0,1,20);
E = 20;
for kk = 1:20
% MATERIAL PROPERTIES
%       E     nu   hypothesis (plane stress=1, plane strain=2)
G = [  E    0.3    1  ];


%
% BOUNDARY CONDITIONS
%
% Prescribed displacements
%(node, dof , prescribed value)
C = [];
for ii = 1:length(X)
     if X(ii,1) > 0.4999
        if X(ii,2) > L/4
            C = [C; ii 1 0.0];
        end
     end
     if X(ii,2) < 0.00001
        C = [C; ii 2 0.0];
        if X(ii,1) > 0.4999
            nslit = ii;
        end
     end
end
        

% Nodal loads
%    node    load
P = [];
for ii = 1:length(X)
    if X(ii,1) < 0.0000001
        if X(ii,2) > 0.4999 || X(ii,2) < 0.0001 
            P = [P; ii -0.5*t(kk)/ely 0.0];
        else
            P = [P; ii -t(kk)/ely 0.0];
        end
    end
end
        

%
% SOLVE: equilbirum
%
% Initialization
%   K: stiffness matrix
%   p: external load vector
%   q: internal load vector 
[K,p,q] = init(rows(X),dof);

% Build global stiffness matrix
K = kq4e(K,T,X,G);

% Applied loads
p = setload(p,P); 

% Impose prescribed values
Penalty=false;

if Penalty
    [K,p] = setbc(K,p,C,dof);
else
    cdof=(C(:,1)-1)*dof+C(:,2);
    p=p-K(:,cdof)*C(:,3);
    K(cdof,:)=[];
    p(cdof)=[];
    K(:,cdof)=[];
end

% Solve system
%  u: displacement vector
u = K\p;
if ~Penalty
    ndof=1:length(q);
    ndof(cdof)=[];
    uf=u;
    u=0*q;
    u(ndof)=uf;
    u(cdof)=C(:,3);
end
U = reshape(u,cols(X),rows(X))';

d(kk) = U(nslit,1);
end
figure()
plot(t,-2*d, 'LineWidth', 1.5)
xlabel("Traction")
ylabel("Wound opening")
xlabel('$\left| \textbf{\textit{t}} \right|$', 'Interpreter', 'latex')
ylabel('$\textbf{\textit{d}}$', 'Interpreter', 'latex')
% title("Wound opening as a function of the boundary traction")

%% Variable wound length

gamma = linspace(0,1,50);
E = 20;

Gs = [  10    0.3    1  ; 20    0.3     1; 40   0.3     1;40   0.3     1];
ts = [1, 1, 0.5, 2];
lineStyles = {'-', '-', '-','--'};
colors = {[0, 0.4470, 0.7410], [0.9290, 0.6940, 0.1250], [0.8500, 0.3250, 0.0980], [0, 0,0]};


figure()
for index = 1:4
    for kk = 1:50
    % MATERIAL PROPERTIES
    %       E     nu   hypothesis (plane stress=1, plane strain=2)
    G = Gs(index,:);


    %
    % BOUNDARY CONDITIONS
    %
    % Prescribed displacements
    %(node, dof , prescribed value)
    C = [];
    for ii = 1:length(X)
         if X(ii,1) > 0.4999
            if X(ii,2) > L/2*gamma(kk)-0.00001
                C = [C; ii 1 0.0];
            end
         end
         if X(ii,2) < 0.00001
            C = [C; ii 2 0.0];
            if X(ii,1) > 0.4999
                nslit = ii;
            end
         end
    end


    % Nodal loads
    %    node    load
    P = [];
    for ii = 1:length(X)
        if X(ii,1) < 0.0000001
            if X(ii,2) > 0.4999 || X(ii,2) < 0.0001     
                P = [P; ii -0.5*ts(index)/ely 0.0];
            else
                P = [P; ii -ts(index)/ely 0.0];
            end
        end
    end

    % SOLVE: equilbirum
    %
    % Initialization
    %   K: stiffness matrix
    %   p: external load vector
    %   q: internal load vector 
    [K,p,q] = init(rows(X),dof);

    % Build global stiffness matrix
    K = kq4e(K,T,X,G);

    % Applied loads
    p = setload(p,P); 

    % Impose prescribed values
    Penalty=false;

    if Penalty
        [K,p] = setbc(K,p,C,dof);
    else
        cdof=(C(:,1)-1)*dof+C(:,2);
        p=p-K(:,cdof)*C(:,3);
        K(cdof,:)=[];
        p(cdof)=[];
        K(:,cdof)=[];
    end

    % Solve system
    %  u: displacement vector
    u = K\p;
    if ~Penalty
        ndof=1:length(q);
        ndof(cdof)=[];
        uf=u;
        u=0*q;
        u(ndof)=uf;
        u(cdof)=C(:,3);
    end
    U = reshape(u,cols(X),rows(X))';
    d(kk) = U(nslit,1);
    end
    plot(gamma, -2*d, lineStyles{index}, 'LineWidth', 1.5, 'color',colors{index});
    hold on
end
xlabel('$\textbf{\textit{l}}$', 'Interpreter', 'latex')
ylabel('$\textbf{\textit{d}}$', 'Interpreter', 'latex')
legend({'$E = 10$, $|t| = 1$','$E = 20$, $|t| = 1$','$E = 40$, $|t| = 0.5$','$E = 40$, $|t| = 2$'}, 'Interpreter', 'latex', 'FontSize', 12)
% title("Wound opening as a function of the wound height")

%%
% POSTPROCESS: 
%  
%  S: stresses
%  E: strains
[q1,S1,E1] = qq4e(q,T,X,G,u);

% Graphical representation of the initial and deformed representation
%Y = ones(size(X))./2;
Y = [ones([length(X),1]),zeros([length(X),1])];
figure(4)
plotelem(T,X,'b-')
hold on
U = reshape(u,cols(X),rows(X))';
plotelem(T,X+U,'g--')  
hold on

plotelem(T,Y-X,'b-')
hold on
plotelem(T,Y-X-U,'g--')  
hold on

U(:,1) = -1*U(:,1);
X(:,1) = -1*X(:,1);
plotelem(T,Y+X,'b-')
hold on
plotelem(T,Y+X+U,'g--')  
hold on 

plotelem(T,-X,'b-')
hold on
plotelem(T,-X-U,'g--')  
hold off


% Compute reactions
RH1 = reaction(q1,C,dof,1);
RV1 = reaction(q1,C,dof,2);

