clc
clear global
%%
N = 70; % number of discretisation points in each direction
global ffunc gfunc total_dof No_Els
ffunc = inline('-16*x.*(1-x).*y.*(1-y)','x','y'); % right hand side f
gfunc = inline('-0.008*(1+2*x+2*y)','x','y');     % obstacle function g
maxit = 1000; % maximum number of iterations

% Number of elements
No_Els=(N-1)^2;
disp(['Total Number of elements=' num2str(No_Els)]);
% total degree of freedom
total_dof=N^2;
disp(['Total DOF=' num2str(total_dof)]);
% element connectivity table
global ele_table
ele_table=elementConnectivity(No_Els,N);

h = 1/(N-1);
global X Y;
[X,Y] = meshgrid(0:h:1); % discretisation of unit square

% Assemble right hand side
ff = h*h*(ffunc(X,Y)); 
ft = ff';
f = ft(:); 

ff2 = ffunc(X,Y); 
ft2 = ff2';
f2 = ft2(:); 

% Assemble finite element matrix
A = gallery('poisson',N); % system matrix
[K,M,F]=assemble(h,N,No_Els,total_dof);
global diagM;
diagM=diag(M);


% Discrete obstacle values
gh = gfunc(X,Y);
gt = gh';
g = gt(:); 

% Insert Dirichlet datas on the whole boundary
[A,f]=dirchletConstraints(A,f,N);
[K,F]=dirchletConstraints(K,F,N);
%% Active set method
solution=activeSet(K,F,maxit);
% u=activeset(A,f,g,maxit,N);
% u2=activeset(K,f,g,maxit,N);
%% solution display
% display solution
show(X,Y,solution,'Active set New solution')
% show(X,Y,u,'Active set with Central Difference')
% show(X,Y,u2,'Active set with FEM')

function solution=activeSet(A,b,maxit)
    old_active_set=[];
    solution=A\b;
    for i=1:maxit
        disp(['Newton Iteration: ' num2str(i)]);
        [A_new,b_new,solution,contact_force,active_set,residual]=updateSolution(A,solution,b);
        disp(['Error: ',num2str(residual)]);
        if length(old_active_set)==length(active_set)
            if all(old_active_set==active_set)
                if isempty(active_set)
                    disp('Empty active set');
                else
                    disp('Solution is found!');
                end
                break;
            end
        end
        old_active_set=active_set;
    end
end
function [A,b,solution,contact_force,active_set,residual]=updateSolution(A,solution,b)
    global gfunc diagM ele_table total_dof No_Els
    pen_par=100;
    lambda=A*solution-b;
    residual=norm(lambda);
    contact_force=lambda;
    contact_force=-contact_force.*diagM;
    active_set=[];
    dof_touched=zeros(total_dof,1);
    for e=1:No_Els
        map = ele_table(e,:);
        for dof_index=map
            if dof_touched(dof_index)==0
               dof_touched(dof_index)=1;
            else
                continue;
            end
            [x_g,y_g]=node2pos(dof_index);
            obstacle_value=gfunc(x_g,y_g);
            solution_value=solution(dof_index);
            active_cond=lambda(dof_index) + ...
            pen_par*diagM(dof_index) * ...
            (solution_value-obstacle_value)<0;
            if active_cond
                active_set=[active_set;dof_index];
                solution(dof_index) = obstacle_value;
                A(dof_index,:)=0;A(dof_index,dof_index)=1;
                b(dof_index)=obstacle_value;
                lambda(dof_index) = 0;
            end

        end
    end
end
function N = BF(xi,eta)
    %This function returns the value of the four linear basis functions
    %at point xi,eta in reference coordinates
    N1 = 0.25*(1-xi).*(1-eta);
    N2 = 0.25*(1+xi).*(1-eta);
    N3 = 0.25*(1+xi).*(1+eta);
    N4 = 0.25*(1-xi).*(1+eta);
    
    N = [N1;N2;N3;N4];
end


function J=jacob(zeta,eta,X,Y)
    x1=X(1);x2=X(2);x3=X(3);x4=X(4); 
    y1=Y(1);y2=Y(2);y3=Y(3);y4=Y(4);
    J=zeros(2);
    J(1,1)= - 0.25*(1-eta)*x1 ...
          + 0.25*(1-eta)*x2 ...
        + 0.25*(1+eta)*x3 ...
        - 0.25*(1+eta)*x4;
    J(1,2)=-0.25*(1-eta)*y1 ...
        + 0.25*(1-eta)*y2 ...
        + 0.25*(1+eta)*y3 ...
        -0.25*(1+eta)*y4;
    J(2,1)=-0.25*(1-zeta)*x1- ...
        0.25*(1+zeta)*x2+ ...
        0.25*(1+zeta)*x3 ...
        +0.25*(1-zeta)*x4;
    J(2,2)=-0.25*(1-zeta)*y1- ...
        0.25*(1+zeta)*y2+ ...
        0.25*(1+zeta)*y3 ...
        +0.25*(1-zeta)*y4;
    
end

function [x,y]=node2pos(node)
    global X Y;
    X2=X'; Y2=Y';
    X3=X2(:); Y3=Y2(:);
    x=X3(node); y=Y3(node);
end

function dN = grad_BF(xi,eta)
    %This function returns the gradient of the four linear basis functions
    %at point xi,eta in reference coordinates
    dN_xi1 = -0.25*(1-eta);
    dN_xi2 =  0.25*(1-eta);
    dN_xi3 =  0.25*(1+eta);
    dN_xi4 = -0.25*(1+eta);

    dN_eta1 = -0.25*(1-xi);
    dN_eta2 = -0.25*(1+xi);
    dN_eta3 =  0.25*(1+xi);
    dN_eta4 =  0.25*(1-xi);

    dN = [dN_xi1,dN_eta1; dN_xi2,dN_eta2; dN_xi3,dN_eta3; dN_xi4,dN_eta4];
end
function fe = Elementalforce(nodes,h)
    % Coordinates and weights of Gaussian quadrature points in 2x2 reference frame
%     GPlocs   = [-5.773502691896258e-01, 5.773502691896258e-01];
    GPweight = [1,1,1,1];
    gp_loc=[-5.773502691896258e-01 5.773502691896258e-01 5.773502691896258e-01 -5.773502691896258e-01;
            -5.773502691896258e-01 -5.773502691896258e-01 5.773502691896258e-01 5.773502691896258e-01];
    gp_loc=gp_loc';    
    global ffunc
    fe=zeros(4,1);
    xiTox=zeros(4,1);etaToy=zeros(4,1);
    for n=1:length(nodes)
        [x,y]=node2pos(nodes(n));
        xiTox(n)=x;etaToy(n)=y;
    end
    J = h/2; %Jacobian of square element

    
    %Loop over quadrature points
    pt=0;
    k=1;
    for iGP=1:length(gp_loc) %Gauss points in x
        N  = BF(gp_loc(iGP,1),gp_loc(iGP,2));
        pt= pt+[xiTox(k)*N,etaToy(k)*N];
        k=k+1;
    end
    
    
    for iGP=1:length(gp_loc) %Gauss points in x
        N  = BF(gp_loc(iGP,1),gp_loc(iGP,2));
        % Integrate by summing the weighted contribution
        fe = fe +ffunc(pt(iGP,1),pt(iGP,2))* N*J * GPweight(iGP)*J * GPweight(iGP);
        % Integrate by summing the weighted contribution
        
    end
%     fe(3:4)=fe(4:-1:3);
end
function Me = ElementalMass(h)
    J=h/2;
    % Trapozoidal rule
    qp_loc=[-1 -1;
             1 -1;
             1  1;
            -1  1];
    QPweight=ones(4,1);
    Me=zeros(4);
    for iGP=1:length(qp_loc) %Gauss points in x
        N  = BF(qp_loc(iGP,1),qp_loc(iGP,2));
        % Integrate by summing the weighted contribution
        Me = Me + N*N'*J*J * QPweight(iGP)*QPweight(iGP);
    end
end
function Ke = ElementalStiffness(h)
    %Compute the element stiffness matrix of a single square element with
    %width h, and four linear nodal basis functions
    
    J = h/2; %Jacobian of square element

    % Coordinates and weights of Gaussian quadrature points in 2x2 reference frame
    GPlocs   = [-5.773502691896258e-01, 5.773502691896258e-01];
    GPweight = [1,1];
    
    %Set element stiffness matrix to zeros
    Ke = zeros(4);
    %Loop over quadrature points
    for iGP=1:length(GPlocs) %Gauss points in y
        for jGP=1:length(GPlocs) %Gauss points in x
            % Compute grad-N matrix in global coordinates
            B  = grad_BF(GPlocs(iGP),GPlocs(jGP)) / J;
            % Integrate by summing the weighted contribution
            Ke = Ke + B*B'*J*J * GPweight(iGP)*GPweight(jGP);
            N  = BF(GPlocs(iGP),GPlocs(jGP));
            
        end
    end
%     Ke=swapRows(Ke);
    
end
function K=swapRows(Ke)
    K=Ke;
    K(1:4,3:4)=Ke(1:4,4:-1:3);
    K(3:4,1:4)=K(4:-1:3,1:4);
    
end
function [K,M,F]=assemble(h,N,No_Els,total_DOF)
    % Elemental stiffness matrix
    Ke = ElementalStiffness(h);
    % Elemental mass matrix
    Me = ElementalMass(h);    

    %Initialize stiffness matrix
    K = zeros( total_DOF,total_DOF );
    M = zeros( total_DOF,total_DOF );
    F = zeros( total_DOF,1);
    ele_table=elementConnectivity(No_Els,N);

    %Assemble global stiffness matrix
    for e=1:No_Els
            %Map and add the element stiffness matrix to the global stiffness matrix
            map = ele_table(e,:);
            fe = Elementalforce(map,h);
            for i=1:4
                for j=1:4
                    K( map(i),map(j) ) = K( map(i),map(j) ) + Ke(i,j);
                    M( map(i),map(j) ) = M( map(i),map(j) ) + Me(i,j);
                end
                F(map(i))=F(map(i))+fe(i);
            end
    end
      
end

function [K,f]=dirchletConstraints(K,f,N)
    % Insert Dirichlet datas on the whole boundary
    for i=1:N
        % Top 
        K(i,:)   = 0; K(i,i)     = 1; f(i)   = 0;
        % Bottom
        ind = i+(N-1)*N;
        K(ind,:) = 0; K(ind,ind) = 1; f(ind) = 0;
        % Left
        ind = (i-1)*N+1;
        K(ind,:) = 0; K(ind,ind) = 1; f(ind) = 0;
        % Right
        ind = i*N;
        K(ind,:) = 0; K(ind,ind) = 1; f(ind) = 0;
    end
end

function u=activeset(A,f,g,maxit,N)
    % initial values
    uold = zeros(N*N,1);
    residual = zeros(N*N,1);
    u = uold;
    it = 0;
    %% Active set Method
    residual = zeros(N*N,1);
    active_workingSet = ones(N*N,1); % active indeces
    while (it < maxit)
        it = it+1;
        u0 = u;
        active_workingSet_0 = active_workingSet;
        % update active indeces
        active_workingSet = g-u > 0 ; % where u is crossing g (> since both are -ve)
        % stopping criterion
        if (active_workingSet == active_workingSet_0) & (it > 1)
            break;
        end
        Mat = A;
        rhs = f;
        % Find indices which are non-zero
        id = find(active_workingSet);
        length(id)
        for i = 1:length(id)
            Mat(id(i),:) = 0;
            Mat(id(i),id(i)) = 1;
            rhs(id(i)) = g(id(i)); 
        end
        % step 5
        u = Mat\rhs;
        residual = f-A*u;
        % assuming omega = 1
%         norm(u-u0)
        variation=residual;
        u= u + variation;
        
    end

end

function u=activeset2(A,f,g,maxit,N)
    % initial values
    uold = zeros(N*N,1);
    residual = zeros(N*N,1);
    u = uold;
    it = 0;
    %% Active set Method
    residual = zeros(N*N,1);
    active_workingSet = ones(N*N,1); % active indeces
    while (it < maxit)
        it = it+1;
        u0 = u;
        active_workingSet_0 = active_workingSet;
        % update active indeces
        active_workingSet = residual + g - u > 0; 
        % stopping criterion
        if (active_workingSet == active_workingSet_0) & (it > 1)
            break;
        end
        [A,f]=solveLinear(A,f)
        % step 5
        u = A\f;
        residual = A*u-f;
        % assuming omega = 1
        omega=1;
        u=u0+omega*u;
        norm(u-u0)
%         u=u+u0;
        
    end

end
function [Mat,rhs]=solveLinear(A,f)
    Mat = A;
    rhs = f;
    % Find indices which are non-zero
    id = find(active_workingSet);
    for i = 1:length(id)
        Mat(id(i),:) = 0;
        Mat(id(i),id(i)) = 1;
        rhs(id(i)) = g(id(i)); 
    end
end