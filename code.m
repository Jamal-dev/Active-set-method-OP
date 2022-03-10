clc
clear global;clear all;
%%
N = 65; % number of discretisation points in each direction
global ffunc active_set_nodes
active_set_nodes=[];
% rhs forcing
% ffunc = inline('-16*x.*(1-x).*y.*(1-y)','x','y'); % right hand side f
ffunc = @(x,y) -1;
maxit = 1000; % maximum number of iterations

% prepare system
[h,No_Els,total_dof,X,Y]=prepareSystem(N);
% Assemble finite element matrix
[K,M,F]=assemble(h,N,No_Els,total_dof);
global diagM;
diagM=diag(M);

% Insert Dirichlet datas on the whole boundary
[K,F]=dirchletConstraints(K,F,N);
%% Active set method
solution=activeSet(K,F,maxit);

%% solution display
% display solution
show(X,Y,solution,'Active set solution')


function val=obstacle_fun(x,y)
%     val=-0.008*(1+2*x+2*y);
    val = -.01; val=val*ones(size(x));
%     val=obs(x,y);
end
function [h,No_Els,total_dof,X,Y]=prepareSystem(N)

    global total_dof No_Els
    % Number of elements
    No_Els=(N-1)^2;
    disp(['Total Number of elements=' num2str(No_Els)]);
    % total degree of freedom
    total_dof=N^2;
    disp(['Total DOF=' num2str(total_dof)]);
    % element connectivity table
    global ele_table
    ele_table=elementConnectivity(No_Els,N);
    
    % step size
    h = 1/(N-1);
    global X Y X_vec Y_vec;
    [X,Y] = meshgrid(0:h:1); % discretisation of unit square
    X2=X'; Y2=Y';
    X_vec=X2(:); Y_vec=Y2(:);

end
function solution=activeSet(A,b,maxit)
    old_active_set=[];
    solution=A\b;
    global inactive_set
    for i=1:maxit
        fprintf('\n');
        disp(['Newton Iteration: ' num2str(i)]);
        [A_new,b_new,solution,active_set,residual]=updateSolution(solution,A,b);
        solution(inactive_set)=A_new(inactive_set,inactive_set)\b_new(inactive_set);
        disp(['Error: ',num2str(residual)]);
        % Condition if previous active set equal to current
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
function [A,b,solution,active_set,residual]=updateSolution(solution,A,b)
    global  diagM ele_table total_dof No_Els active_set_nodes inactive_set
    pen_par=100;
    lambda=A*solution-b;
    active_set=[];
    dof_touched=zeros(total_dof,1);
    for e=1:No_Els
        map = ele_table(e,:);
        for dof_index=map
            if dof_touched(dof_index)==0
               dof_touched(dof_index)=1;
            else
                % If dof is already touched then just continue
                continue;
            end
            [x_g,y_g]=node2pos(dof_index);
            obstacle_value=obstacle_fun(x_g,y_g);
            solution_value=solution(dof_index);
            active_cond=lambda(dof_index) + ...
            pen_par*diagM(dof_index) * ...
            (solution_value-obstacle_value)<0;
            if active_cond
                active_set=[active_set;dof_index];
                active_set_nodes=[active_set_nodes;dof_index];
                
%                 A(dof_index,:)=0;
%                 A(dof_index,dof_index)=1;
%                 b(dof_index)=obstacle_value;
                solution(dof_index) = obstacle_value;
                lambda(dof_index) = 0;
            end

        end
    end
    active_set_nodes=unique(active_set_nodes);
    inactive_set=setdiff(1:total_dof,active_set_nodes)';
    [A,b]=dirchlet_constraint(inactive_set,active_set,A,b);
    residual=norm(lambda(inactive_set));
end
function [A,b]=dirchlet_constraint(int_dofs,active_set,A,b)
    for k=1:length(active_set)
        dof=active_set(k);
        A(dof,:)=0;
        A(dof,dof)=1;
        [x,y]=node2pos(dof);
        b(dof)=obstacle_fun(x,y);
        b(int_dofs) = b(int_dofs) - A(int_dofs, dof)*obstacle_fun(x,y);
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



function [x,y]=node2pos(node)
    global X_vec Y_vec;
    
    x=X_vec(node); y=Y_vec(node);
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

function table=elementConnectivity(N_e,N)
    node_pos=reshape(1:N^2,N,N)';
    table=zeros(N_e,4);
    k=0;l=0;
    for e=1:N_e
    ele=zeros(2,2);
    for i=1:2
        for j=1:2
            ele(i,j)=node_pos(i+k,j+l);
        end
    end

    if rem(e,sqrt(N_e))~=0
        l=l+1;
    else
        k=k+1;
        l=0;
    end
    table(e,:)=[ele(1,1), ele(1,2),ele(2,2),ele(2,1)];
    end
end

function show(X,Y,u,tit)
    %% solution display
    % display solution
    N=length(X);
    for i=1:N
        U(i,1:N) = u((i-1)*N+1:i*N);
    end

    figure
    surf(X,Y,U);
    title(tit)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
end

function val=obs(x,y)
    val=zeros(size(x));
    r1=0.2;r2=0.1;c=0.5;
    h1=-0.02; h2=-0.01;
    for k=1:numel(x)
        pt=[x(k) y(k)];
        val(k)=cond(pt,c,r1,h1);
        val(k)=cond(pt,c,r2,h2,val(k));
        
    end
        
end

function val=cond(pt,c,r,v,elif)
    if nargin==4
        elif=true;
    else
        val=elif;
        elif=false;
    end
    x=pt(1); y=pt(2);
    a=c(1);
    if length(c)==1
        b=a;
    elseif length(c)==2
        b=c(2);
    end
    if ( ( x-a)^2  + ( y -b )^2)<=r^2
            val=v;
    else
        if elif
            val=-1;
        end
    end
end