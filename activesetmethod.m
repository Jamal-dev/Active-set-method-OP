function u=activesetmethod(A,f,g,maxit,N)

% initial values
residual = zeros(N*N,1);
u = zeros(N*N,1);
it = 0;
active_set = ones(N*N,1); % active indeces
    
while (it <= maxit)
    it = it + 1;
    u_old = u;
    active_set_old = active_set;
        
    
    active_set = -residual+ g - u > 0;
    % stopping criterion
    if all(active_set_old == active_set) 
        break;
    end
        
    G = A;
    F = f; % this will be f as per eq 11 in paper
    % Find indices which are non-zero
    id = find(active_set);
    length(id)
    for i = 1:length(id)
        G(id(i),:) = 0;
        G(id(i),id(i)) = 1;
        F(id(i)) = g(id(i));
    end
    
    u = G\F;
    
%     u = u_old + u; % needed only for G\F if own algorithim is made
    residual = f-A*u;
%     disp(['residual=' num2str(norm(residual))])
end
if it==maxit
    disp('Max itreration reached. not converged')
end
end