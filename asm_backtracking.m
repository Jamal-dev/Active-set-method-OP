function u = asm_backtracking(A,f,g,maxit,N)
    %initial values
    it = 0;
    tol = 1e-15;
    u = zeros(N*N,1);
    u = g;
%     delta_u = zeros(N*N,1);
    active_set = ones(N*N,1);
    
    while (it < maxit)
       F = -A*u+f;
       it = it + 1;
       lambda = A*u - f;
       active_set0 = active_set;
       
       active_set = lambda + g - u > 0;
       inactive_set = lambda + g - u <= 0;
       
       id = find(active_set);
       id_i = find(inactive_set);
       
       residual = norm(lambda(id_i),2);
       
       if (active_set == active_set0) & (residual < tol)
            break;
       end
        
       Mat = A;
       rhs = F;
       for i = 1:length(id)
           Mat(id(i),:) = 0;
           Mat(id(i),id(i)) = 1;
           rhs(id(i)) = 0;
       end
       
       delta_u = Mat\rhs;
       u0 = u;
       alpha = 1;
       u = alpha*delta_u + u0;
       new_lambda = A*u - f;
       new_inactive_set = lambda + g - u <= 0;
       new_id_i = find(new_inactive_set);
       new_residual = norm(new_lambda(new_id_i),2);
       while (new_residual >= residual)
          alpha = 0.9*alpha;
          u = alpha*delta_u + u0;
          new_lambda = A*u - f;
          new_inactive_set = lambda + g - u <= 0;
          new_id_i = find(new_inactive_set);
          new_residual = norm(new_lambda(new_id_i),2);
       end
    end
end