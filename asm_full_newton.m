function u = asm_full_newton(A,f,g,maxit,N)
    %initial values
    it = 0;
    tol = 1e-15;
%     u = zeros(N*N,1);
    u = g;
    active_set = ones(N*N,1);
    
    while (it < maxit)
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
        
       Mat = A;  % taking here fresh A and f
       rhs = f;
       for i = 1:length(id)
           Mat(id(i),:) = 0;
           Mat(id(i),id(i)) = 1;
           rhs(id(i)) = g(id(i));
       end
       
       u = Mat\rhs;
    end
end