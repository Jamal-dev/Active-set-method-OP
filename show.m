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