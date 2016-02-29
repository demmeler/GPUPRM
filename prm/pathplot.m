function [] = pathplot( dhparams, ndof, polys, N, Q )
    
    for i=1:N
        polys{i}.K=convhulln(polys{i}.vertices);
    end
    
    for t=1:size(Q,2)
        clf;
        stateplot(dhparams,ndof,polys,N,Q(:,t));
        pause(0.05);
    end
    
end