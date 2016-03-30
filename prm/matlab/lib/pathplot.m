function [] = pathplot( dhparams, ndof, polys, N, Q, dt, infinite)
    
    for i=1:N
        polys{i}.K=convhulln(polys{i}.vertices);
    end
    
    colors=repmat([1,1,0.5; 0.5,0.5,1],ndof+1,1);
    
    l=size(Q,2);
    while 1
    for t=[1:l,l:-1:1]
        lim=[xlim;ylim;zlim];
        
        clf;
        stateplot(dhparams,ndof,polys,N,Q(:,t),colors);
        
        limneu=[xlim;ylim;zlim];
        limneu=[min(limneu(:,1),lim(:,1)), max(limneu(:,2),lim(:,2))];
        xlim(limneu(1,:));
        ylim(limneu(2,:));
        zlim(limneu(3,:));
        
        pause(dt);
    end
        if ~infinite
            break;
        end
    end
    
    hold off
    
end