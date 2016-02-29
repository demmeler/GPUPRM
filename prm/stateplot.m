function [] = stateplot( dhparams, ndof, polys, N, q )
% Q= [q1, ...., ql] describing a path

    hold on
    
    T=calctrafos(dhparams,ndof,q);
    
    
    for i=1:N
        V=[polys{i}.vertices,ones(size(polys{i}.vertices,1),1)];
        X=V*T{polys{i}.sys+1}';
        X=X(:,1:3);
        patch('Faces',polys{i}.K,'Vertices',X,'FaceColor',[1,1,0.5],'FaceAlpha',1.0);
    end
    
    axis equal
    view(3);
    
    hold off

end