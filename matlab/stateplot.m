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
    
    x=zeros(1,ndof);
    y=zeros(1,ndof);
    z=zeros(1,ndof);
    xr=zeros(1,ndof);
    yr=zeros(1,ndof);
    zr=zeros(1,ndof);
    
    
    for i=1:(ndof+1)
        Ti=T{i};
        x(i)=Ti(1,4);
        y(i)=Ti(2,4);
        z(i)=Ti(3,4);
        xr(i)=Ti(1,3);
        yr(i)=Ti(2,3);
        zr(i)=Ti(3,3);
    end
    
    
    axis equal
    view(3);
    
    lim=[xlim;ylim;zlim];
    plot3(x,y,z,'x');
    
    labels = cellstr( num2str((0:ndof)') );
    text(x,y,z,labels, 'HorizontalAlignment', 'right');
    for i=1:(ndof+1)
        l=1000;
        plot3(x(:,i)+[1,-1]*l*xr(:,i),y(:,i)+[1,-1]*l*yr(:,i),z(:,i)+[1,-1]*l*zr(:,i),'-.');
    end
    xlim(lim(1,:));
    ylim(lim(2,:));
    zlim(lim(3,:));
    
    
    hold off

end