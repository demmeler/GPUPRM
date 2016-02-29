
Node = [0 0 0; 
        1 0 0; 
        1 1 0; 
        0 1 0; 
        0 0 1; 
        1 0 1; 
        1 1 1; 
        0 1 1;
        2 0 0;
        0 -1 0];
Elem = cell(1); Elem{1} = 1:10;
figure
for elm = 1:size(Elem,1)    
    X = Node(Elem{elm},:); K = convhulln(X); hold on;    
    patch('Faces',K,'Vertices',X,'FaceColor',rand(1,3),'FaceAlpha',1.0);
end
view(3); 
%grid off;
axis equal; 
cameramenu; 
%axis off;