clear all
close all

Node = [0 0 0; 
        1 0 0; 
        0 1 0; 
        1 1 0; 
        0 0 1; 
        1 0 1; 
        0 1 1; 
        1 1 1];
Elem = cell(1); Elem{1} = 1:size(Node,1); %Elem{2} = 6:10;
figure
hold on;
for elm = 1:size(Elem,2)    
    X = Node(Elem{elm},:); K = convhulln(X);    
    patch('Faces',K,'Vertices',X,'FaceColor',rand(1,3),'FaceAlpha',1.0);
end
view(3); 
%grid off;
axis equal; 
%cameramenu; 
%axis off;

