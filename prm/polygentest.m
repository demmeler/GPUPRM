clear all;
close all;

wurfel =[0 0 0; 
        1 0 0; 
        1 1 0; 
        0 1 0; 
        0 0 1; 
        1 0 1; 
        1 1 1; 
        0 1 1
        ];
simplex=[0 0 0; 
        1 0 0;  
        0 1 0; 
        0 0 1;
        ];
test=   [0 0 0; 
        1 0 0;  
        0 1 0; 
        0 0 1;
        1 1 0
        ];
    
X=simplex;

[dsp cnt dest]=polygen(X)


K = convhulln(X);

if 1
patch('Faces',K,'Vertices',X,'FaceColor',rand(1,3),'FaceAlpha',1.0);
view(3);
axis equal;
end


E=[K(:,1),K(:,2);K(:,2),K(:,3);K(:,3),K(:,1)];
if 0
mask1=E(:,1)<E(:,2);
mask2=~mask1;
Efrom=mask1.*E(:,1)+mask2.*E(:,2);
Eto=mask2.*E(:,1)+mask1.*E(:,2);
E=[Efrom,Eto;Eto,Efrom];
%remove duplicates...
end
E=sortrows(E);


Efrom=E(:,1);
Eto=E(:,2);

ele=1:max(Efrom);

[A,B]=meshgrid(ele,Efrom);
cnt=sum(A==B,1);
dsp=cumsum(cnt)-cnt(1);
dest=Eto';







