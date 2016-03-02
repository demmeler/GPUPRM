clear all;
close all;

% Polytopes

wurfel =[0 0 0; 
        1 0 0; 
        1 1 0; 
        0 1 0; 
        0 0 1; 
        1 0 1; 
        1 1 1; 
        0 1 1
        ];
wurfeld=@(x,y,z)[wurfel(:,1)+x,wurfel(:,2)+y,wurfel(:,3)+z];

simplex=[0 0 0; 
        1 0 0;  
        0 1 0; 
        0 0 1;
        ];
body1 =[0 0 0; 
        1 0 0; 
        1 1 0; 
        0 2 0; 
        0 0 1; 
        1 0 1; 
        1 1 1; 
        0 1 1
        ];
body2 = [wurfel*0.5;
         ];


% DH params

ndof=2;
if 0
    dhparams.a=[1.1, 0];
    dhparams.alpha=[pi/2, 0];
    dhparams.q=[0, 0];
    dhparams.d=[0, -1.0];
    dhparams.types=[1, 0];
    
    P1.vertices=wurfel;
    P1.sys=0;
    P3.vertices=wurfel;
    P3.sys=2;
    polys={P1,P3};
    pairs=[1,2];
else
    dhparams.a=[0, 4.5];
    dhparams.alpha=[0, 0];
    dhparams.q=[0, 0];
    dhparams.d=[0, 0];
    dhparams.types=[0, 0];
    
    P1.vertices=wurfeld(2.0,0,0);
    P1.sys=0;
    P2.vertices=wurfeld(4.5,0,0);
    P2.sys=1;
    P3.vertices=wurfeld(1.7,0,0);
    P3.sys=2;
    polys={P1,P2,P3};
    pairs=[1,2;1,3;2,3];
end


N=length(polys);



%%% write files %%%

configpath='config1';
configwrite(configpath,dhparams,ndof,polys,N, pairs);


pathplot(dhparams,ndof,polys,N,[0;0], 0.03, false);





