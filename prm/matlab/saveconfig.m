clear all;
%close all;

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

quader=@(x,y,z,lx,ly,lz)repmat([x,y,z],8,1)+wurfel*diag([lx,ly,lz]);
quaderm=@(mx,my,mz,lx,ly,lz) repmat([mx-lx/2,my-ly/2,mz-lz/2],8,1) ...
                             +wurfel*diag([lx,ly,lz]);


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

addpath('lib');

if 0
    ndof=2;
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
end
if 1
    ndof=2;
    dhparams.a=[0, 4.5];
    dhparams.alpha=[0, 0];
    dhparams.q=[0, 0];
    dhparams.d=[0, 0];
    dhparams.types=[0, 0];
    
    P1.vertices=quaderm(2,0,0,1,1,1);
    P1.sys=0;
    P2.vertices=quaderm(4.5,0,0,1,1,1);
    P2.sys=1;
    P3.vertices=quaderm(2,0,0,1,1,1);
    P3.sys=2;
    polys={P1,P2,P3};
    pairs=[1,2;1,3;2,3];
end
if 0
    ndof=2;
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
if 1
    ndof=3;
    dhparams.a=[0,0,0];
    dhparams.alpha=[0,-pi/2,-pi/2];
    dhparams.q=[0,0,0];
    dhparams.d=[4,3,2.25];
    dhparams.types=[0,0,0];
    
    P0.vertices=quaderm(0,0,-0.25,8,8,0.5);
    P0.sys=0;
    P1.vertices=quaderm(0,0,1.75,1,1,3.5);
    P1.sys=0;
    P2.vertices=quaderm(0,1,0,1,3.0,1);
    P2.sys=1;
    P3.vertices=quaderm(0,1.125-0.5,0,1,2.25,1);
    P3.sys=2;
    P4.vertices=quaderm(1,0,0,3,1,1);
    P4.sys=3;
    
    Pe1.vertices=quader(-3,2,0,1,1,2);
    Pe1.sys=0;
    
    polys={P0,P1,P2,P3,P4,Pe1};
    pairs=[1,3;1,4;1,5;2,4;2,5;3,5;6,4;6,5];
    
end



N=length(polys);



%%% write files %%%

configpath='../config1';
configwrite(configpath,dhparams,ndof,polys,N, pairs);


pathplot(dhparams,ndof,polys,N,zeros(ndof,1)+0.2, 0.03, false);





