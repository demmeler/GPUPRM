addpath('lib');
addpath('geodata');

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
if 0
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
    
    Pe1.vertices=quader(-3,2,0,1,1,3);
    Pe1.sys=0;
    
    Pe2.vertices=quader(-3,-3,0,1,1,2);
    Pe2.sys=0;
    
    Pe3.vertices=quader(-6,0,6,8,1,1);
    Pe3.sys=0;
    
    polys={P0,P1,P2,P3,P4,Pe1,Pe2,Pe3};
    pairs=[1,3;1,4;1,5;2,4;2,5;3,5;6,4;6,5;7,4;7,5;8,3;8,4;8,5];
    
end
if 1
    ndof=4;
    dhparams.a=[0,0,0,2];
    dhparams.alpha=[0,-pi/2,-pi/2,pi/2];
    dhparams.q=[0,0,0,0];
    dhparams.d=[4,3,2.25,1];
    dhparams.types=[0,0,0,0];
    
    P0.vertices=quaderm(0,0,-0.25,12,12,0.5);
    P0.sys=0;
    P1.vertices=quaderm(0,0,1.75,1,1,3.5);
    P1.sys=0;
    P2.vertices=quaderm(0,1,0,1,3.0,1);
    P2.sys=1;
    P3.vertices=quaderm(0,1.125-0.5,0,1,2.25,1);
    P3.sys=2;
    P4.vertices=quaderm(1,0,0,3,1,1);
    P4.sys=3;
    P5.vertices=quaderm(1,0,0.3,3,1,1);
    P5.sys=4;
    P6.vertices=quaderm(0,0,-0.5+0.15,1,1,0.3);
    P6.sys=4;
    
    Pe1.vertices=quader(0,3.7,0,1,1,3);
    Pe1.sys=0;
    
    Pe2.vertices=quader(-3,-3,0,1,1,2);
    Pe2.sys=0;
    
    Pe3.vertices=quader(-6,0.6,5,8,1,1);
    Pe3.sys=0;
    
    polys={P0,P1,P2,P3,P4,Pe1,Pe2,Pe3,P5,P6};
    pairs=[1,3;1,4;1,5;
           2,4;2,5;
           3,5;
           6,4;6,5;
           7,4;7,5;
           8,2;8,3;8,4;8,5;
           9,1;9,2;9,3;9,4;9,6;9,7;9,8;
           10,1;10,2;10,3;10,4;10,6;10,7;10,8
           ];
       
     %pairs=[8,2];
     
     mins=[-3,-3,-3,-3];
     maxs=[4.5, 4.5, 4.5, 4.5];
     
    
end
if 0
    ndof=4;
    dhparams.a=[0,0,0,2];
    dhparams.alpha=[0,-pi/2,-pi/2,pi/2];
    dhparams.q=[0,0,0,0];
    dhparams.d=[4,3,2.25,1];
    dhparams.types=[0,0,0,0];
    
    P0.vertices=quaderm(0,0,-0.25,12,12,0.5);
    P0.sys=0;
    P1.vertices=quaderm(0,0,1.75,1,1,3.5);
    P1.sys=0;
    P2.vertices=quaderm(0,1,0,1,3.0,1);
    P2.sys=1;
    P3.vertices=quaderm(0,1.125-0.5,0,1,2.25,1);
    P3.sys=2;
    P4.vertices=quaderm(1,0,0,3,1,1);
    P4.sys=3;
    P5.vertices=quaderm(1,0,0.3,3,1,1);
    P5.sys=4;
    P6.vertices=quaderm(0,0,-0.5+0.15,1,1,0.3);
    P6.sys=4;
    
    Pe1.vertices=quader(0,3.7,0,1,1,3);
    Pe1.sys=0;
    
    Pe2.vertices=quader(-3,-3,0,1,1,2);
    Pe2.sys=0;
    
    Pe3.vertices=quader(-6,0.6,5,8,1,1);
    Pe3.sys=0;
    
    
    Pe4.vertices=movevertices( readstl('geodata/testobj.stl') , [1,0,0]);
    Pe4.sys=0;
    
    
    polys={P0,P1,P2,P3,P4,Pe1,Pe2,Pe3,P5,P6,Pe4};
    pairs=[1,3;1,4;1,5;
           2,4;2,5;
           3,5;
           6,4;6,5;
           7,4;7,5;
           8,2;8,3;8,4;8,5;
           9,1;9,2;9,3;9,4;9,6;9,7;9,8;
           10,1;10,2;10,3;10,4;10,6;10,7;10,8
           ];
       
     %pairs=[8,2];
     
     mins=[-3,-3,-3,-3];
     maxs=[4.5, 4.5, 4.5, 4.5];
     
    
end


N=length(polys);



%%% write files %%%

configpath='../config1';
configwrite(configpath,dhparams,ndof,polys,N, pairs, mins, maxs);


pathplot(dhparams,ndof,polys,N,zeros(ndof,1)+0.0, 0.03, false);





