addpath('lib');

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


% DH params

if 1
    ndof=4;
    dhparams.a=[0,0,0,2];
    dhparams.alpha=[0,-pi/2,-pi/2,pi/2];
    dhparams.q=[0,0,0,0];
    dhparams.d=[4,3,2.25,1];
    dhparams.types=[0,0,0,0];
    
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
    
    Pe0.vertices=quaderm(0,0,-0.25,12,12,0.5);
    Pe0.sys=0;
    Pe1.vertices=quader(0,3.7,0,1,1,3);
    Pe1.sys=0;
    Pe2.vertices=quader(-3,-3,0,1,1,2);
    Pe2.sys=0;
    Pe3.vertices=quader(-6,0.6,5,8,1,1);
    Pe3.sys=0;
    Pe4.vertices=movevertices( readstl('geodata/testobj2.stl') , [1,0,0]);
    Pe4.sys=0;
    
    polysbase={P1};
    Nbase=length(polysbase);
    polysself={P2,P3,P4,P5,P6};
    Nself=length(polysself);
    polysenv={Pe0,Pe1,Pe2,Pe3,Pe4};
    Nenv=length(polysenv);
    polys=[polysenv, polysbase, polysself];
    
    pairsenv=[repelem((1:Nenv)',Nself,1),repmat(Nenv+Nbase+(1:Nself)',Nenv,1)];
    pairsbase=Nenv+[repelem((1:Nbase)',Nself-1,1),repmat(Nbase+(2:Nself)',Nbase,1)];
    
    pairsself=Nenv+Nbase+[2,4;2,5;2,6;
                          3,5;3,6]-1;
                
    
    pairs=[pairsenv;pairsbase;pairsself];
    
    mins=[-3,-3,-3,-3];
    maxs=[4.5, 4.5, 4.5, 4.5];
    
end


N=length(polys);



%%% write files %%%

configpath='../config1';
configwrite(configpath,dhparams,ndof,polys,N,pairs, mins, maxs);

% check if evering stored correctly
if 0
    clear all;
    [dhparams, ndof, polys, N, pairs, mins, maxs] = configread( '../config1' );
end

pathplot(dhparams,ndof,polys,N,zeros(ndof,1)+0.0, 0.03, false);





