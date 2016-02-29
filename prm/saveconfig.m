% DH params

ndof=2;

dhparams.a=[1.1, 0];
dhparams.alpha=[0, 0];
dhparams.q=[0, 0];
dhparams.d=[0, 0];
dhparams.types=[1, 0];

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
simplex=[0 0 0; 
        1 0 0;  
        0 1 0; 
        0 0 1;
        ];

P1.vertices=wurfel;
P1.sys=0;
P2.vertices=wurfel;
P2.sys=2;

polys={P1,P2,P2};
N=3;


%%% write files %%%

configpath='config1';
configwrite(configpath,dhparams,ndof,polys,N);