% DH params

ndof=2;

dhparams.a=[1.1, 0];
dhparams.alpha=[pi/2, 0];
dhparams.q=[0, 0];
dhparams.d=[0, -1.0];
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
body1 =[0 0 0; 
        1 0 0; 
        1 1 0; 
        0 2 0; 
        0 0 1; 
        1 0 1; 
        1 1 1; 
        0 1 1
        ];


P1.vertices=wurfel;
P1.sys=0;
P3.vertices=wurfel;
P3.sys=2;

polys={P1,P3};
N=length(polys);


%%% write files %%%

configpath='config1';
configwrite(configpath,dhparams,ndof,polys,N);


pathplot(dhparams,ndof,polys,N,zeros(ndof,1), 0.03, false);