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

polys={P1,P2};
N=2;


configpath='config1';
listpath=[configpath '/polys'];
dhpath=[configpath '/dh'];

polylistwrite(listpath, polys, N);
dhwrite(dhpath, a, alpha, q, d, ndof);