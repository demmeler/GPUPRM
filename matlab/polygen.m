function [ dsp, cnt, dest ] = polygen( X )
% generate polytope edges in crs
   
K = convhulln(X);
E=[K(:,1),K(:,2);K(:,2),K(:,3);K(:,3),K(:,1)];
E=sortrows(E);
    
Efrom=E(:,1);
Eto=E(:,2);

ele=1:max(Efrom);

[A,B]=meshgrid(ele,Efrom);
cnt=sum(A==B,1);
dsp=[0,cumsum(cnt)];
dsp=dsp(1:length(cnt));
dest=Eto';

end