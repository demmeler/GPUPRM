clear all

if 0
    A=imread('labyrinth3.jpg');
    A=A(:,:,1);

    [h,b]=size(A);
end

Ql=binread('prmoutput/graphl/qstoragec.bin','float');
Qr=binread('prmoutput/graphr/qstoragec.bin','float');
Ql=reshape(Ql,2,[]);
Qr=reshape(Qr,2,[]);
if 0
Ql=[h,0;
    0,h]*Ql;
Qr=[h,0;
    0,h]*Qr;
end

Elfrom=binread('prmoutput/graphl/edgesfromc.bin','int')'+1;
Elto=binread('prmoutput/graphl/edgestoc.bin','int')'+1;
index=Elfrom<Elto;
Elfrom=Elfrom(index);
Elto=Elto(index);

Qlfrom=Ql(:,Elfrom);
Qlto=Ql(:,Elto);

Erfrom=binread('prmoutput/graphr/edgesfromc.bin','int')'+1;
Erto=binread('prmoutput/graphr/edgestoc.bin','int')'+1;
index=Erfrom<Erto;
Erfrom=Erfrom(index);
Erto=Erto(index);
Qrfrom=Qr(:,Erfrom);
Qrto=Qr(:,Erto);

conl=binread('prmoutput/graphl/endc.bin','int')+1;
conr=binread('prmoutput/graphr/startc.bin','int')+1;

con=binread('prmoutput/connection.bin','int');

pathl=binread('prmoutput/graphl/pathc.bin','int')'+1;
pathr=binread('prmoutput/graphr/pathc.bin','int')'+1;

Qpathl=Ql(:,pathl);
Qpathr=Qr(:,pathr);

close all;
pause(0.1);

%image(A);
hold on;

Gl=graph(Elfrom,Elto);
pl=plot(Gl,'-b','XData',Ql(1,:),'YData',Ql(2,:),'NodeLabel',{});
pl.NodeColor='y';

Gr=graph(Erfrom,Erto);
pr=plot(Gr,'-b','XData',Qr(1,:),'YData',Qr(2,:),'NodeLabel',{});
pr.NodeColor='g';

plot(Qpathl(1,:),Qpathl(2,:),'r-');
plot(Qpathr(1,:),Qpathr(2,:),'r-');

X=[Ql(1,conl),Qr(1,conr)];
Y=[Ql(2,conl),Qr(2,conr)];
plot(X,Y,'r-');

%plot(Ql(1,:),Ql(2,:),'r.');
%plot(Qr(1,:),Qr(2,:),'g.');
