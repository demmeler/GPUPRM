clear all

if 1
    A=imread('labyrinth2.jpg');
    A=A(:,:,1);

    [h,b]=size(A);
end

Ql=binread('prmoutput/graphl/qstoragec.bin','float');
Qr=binread('prmoutput/graphr/qstoragec.bin','float');
Ql=reshape(Ql,2,[]);
Qr=reshape(Qr,2,[]);
Ql=[h,0;
    0,h]*Ql;
Qr=[h,0;
    0,h]*Qr;



Elfrom=binread('prmoutput/graphl/edgesfromc.bin','int')'+1;
Elto=binread('prmoutput/graphl/edgestoc.bin','int')'+1;
Qlfrom=Ql(:,Elfrom);
Qlto=Ql(:,Elto);

Erfrom=binread('prmoutput/graphr/edgesfromc.bin','int')'+1;
Erto=binread('prmoutput/graphr/edgestoc.bin','int')'+1;
Qrfrom=Qr(:,Erfrom);
Qrto=Qr(:,Erto);


conl=binread('prmoutput/graphl/endc.bin','int')+1;
conr=binread('prmoutput/graphr/endc.bin','int')+1;


image(A);
hold on;

for i=1:length(Elfrom)
    X=[Qlfrom(1,i),Qlto(1,i)];
    Y=[Qlfrom(2,i),Qlto(2,i)];
    plot(X,Y,'w-');
end
for i=1:length(Erfrom)
    X=[Qrfrom(1,i),Qrto(1,i)];
    Y=[Qrfrom(2,i),Qrto(2,i)];
    plot(X,Y,'w-');
end

X=[Ql(1,conl),Qr(1,conr)];
Y=[Ql(2,conl),Qr(2,conr)];
plot(X,Y,'w--');

plot(Ql(1,:),Ql(2,:),'r.');
plot(Qr(1,:),Qr(2,:),'g.');

