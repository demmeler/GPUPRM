
if 0
    A=imread('labyrinth.jpg');
    A=A(:,:,1);

    [h,b]=size(A);
end

if 1
    A=zeros(3,4);
    A(1,1)=1;
    A(1,2)=1;
    A(2,3)=1;
    
    h=3;
    b=4;
end

if 1
binwrite('array.bin',A','int');

B=reshape(binread('array.bin','int'),b,[])';
end

