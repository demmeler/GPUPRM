function [] = polylistwrite( path, polys, N )
% generate polytope edges in crs and store data in folder path
% uses polys{i}.vertices
%      polys{i}.sys

    listpath=path;%[path '/polylist'];
    
    sys=zeros(1,N);
    binwrite([listpath '/N.bin'],N,'int');
    for i=1:N
       sys(i)=polys{i}.sys;
       polywrite( [listpath '/poly' num2str(i-1)] ,polys{i}.vertices );
    end
    
    binwrite([listpath '/sys.bin'],sys,'int');
    
end