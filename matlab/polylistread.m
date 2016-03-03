function [ polys, N ] = polylistread( path )
% generate polytope edges in crs and store data in folder path
% produces polys{i}.vertices
%          polys{i}.sys

    listpath=path;
    
    sys=binread([listpath '/sys.bin'],'int');
    N=binread([listpath '/N.bin'],'int');
    
    polys=cell(1);
    for i=1:N
       polys{i}=polyread( [listpath '/poly' num2str(i-1)] );
       polys{i}.sys=sys(i);
    end
    
    
end