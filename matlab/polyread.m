function [ poly ] = polyread( path )
% generate polytope edges in crs and store data in folder path
   
    poly.vertices=binread([path '/vertices.bin'],'float'); %array of structs
    poly.vertices=reshape(poly.vertices,3,[])';
    [poly.dsp,poly.cnt,poly.dest]=polygen(poly.vertices);
    
end