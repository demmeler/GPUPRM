function [ cnt, dsp, dest ] = polywrite( path, vertices )
% generate polytope edges in crs and store data in folder path
   
    n=size(vertices,1);
    [dsp,cnt,dest]=polygen(vertices);

    binwrite([path '/vertices.bin'],vertices,'float');
    binwrite([path '/dsp.bin'],dsp,'int');
    binwrite([path '/cnt.bin'],cnt,'int');
    binwrite([path '/dest.bin'],dest,'int');
    binwrite([path '/size.bin'],[n length(dest)],'int');
    
end