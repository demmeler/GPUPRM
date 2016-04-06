function [ cnt, dsp, dest ] = polywrite( path, vertices )
% generate polytope edges in crs and store data in folder path
   
    [dsp,cnt,dest,verticesneu]=polygen2(vertices);
    n=size(verticesneu,1);
    
    mkdir(path);
    binwrite([path '/vertices.bin'],verticesneu','float'); %array of structs
    binwrite([path '/dsp.bin'],dsp,'int');
    binwrite([path '/cnt.bin'],cnt,'int');
    binwrite([path '/dest.bin'],dest-1,'int');
    binwrite([path '/size.bin'],[n length(dest)],'int');
    
end