function [] = configwrite( configpath, dhparams, ndof, polys, N )
% write configuration to file structure
    
    listpath=[configpath '/polys'];
    dhpath=[configpath '/dh'];
    
    mkdir(configpath)
    binwrite([configpath '/ndof.bin'],ndof,'int');
    polylistwrite(listpath, polys, N);
    dhwrite(dhpath, dhparams);  
    
end