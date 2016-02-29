function [dhparams, ndof, polys, N] = configread( configpath )
% write configuration to file structure
    
    listpath=[configpath '/polys'];
    dhpath=[configpath '/dh'];
    
    ndof=binread([configpath '/ndof.bin'],'int');
    [ polys, N]=polylistread(listpath);
    dhparams=dhread(dhpath);
    
end