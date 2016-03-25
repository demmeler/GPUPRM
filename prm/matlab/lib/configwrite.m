function [] = configwrite( configpath, dhparams, ndof, polys, N, pairs, mins, maxs )
% write configuration to file structure
% pairs: [from,to]
    
    listpath=[configpath '/polys'];
    dhpath=[configpath '/dh'];
    pairspath=[configpath '/pairs'];
    
    mkdir(configpath)
    binwrite([configpath '/ndof.bin'],ndof,'int');
    binwrite([configpath '/mins.bin'],mins,'float');
    binwrite([configpath '/maxs.bin'],maxs,'float');
    
    polylistwrite(listpath, polys, N);
    dhwrite(dhpath, dhparams);
    
    mkdir(pairspath);
    from=pairs(:,1);
    to=pairs(:,2);
    M=length(from);
    binwrite([pairspath '/M.bin'], M, 'int');
    binwrite([pairspath '/from.bin'], from-1, 'int');
    binwrite([pairspath '/to.bin'], to-1, 'int');
    
end