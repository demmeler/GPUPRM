function [dhparams, ndof, polys, N, pairs, mins, maxs] = configread( configpath )
% write configuration to file structure


    listpath=[configpath '/polys'];
    dhpath=[configpath '/dh'];
    pairspath=[configpath '/pairs'];
    
    ndof=binread([configpath '/ndof.bin'],'int');
    mins=binread([configpath '/mins.bin'],'float');
    maxs=binread([configpath '/maxs.bin'],'float');
    [ polys, N]=polylistread(listpath);
    dhparams=dhread(dhpath);
    
    from=binread([pairspath '/from.bin'], 'int')+1;
    to=binread([pairspath '/to.bin'], 'int')+1;
    pairs=[reshape(from,[],1),reshape(to,[],1)];
    
end