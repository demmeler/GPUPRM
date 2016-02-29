function [] = dhwrite( path, a, alpha, q, d, ndof )
% write DH Parameters to path

    binwrite([path '/ndof.bin'], ndof, 'int');
    binwrite([path '/a.bin'], a, 'float');
    binwrite([path '/alpha.bin'], alpha, 'float');
    binwrite([path '/q.bin'], q, 'float');
    binwrite([path '/d.bin'], d, 'float');
    
end