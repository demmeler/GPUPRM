function [] = dhwrite( path, params )
% write DH Parameters to path
% using:
% params.a
% params.alpha
% params.q
% params.d
% params.types: rotational = 0
%               prismatic  = 1

    mkdir(path);
    binwrite([path '/a.bin'], params.a, 'float');
    binwrite([path '/alpha.bin'], params.alpha, 'float');
    binwrite([path '/q.bin'], params.q, 'float');
    binwrite([path '/d.bin'], params.d, 'float');
    binwrite([path '/types.bin'], params.types, 'int');
    
end