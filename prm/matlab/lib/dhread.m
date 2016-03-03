function [ params ] = dhread( path )
% read DH Parameters from path
% writing:
% params.a
% params.alpha
% params.q
% params.d
% params.types: rotational = 0
%               prismatic  = 1

    params.a=binread([path '/a.bin'], 'float');
    params.alpha=binread([path '/alpha.bin'], 'float');
    params.q=binread([path '/q.bin'], 'float');
    params.d=binread([path '/d.bin'], 'float');
    params.types=binread([path '/types.bin'], 'int');
    
end