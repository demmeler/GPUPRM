addpath('lib');

clear all;
close all;

[dhparams, ndof, polys, N]=configread('../config1');


prmoutput=prmoutputread('../prmoutput',ndof);

l=size(prmoutput.Qpath,2);

part=0:0.1:0.9;
seg=3*part.^2-2*part.^3;
t=repmat(seg, 1, l)+repelem(1:l,length(seg));

Qpathcont=interp1q((1:l)',prmoutput.Qpath',(1:0.2:l)')';
%Qpathcont=interp1q((1:l)',prmoutput.Qpath',t')';



pathplot(dhparams,ndof,polys,N,Qpathcont, 0.01, true);


