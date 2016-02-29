clear all;
close all;

[dhparams, ndof, polys, N]=configread('config1');


prmoutput=prmoutputread('prmoutput',ndof);

l=size(prmoutput.Qpath,2);

Qpathcont=interp1q((1:l)',prmoutput.Qpath',(1:0.2:l)')';


pathplot(dhparams,ndof,polys,N,Qpathcont, 0.01, true);


