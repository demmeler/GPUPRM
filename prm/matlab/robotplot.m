addpath('lib');

clear all;
close all;

[dhparams, ndof, polys, N, pairs, mins, maxs]=configread('../config1');


prmoutput=prmoutputread('../prmoutput',ndof);

l=size(prmoutput.Qpath,2);

part=0:0.05:0.9;
seg=3*part.^2-2*part.^3;
t=repmat(seg, 1, l)+repelem(1:l,length(seg));

Qpathcont=interp1q((1:l)',prmoutput.Qpath',(1:0.4:l)')';
%Qpathcont=interp1q((1:l)',prmoutput.Qpath',t')';
%Qpathcont=interp1q((1:l)',prmoutput.Qpath',[1;18;20;21;25;l])';


%pathplot(dhparams,ndof,polys,N,Qpathcont, 0.06, true);

pathplot(dhparams,ndof,polys,N, ...
    [prmoutput.Qpath(:,1),prmoutput.Qpath(:,round(l/4)),...
    prmoutput.Qpath(:,round(l*0.7)),prmoutput.Qpath(:,l)], ...
    0.06, false);



