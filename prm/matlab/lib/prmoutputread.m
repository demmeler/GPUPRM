function [prmoutput] = prmoutputread( prmpath, ndof )
% write configuration to file structure
    
    Ql=binread([prmpath '/graphl/qstoragec.bin'],'float');
    Qr=binread([prmpath '/graphr/qstoragec.bin'],'float');
    prmoutput.Ql=reshape(Ql,ndof,[]);
    prmoutput.Qr=reshape(Qr,ndof,[]);

    
    Elfrom=binread([prmpath '/graphl/edgesfromc.bin'],'int')'+1;
    Elto=binread([prmpath '/graphl/edgestoc.bin'],'int')'+1;
    index=Elfrom<Elto;
    prmoutput.Elfrom=Elfrom(index);
    prmoutput.Elto=Elto(index);

    prmoutput.Qlfrom=prmoutput.Ql(:,prmoutput.Elfrom);
    prmoutput.Qlto=prmoutput.Ql(:,prmoutput.Elto);

    Erfrom=binread([prmpath '/graphr/edgesfromc.bin'],'int')'+1;
    Erto=binread([prmpath '/graphr/edgestoc.bin'],'int')'+1;
    index=Erfrom<Erto;
    prmoutput.Erfrom=Erfrom(index);
    prmoutput.Erto=Erto(index);
    
    prmoutput.Qrfrom=prmoutput.Qr(:,prmoutput.Erfrom);
    prmoutput.Qrto=prmoutput.Qr(:,prmoutput.Erto);

    
    prmoutput.conl=binread([prmpath '/graphl/endc.bin'],'int')+1;
    prmoutput.conr=binread([prmpath '/graphr/startc.bin'],'int')+1;

    prmoutput.connected=binread([prmpath '/connection.bin'],'int');

    prmoutput.pathl=binread([prmpath '/graphl/pathc.bin'],'int')'+1;
    prmoutput.pathr=binread([prmpath '/graphr/pathc.bin'],'int')'+1;

    prmoutput.Qpathl=prmoutput.Ql(:,prmoutput.pathl);
    prmoutput.Qpathr=prmoutput.Qr(:,prmoutput.pathr);
    
    prmoutput.Qpath=[prmoutput.Qpathr, prmoutput.Qpathl];
    
end