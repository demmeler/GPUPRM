function [ res ] = binread( path, type )
% read a file

    bufsz=20000000;

    fd2= fopen(path,'r');
    res = fread(fd2, bufsz,type);
    fclose(fd2);

end

