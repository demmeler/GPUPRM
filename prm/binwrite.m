function [ res ] = binwrite( path, data, type )
% write a file

    fd2= fopen(path,'w');
    res = fwrite(fd2, data,type);
    fclose(fd2);
    
end