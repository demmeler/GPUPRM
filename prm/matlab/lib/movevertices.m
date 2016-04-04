function V=movevertices(Vin,dir)
    
    if length(dir)~=3
        disp('Error in movecertices: length(dir)!=3');
        return;
    end
    
    if size(dir,1)==3
        dir=dir';
    end
    
    V=Vin+repmat(dir,size(Vin,1),1);
    
end