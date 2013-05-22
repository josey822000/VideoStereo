function D = ReadJiaya( filename )
%READJIAYA Summary of this function goes here
%   Detailed explanation goes here
    File = fopen(filename,'r');
    if(~File)
        sprintf('cant open file');
    else
        sz = [576 352];
        D = fscanf(File,'%lg',prod(sz));
        D= reshape(D,sz)';
        
    end

end

