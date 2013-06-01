function D = ReadJiaya( filename, sz )
%READJIAYA Summary of this function goes here
%   Detailed explanation goes here
    File = fopen(filename,'r');
    if(~File)
        sprintf('cant open file');
    else
        D = fscanf(File,'%lg',prod(sz));
        D= reshape(D,sz(2),sz(1))';
        
    end

end

