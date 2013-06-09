function [ output_args ] = savetoFile( InputDir )
%SAVETOFIG Summary of this function goes here
%   Detailed explanation goes here
    mkdir(fullfile(InputDir,'TXT'));
    for i = 1:105
        tmp = load(fullfile(InputDir,[InputDir '_' num2str(i)]));
        
        fout = fopen(fullfile(InputDir,'TXT',sprintf('%03d.txt',i-1)),'w');
        D = tmp.finalD;
        D = D';
        
        fprintf(fout,'%.9f ',D(:));
        fclose(fout);
        
    end

end

