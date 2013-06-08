function [ output_args ] = savetoFile( InputDir, OutputDir )
%SAVETOFIG Summary of this function goes here
%   Detailed explanation goes here
    mkdir(fullfile(InputDir,'TXT'));
    for i = 15:105
        tmp = load(fullfile(InputDir,[InputDir '_' num2str(i)]));
        fout = fopen(fullfile(InputDir,'TXT',sprintf('%03d.txt',num2str(i-1))),'w');
        D = tmp.finalD;
        D = D';
        for d = 1:numel(D)
            fprintf(fout,'%.9f ',D(d));
        end
        fclose(fout);
        
    end

end

