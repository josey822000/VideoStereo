for i=15:105
    a = sprintf('Angkor/TXT/%d.txt',i-1)
    b = sprintf('Angkor/TXT/%03d.txt',i-1);
    movefile(a,b);
    
end