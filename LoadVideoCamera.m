function KRT = LoadVideoCamera( fileName )
%LOADVIDEOCAMERA Summary of this function goes here
%   Detailed explanation goes here
    File = fopen(fileName,'r');
    if(~File)
        disp('cant open file');
    else
        disp('open file');
        fgets(File);
        fgets(File);
        fgets(File);
        fgets(File);
        frame_num = fscanf(File,'%d',1);
        KRT = cell(frame_num,1);
        K = zeros(3,3);
        R = zeros(3,3);
        T = zeros(3,1);
        for f=1:frame_num
            K = reshape(fscanf(File,'%lg',9),size(K))';
            R = reshape(fscanf(File,'%lg',9),size(R))';
            T = fscanf(File,'%lg',3);
            KRT{f}.K = K;
            KRT{f}.R = R;
            KRT{f}.T = T;
        end
        save([fileName(1:end-4) '_KRT'],'KRT');
        fclose(File);
    end
    
    

end

