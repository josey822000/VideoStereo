function WarpedInfo = WarpAll( Key, KeyP)
%OBJWARP Summary of this function goes here
%   Detailed explanation goes here
    WarpedInfo = cell(size(Key,1),size(Key,1));
    sz = size(Key{1}.D);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    pts = ones(sz(1)*sz(2), 3);
    pts(:,1) = X(:);
    pts(:,2) = Y(:);
    for i=1:size(Key,1)
        Kf = KeyP.K(:,:,i);
        Rf = KeyP.R(:,:,i);
        Tf = KeyP.T(:,i);
        
        for j = 1:size(Key,1)
            if i==j
                WarpedInfo{i,j}.F = Key{i}.Model.F;
                WarpedInfo{i,j}.segMap = Key{i}.Model.segMap;
                WarpedInfo{i,j}.D = Key{i}.D;
                continue;
            end
            K2 = KeyP.K(:,:,j);
            R2 = KeyP.R(:,:,j);
            T2 = KeyP.T(:,j);
            tmp_mat = K2 * R2' * Rf;
            tmp_vec = K2 * R2' * ( Tf - T2);
            tmp_vec = repmat(tmp_vec, [1 prod(sz)]);

            epl_pts = ( tmp_mat * (Kf\(pts')) + repmat(Key{i}.D(:)',[3 1]) .* tmp_vec )';
            epl_pts = epl_pts ./ repmat(epl_pts(:,3),[1 3]);
            [WarpedInfo{i,j}.F WarpedInfo{i,j}.segMap WarpedInfo{i,j}.D] = GetWarpNoFillHole(Key{i}.Model.F,Key{i}.Model.segMap,Key{i}.D,epl_pts(:,1)-0.5,epl_pts(:,2)-0.5);
            
        end
    end
end

