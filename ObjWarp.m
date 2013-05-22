function otherView = ObjWarp( Dproposal,F,segMap, P)
%OBJWARP Summary of this function goes here
%   Detailed explanation goes here
    Kf = P.K(:,:,1);
    Rf = P.R(:,:,1);
    Tf = P.T(:,1);
    sz = size(Dproposal);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    pts = ones(sz(1)*sz(2), 3);
    pts(:,1) = X(:);
    pts(:,2) = Y(:);
    
    % For each image...
    otherView = cell(size(P.K,3)-1,1);
    
    for a = 2:size(P.K,3)
        K2 = P.K(:,:,a);
        R2 = P.R(:,:,a);
        T2 = P.T(:,a);
        tmp_mat = K2 * R2' * Rf;
        tmp_vec = K2 * R2' * ( Tf - T2);
        tmp_vec = repmat(tmp_vec, [1 prod(sz)]);
        
        epl_pts = ( tmp_mat * (Kf\(pts')) + repmat(Dproposal(:)',[3 1]) .* tmp_vec )';
        epl_pts = epl_pts ./ repmat(epl_pts(:,3),[1 3]);
        [otherView{a-1}.F otherView{a-1}.segMap otherView{a-1}.Disparity] = GetWarp(F,segMap,Dproposal,epl_pts(:,1)-0.5,epl_pts(:,2)-0.5);
    end
end

