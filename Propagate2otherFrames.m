function otherView = Propagate2otherFrames( D,ObjModel, images, P)
%PROPAGATE2OTHERFRAMES Summary of this function goes here
%   Detailed explanation goes here
    sz = size(D);
    Kf = P.K(:,:,1);
    Rf = P.R(:,:,1);
    Tf = P.T(:,1);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    pts = ones(sz(1)*sz(2), 3);
    pts(:,1) = X(:);
    pts(:,2) = Y(:);
    P = vals.P;
    Kf = P.K(:,:,1);
    Rf = P.R(:,:,1);
    Tf = P.T(:,1);

    segment_params = [1 1.5 10 100];
    info.segments = zeros(sz(1), sz(2), nMaps, 'uint32');
    sp = segment_params;
    otherView = cell(size(P.K,3)-1,1);
    for a = 1:size(P.K,3)-1
        % Segment the image using mean shift
        otherView{a}.segments = vgg_segment_ms(images{a}, sp(1), sp(2), sp(3));
    end
    % For each image...
    
    for a = 1:size(P.K,3)-1
        % Project the points
        K2 = P.K(:,:,a+1);
        R2 = P.R(:,:,a+1);
        T2 = P.T(:,a+1);
        tmp_mat = K2 * R2' * Rf;
        tmp_vec = K2 * R2' * ( Tf - T2);
        tmp_vec = repmat(tmp_vec, [1 prod(sz)]);
        epl_pts = ( tmp_mat * (Kf\(pts(:,:)')) + repmat(D(:)',[3 1]) .* tmp_vec )';
        [otherView{a}.F otherView{a}.segMap otherView{a}.Disparity] = GetWarpNoFillHole(F,segMap,D,epl_pts(:,1),epl_pts(:,2));
        Hole = otherView{a}.segMap == 0;
        HoleSeg = unique(otherView{a}.segments(Hole));
        img = images{a+1};
        % predict the hole samples in seg
        for s = 1:max(HoleSeg)
            seg = otherView{a}.segments == HoleSeg(HoleSeg(s));
            % get ObjId on segment
            ObjInSeg = unique(otherView{a}.segMap(seg));
            % get hole position
            seg(seg>0) = otherView{a}.segments(seg) == 0;
            
            samples = img(seg,:);
            probs = zeros(numel(ObjInSeg),1);
            for o =1:max(ObjInSeg)
                prob = PredictGMM(samples',ObjModel.GMM_Name{ObjInSeg(o)});
                probs(o) = mean(prob);
            end
            [v id] = max(probs);
            otherView{a}.segMap(seg) = id;
        end
    end
    
end

