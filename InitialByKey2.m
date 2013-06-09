function [ output_args ] = InitialByKey2( Key, vals, options)
%INITIALALL Summary of this function goes here
%   Detailed explanation goes here
    P = vals.P;
    KeyDir = 'key/';
   
	
	keyFrameNum = numel(options.KeyFrame);
   
    % rec key frame 
    ObjNum = 0;
    FNum = 0;

    sz = size(Key{1}.D);
	for i=1:size(P.K,3)
% 		if ~isempty(find(useKey == i))
% 			continue;
% 		end
		% warp to now frame
		
		% hole stay -1 label and assign big cost(no information from this)
        
		
        [X Y] = meshgrid(0:sz(2)-1, 0:sz(1)-1);
        pts = ones(sz(1)*sz(2), 3);
        pts(:,1) = X(:);
        pts(:,2) = Y(:);


        D= zeros(size(Key{1}.D,1),size(Key{1}.D,2),keyFrameNum);
        Model = cell(2,1);
        % use 2 keyframe, one warp to current view as reference, other one warp
        % to matching frame
        for k=1:keyFrameNum
            % hole stay -1 label and assign big cost(no information from this)
            % Kf: keyframe1, K2: current
            Kf = vals.P.K(:,:,options.KeyFrame(k));
            Rf = vals.P.R(:,:,options.KeyFrame(k));
            Tf = vals.P.T(:,options.KeyFrame(k));
            K2 = vals.P.K(:,:,i);
            R2 = vals.P.R(:,:,i);
            T2 = vals.P.T(:,i);

            tmp_mat = K2 * R2' * Rf;
            tmp_vec = K2 * R2' * ( Tf - T2);
            tmp_vec = repmat(tmp_vec, [1 prod(sz)]);
            epl_pts = ( tmp_mat * (Kf\(pts')) + repmat(Key{k}.D(:)',[3 1]) .* tmp_vec )';
            epl_pts = epl_pts ./ repmat(epl_pts(:,3),[1 3]);
            epl_pts(:,1:2) = epl_pts(:,1:2) + 1;
            [Model{k}.F  Model{k}.segMap  D(:,:,k)] = GetWarpNoFillHole(Key{k}.F,Key{k}.segMap,Key{k}.D,epl_pts(:,1),epl_pts(:,2));
            figure(2*k); imshow(Key{k}.D/vals.d_step);
            figure(2*k+1); imshow(D(:,:,k)/vals.d_step);
        end
        finalD = D(:,:,2);
        Holemap = double(Model{2}.segMap ~= 0);
        for k = 2:keyFrameNum
            kD = D(:,:,k);
            newHoleMap = (Model{k}.segMap ~= 0);
            finalD(newHoleMap) = finalD(newHoleMap) + kD(newHoleMap);
            Holemap = Holemap + double(newHoleMap);
        end
        finalD = (finalD + 1.e-10) ./ Holemap;
        finalD(finalD == Inf) = 0;
        figure(10); imshow(finalD/vals.d_step);
        delete([options.sequence '_' num2str(i) '.mat']);
        save([options.sequence '_' num2str(i) '.mat'],'finalD');	
	end
end

