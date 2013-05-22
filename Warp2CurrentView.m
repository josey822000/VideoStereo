function [D Model] = Warp2CurrentView(Key,P,curf,step)
%WARP2CURRENTVIEW Summary of this function goes here
%   Detailed explanation goes here
	sz = size(Key{1}.D);
	[X Y] = meshgrid(1:sz(2), 1:sz(1));
	pts = ones(sz(1)*sz(2), 3);
	pts(:,1) = X(:);
	pts(:,2) = Y(:);
	
	
	D= zeros(size(Key{1}.D,1),size(Key{1}.D,2),2);
	Model = cell(2,1);
    % use 2 keyframe, one warp to current view as reference, other one warp
    % to matching frame
    ka = [1 2 ; 2 1];
    for k=1:2
		% hole stay -1 label and assign big cost(no information from this)
        % Kf: keyframe1, K2: current
        Kf = P.K(:,:,ka(k,1)+1);
        Rf = P.R(:,:,ka(k,1)+1);
        Tf = P.T(:,ka(k,1)+1);
		Model{k}.frame = k*step;
		K2 = P.K(:,:,1);
		R2 = P.R(:,:,1);
		T2 = P.T(:,1);
		
		tmp_mat = K2 * R2' * Rf;
		tmp_vec = K2 * R2' * ( Tf - T2);
		tmp_vec = repmat(tmp_vec, [1 prod(sz)]);
		Model{k}.plane = ((Key{k}.Model.plane * Kf * Rf' * R2)/K2)./repmat(1 + Key{k}.Model.plane * Kf * Rf' * (T2-Tf),[1 3]);
        Model{k}.table = Key{k}.Model.table;
		epl_pts = ( tmp_mat * (Kf\(pts')) + repmat(Key{ka(k,1)}.D(:)',[3 1]) .* tmp_vec )';
		epl_pts = epl_pts ./ repmat(epl_pts(:,3),[1 3]);
		[Model{k}.F  Model{k}.segMap  D(:,:,k)] = GetWarpNoFillHole(Key{ka(k,1)}.Model.F,Key{ka(k,1)}.Model.segMap,Key{ka(k,1)}.D,epl_pts(:,1)-0.5,epl_pts(:,2)-0.5);
        % test
        Objdepth = zeros(sz);
        ObjIds = unique(Model{k}.segMap);
        if ObjIds(1) ==0
           ObjIds = ObjIds(2:end); 
        end
        for i= 1:numel(ObjIds)
            s = ObjIds(i);
            M = Model{k}.segMap==s;
            if nnz(M) ==0
                continue;
            end
            N = Model{k}.plane(Model{k}.table(s),:)';
            Dtmp = -(X(M) * N(1) + Y(M) * N(2) + N(3));
            Objdepth(M) = Dtmp;
        end
%         figure(8); imshow(Objdepth/0.0087);
        %test
% 		figure(2);imshow(D(:,:,k)/0.0087);
%         figure(3);sc(Model{k}.F,'rand');
%         figure(4);sc(Model{k}.segMap,'rand');
        % Kf: keyframe1, K2: keyframe2
		mpP.K = Kf;
		mpP.R = Rf;
		mpP.T = Tf;
		mpP.K(:,:,2) = P.K(:,:,ka(k,2)+1);
		mpP.R(:,:,2) = P.R(:,:,ka(k,2)+1);
		mpP.T(:,2) = P.T(:,ka(k,2)+1);
        % use k's D and warp k's info to k2 position
		Model{k}.otherView = ObjWarp(Key{k}.D,Key{k}.Model.F,Key{k}.Model.segMap,mpP);
% 		figure(5);imshow(Model{k}.otherView{1}.Disparity/0.0087);
%         figure(6);sc(Model{k}.otherView{1}.F,'rand');
%         figure(7);sc(Model{k}.otherView{1}.segMap,'rand');
		Model{k}.parallax = Key{k}.Model.parallax;
		Model{k}.GMM_Name = Key{k}.Model.GMM_Name;
        Model{k}.Key = Key{k}.Key;
    end
    % model1:  k1 warp to curf as reference ,K2 warp to K1 as matching
    % model1:  k2 warp to curf as reference ,K1 warp to K2 as matching
    tmp = Model{1}.otherView{1};
    other1.segMap = Key{1}.Model.segMap;
    other1.F = Key{1}.Model.F;
    other1.Disparity = Key{1}.D;
    Model{1}.otherView = {other1 tmp};
    
    tmp = Model{2}.otherView{1};
    other2.segMap = Key{2}.Model.segMap;
    other2.F = Key{2}.Model.F;
    other2.Disparity = Key{2}.D;
    Model{2}.otherView = {tmp other2};
    % here should fill the hole with neighbor object 
end

