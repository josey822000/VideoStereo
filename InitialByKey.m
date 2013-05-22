function [ output_args ] = InitialByKey( vals, options)
%INITIALALL Summary of this function goes here
%   Detailed explanation goes here
    P = vals.P;
    KeyDir = 'keyframe';
    options.KeyDir = KeyDir;
	step = 15;
	DirName = dir(KeyDir);
	Key = cell(numel(DirName)-2,1);
	keyFrameNum = numel(DirName)-2;
    % rec key frame 
    ObjNum = 0;
    FNum = 0;
	for i= 3:numel(DirName)
		keypath = DirName(i).name;
        fprintf('%s',keypath);
		tmp = load(fullfile(KeyDir,keypath,'Final_D'));
		Key{i-2}.D = tmp.D;
		tmp = load(fullfile(KeyDir,keypath,'Final_Model'));
        Key{i-2}.Key = 'KeyGMM';
		Key{i-2}.Model = tmp.ObjModel;
        oldObjId = 1:max(tmp.ObjModel.segMap(:));
        oldObjId(unique(tmp.ObjModel.segMap)) = 1:numel(unique(tmp.ObjModel.segMap));
        Key{i-2}.Model.segMap = oldObjId(Key{i-2}.Model.segMap);
        Key{i-2}.Model.segMap = Key{i-2}.Model.segMap + ObjNum;
        oldObjId = 1:max(tmp.ObjModel.F(:));
        oldObjId(unique(tmp.ObjModel.F)) = 1:numel(unique(tmp.ObjModel.F));
        Key{i-2}.Model.F = oldObjId(Key{i-2}.Model.F);
        Key{i-2}.Model.F = Key{i-2}.Model.F + FNum;
        ObjNum = max(Key{i-2}.Model.segMap(:));
        FNum = max(Key{i-2}.Model.F(:));
	end
	% load P
    useKey = [1 15 30 45 60 75 90 105];
    for i=1:keyFrameNum
        table = zeros(ObjNum,1);
        table(unique(Key{i}.Model.segMap)) = 1:numel(unique(Key{i}.Model.segMap));
        Key{i}.Model.table = table;
    end
    
	for i=1:size(P.K,3)
% 		if ~isempty(find(useKey == i))
% 			continue;
% 		end
		% warp to now frame
		
		% hole stay -1 label and assign big cost(no information from this)
        tmpvals = vals;
        tmpvals.I = vals.I(useKey);
        tmpvals.R = repmat(reshape(single(vals.I{i}), [], 3), [2 1]);
        tmpvals.P.K = P.K(:,:,[i useKey]);
        tmpvals.P.R = P.R(:,:,[i useKey]);
        tmpvals.P.T = P.T(:,[i useKey]);
        tmpvals.useKey = useKey;
		sz = size(Key{1}.D);
        [X Y] = meshgrid(1:sz(2), 1:sz(1));
        pts = ones(sz(1)*sz(2), 3);
        pts(:,1) = X(:);
        pts(:,2) = Y(:);


        D= zeros(size(Key{1}.D,1),size(Key{1}.D,2),keyFrameNum);
        Model = cell(2,1);
        % use 2 keyframe, one warp to current view as reference, other one warp
        % to matching frame
        for k=1:size(tmpvals.P.K,3)-1
            % hole stay -1 label and assign big cost(no information from this)
            % Kf: keyframe1, K2: current
            Kf = tmpvals.P.K(:,:,k+1);
            Rf = tmpvals.P.R(:,:,k+1);
            Tf = tmpvals.P.T(:,k+1);
            Model{k}.frame = k*step;
            K2 = tmpvals.P.K(:,:,1);
            R2 = tmpvals.P.R(:,:,1);
            T2 = tmpvals.P.T(:,1);

            tmp_mat = K2 * R2' * Rf;
            tmp_vec = K2 * R2' * ( Tf - T2);
            tmp_vec = repmat(tmp_vec, [1 prod(sz)]);
            Model{k}.plane = ((Key{k}.Model.plane * Kf * Rf' * R2)/K2)./repmat(1 + Key{k}.Model.plane * Kf * Rf' * (T2-Tf),[1 3]);
            Model{k}.table = Key{k}.Model.table;
            epl_pts = ( tmp_mat * (Kf\(pts')) + repmat(Key{k}.D(:)',[3 1]) .* tmp_vec )';
            epl_pts = epl_pts ./ repmat(epl_pts(:,3),[1 3]);
            [Model{k}.F  Model{k}.segMap  D(:,:,k)] = GetWarpNoFillHole(Key{k}.Model.F,Key{k}.Model.segMap,Key{k}.D,epl_pts(:,1)-0.5,epl_pts(:,2)-0.5);
            figure(k+1); imshow(D(:,:,k)/0.0087);
        end
        finalD = D(:,:,1);
        Holemap = double(Model{1}.segMap ~= 0);
        for k = 2:keyFrameNum
            kD = D(:,:,k);
            newHoleMap = (Model{k}.segMap ~= 0);
            finalD(newHoleMap) = finalD(newHoleMap) + kD(newHoleMap);
            Holemap = Holemap + double(newHoleMap);
        end
        finalD = (finalD + 1.e-10) ./ Holemap;
        figure(10); imshow(finalD/0.0087);
        delete([options.sequence '_' num2str(i) '.mat']);
        save([options.sequence '_' num2str(i) '.mat'],'finalD');	
	end
end

