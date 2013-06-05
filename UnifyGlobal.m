function Key = UnifyGlobal( Key,P , ref, FThr, ObjThr)
%UNIFYUSINGTHR Summary of this function goes here
%   Detailed explanation goes here
    ObjOffset = zeros(size(Key,1)+1,1);
    FOffset = zeros(size(Key,1)+1,1);
    ONum = 0;
    FNum = 0;
    for k=1:size(Key,1)
        ObjOffset(k) = ONum;
        FOffset(k) = FNum;
        Key{k}.segMap = Key{k}.segMap + ONum;
        Key{k}.F = Key{k}.F + FNum;
        ONum = max(Key{k}.segMap(:));
        FNum = max(Key{k}.F(:));
    end
    ObjOffset(end,1) = ONum;
    FOffset(end,1) = FNum;
    disp(['Object Total Number:' num2str(ONum)]);
    disp(['Depth Total Number:' num2str(FNum)]);
    
    OTab = 1:ONum;
    FTab = 1:FNum;
    
    sz = size(Key{1}.D);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    pts = ones(sz(1)*sz(2), 3);
	pts(:,1) = X(:);
	pts(:,2) = Y(:);
    
    %% project plane to ref frame
    
    % object
    ObPln = zeros(0,3);
    ObjTable = 1:ONum;
    % depth plane
    plane = zeros(0,3);
    FTable = 1:FNum;
    DispFromObj = zeros(sz);
    Oid = unique(Key{ref}.segMap);
    % two neighbor combine
    
    % combine
    WarpInfo = cell(size(Key));
    % init ObjPln 
    ObPln = Key{1}.ObPln;
    plane = Key{1}.plane;
    for k=1:size(Key,1)-1
        % Kf: trans from keyframe, K2: trans to
        Kf = P.K(:,:,k+1);
        Rf = P.R(:,:,k+1);
        Tf = P.T(:,k);        
        K2 = P.K(:,:,k);
        R2 = P.R(:,:,k);
        T2 = P.T(:,k);
 
        tmp_mat = K2 * R2' * Rf;
		tmp_vec = K2 * R2' * ( Tf - T2);
		tmp_vec = repmat(tmp_vec, [1 prod(sz)]);
        
        ObjPln = [ObjPln ; ((Key{k+1}.ObPln * Kf * Rf' * R2)/K2)./repmat(1 + Key{k}.ObPln * Kf * Rf' * (T2-Tf),[1 3])];
        plane = [plane ; ((Key{k+1}.plane * Kf * Rf' * R2)/K2)./repmat(1 + Key{k}.plane * Kf * Rf' * (T2-Tf),[1 3])];
        
        % warp info to next frame
        %{
        for p = 1:sz(3)
            epl_pts = ( tmp_mat * (Kf\(pts')) + repmat(Key{k}.D(:,:,p)',[3 1]) .* tmp_vec )';
            epl_pts = epl_pts ./ repmat(epl_pts(:,3),[1 3]);
            [WarpInfo{k}.F WarpInfo{k}.segMap] = GetWarp(Key{k}.F(:,:,p),Key{k}.segMap(:,:,p),Key{k}.D(:,:,p),epl_pts(:,1)-0.5,epl_pts(:,2)-0.5);
            % recover D on ref frame using converted plane parameter
            DObj = zeros(sz);
            DF = zeros(sz);
            
        end
        %}
        % from depth plane
        K1_FNum = numel(unique(Key{k}.F));
        K2_FNum = numel(unique(Key{k+1}.F));
        Fid = unique(Key{k+1}.F);
        disp(['origin F:' num2str(K1_FNum)]);
        % calc difference between planes
        for i=1:numel(Fid)
            N = plane(  Fid(i),:);
            DF = -( X * N(1) + Y * N(2) + N(3));
            samePln = [];
            DiffOfKandRef = abs(repmat(DF,[1 1 sz(3)]) - Key{k}.D);
            % find not same set
            samePln = unique(Key{k}.F(DiffOfKandRef >= FThr));
            
            % find diff set
            samePln = setdiff(FOffset(k)+(1:K1_FNum),samePln);
            FTable(samePln) = Fid(i);
            for k2=1:k
                Key{k2}.F(ismember(Key{k2}.F,samePln)) = Fid(i);
            end
        end
        disp(['Unify F:' num2str(numel(unique(Key{k}.F)))]);
        
        % from object plane
        K1_ONum = numel(unique(Key{k}.segMap));
        K2_ONum = numel(unique(Key{k+1}.segMap));
        Oid = unique(Key{k+1}.segMap);
        
        disp(['origin Obj:' num2str(K1_ONum)]);
        % calc difference between planes
        for i=1:numel(Oid)
            N = ObjPln(  Oid(i),:);
            DF = -( X * N(1) + Y * N(2) + N(3));
            samePln = [];
            DiffOfKandRef = abs(repmat(DF,[1 1 sz(3)]) - Key{k}.D);
            % find not same set
            samePln = unique(Key{k}.segMap(DiffOfKandRef >= ObjThr));
            
            % find diff set
            samePln = setdiff(ObjOffset(k)+(1:K1_ONum),samePln);
            ObjTable(samePln) = Oid(i);
            for k2=1:k
                Key{k2}.segMap(ismember(Key{k2}.segMap,samePln)) = Oid(i);
            end
        end
        disp(['Unify Obj:' num2str(numel(unique(Key{k}.segMap)))]);
        
    end
    % according the sort, adjust plane sort
    for k=1:size(Key,1)
        [x y] = sort(FTable((FOffset(k)+1):FOffset(k+1)));
        Key{k}.plane = Key{k}.plane(y,:);
        [x y] = sort(ObjTable((ObjOffset(k)+1):ObjOffset(k+1)));
        Key{k}.ObjPln = Key{k}.ObjPln(y,:);
    end
end

function [D ObjModel] = RefitObj(D,ObjModel,R,imout,d_step,SEI,P)
    % Switch off annoying warnings
    warning_state = warning('query', 'all');
    warning off MATLAB:divideByZero
    warning off MATLAB:singularMatrix
    warning off MATLAB:nearlySingularMatrix
    warning off MATLAB:illConditionedMatrix
    warning off MATLAB:rankDeficientMatrix

    ObjIds = unique(ObjModel.segMap);
    if ObjIds(1) == 0
        ObjIds = ObjIds(2:end);
    end
    table = 1:max(ObjModel.segMap(:));
    table(ObjIds) = 1:numel(ObjIds);
    ObjModel.segMap(ObjModel.segMap>0) = table(ObjModel.segMap(ObjModel.segMap>0));
    
    segIds = unique(ObjModel.F);
    if segIds(1) == 0
        segIds = segIds(2:end);
    end
	table = 1:max(ObjModel.F(:));
    table(segIds) = 1:numel(segIds);
    ObjModel.F(ObjModel.F>0) = table(ObjModel.F(ObjModel.F>0));
    % every ObjSeg must have more than 5 pixels
    segMap = ObjModel.segMap;
    tmp = int32(segMap(SEI));
    DiffObjIdx = repmat(any(diff(int32(tmp)),1), [2 1]);
    Neigh = reshape(tmp(DiffObjIdx),2,[]);
    ObjIds = unique(segMap);
    if ObjIds(1) == 0
        ObjIds = ObjIds(2:end);
    end
    for s=1:numel(ObjIds)
       ObjId = ObjIds(s);
       if numel(find(segMap == ObjId)) < 5
           idx = find(Neigh(1,:)==ObjId & Neigh(2,:) ~= ObjId & Neigh(2,:) ~= 0);
           if ~isempty(idx)
               segMap(segMap == ObjId) = Neigh(2,idx(1));
               Neigh(Neigh==ObjId) = Neigh(2,idx(1));
           elseif ~isempty(find(Neigh(2,:)==ObjId & Neigh(1,:) ~= ObjId & Neigh(1,:) ~= 0))
               idx = find(Neigh(2,:)==ObjId & Neigh(1,:) ~= ObjId);
               segMap(segMap == ObjId) = Neigh(1,idx(1));
               Neigh(Neigh==ObjId) = Neigh(1,idx(1));
           else
               disp('filter too few pixels fail');
           end
       end
    end
    ObjIds = unique(segMap);
    if ObjIds(1) == 0
        ObjIds = ObjIds(2:end);
    end
    table = 1:max(segMap(:));
    table(ObjIds) = 1:numel(ObjIds);
    segMap(segMap>0) = table(segMap(segMap>0));
    % segMap
    
    ObjModel.segMap = segMap;
    [ObjModel.plane ObjModel.parallax ObjModel.GMM_Name] = calcObjModel(d_step,D,ObjModel,R,imout,2);
    
    ObjModel.otherView = ObjWarp(D,ObjModel.F,ObjModel.segMap,P);
    disp(numel(unique(ObjModel.segMap)));
    
    % Reset warnings
    warning(warning_state);
return
end