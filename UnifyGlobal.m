function Key = UnifyGlobal( Key,P , ref, threshold)
%UNIFYUSINGTHR Summary of this function goes here
%   Detailed explanation goes here
    ONum = 0;
    FNum = 0;
    for k=1:size(Key,1)
        ONum = ONum + max(Key{k}.segMap(:));
        FNum = FNum + max(Key{k}.F(:));
    end
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
    K2 = P.K(:,:,ref);
    R2 = P.R(:,:,ref);
    T2 = P.T(:,ref);
    % object
    plane = zeros(0,3);
    GMM_Name = [];
    parallax = zeros(0,size(Key{ref}.Model.parallax,2));
    ObjTable = 1:ONum;
    % depth plane
    depthpln = zeros(0,3);
    FTable = 1:FNum;
    DispFromObj = zeros(sz);
    Oid = unique(Key{ref}.Model.segMap);
    for i = 1:numel(Oid)
        M = Key{ref}.Model.segMap == Oid(i);
        N = Key{ref}.Model.plane(  Key{ref}.Model.table(Oid(i)),:);
        DispFromObj(M) = -( X(M) * N(1) + Y(M) * N(2) + N(3));
    end
    % combine
    WarpInfo = cell(size(Key));
    for k=1:size(P.K,3)
        % Kf: trans from keyframe, K2: trans to
        if k == ref
            plane = [plane ; Key{ref}.Model.plane];
            GMM_Name = [GMM_Name ; Key{ref}.Model.GMM_Name];
            parallax = [parallax ; Key{ref}.Model.parallax];
            depthpln = [depthpln ; Key{ref}.Model.depthpln];
            continue;
        end
        Kf = P.K(:,:,k);
        Rf = P.R(:,:,k);
        Tf = P.T(:,k);
        
        
        tmp_mat = K2 * R2' * Rf;
		tmp_vec = K2 * R2' * ( Tf - T2);
		tmp_vec = repmat(tmp_vec, [1 prod(sz)]);
        epl_pts = ( tmp_mat * (Kf\(pts')) + repmat(Key{k}.D(:)',[3 1]) .* tmp_vec )';
		epl_pts = epl_pts ./ repmat(epl_pts(:,3),[1 3]);
		[WarpInfo{k}.F WarpInfo{k}.segMap] = GetWarp(Key{k}.Model.F,Key{k}.Model.segMap,Key{k}.D,epl_pts(:,1)-0.5,epl_pts(:,2)-0.5);
        % recover D on ref frame using converted plane parameter
        DObj = zeros(sz);
        DF = zeros(sz);
        plane = [plane ; ((Key{k}.Model.plane * Kf * Rf' * R2)/K2)./repmat(1 + Key{k}.Model.plane * Kf * Rf' * (T2-Tf),[1 3])];
        depthpln = [depthpln ; ((Key{k}.Model.depthpln * Kf * Rf' * R2)/K2)./repmat(1 + Key{k}.Model.depthpln * Kf * Rf' * (T2-Tf),[1 3])];
        
        % from depth plane
        Fid = unique(WarpInfo{k}.F);
        Fid = Fid(Fid>1);
        for i=1:numel(Fid)
            M = WarpInfo{k}.F == Fid(i);
            N = depthpln(  Fid(i),:);
            DF(M) = -( X(M) * N(1) + Y(M) * N(2) + N(3));
        end
        WarpInfo{k}.D = DF;
        DiffOfKandRef = abs(DF - Key{ref}.D);
        for i=1:numel(Fid)
            M = WarpInfo{k}.F == Fid(i);
            refID = unique(Key{ref}.Model.F(M));
            meanVal = ones(numel(refID),1);
            for j=1:numel(refID)
                boolMap = M & (Key{ref}.Model.F == refID(j));
                if mean(DiffOfKandRef(boolMap)) < threshold
                    meanVal(j) = 1/nnz(boolMap);
                end
            end
            [val idx] = min(meanVal);
            if val< threshold
                FTable(Fid(i)) = refID(idx);
%                 fprintf('%d -> %d \n',Fid(i),refID(idx));
            end
        end
        % from object plane
        Oid = unique(WarpInfo{k}.segMap);
        Oid = Oid(Oid>1);
        for i=1:numel(Oid)
            M = WarpInfo{k}.segMap == Oid(i);
            N = plane(  Oid(i),:);
            DObj(M) = -( X(M) * N(1) + Y(M) * N(2) + N(3));
        end
        DiffOfKandRef = abs(DObj - DispFromObj);
        for i=1:numel(Oid)
            M = WarpInfo{k}.segMap == Oid(i);
            refID = unique(Key{ref}.Model.segMap(M));
            for j=1:numel(refID)
                if mean(DiffOfKandRef(M & (Key{ref}.Model.segMap == refID(j)))) < threshold
                    ObjTable(Oid(i)) = refID(j);
%                     fprintf('%d -> %d \n',Oid(i),refID(j));
                end
            end
        end
        % rebuild F and segMap
        Key{k}.Model.F = FTable(Key{k}.Model.F);
        Key{k}.Model.segMap = ObjTable(Key{k}.Model.segMap);
        Key{k}.Model.table(ObjTable(Key{k}.Model.table>0)) = Key{k}.Model.table(Key{k}.Model.table>0);
        
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