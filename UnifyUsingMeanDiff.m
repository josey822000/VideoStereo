function Key = UnifyUsingMeanDiff( Key,P , ref, threshold)
%UNIFYUSINGTHR Summary of this function goes here
%   Detailed explanation goes here
    ONum = max(Key{end}.Model.segMap(:));
    FNum = max(Key{end}.Model.F(:));
    Table = 1:ONum;
    
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
                fprintf('%d -> %d \n',Fid(i),refID(idx));
            end
        end
        figure(2); imshow(DF/0.0053);
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
                    fprintf('%d -> %d \n',Oid(i),refID(j));
                end
            end
        end
        figure(3); imshow(DObj/0.0053);
        % rebuild F and segMap
        Key{k}.Model.F = FTable(Key{k}.Model.F);
        Key{k}.Model.segMap = ObjTable(Key{k}.Model.segMap);
        figure(3*k+2); imshow(DiffOfKandRef/threshold);
        figure(3*k+3); sc([Key{k}.Model.F Key{ref}.Model.F],'rand');
        figure(3*k+4); sc([Key{k}.Model.segMap Key{ref}.Model.segMap],'rand');
        pause;
    end
    %% filter Obj less than Threshold
    sameTable = false(ONum,ONum);
    channelTable = false(ONum,ONum);
    aThr = 1.e6;
    bThr = 1.e6;
    cThr = 1.e4;
    plane(:,1) = floor(plane(:,1)*aThr)/aThr;
    plane(:,2) = floor(plane(:,2)*bThr)/bThr;
    plane(:,3) = floor(plane(:,3)*cThr)/cThr;
    
    %% looking for the same object plane
    [T M] = sortrows(plane(:,1));
    % find a within aThr
    T = abs(diff(T));
    T = find(T == 0);
    % if the index is continue than this plane is inside a Thr
    diffT = diff(T);
    diffT = find(diffT > 1);
    diffT(end+1) = numel(T);
    start = 0;
    for s = 1:numel(diffT)
        a = M([T(start+1:diffT(s)); T(diffT(s))+1] );
        c = reshape( repmat(a', numel(a), 1), numel(a) * numel(a), 1 );
        d = repmat(a(:), length(a), 1);
        sameTable((c-1)*double(ONum) + d) = 1;
        start = diffT(s);
    end
    fprintf('samePlane at channel 1:%d\n',numel(T));
    [T M] = sortrows(plane(:,2));
    T = abs(diff(T));
    T = find(T == 0);
    diffT = diff(T);
    diffT = find(diffT > 1);
    diffT(end+1) = numel(T);
    start = 0;
    for s = 1:numel(diffT)
        a = M([T(start+1:diffT(s)); T(diffT(s))+1] );
        c = reshape( repmat(a', numel(a), 1), numel(a) * numel(a), 1 );
        d = repmat(a(:), length(a), 1);
        channelTable((c-1)*double(ONum) + d) = 1;
        start = diffT(s);
    end
    fprintf('samePlane at channel 2:%d\n',numel(T));
    sameTable = sameTable & channelTable;
    channelTable = false(ONum,ONum);
    [T M] = sortrows(plane(:,3));
    T = abs(diff(T));
    T = find(T == 0);
    diffT = diff(T);
    diffT = find(diffT > 1);
    diffT(end+1) = numel(T);
    start = 0;
    for s = 1:numel(diffT)
        a = M([T(start+1:diffT(s)); T(diffT(s))+1] );
        c = reshape( repmat(a', numel(a), 1), numel(a) * numel(a), 1 );
        d = repmat(a(:), length(a), 1);
        channelTable((c-1)*double(ONum) + d) = 1;
        start = diffT(s);
    end
    fprintf('samePlane at channel 3:%d\n',numel(T));
    sameTable = sameTable & channelTable;
    sameTable = sameTable | sameTable';
    
    ObjectTable = 1:ONum;
    [X Y] = find(sameTable);
    unqY = unique(Y);
    for o = 1: numel(unqY)
        ObjectTable(X(Y == unqY(o))) = min(X(Y == unqY(o)));
    end
    fprintf('Total plane:%d\n',numel(unique(ObjectTable)));
    for k=1:size(P.K,3)
        Key{k}.Model.table(ObjectTable(Key{k}.Model.table>0)) = Key{k}.Model.table(Key{k}.Model.table > 0);
        Key{k}.Model.segMap = ObjectTable(Key{k}.Model.segMap);
        Key{k}.Model.GMM_Name = GMM_Name(ObjectTable(unique(Key{k}.Model.segMap)));
        Key{k}.Model.parallax = parallax(ObjectTable(unique(Key{k}.Model.segMap)),:);
%         Key{k}.Model.plane = plane( ObjectTable(unique(Key{k}.Model.segMap)),:);
    end
    %% filter Depth less than Threshold
    sameTable = false(FNum,FNum);
    
    
    depthpln(:,1) = floor(depthpln(:,1)*aThr)/aThr;
    depthpln(:,2) = floor(depthpln(:,2)*bThr)/bThr;
    depthpln(:,3) = floor(depthpln(:,3)*cThr)/cThr;
    
    %% looking for the same F plane
    [T M] = sortrows(depthpln(:,1));
    % find a within aThr
    T = abs(diff(T));
    T = find(T == 0);
    % if the index is continue than this plane is inside a Thr
    diffT = diff(T);
    diffT = find(diffT > 1);
    diffT(end+1) = numel(T);
    start = 0;
    for s = 1:numel(diffT)
        a = M([T(start+1:diffT(s)); T(diffT(s))+1] );
        c = reshape( repmat(a', numel(a), 1), numel(a) * numel(a), 1 );
        d = repmat(a(:), length(a), 1);
        sameTable((c-1)*double(FNum) + d) = 1;
        start = diffT(s);
    end
    fprintf('samePlane:%d\n',numel(T));
    channelTable = false(FNum,FNum);
    [T M] = sortrows(depthpln(:,2));
    T = abs(diff(T));
    T = find(T == 0);
    diffT = diff(T);
    diffT = find(diffT > 1);
    diffT(end+1) = numel(T);
    start = 0;
    for s = 1:numel(diffT)
        a = M([T(start+1:diffT(s)); T(diffT(s))+1] );
        c = reshape( repmat(a', numel(a), 1), numel(a) * numel(a), 1 );
        d = repmat(a(:), length(a), 1);
        channelTable((c-1)*double(FNum) + d) = 1;
        start = diffT(s);
    end
    fprintf('samePlane:%d\n',numel(T));
    sameTable = sameTable & channelTable;
    channelTable = false(FNum,FNum);
    [T M] = sortrows(depthpln(:,3));
    T = abs(diff(T));
    T = find(T == 0);
    diffT = diff(T);
    diffT = find(diffT > 1);
    diffT(end+1) = numel(T);
    start = 0;
    for s = 1:numel(diffT)
        a = M([T(start+1:diffT(s)); T(diffT(s))+1] );
        c = reshape( repmat(a', numel(a), 1), numel(a) * numel(a), 1 );
        d = repmat(a(:), length(a), 1);
        channelTable((c-1)*double(FNum) + d) = 1;
        start = diffT(s);
    end
    fprintf('samePlane:%d\n',numel(T));
    sameTable = sameTable & channelTable;
    sameTable = sameTable | sameTable';
    
    FTable = 1:FNum;
    [X Y] = find(sameTable);
    unqY = unique(Y);
    for f = 1: numel(unqY)
        FTable(X(Y == unqY(f))) = min(X(Y == unqY(f)));
    end
    fprintf('samePlane:%d\n',numel(unique(FTable)));

    for k=1:size(P.K,3)
        Key{k}.Model.F = FTable(Key{k}.Model.F);
    end
    for k=1:size(P.K,3)
        tmp = Key{k};
        save(['UnifiedModel[' num2str(k) ']'],'tmp');
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