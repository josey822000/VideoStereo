function [plane F] = UnifyLocal( D, pln, F, threshold)
%UNIFYUSINGTHR Summary of this function goes here
%   Detailed explanation goes here
    
    
    sz = size(F);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    pts = ones(sz(1)*sz(2), 3);
	pts(:,1) = X(:);
	pts(:,2) = Y(:);
    
    

    DF = zeros(sz(1:2));
    % combine inside k
    Fnum = 0;
    plane = [];
    for k=1:sz(3)
        F(:,:,k) = F(:,:,k) + Fnum;
        plane = [plane ;pln{k}];
        Fnum = max(reshape(F(:,:,k),1,[]));
    end
    
    % from depth plane
    Fid = unique(F);
    Ftable = 1:numel(Fid);
    
    % calc difference between planes
    for i=1:numel(Fid)
        if Ftable(i) ~= i
            continue;
        end
        N = plane(  Fid(i),:);
        DF = -( X * N(1) + Y * N(2) + N(3));
        samePln = [];
        for k=1:sz(3)
            Fk = F(:,:,k);
            Dk = D(:,:,k);
            DiffOfKandRef = abs(DF - Dk);
            % find not same set
            samePln = [samePln ;unique(Fk(DiffOfKandRef >= threshold))];
        end
        % find diff set
        samePln = setdiff(1:numel(Fid),samePln);
        samePln = samePln(samePln>i);
        Ftable(samePln) = i;
        F(ismember(F,samePln)) = i;       
    end
    Ftable(unique(Ftable)) = 1:numel(unique(Ftable));
    plane = plane(unique(Ftable),:);
    disp([num2str(k) 'map:' num2str(numel(Fid)) '->' num2str(numel(unique(Ftable)))]);
    F = Ftable(F);
    
    % combine between k
    
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