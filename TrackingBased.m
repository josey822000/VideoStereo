function Model = TrackingBased( Key, KeyP, Pos , Ref ,Winsz, WarpedInfo)
%TRACKINGBASED Summary of this function goes here
%   Key: model of frames in the tube
%   KeyP: projection matrix of frames in the tube
%   Pos: position(y,x) of window
%   Ref: reference frame
%   WinSz: window size
%   Detailed explanation goes here
    % take out window
    ObjTube = zeros([Winsz size(Key)],class(Key{1}.Model.segMap));
    DispTube = zeros([Winsz size(Key)]);
    sz = size(Key{1}.D);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    ptsAll = ones(sz(1)*sz(2), 3);
    ptsAll(:,1) = X(:);
    ptsAll(:,2) = Y(:);
    pts = [Pos 1]';
    % get tube from a ref frame warp by disparity
    WindowModel = cell(size(Key),1);
    Kf = KeyP.K(:,:,Ref);
    Rf = KeyP.R(:,:,Ref);
    Tf = KeyP.T(:,Ref);
    ObjCollect = [];
    for i=1:size(Key)
        if i == Ref
            epl_pts = [Pos 1]';
        end
        K2 = KeyP.K(:,:,i);
        R2 = KeyP.R(:,:,i);
        T2 = KeyP.T(:,i);

        tmp_mat = K2 * R2' * Rf;
        tmp_vec = K2 * R2' * ( Tf - T2);
        
        WindowModel{i}.plane = ((Key{i}.Model.plane * Kf * Rf' * R2)/K2)./repmat(1 + Key{i}.Model.plane * Kf * Rf' * (T2-Tf),[1 3]);
        WindowModel{i}.table = Key{i}.Model.table;
        % one point
        epl_pts = ( tmp_mat * (Kf\(pts)) + Key{Ref}.D(sub2ind(size(Key{Ref}.D),Pos(2),Pos(1))) .* tmp_vec )';
        epl_pts = epl_pts ./ repmat(epl_pts(:,3),[1 3]);
        % save the window point
        WindowModel{i}.pts = epl_pts(1:2);
        figure(i+1); imshow(Key{i}.D/0.0053); hold on; DrawCircle(epl_pts(1), epl_pts(2), 3, 5, 'r');
        % object
        WindowModel{i}.segMap = Key{i}.Model.segMap(epl_pts(2):epl_pts(2)+Winsz(2), epl_pts(1):epl_pts(1)+Winsz(1));
        ObjCollect = [ObjCollect ;unique(WindowModel{i}.segMap)];
        ObjId = unique(WindowModel{i}.segMap);
        WindowModel{i}.table = 1:numel(Key{i}.Model.table);
        WindowModel{i}.table(ObjId) = 1:numel(ObjId);
        WindowModel{i}.plane = Key{i}.Model.plane(Key{i}.Model.table(ObjId),:);
        WindowModel{i}.parallax = Key{i}.Model.parallax(Key{i}.Model.table(ObjId),:);
        WindowModel{i}.GMM_Name = Key{i}.Model.GMM_Name(Key{i}.Model.table(ObjId));
        % f
        WindowModel{i}.F = Key{i}.Model.F(Pos(2):Pos(2)+Winsz(2), Pos(1):Pos(1)+Winsz(1));
        % D
        WindowModel{i}.D = Key{i}.D(Pos(2):Pos(2)+Winsz(2),Pos(1):Pos(1)+Winsz(1));
        
        
    end
    
    % Model -> for each object 
    %       Frame -> for each frame in Tube
    %             F, segMap , table , GMM, parallax
    Model = cell(numel(unique(ObjTube)),1);
    for i=1:size(numel(unique(ObjTube)))
        Model{i} = cell(size(Key),1);
    end
    for i=1:size(Key)
        ObjId = unique(WindowModel{i}.segMap);
        
        % warp segMap warp F of this object
        for view = 1:size(Key)
            % get warped info F,segMap,D from WarpedInfo
            WarpedF = WarpedInfo{i,view}.F(WindowModel{view}.pts(2):WindowModel{view}.pts(2)+Winsz(2), WindowModel{view}.pts(1):WindowModel{view}.pts(1)+Winsz(1));
            WarpedObj = WarpedInfo{i,view}.segMap(WindowModel{view}.pts(2):WindowModel{view}.pts(2)+Winsz(2), WindowModel{view}.pts(1):WindowModel{view}.pts(1)+Winsz(1));
            WarpedD = WarpedInfo{i,view}.D(WindowModel{view}.pts(2):WindowModel{view}.pts(2)+Winsz(2), WindowModel{view}.pts(1):WindowModel{view}.pts(1)+Winsz(1));
            for j=1:numel(ObjId)
                if view == i 
                    Model{j}{i}.D = WindowModel{i}.D;
                    Model{j}{i}.F = WindowModel{i}.F;
                    Model{j}{i}.segMap = WindowModel{i}.segMap;
                    Model{j}{i}.GMM_Name = WindowModel{i}.GMM_Name;
                    Model{j}{i}.parallax = WindowModel{i}.parallax;
                    Model{j}{i}.table = WindowModel{i}.table;
                    Model{j}{i}.plane = WindowModel{i}.plane;
                    continue;
                end
                % expand object and convert in tube view
                Map = WarpedObj == ObjId(j);
                Model{j}{view}.D(Map) = WarpedD(Map);
                Model{j}{view}.D(~Map) = WindowModel{view}.D(~Map);
                Model{j}{view}.F(Map) = WarpedF(Map);
                Model{j}{view}.F(~Map) = WindowModel{view}.F(~Map);
                warpObjId = unique(WarpedObj(Map));
                Model{j}{view}.segMap(Map) = WarpedObj(Map);
                viewObjId = unique(WindowModel{view}.segMap(~Map));
                Model{j}{view}.segMap(~Map) = WindowModel{view}.segMap(~Map);
                Model{j}{view}.table = zeros(numel(WindowModel{view}.table),1);
                % put warp index then view index
                Model{j}{view}.table(warpObjId) = 1:numel(warpObjId);
                Model{j}{view}.table(viewObjId) = numel(warpObjId)+(1:numel(viewObjId));
                Model{j}{view}.GMM_Name = WindowModel{i}.GMM_Name(WindowModel{i}.table(warpObjId));
                Model{j}{view}.GMM_Name = [Model{j}{view}.GMM_Name; WindowModel{view}.GMM_Name(WindowModel{view}.table(viewObjId) )];
                Model{j}{view}.parallax = WindowModel{i}.parallax(WindowModel{i}.table(warpObjId),:);
                Model{j}{view}.parallax = [Model{j}{view}.parallax ; WindowModel{view}.parallax(WindowModel{view}.table(viewObjId),:)];
                Model{j}{view}.plane = WindowModel{i}.plane(WindowModel{i}.table(warpObjId),:);
                Model{j}{view}.plane = [Model{j}{view}.plane ; WindowModel{view}.plane(WindowModel{view}.table(viewObjId),:)];

            end
        end
    end
    % extend the object plane to other proposals and fuse them
    
    
end
function [FMap segMap DMap] = TubeWarp(src,dst,srcModel,srcD,P)
    Kf = P.K(:,:,src);
    Rf = P.R(:,:,src);
    Tf = P.T(:,src);
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
