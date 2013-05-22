function [ output_args ] = BackTrackBasedRefine( P,Frame,Model , SEI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     from middle frame,select a point and refine neighbor
    sz = size(Model{middle}.segMap);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    pts = ones(sz(1)*sz(2), 3);
    pts(:,1) = X(:);
    pts(:,2) = Y(:);

    %-----video disparity calculation ----%
    epl_pts = zeros(prod(sz), 3);

    Kf = P.K(:,:,1);
    Rf = P.R(:,:,1);
    Tf = P.T(:,1);
    
    DiffPos = any(diff(int32(Model{middle}.segMap(SEI))),1);
    NeighborId = SEI(1,DiffPos);
    RefinedMap = zeros(sz,'logical');
    Window_Center = ind2sub(sz,NeighborId);
    Window_Center = Window_Center(Window_Center(:,1)>30 & Window_Center(:,1)<sz(2)-30 & Window_Center(:,2)>30 & Window_Center(:,2)>sz(1)-30);
    for i= 1:size(NeighborId,2)
%       filter   overlap points
        center = Window_Center(i,:);
        Window_Center = Window_Center(~(Window_Center(:,1)>center(1)-30 & Window_Center(:,1)<center(1)+30 & Window_Center(:,2)>center(2)-30 & Window_Center(:,2)<center(2)+30));
        % create graph
        sites = 30*30*size(Frame,3);
        ObjNum = numel(unique(Model{middle}.segMap));
        h = GCO_Create(sites,ObjNum);
        data = zeros(sites,ObjNum,'int32');
        smooth = zeros(sites,sites,'int32');
        for f = 1:max(Frame)
            
    %         find corresponding points at other view

    %       open a window on corresponding position

    %       set data cost ()

    %       set smooth cost ()

    %       solve

    %       fill label into windows

        end
        GCO_SetDataCost(h,data);
        GCO_SetNeighbors(h,Smooth);
        GCO_Expansion(h);
        Label = GCO_GetLabeling(h);
        GCO_Delete(h);
    end
end

