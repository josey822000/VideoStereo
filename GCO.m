function [segMap ObjPln] = GCO( Key, img, vals )
%GCO Summary of this function goes here
%   Detailed explanation goes here
    
    Pln = Key.plane;
    sz = size(Key.F);
    [X Y] = meshgrid(1:sz(2),1:sz(1));
    segMap = zeros(sz,'int32');
    ObjPln = cell(sz(3),1);
    img = reshape(img,[],3);
    ColorDiv = [20 25 30 35 40 45 50 55];
    DepthMinus = 4;
    % for proposal, do depth segmentation
    for i = 1:sz(3)
        %update segMap(object map) and plane to global numbers
        segment = Key.F(:,:,i);
        plane = Pln(unique(segment),:);
        Ftable = 1:max(segment);
        Ftable(unique(segment)) = 1:numel(unique(segment));
        segment = Ftable(segment);
        disp(['originseg:' num2str(numel(unique(segment)))]);
        figure(2);
        sc(segment,'rand');
        figure(6);imshow(Key.D(:,:,i)/0.0053);
        segNum = max(segment(:));
        
        % set data term
        data = zeros(segNum,segNum,'int32');
        segCnt = histc(segment(:),1:segNum)';
        ColorList = zeros(segNum,3);
        recoverImg = zeros([prod(sz(1:2)) 3]);
        recoverD = zeros(sz(1:2));
        for sid = 1:segNum
            N = plane(sid,:)';
            M = segment(:) == sid;
            recoverD(M) = -(X(M) * N(1) + Y(M) * N(2) + N(3));
        end
        figure(4);
        imshow(recoverD/0.0053);
        for sid = 1:segNum %row
            N = plane(sid,:)';
            planeD = -(X * N(1) + Y * N(2) + N(3));
            planeD(planeD < vals.d_min) = -vals.d_step;
            planeD(planeD > vals.d_min+vals.d_step) = 2*vals.d_step;
            planeDiff = abs(planeD(:)-reshape(Key.D(:,:,i),[],1))/vals.step;
            data(sid,:) = accumarray(reshape(segment,[],1),planeDiff)'./segCnt;
            ColorList(sid,:) = sum(img(segment(:) == sid,:));
        end
        
        % color
        ColorList = ColorList./repmat(segCnt(:),[1 3]);
        
        for sid=1:segNum
            recoverImg(segment(:)==sid,:) = repmat(ColorList(sid,:),[segCnt(sid) 1]);
        end
        recoverImg = reshape(recoverImg,[sz(1:2) 3]);
        figure(5);
        imshow(recoverImg);
        ColorList = repmat(ColorList,[segNum 1]) - reshape(repmat(ColorList',[segNum 1]),3,[])';
        ColorList = sum(ColorList.^2,2) ;
        clear planeD planeDiff mapp
        disp(['data0:' num2str(nnz(data==0))]);
        tmp = int32(segment(vals.SEI));
        % depth diff on boundary
        DiffObjIdx = repmat(any(diff(int32(segment(vals.SEI))),1), [2 1]);
        Neigh = reshape(tmp(DiffObjIdx),2,[]);
        
        %
        clear DiffObjIdx meanColor4eachSeg
        List = (Neigh(1,:)-1)*int32(segNum)+Neigh(2,:);
        List = [List (Neigh(2,:)-1)*int32(segNum)+Neigh(1,:)];

        % boundary length
%         pixNumInSeg = histc(segment(:),1:segNum);
%         pixNumInSeg = repmat(pixNumInSeg,[1 segNum]) + repmat(pixNumInSeg',[segNum 1]);
%         Smooth = zeros(1,segNum*segNum, 'int32');
        Neigh = histc(List,1:segNum*segNum);
        
        
%         Smooth = 1+(sqrt(ColorList) + 1e-3)/ColorDiv(t);
%         disp(['mean:' num2str(mean2(ColorList(Neigh>0)))]);
%         Smooth = reshape(Smooth,[segNum segNum]);
%         Smooth(logical(eye(segNum,segNum))) = 0;
        clear tmp List ColorD SmoothNormalize

        h = GCO_Create(segNum,segNum);
        GCO_SetDataCost(h,data);
        GCO_SetNeighbors(h,double(reshape(Neigh,segNum,segNum)>0));
%             GCO_SetSmoothCost(h,Smooth);
        GCO_Expansion(h);
        Label = GCO_GetLabeling(h);
        [E D S] = GCO_ComputeEnergy(h);
        disp(['E:' num2str(E) '  D:' num2str(D) '  S:' num2str(S)]);
        GCO_Delete(h);
        % after relabel
        uq = unique(Label);
        NewSegNum = numel(uq);
        tmpSegment = Label(segment);
        %---- store segment plane for each object (object plane) for later use 
        % info.obj_pln{i} = plane{i}(segment, :);
        %----------------------------------------------%
        disp(['newseg:' num2str(NewSegNum)]);
        tmpL = 1:segNum;
        tmpL(uq) = 1:NewSegNum;
        tmpSegment = tmpL(tmpSegment);
        segMap(:,:,i) = segment;
        Objpln{i} = FitPlane(segment,Key.D(:,:,i));
        figure(3);
        sc(tmpSegment,'rand');
        disp(['solving GCO for proposal' num2str(i)]);
        
        clear segment tmpL
        clear data smooth
    end
    % update table
    clear R tmp
end

