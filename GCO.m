function [segMap ObjPln] = GCO( Key, vals )
%GCO Summary of this function goes here
%   Detailed explanation goes here
    
    plane = Key.plane;
    sz = size(Key.F);
    [X Y] = meshgrid(1:sz(2),1:sz(1));
    segMap = zeros(sz,'int32');
    ObjPln = cell(sz(3),1);
    
    % for proposal, do depth segmentation
    for i = 1:sz(3)
        %update segMap(object map) and plane to global numbers
        segment = Key.F(:,:,i);
        segment = segment - (min(segment(:))-1);
        figure(2);
        sc(segment,'rand');
        figure(6);imshow(Key.D(:,:,i)/0.0053);
        segNum = max(segment(:));
        h = GCO_Create(segNum,segNum);
        % set data term
        data = zeros(segNum,segNum,'int32');
        for sid = 1:segNum %row
            N = plane(sid,:)';
            planeD = -(X * N(1) + Y * N(2) + N(3));
            planeD(planeD < vals.d_min) = -vals.d_step;
            planeD(planeD > vals.d_min+vals.d_step) = 2*vals.d_step;
            planeDiff = abs(planeD(:)-reshape(Key.D(:,:,i),[],1))/vals.step;
            data(sid,:) = accumarray(reshape(segment,[],1),planeDiff)';
        end
        clear planeD planeDiff mapp
        disp(['data0:' num2str(nnz(data==0))]);
        tmp = int32(segment(vals.SEI));
        % depth diff on boundary
        DiffObjIdx = repmat(any(diff(int32(segment(vals.SEI))),1), [2 1]);
        Neigh = reshape(tmp(DiffObjIdx),2,[]);
%       % depth difference on smooth term
%       tmpD = Dproposals(:,:,i);
%       tmpD = tmpD(vals.SEI);
%       tmpD = abs(diff(reshape(tmpD(DiffObjIdx),2,[])));
%       tmpD = repmat(tmpD,[1 2]);
%       % color difference on smooth term
%             meanColor4eachSeg = accumarray(reshape(segment,[],1),ColorD(:,1));
%             meanColor4eachSeg(:,2) = accumarray(reshape(segment,[],1),ColorD(:,2));
%             meanColor4eachSeg(:,3) = accumarray(reshape(segment,[],1),ColorD(:,3));
%             ColorD = meanColor4eachSeg ./ repmat(histc(segment(:),1:max(segment(:))),[1 3]);
        clear DiffObjIdx meanColor4eachSeg
        List = (Neigh(1,:)-1)*int32(segNum)+Neigh(2,:);
        List = [List (Neigh(2,:)-1)*int32(segNum)+Neigh(1,:)];

%             ColorMap = zeros(segNum,segNum,'int32');
%             ColorMap((Neigh(1,:)-1)*int32(segNum)+Neigh(2,:)) = sum((ColorD(Neigh(1,:),:) - ColorD(Neigh(2,:),:)).^2,2);
%             ColorMap((Neigh(2,:)-1)*int32(segNum)+Neigh(1,:)) = sum((ColorD(Neigh(2,:),:) - ColorD(Neigh(1,:),:)).^2,2);
%             clear ColorD
%             ColorMap = ColorMap/(30*3);
        % boundary length
        pixNumInSeg = histc(segment(:),1:segNum);
        pixNumInSeg = repmat(pixNumInSeg,[1 segNum]) + repmat(pixNumInSeg',[segNum 1]);
        Smooth = zeros(1,segNum*segNum, 'int32');
        Smooth = histc(List,1:segNum*segNum);
%         pixNumInSeg = pixNumInSeg(:) .* (Smooth(:)>0);
%         pixNumInSeg = 9000/pixNumInSeg(:);
%         Smooth = ceil(Smooth*0.025);
% %                     Smooth = ceil(Smooth.*pixNumInSeg);
        Smooth = reshape(Smooth,[segNum segNum]);
        Smooth(logical(eye(segNum,segNum))) = 0;
        Smooth(isnan(Smooth))=0;
        Smooth = double(Smooth>0);
        Smooth = zeros(segNum,segNum);
        clear tmp Neigh List ColorD SmoothNormalize
        
        GCO_SetDataCost(h,data);
        GCO_SetNeighbors(h,Smooth);
        GCO_Expansion(h);
        Label = GCO_GetLabeling(h);
        GCO_Delete(h);
        % after relabel
        uq = unique(Label);
        NewSegNum = numel(uq);
        segment = Label(segment);
        %---- store segment plane for each object (object plane) for later use 
        % info.obj_pln{i} = plane{i}(segment, :);
        %----------------------------------------------%

        tmpL = 1:segNum;
        tmpL(uq) = 1:NewSegNum;
        segment = tmpL(segment);
        segMap(:,:,i) = segment;
        Objpln{i} = plane(uq,:);
        figure(3);
        sc(segment,'rand');
        disp(['solving GCO for proposal' num2str(i)]);
        
        clear segment tmpL
        clear data smooth
        pause;
    end
    % update table
    clear R tmp
end

