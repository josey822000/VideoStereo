function [ExpD ExpModel] = FillHole(DInit,ModelInit,SEI)
%FILLHOLE Summary of this function goes here
%   Detailed explanation goes here
    ExpD = zeros(size(DInit,1),size(DInit,2),0);
    sz = size(DInit);
    ExpModel = cell(0,0);
    for i = 1: size(DInit,3)
        D = DInit(:,:,i);
        Model = ModelInit{i};
        tmpModel = Model;
        tmp = Model.segMap(SEI);
        DiffObjIdx = repmat(any(diff(int32(tmp)),1), [2 1]);
        Neigh = reshape(tmp(DiffObjIdx),2,[]);
        Neigh = Neigh(:,Neigh(1,:) == 0 | Neigh(2,:) == 0);
        % neighbor with occlude
        Neigh = unique(Neigh);
        Neigh = Neigh(2:end);
        M = Model.segMap == 0;
        [X Y] = ind2sub(sz(1:2),find(M));
        InsideD = zeros(size(DInit,1),size(DInit,2),numel(Neigh));
        InsideModel = cell(numel(Neigh),1);
        for s = 1:numel(Neigh)
            % expand the obj index to hole
            plane = Model.plane(Model.table(Neigh(s)),:);
            Dtmp = -( X * plane(1) + Y * plane(2) + plane(3) );
            D(M) = Dtmp;
            InsideD(:,:,s) = D;
            tmpModel.segMap(M) = Neigh(s);
            for o = 1:size(tmpModel.otherView,1)
                M2 = tmpModel.otherView{o}.segMap == 0;
                tmpModel.otherView{o}.segMap(M2) = Neigh(s);
                tmpModel.otherView{o}.F(M2) = Neigh(s);
            end
            InsideModel{s} = tmpModel;
        end
        fprintf('\n%d',numel(Neigh));
        ExpD = cat(3,ExpD,InsideD);
        ExpModel = cat(1,ExpModel,InsideModel);
    end

end

