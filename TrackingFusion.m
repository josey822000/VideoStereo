function [D info ObjModel] = TrackingFusion(vals, Model, options)
% $Id: ojw_stereo_optim.m,v 1.2 2008/11/17 11:27:35 ojw Exp $

% Create initial arrays
info.map = zeros(vals.sz, 'uint16');
max_iters = options.max_iters + options.average_over + 1;
info.energy = zeros(max_iters, numel(vals.improve));
info.numbers = zeros(4, max_iters, numel(vals.improve), 'uint32'); % Number [updated; unlabelled; independent regions]
info.timings = zeros(4, max_iters, numel(vals.improve)); % Cumulative timings of [proposal; data; smoothness; optimization; finish]

% Initialise depth map

nowModel = Model{1};
options.converge = options.converge * 0.01 * options.average_over;
iter = options.average_over + 2;
info.energy(1:options.average_over+1,end) = realmax;
info.energy(iter,end) = realmax / 1e20;
% R = vals.R(1:end/2,:);
% R = uint8(reshape(permute(R,[1 3 2]),size(vals.I{1})));
while (1-(info.energy(iter,end)/info.energy(iter-options.average_over,end)) > options.converge) && iter < max_iters
    iter = iter + 1;
    % Set the new (proposal) depth map
    t_start = cputime; % Start timing
    n = iter-(options.average_over+1);
    n = mod(n-1, size(Model, 1))+1;
    newModel = Model{n};
    
    
    %Dnew((Dnew>vals.d_min+vals.d_step)) = vals.d_min+vals.d_step;
    
    info.timings(1,iter,:) = cputime - t_start; % Time proposal generation

    % Fuse the depths
%     [M stats info.energy(iter,:) V] = obj_fuse_depths(D, Dnew, vals, ObjModel,Model.segpln_Obj{n},options);
    
    [M stats info.energy(iter,:) V] = Tube_fuse_depths(nowModel, newModel, vals, options);
        
    try
        a = 1;
    catch
        if iter == (options.average_over + 2)
            rethrow(lasterror);
        end
        % Fusing failed - perhaps user bailed early.
        % Output current state.
        info.error = lasterr;
        disp(info.error);
        iter = iter - 1;
        break;
    end
    fprintf(['iter: ' num2str(n)]);
    
    
    nowModel = RefreshObj(M,nowModel,newModel,R,iter,vals.disps,vals.SEI,P);
    
    % 
    info.map(M) = iter;
    info.timings(2:4,iter,:) = stats.timings + info.timings(1,iter,end);
    info.numbers(:,iter,:) = stats.numbers;
    
    % Record progress of disparity map
    save_progress(options.save_name, 'info');
end
[D ObjModel] = RefitObj(nowModel,R,options.imout,vals.disps,vals.SEI,P);
if nargout > 1
    % Save stats
    info.energy = info.energy(options.average_over+2:iter,:);
    info.numbers = info.numbers(:,options.average_over+2:iter,:);
    info.timings = info.timings(:,options.average_over+2:iter,:);
    % Save visibility & mapping
    info.vis = reshape(V, vals.sz(1), vals.sz(2), []);
    info.map = info.map - uint16(options.average_over + 1);
end
return

function [D ObjModel] = RefreshObj(M,nowModel,newModel,R,iter,disps,SEI,P)
    % Switch off annoying warnings
    table = zeros(numel(ObjModel.table),1,'uint32');
    D(M) = Dnew(M);
    oldNum = numel(unique(ObjModel.segMap(~M)));
    oldObjId = ObjModel.table(unique(ObjModel.segMap(~M)));
    table(unique(ObjModel.segMap(~M))) = 1:numel(oldObjId);
    newObjId = DnewModel.table(unique(DnewModel.segMap(M)));
    table(unique(DnewModel.segMap(M))) = oldNum+uint32((1:numel(newObjId)));
    
    
    ObjModel.GMM_Name = ObjModel.GMM_Name(ObjModel.table(unique(ObjModel.segMap(~M))));
    ObjModel.GMM_Name = [ObjModel.GMM_Name; DnewModel.GMM_Name(DnewModel.table(unique(DnewModel.segMap(M))))];
    
    ObjModel.parallax = ObjModel.parallax(ObjModel.table(unique(ObjModel.segMap(~M))),:);
    ObjModel.parallax = [ObjModel.parallax; DnewModel.parallax(DnewModel.table(unique(DnewModel.segMap(M))),:)];
    
    ObjModel.plane = ObjModel.plane(ObjModel.table(unique(ObjModel.segMap(~M))),:);
    ObjModel.plane = [ObjModel.plane; DnewModel.plane(DnewModel.table(unique(DnewModel.segMap(M))),:)];
    
    ObjModel.table = table;
    
    ObjModel.segMap(M) = DnewModel.segMap(M);
    ObjModel.F(M) = DnewModel.F(M);
	
    ObjModel.otherView = ObjWarp(D,ObjModel.F,ObjModel.segMap,P);
    disp(numel(unique(ObjModel.segMap)));
    
    % Reset warnings
return
% problem remain here about index
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


function [D ObjModel] = ExpanRefreshObj(M,D,Dnew,ObjModel,DnewModel,R,iter,disps,SEI,P)
    % Switch off annoying warnings
    table = zeros(numel(ObjModel.table),1,'uint32');
    D(M) = Dnew(M);
    ObjModel.segMap(M) = DnewModel.segMap(M);
    ObjModel.F(M) = DnewModel.F(M);
    % refresh object use old table
    ObjIds = unique(ObjModel.segMap(~M));
    if ObjIds(1) ==0
        ObjIds = ObjIds(2:end);
    end
    oldNum = numel(ObjIds);
    table(ObjIds) = 1:oldNum;
    oldObjPos = ObjModel.table(ObjIds);
    ObjModel.GMM_Name = ObjModel.GMM_Name(oldObjPos);
    ObjModel.parallax = ObjModel.parallax(oldObjPos,:);
    ObjModel.plane = ObjModel.plane(oldObjPos,:);
    newObjIds = unique(DnewModel.segMap(M));
    if numel(newObjIds) ~= 0
        if newObjIds(1) ==0
            newObjIds = newObjIds(2:end);
        end
        newObjPos = DnewModel.table(newObjIds);
        table(newObjIds) = oldNum+uint32((1:numel(newObjIds)));
        ObjModel.GMM_Name = [ObjModel.GMM_Name; DnewModel.GMM_Name(newObjPos)];
        ObjModel.parallax = [ObjModel.parallax; DnewModel.parallax(newObjPos,:)];
        ObjModel.plane = [ObjModel.plane; DnewModel.plane(newObjPos,:)];
    end
    ObjModel.table = table;
    ObjModel.otherView = ObjWarp(D,ObjModel.F,ObjModel.segMap,P);
    disp(numel(unique(ObjModel.segMap)));
    % Reset warnings
return


