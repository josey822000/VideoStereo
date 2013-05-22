function [ output_args ] = InitialAll( vals, options)
%INITIALALL Summary of this function goes here
%   Detailed explanation goes here
    P = vals.P;
    KeyDir = 'keyframe';
    options.KeyDir = KeyDir;
	step = 15;
	DirName = dir(KeyDir);
	Key = cell(numel(DirName)-2,1);
	keyFrameNum = numel(DirName)-2;
    % rec key frame 
    ObjNum = 0;
    FNum = 0;
	for i= 3:numel(DirName)
		keypath = DirName(i).name;
        fprintf('%s',keypath);
		tmp = load(fullfile(KeyDir,keypath,'Final_D'));
		Key{i-2}.D = tmp.D;
		tmp = load(fullfile(KeyDir,keypath,'Final_Model'));
        Key{i-2}.Key = 'KeyGMM';
		Key{i-2}.Model = tmp.ObjModel;
        oldObjId = 1:max(tmp.ObjModel.segMap(:));
        oldObjId(unique(tmp.ObjModel.segMap)) = 1:numel(unique(tmp.ObjModel.segMap));
        Key{i-2}.Model.segMap = oldObjId(Key{i-2}.Model.segMap);
        Key{i-2}.Model.segMap = Key{i-2}.Model.segMap + ObjNum;
        oldObjId = 1:max(tmp.ObjModel.F(:));
        oldObjId(unique(tmp.ObjModel.F)) = 1:numel(unique(tmp.ObjModel.F));
        Key{i-2}.Model.F = oldObjId(Key{i-2}.Model.F);
        Key{i-2}.Model.F = Key{i-2}.Model.F + FNum;
        ObjNum = max(Key{i-2}.Model.segMap(:));
        FNum = max(Key{i-2}.Model.F(:));
	end
	% load P
    for i=1:keyFrameNum
        table = zeros(ObjNum,1);
        table(unique(Key{i}.Model.segMap)) = 1:numel(unique(Key{i}.Model.segMap));
        Key{i}.Model.table = table;
    end
	useKey = [1:size(P.K,3)];
	for i=16:size(P.K,3)
		if mod(i,15) == 0
			continue;
		end
		% warp to now frame
		
		% hole stay -1 label and assign big cost(no information from this)
        key1 = floor(i/step)*step;
        key2 = ceil(i/step)*step;
        if key1 ==0
            key1 = 1;
            
        end
        
        useKey = [key1 key2];
        tmpvals = vals;
        tmpvals.I = vals.I(useKey);
        tmpvals.R = repmat(reshape(single(vals.I{i}), [], 3), [2 1]);
        tmpvals.P.K = P.K(:,:,[i useKey]);
        tmpvals.P.R = P.R(:,:,[i useKey]);
        tmpvals.P.T = P.T(:,[i useKey]);
        tmpvals.useKey = useKey;
        tmpKey = Key((useKey/step)+1);
        fprintf('Key:%d %d',useKey(1)/step+1,useKey(2)/step+1);
		[D Model] = Warp2CurrentView(tmpKey,tmpvals.P,i,step);
        finalD = D(:,:,1);
        Holemap = double(Model{1}.segMap ~= 0);
        kD = D(:,:,2);
        newHoleMap = (Model{2}.segMap ~= 0);
        finalD(newHoleMap) = finalD(newHoleMap) + kD(newHoleMap);
        Holemap = Holemap + double(newHoleMap);
        finalD = (finalD + 1.e-10) ./ Holemap;
        delete([options.sequence '_' num2str(i-1) '.mat']);
        figure(10); imshow(finalD/0.0087);
        save([options.sequence '_' num2str(i-1) '.mat'],'finalD');
        %
        fpath = [num2str(f) '_InitAll'] ;
        mkdir(fpath);
        cd(fpath);
        mkdir('GMM');

        useKey = [key1 key2];
        tmpvals = vals;
        tmpvals.I = images(useKey);
        tmpvals.R = repmat(reshape(single(images{f}), [], 3), [2 1]);
        tmpvals.P.K = P.K(:,:,[f useKey]);
        tmpvals.P.R = P.R(:,:,[f useKey]);
        tmpvals.P.T = P.T(:,[f useKey]);
        tmpvals.useKey = useKey;
        tmpKey = Key(useKey/step);
        tmpimages = images([f useKey]);
        tmpoptions = options;
        tmpoptions.imout = f;
        if ~exist('Dproposals.mat');
            [Dproposals info.segpln_gen] = ojw_segpln(tmpimages, tmpvals.P, disps, images{f}, tmpoptions);
            tmp = info.segpln_gen;
            delete('Dproposals.mat','segpln_gen.mat');
            save('Dproposals','Dproposals');
            save('segpln_gen','tmp');
            clear tmp
        else
        % load segpln_gen & Dproposals
            tmp = load('segpln_gen');
            info.segpln_gen = tmp.tmp;
            Dproposals = load('Dproposals');
            Dproposals = Dproposals.Dproposals;
            clear tmp
        end
        % Truncate extreme values in Dproposal
        % if Truncate the cost may not right

        % object depth segmentation
        info.segpln_Obj = cell(size(info.segpln_gen.segments,3),1);
        [X Y] = meshgrid(1:sz(2),1:sz(1));
        plane = info.segpln_gen.plane;

        % load segpln_Obj

        FsegNum = 0;
        OsegNum = 0;


        if ~exist('segpln_Obj.mat')
            for i = 1:size(Dproposals,3)
                %update segMap(object map) and plane to global numbers
                if ~exist(['Obj[' num2str(i) '].mat'])
                    % try gco

                    segment = info.segpln_gen.segments(:,:,i);
                    figure(2); 
                    sc(segment,'rand');
                    info.segpln_Obj{i}.F = segment;
                    info.segpln_Obj{i}.depthPlane = info.segpln_gen.plane{i};
                    %info.segpln_Obj{i}.F = label_table(segment);

                    segNum = max(max(segment));
                    h = GCO_Create(segNum,segNum);

                    % set data term
                    data = zeros(segNum,segNum, 'int32');

                    for sid = 1:segNum %row
                        N = plane{i}(sid,:)';

                        planeD = -(X * N(1) + Y * N(2) + N(3));
                        planeD(planeD < vals.d_min) = -2*vals.d_step;
                        planeD(planeD > vals.d_min+vals.d_step) = 2*vals.d_step;
                        % planeDiff: 1*pix
                        planeDiff = abs(planeD(:)-reshape(Dproposals(:,:,i),[],1))/vals.d_step;
                        data(sid,:) = accumarray(reshape(segment,[],1),planeDiff)';
                    end
                    %data = data/max(max(data));
                    clear planeD planeDiff mapp
                    tmp = int32(segment(vals.SEI));
                    % depth diff on boundary

                    DiffObjIdx = repmat(any(diff(int32(segment(vals.SEI))),1), [2 1]);
                    Neigh = reshape(tmp(DiffObjIdx),2,[]);

                    clear DiffObjIdx meanColor4eachSeg
                    List = (Neigh(1,:)-1)*int32(segNum)+Neigh(2,:);
                    List = [List (Neigh(2,:)-1)*int32(segNum)+Neigh(1,:)];


                    % boundary length
                    pixNumInSeg = histc(segment(:),1:segNum);
                    pixNumInSeg = repmat(pixNumInSeg,[1 segNum]) + repmat(pixNumInSeg',[segNum 1]);
                    Smooth = zeros(1,segNum*segNum, 'int32');
                    Smooth = histc(List,1:segNum*segNum);
                    pixNumInSeg = pixNumInSeg(:) .* (Smooth(:)>0);
                    pixNumInSeg = 9000/pixNumInSeg(:);
                    Smooth = Smooth*0.5;
                    Smooth = ceil(Smooth.*pixNumInSeg);
                    Smooth = reshape(Smooth,[segNum segNum]);
                    Smooth(isnan(Smooth))=0;


                    clear tmp Neigh List ColorD SmoothNormalize
                    GCO_SetDataCost(h,data);
                    GCO_SetNeighbors(h,Smooth);
                    GCO_Expansion(h);
                    Label = GCO_GetLabeling(h);
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
                    info.segpln_Obj{i}.segMap = segment;
                    info.segpln_Obj{i}.plane = plane{i}(uq,:);
                    figure(3);
                    sc(segment,'rand');
                    fprintf('solving GCO for proposal %d\n', i);

                    clear segment tmpL
                    clear data smooth
                    GCO_Delete(h);

                    info.segpln_Obj{i}.F = info.segpln_Obj{i}.F + FsegNum;
                    info.segpln_Obj{i}.segMap = info.segpln_Obj{i}.segMap + OsegNum;
                    FsegNum = max(max(info.segpln_Obj{i}.F));
                    OsegNum = max(max(info.segpln_Obj{i}.segMap));

                    %warp to other view
                    info.segpln_Obj{i}.otherView = ObjWarp(Dproposals(:,:,i),info.segpln_Obj{i}.F,info.segpln_Obj{i}.segMap,tmpvals.P);           
                    %calc Model            
                    ExpFuse = 0;
                    [info.segpln_Obj{i}.plane info.segpln_Obj{i}.parallax info.segpln_Obj{i}.GMM_Name] = calcObjModel(disps,Dproposals(:,:,i),info.segpln_Obj{i},images{f},f,ExpFuse);
                    delete(['Obj[' num2str(i) ']']);
                    tmp = info.segpln_Obj{i};
                    save(['Obj[' num2str(i) ']'],'tmp');
                else
                    tmp = load(['Obj[' num2str(i) ']']);
                    info.segpln_Obj{i} = tmp.tmp;
                    clear tmp;
                end

            end
            % update table
            OsegNum = max(info.segpln_Obj{14}.segMap(:));
            for i = 1:size(Dproposals,3)
                table = zeros(OsegNum,1,'uint32');
                table(unique(info.segpln_Obj{i}.segMap)) = 1:numel(unique(info.segpln_Obj{i}.segMap));
                info.segpln_Obj{i}.table = table;
            end

            % started to fuse proposals
            tmp = info.segpln_Obj;
            delete('segpln_Obj.mat');
            save('segpln_Obj','tmp');
            clear R tmp
        else  % load segplnObj
            tmp = load('segpln_Obj');
            info.segpln_Obj = tmp.tmp;
            clear tmp
        end
        Expan = false;
        Dproposals = Dproposals(:,:,1:7);
        tmpvals.visibility = false;
        %
        Expan = true;
		tmpvals.visibility = false;
        
        options.max_iters = 5;
        options.average_over = 5;
        options.imout = i;
        % here only two proposals used, fucking lack information, what is
        % the pros. to mark the black region? black region is occluded, the
        % info can be recoverd by another proposal
        close all;
		[D info.InitialEach ObjModel] = obj_fuse_proposals(tmpvals, D, options, Model, Expan);
		delete(['Frame_' num2str(i) '_info.mat'],['Final_' num2str(i) '_D.mat'],['Final_' num2str(i) '_Model.mat']);
		save(['Frame_' num2str(i) '_info.mat'],'info');
		save(['Final_' num2str(i) '_D.mat'],'D');
		save(['Final_' num2str(i) '_Model.mat'],'ObjModel');
			
	end
end

