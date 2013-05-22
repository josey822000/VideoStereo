function [D info] = ojw_stereo(images, P, disps, sz, options)

% $Id: ojw_stereo.m,v 1.5 2008/11/17 11:27:35 ojw Exp $

% Crude check for a reference image
%if max(abs(P.P([1:6 9])-[1 0 0 0 1 0 1])) > 1e-12
%    error('First image must be reference image');
%end

% Initialize data arrays
%R = images{1}(round(P.P(8))+(1:sz(1)),round(P.P(7))+(1:sz(2)),:);
R = images{1};
vals.I = images(2:end);
vals.P = P;
vals.disps = disps;
vals.sz = sz;
colors = size(R, 3);
num_in = numel(images);
Rorig = uint8(R);
if colors == 1
    Rorig = repmat(Rorig, [1 1 3]);
end
vals.R = repmat(reshape(single(R), [], colors), [2 1]);
vals.d_min = disps(end);
vals.d_step = disps(1) - vals.d_min;
vals.ndisps = numel(disps);

T = reshape(uint32(1:prod(sz)), sz);
if options.planar
    % Use 2nd order smoothness prior
    SEI = [reshape(T(1:end-2,:), 1, []) reshape(T(:,1:end-2), 1, []); ...
           reshape(T(2:end-1,:), 1, []) reshape(T(:,2:end-1), 1, []); ...
           reshape(T(3:end,:), 1, []) reshape(T(:,3:end), 1, [])];
    if options.connect == 8
        SEI = [SEI [reshape(T(1:end-2,1:end-2), 1, []) reshape(T(3:end,1:end-2), 1, []); ...
                    reshape(T(2:end-1,2:end-1), 1, []) reshape(T(2:end-1,2:end-1), 1, []); ...
                    reshape(T(3:end,3:end), 1, []) reshape(T(1:end-2,3:end), 1, [])]];
    end
else
    % Use 1st order smoothness prior
    SEI = [reshape(T(1:end-1,:), 1, []) reshape(T(:,1:end-1), 1, []); ...
           reshape(T(2:end,:), 1, []) reshape(T(:,2:end), 1, [])];
    if options.connect == 8
        SEI = [SEI [reshape(T(1:end-1,1:end-1), 1, []) reshape(T(2:end,1:end-1), 1, []); ...
                    reshape(T(2:end,2:end), 1, []) reshape(T(1:end-1,2:end), 1, [])]];
    end
end
clear T

% Initialise display
vals.show_output = options.show_output;
if vals.show_output
    vals.show_output = gcf;
    set(0, 'CurrentFigure', vals.show_output);
    subplot('Position', [0 0.5 1/3 0.5]);
    sc(R, [0 255]);
end

% Segment the image using mean shift
info.segment = vgg_segment_ms(Rorig, options.seg_params(1), options.seg_params(2), options.seg_params(3));
% Find smoothness edges which don't cross segmentation boundaries
EW = reshape(~any(diff(int32(info.segment(SEI))), 1), 1, []);
EW = EW * options.lambda_h + ~EW * options.lambda_l;
EW = EW * (num_in / ((options.connect==8) + 1));
EW = reshape(repmat(EW, [4*(1+(options.planar~=0)) 1]), [], 1);

% Set up values for ibr_fuse_depths
vals.visibility = (options.visibility ~= 0) * 1e4;
vals.improve = options.improve;
vals.contract = options.contract;
vals.independent = options.independent;
vals.compress_graph = options.compress_graph;

% Set up our robust kernels
% vals.ephoto = @(F) log(2) - log(exp(sum(F .^ 2, 2)*(-1/(options.col_thresh*colors)))+1);
vals.ephoto = @(F) sqrt(sum(F .^ 2, 2)/3);
switch options.smoothness_kernel
    case 1
        vals.esmooth = @(F) EW .* min(abs(F), options.disp_thresh);
    case 2
        EW = EW / options.disp_thresh;
        vals.esmooth = @(F) EW .* min(F.^2, options.disp_thresh^2);
    otherwise
        error('Unknown smoothness kernel specified');
end
vals.occl_val = options.occl_const + log(2);
vals.SEI = SEI;
clear T SEI EW Rorig

if nargout > 1
    % Save parameters
    % disp_thresh ??????
    info.params.disp_thresh = options.disp_thresh;
    info.params.col_thresh = options.col_thresh;
    info.params.occl_const = options.occl_const;
    info.params.lambda_l = options.lambda_l;
    info.params.lambda_h = options.lambda_h;
end
sigmaColor = 10;
if isnumeric(options.proposal_method) && size(options.proposal_method, 1) == 1
    switch options.act
    case 'InitialKey'
        mkdir(num2str(options.imout));
        cd(num2str(options.imout));
        mkdir GMM
        % Use the proposal methods:
        % SegPln (prototypical segment-based stereo proposals)
        if ~exist('Dproposals.mat');
            [Dproposals info.segpln_gen] = ojw_segpln(images, P, disps, R, options);
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
    %     Dproposals(Dproposals < vals.d_min) = vals.d_min;
    %     Dproposals(Dproposals > vals.d_min+vals.d_step) = vals.d_min+vals.d_step;

        % object depth segmentation
        info.segpln_Obj = cell(size(info.segpln_gen.segments,3),1);
        [X Y] = meshgrid(1:sz(2),1:sz(1));
        plane = info.segpln_gen.plane;

        % load segpln_Obj
    %   tmp = load('segplnRight6_Obj');
    %   info.segpln_Obj = tmp.tmp;
    %   clear tmp
    %
        FsegNum = 0;
        OsegNum = 0;


        if ~exist('segpln_Obj.mat')
            Kf = P.K(:,:,1);
            Rf = P.R(:,:,1);
            Tf = P.T(:,1);
            for i = 1:size(Dproposals,3)
                %update segMap(object map) and plane to global numbers
                if ~exist(['Obj[' num2str(i) '].mat'])
                    % try gco

                    segment = info.segpln_gen.segments(:,:,i);
                    figure(2); sc(segment,'rand');
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
                        planeDiff = abs(planeD(:)-reshape(Dproposals(:,:,i),[],1))/vals.d_step;
                        data(sid,:) = accumarray(reshape(segment,[],1),planeDiff)';
                    end
                    clear planeD planeDiff mapp
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
                    pixNumInSeg = pixNumInSeg(:) .* (Smooth(:)>0);
                    pixNumInSeg = 9000/pixNumInSeg(:);
                    Smooth = ceil(Smooth*0.025);
%                     Smooth = ceil(Smooth.*pixNumInSeg);
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
                    info.segpln_Obj{i}.otherView = ObjWarp(Dproposals(:,:,i),info.segpln_Obj{i}.F,info.segpln_Obj{i}.segMap,P);           
                    %calc Model            
                    ExpFuse = 0;
                    [info.segpln_Obj{i}.plane info.segpln_Obj{i}.parallax info.segpln_Obj{i}.GMM_Name] = calcObjModel(disps,Dproposals(:,:,i),info.segpln_Obj{i},images{1},i,ExpFuse);
                    delete(['Obj[' num2str(i) ']']);
                    tmp = info.segpln_Obj{i};
                    save(['Obj[' num2str(i) ']'],'tmp');
                else
                    tmp = load(['Obj[' num2str(i) ']']);
                    info.segpln_Obj{i} = tmp.tmp;
                    clear tmp;
%                     %warp to other view
%                     info.segpln_Obj{i}.otherView = ObjWarp(Dproposals(:,:,i),info.segpln_Obj{i}.F,info.segpln_Obj{i}.segMap,P);           
%                     %calc Model            
%                     ExpFuse = 0;
%                     [info.segpln_Obj{i}.plane info.segpln_Obj{i}.parallax info.segpln_Obj{i}.GMM_Name] = calcObjModel(disps,Dproposals(:,:,i),info.segpln_Obj{i},images{1},i,ExpFuse,P);
%                     delete(['Obj[' num2str(i) '].mat']);
%                     tmp = info.segpln_Obj{i};
%                     save(['Obj[' num2str(i) ']'],'tmp');
%                     clear tmp;
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
        vals.visibility = false;
        [D info.segpln_optim ObjModel] = obj_fuse_proposals(vals, Dproposals, options, info.segpln_Obj, Expan);
        delete('Final_info.mat','Final_D.mat','Final_Model.mat');
        save('Final_info','info');
        save('Final_D','D');
        save('Final_Model','ObjModel');
        figure(13);imshow(D);
        clear Dproposals    
        cd ..
    case 'InitialAll'
        vals.I = images;
        InitialAll(vals,options);
    case 'InitialByKey'
        vals.I = images;
%         InitialByKey(vals,options);
        step = 15;
        Key = LoadKeyFrame();
        KeyP.K = P.K(:,:,options.KeyFrame);
        KeyP.R = P.R(:,:,options.KeyFrame);
        KeyP.T = P.T(:,options.KeyFrame);
        ref = floor(size(KeyP.K,3)/2);
        Key = UnifyUsingMeanDiff(Key,KeyP,ref,disps(1)/300);
        % Use the proposal methods:
        % SegPln (prototypical segment-based stereo proposals)
        for f=27:size(P.K,3)
            if isempty(find(options.KeyFrame == f, 1))
                continue;
            end
            fpath = num2str(f);
            mkdir(fpath);
            cd(fpath);
            mkdir('GMM');
            key1 = floor(f/step)*step;
            key2 = ceil(f/step)*step;
            if key1 == key2 || key2 ==0
                key2 = key1+step;
            end
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
                        figure(2); sc(segment,'rand');
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
            [D info.segpln_optim ObjModel] = obj_fuse_proposals(tmpvals, Dproposals, tmpoptions, info.segpln_Obj, Expan);
            delete('Final_info.mat','Final_D.mat','Final_Model.mat');
            save('Final_info','info');
            save('Final_D','D');
            save('Final_Model','ObjModel');
            figure(13);imshow(D);
            clear Dproposals    
            cd ..
        end
    case 'FilterOut'
        vals.I = images;
%         InitialByKey(vals,options);
        step = 15;
        Key = LoadKeyFrame();
        KeyP.K = P.K(:,:,options.KeyFrame);
        KeyP.R = P.R(:,:,options.KeyFrame);
        KeyP.T = P.T(:,options.KeyFrame);
        Key = UnifyUsingThr(Key,KeyP);
        % Use the proposal methods:
        FilterKey(vals,options,Key,KeyP);
    end
else
    sprintf('dont go here');
end

return