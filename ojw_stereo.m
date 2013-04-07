function [D info] = ojw_stereo(images, P, disps, sz, options)
%OJW_STEREO  Global stereo with 2nd-order smoothness prior & occlusion model
%
%   [D info] = ojw_stereo(images, P, disps, sz, options)
%
% Generates a disparity (1/depth) map for the first input image. This
% algorithm implements a "global" stereo algorithm, with asymmetrical
% occlusion modelling, using an alpha-expansion graph cuts style approach,
% but with arbitrary disparity proposals.
%
%IN:
%   images - 1xN cell array of input images, the reference image being
%            images{1}.
%   P - 3x4xN array of projection matrices for the input images, relative
%       to the output image.
%   disps - 1xM list of disparities to sample at.
%   sz - 1x2 vector of output image dimensions: [H W].
%   options - a structure containing the following input parameters:
%       col_thresh - scalar noise parameter for data likelihood.
%       occl_const - scalar occlusion cost.
%       disp_thresh - scalar disparity threshold for smoothness prior.
%       smoothness_kernel - index denoting which smoothness kernel to use.
%                           1: truncated linear; 2: truncated quadratic.
%       lambda_l - scalar smoothness prior weight for cliques crossing
%                  segmentation boundaries.
%       lambda_h - scalar smoothness prior weight for cliques not crossing
%                  segmentation boundaries.
%       seg_params - 1x3 vector of parameters for the mean-shift
%                    over-segmentation of the reference image.
%       visibility - boolean indicating whether to employ the geometrical
%                    visbility contstraint.
%       connect - scalar neighbourhood system of the graph, 4 or 8
%                 connected.
%       max_iters - scalar number of iterations to halt after, if
%                   convergence is not achieved first.
%       converge - scalar percentage decrease in energy per iteration at
%                  which optimization stops.
%       average_over - scalar number of iterations to average over when
%                      checking convergence.
%       contract - scalar number of QPBOP iterations to do.
%       improve - scalar indicating which method to use to label unlabelled
%                 nodes. 0: QPBO-F, 1: QPBOI-F, 2: QPBO-R, 3: QPBO-L,
%                 4: QPBOI-R.
%       independent - boolean indicating whether to use independent, or
%                     merely strongly-connected, regions for improve
%                     methods 2 & 4.
%
%OUT:
%   D - HxW disparity map.
%   info - structure containing other outputs from the algorithm.

% $Id: ojw_stereo.m,v 1.5 2008/11/17 11:27:35 ojw Exp $

% Crude check for a reference image
if max(abs(P([1:6 9])-[1 0 0 0 1 0 1])) > 1e-12
    error('First image must be reference image');
end

% Initialize data arrays
R = images{1}(round(P(8))+(1:sz(1)),round(P(7))+(1:sz(2)),:);
vals.I = images(2:end);
vals.P = permute(P(:,:,2:end), [2 1 3]);
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
    info.params.disp_thresh = options.disp_thresh;
    info.params.col_thresh = options.col_thresh;
    info.params.occl_const = options.occl_const;
    info.params.lambda_l = options.lambda_l;
    info.params.lambda_h = options.lambda_h;
end
sigmaColor = 10;
if isnumeric(options.proposal_method) && size(options.proposal_method, 1) == 1
    % Use the proposal methods:
    
                % SegPln (prototypical segment-based stereo proposals)
%                 [Dproposals info.segpln_gen] = ojw_segpln(images, P, disps, R, options);
% 				tmp = info.segpln_gen;
%                 save('Dproposals','Dproposals');
% 				save('segpln_gen','tmp');
% 				clear tmp
                %   load segpln_gen & Dproposals
                tmp = load('segpln_gen');
                info.segpln_gen = tmp.tmp;
                Dproposals = load('Dproposals');
                Dproposals = Dproposals.Dproposals;
                clear tmp
                
                % object depth segmentation
				info.segpln_Obj = cell(size(info.segpln_gen.segments,3),1);
                [X Y] = meshgrid(1:sz(2),1:sz(1));
                plane = info.segpln_gen.plane;
                % load segpln_Obj
                tmp = load('segplnRight7_Obj');
                info.segpln_Obj = tmp.tmp;
                clear tmp
                %
                FsegNum = 0;
                OsegNum = 0;
				for i = 1:size(Dproposals,3)
              
                    % try gco
%                     segment = info.segpln_gen.segments(:,:,i);
%                     info.segpln_Obj{i}.F = segment;
%                     segNum = max(max(segment));
%                     h = GCO_Create(segNum,segNum);
%                     % set data term
% 
%                     data = zeros(segNum,segNum);
%                    
% 
%                     
%                     for sid = 1:segNum %row
%                         planeD = -(X * plane{i}(sid,1) + Y * plane{i}(sid,2) + plane{i}(sid,3));
%                         planeD(planeD<vals.d_min | planeD>vals.d_min+vals.d_step) = 2*vals.d_step;
%                         % planeDiff: 1*pix
%                         planeDiff = reshape(abs(planeD-Dproposals(:,:,i))/vals.d_step,[],1);
%                         data(sid,:) = accumarray(reshape(segment,[],1),planeDiff)';
%                         %data(sid,:) = accumarray(reshape(segment,[],1),planeDiff)';
%                     end
%                     %data = data/max(max(data));
%                     clear planeD planeDiff mapp
%                     tmp = int32(segment(vals.SEI));
%                     % depth diff on boundary
%                     ColorD = double(images{1})/255.;
%                     ColorD = reshape(ColorD,[],3);
%                     DiffObjIdx = repmat(any(diff(int32(segment(vals.SEI))),1), [2 1]);
%                     Neigh = reshape(tmp(DiffObjIdx),2,[]);
% %                     % depth difference on smooth term
% %                     tmpD = Dproposals(:,:,i);
% %                     tmpD = tmpD(vals.SEI);
% %                     tmpD = abs(diff(reshape(tmpD(DiffObjIdx),2,[])));
% %                     tmpD = repmat(tmpD,[1 2]);
% %                     % color difference on smooth term
% %                     Idx = reshape(vals.SEI(DiffObjIdx),2,[]);
% %                     ColorD = 5/exp(sqrt(sum((ColorD(Idx(1,:),:) - ColorD(Idx(2,:),:)).^2,2)))+1;
% %                     ColorD = repmat(ColorD,[1 2]);
%                     clear DiffObjIdx Idx
%                     pixObj = histc(reshape(segment,1,[]),1:segNum)';
%                     List = (Neigh(1,:)-1)*int32(segNum)+Neigh(2,:);
%                     List = [List (Neigh(2,:)-1)*int32(segNum)+Neigh(1,:)];
%                     % boundary length
%                     Smooth = zeros(1,segNum*segNum);
%                     Smooth(List) = [pixObj(Neigh(1,:))+pixObj(Neigh(2,:)) pixObj(Neigh(1,:))+pixObj(Neigh(2,:))];
% %                     SmoothNormalize = [pixOnBoundary' ones(segNum,1)]*[ones(1,segNum) ; pixOnBoundary];
% %                     SmoothNormalize = SmoothNormalize.^2;
% %                     SmoothNormalize(logical(eye(size(SmoothNormalize)))) = 1;
% %                     tmpD = accumarray(List',tmpD');
% %                     tmpD = [tmpD ; zeros(segNum*segNum-size(tmpD,1),1)];
% %                     ColorD = accumarray(List',ColorD');
% %                     ColorD = [ColorD ; zeros(segNum*segNum-size(ColorD,1),1)];
% %                     ColorD = ColorD';
%  					Smooth = histc(List,1:segNum*segNum);
%                     Smooth = Smooth*0.5;
% %                      Smooth = (ColorD .* tmpD)./(Smooth.^2)' ; COLOR& DetphDIFF
% %                     Smooth = ColorD./(Smooth)' ;
% 					Smooth = reshape(Smooth,[segNum segNum]);
%                     Smooth(isnan(Smooth))=0;
%                     %Smooth = Smooth * 0.5 * ;
%                     %Smooth = reshape(Smooth,[segNum segNum])./SmoothNormalize;
%                     %Smooth = Smooth.*ColorD;
%                     %Smooth = max(max(Smooth))-Smooth;
%                     %Smooth(logical(eye(size(Smooth)))) = 0; % avoid NaN
% %                     Smooth = histc(List,1:segNum*segNum);
% % 					Smooth = reshape(Smooth,[segNum segNum])./SmoothNormalize;
%                     clear tmp Neigh List ColorD SmoothNormalize
%                     GCO_SetDataCost(h,data);
% 					GCO_SetNeighbors(h,Smooth);
%                     GCO_Expansion(h);
%                     Label = GCO_GetLabeling(h);
%                     % after relabel
%                     uq = unique(Label);
%                     NewSegNum = numel(uq);
%                     segment = Label(segment);
%                     tmpL = 1:segNum;
%                     tmpL(uq) = 1:NewSegNum;
%                     segment = tmpL(segment);
%                     info.segpln_Obj{i}.segMap = segment;
%                     info.segpln_Obj{i}.plane = plane{i}(uq,:);
%                     save(['OnPaperMap[' num2str(i) ']'],'segment');
%                     clear segment tmpL
%                     GCO_Delete(h);
                    % symmetric - on other view
%                   % warp to other view
                    info.segpln_Obj{i}.F = info.segpln_Obj{i}.F + FsegNum;
                    info.segpln_Obj{i}.segMap = info.segpln_Obj{i}.segMap + OsegNum;
                    FsegNum = max(max(info.segpln_Obj{i}.F));
                    OsegNum = max(max(info.segpln_Obj{i}.segMap));
                    info.segpln_Obj{i}.otherView = ObjWarp(Dproposals(:,:,i),info.segpln_Obj{i}.F,info.segpln_Obj{i}.segMap,P(:,:,2:end));
                    
                    %calc Model
                    
%                     [info.segpln_Obj{i}.plane info.segpln_Obj{i}.parallax info.segpln_Obj{i}.GMM_Name] = calcObjModel(vals.d_step,Dproposals(:,:,i),info.segpln_Obj{i},images{1},i,0);
                end
                
				tmp = info.segpln_Obj;
				save('segplnRight7_Obj','tmp');
                clear R tmp

                
                Dproposals = Dproposals(:,:,1:9);
                info.segpln_Obj = info.segpln_Obj(1:9);
 
                [D info.segpln_optim ObjModel] = obj_fuse_proposals(vals, Dproposals, options, info.segpln_Obj,0);
                delete('Final_info.mat');
                delete('Final_D.mat');
                delete('Final_Model.mat');
                save('Final_D','D');
                save('Final_info','info');
                save('Final_Model','ObjModel');
                tmp = load('Final_info');
                info = tmp.info;
                tmp = load('Final_Model');
                ObjModel = tmp.ObjModel;
                tmp = load('Final_D');
                D = tmp.D;
                clear tmp;
                [ExpanProposals ObjExpansion] = ExpansionProposals(D,ObjModel,P);
                [D info.fuseExp ObjModel] = obj_fuse_proposals(vals, ExpanProposals, options, ObjExpansion,1);
                Propagate2otherFrames(D,ObjModel, images,P(:,:,2:end))
                figure(13);imshow(D);
%                 clear Dproposals
   
    
else
    sprintf('dont go here');
end
return
end

function [ExpanProposals ObjExpansion] = ExpansionProposals(D,ObjModel,P)
ObjIds = unique(ObjModel.segMap);
ObjNum = numel(ObjIds);
pixNumInSeg = histc(ObjModel.F(:),1:max(max(ObjModel.F)));
ExpanProposals = D;
ObjExpansion{1} = ObjModel;
sz = size(D);
[X Y] = meshgrid(1:sz(2),1:sz(1));
WC = ones(sz(1)*sz(2), 3);
WC(:,3) = 1./D(:);
WC(:,1) = X(:).* WC(:,3);
WC(:,2) = Y(:).* WC(:,3);
iter = 1;
rt = 0.1;
for i = 1:ObjNum
    s = ObjIds(i);
    M = ObjModel.segMap == s;
    depthPlaneId = unique(ObjModel.F(M));
    for d = 1:numel(depthPlaneId)
        if pixNumInSeg(depthPlaneId(d)) > 500
            iter = iter+1;
            N = WC(ObjModel.F == depthPlaneId(d),:);
            if size(N, 1) > 3
                % Ransac to weed out outliers
                M_ = rplane(N, rt);
                N = N(M_,:);
            end
            % Find least squares plane from inliers
            N = N \ repmat(-1, [size(N, 1) 1]);
            % DisparityMap
            ExpanProposals(:,:,iter) = -(X * N(1) + Y * N(2) + N(3));
            tmpObj.segMap = ones(size(ObjModel.segMap),'uint32')*s;
            tmpObj.F = ones(size(ObjModel.segMap),'uint32')*depthPlaneId(d);
            tmpObj.GMM_Name{1} = ObjModel.GMM_Name{i};
            tmpObj.parallax = ObjModel.parallax(i,:);
            tmpObj.otherView = ObjWarp(ExpanProposals(:,:,iter),tmpObj.F,tmpObj.segMap,P(:,:,2:end));
            tmpObj.plane = ObjModel.plane(i,:);
            ObjExpansion{iter} = tmpObj;
        end
    end
end


return 
end


% LO-RANSAC functions
function inls = rplane(pts, th)

MAX_SAM = 500;
conf = .95;

len = size(pts, 1);
max_i = 3;
max_sam = MAX_SAM;
no_sam = 0;
div = repmat(-1, [3 1]);
inls = false(len, 1);

while no_sam < max_sam
    no_sam = no_sam + 1;
    sam = randperm(len);
    sam = sam(1:3);

    %%% compute a distance of all points to a plane given by
    %%% pts(:,sam) to dist
    N = pts(sam,:) \ div;
    dist = abs((pts * N) + 1);  % -1 move to 
    v = dist < th;
    no_i = sum(v);

    if max_i < no_i
        % Re-estimate plane and inliers
        N = pts(v,:) \ repmat(-1, [no_i 1]);
        dist = abs((pts * N) + 1);
        v = dist < th;
        
        if sum(v) > sum(inls)
            inls = v;
            max_i = no_i;
            max_sam = min([max_sam,nsamples(sum(inls), len, 3, conf)]);
        end
    end
end
return
end
function SampleCnt = nsamples(ni, ptNum, pf, conf)
q  = prod ([(ni-pf+1) : ni] ./ [(ptNum-pf+1) : ptNum]);

if (1 -q) < eps
    SampleCnt = 1;
else
    SampleCnt  = log(1 - conf) / log(1 - q);
end

if SampleCnt < 1
    SampleCnt = 1;
end
return
end