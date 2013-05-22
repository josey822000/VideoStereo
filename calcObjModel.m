function [ plane parallax GMM_Name] = calcObjModel( disps, F, segpln_Obj, Img ,id, ExpFuse)
%CALCOBJMODEL Summary of this function goes here
%   Detailed explanation goes here
warning_state = warning('query', 'all');
warning off MATLAB:divideByZero
warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix
warning off MATLAB:illConditionedMatrix
warning off MATLAB:rankDeficientMatrix

ObjIds = unique(segpln_Obj.segMap);
segNum = numel(ObjIds);
sz = size(F);
plane = zeros(segNum,3);
parallax = zeros(segNum,ceil((2*numel(disps)+1)/2));
GMM_Name = cell(segNum,1);
% World coordinates for plane fitting
% Kf = P.K(:,:,1);
% Rf = P.R(:,:,1);
% Tf = P.T(:,1);
[X Y] = meshgrid(1:sz(2), 1:sz(1));
WC = zeros(sz(2)*sz(1), 3);
WC(:,3) = 1 ./ F(:);
WC(:,2) = WC(:,3) .* Y(:);
WC(:,1) = WC(:,3) .* X(:);
% WC(:,:) = ( Kf \ WC(:,:)' )';
% tmp_vec = repmat(Tf, [1 prod(sz)]);
% WC(:,:) = (Rf * WC(:,:)' + tmp_vec)';

rt = 0.1;
Img = single(reshape(Img,[],3));
dataCnt = 0;
for i=1:segNum
	% Choose a segment
    s = ObjIds(i);
	M = segpln_Obj.segMap == s;
	N = WC(M,:);
	if size(N, 1) > 3
		% Ransac to weed out outliers
		M_ = rplane(N, rt);
		N = N(M_,:);
    else
        continue;
    end
    dataCnt = dataCnt + 1;
	% Find least squares plane from inliers
	N = N \ repmat(-1, [size(N, 1) 1]);
	plane(dataCnt,:) = N;
	findM = find(M);
	[Y X] = ind2sub(sz, findM);
%     N = [ (N' * Rf(:, 1))  (N' * Rf(:, 2))  (N' * Rf(:, 3)) ] / ( N' * Tf + 1) ; 
%     pts = [X';Y';ones(size(X))'];
%     pts = Kf\ pts;
%     
%     OM = -N*pts;
    OM = -(X * N(1) + Y * N(2) + N(3));
	parallax(dataCnt,:) = hist(F(findM)-OM(:),-disps(1):2*disps(1)/numel(disps):disps(1))/numel(OM);
    switch ExpFuse
        case 2
            GMM_Name{dataCnt} = ['GMM\ExpFuse_' num2str(id) '_' num2str(s) '.yml'];
        case 1
            GMM_Name{dataCnt} = ['GMM\fuse_' num2str(id) '_' num2str(s) '.yml'];
        case 0
            GMM_Name{dataCnt} = ['GMM\' num2str(s) '.yml'];
    end
    Pix = Img(reshape(M,[],1),:);
   
    TrainGMM(Pix',GMM_Name{dataCnt});      % element column major
    
end
% Reset warnings
warning(warning_state);
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