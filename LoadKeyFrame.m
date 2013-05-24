function Key = LoadKeyFrame(  )
%LOADKEYFRAME Summary of this function goes here
%   Detailed explanation goes here
    KeyDir = 'keyframe';
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
        
        Key{i-2}.Model.depthpln = calcPlane(Key{i-2}.D,Key{i-2}.Model.F);
        ObjNum = max(Key{i-2}.Model.segMap(:));
        FNum = max(Key{i-2}.Model.F(:));
	end
	% load P
    for i=1:keyFrameNum
        table = zeros(ObjNum,1);
        table(unique(Key{i}.Model.segMap)) = 1:numel(unique(Key{i}.Model.segMap));
        Key{i}.Model.table = table;
    end
end
function [ plane ] = calcPlane( D, F)
%CALCOBJMODEL Summary of this function goes here
%   Detailed explanation goes here
warning_state = warning('query', 'all');
warning off MATLAB:divideByZero
warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix
warning off MATLAB:illConditionedMatrix
warning off MATLAB:rankDeficientMatrix

sz = size(F);
FIds = unique(F);
segNum = numel(unique(F));
plane = zeros(segNum,3);
% figure(2);imshow(D/0.0087);

[X Y] = meshgrid(1:sz(2), 1:sz(1));
WC = zeros(sz(2)*sz(1), 3);
WC(:,3) = 1 ./ D(:);
WC(:,2) = WC(:,3) .* Y(:);
WC(:,1) = WC(:,3) .* X(:);

rt = 0.1;

tmpD = zeros(size(D));
for i=1:segNum
	% Choose a segment
    s = FIds(i);
	M = F == s;
	N = WC(M,:);
	if size(N, 1) > 3
		% Ransac to weed out outliers
		M_ = rplane(N, rt);
		N = N(M_,:);
    end
    
	% Find least squares plane from inliers
	N = N \ repmat(-1, [size(N, 1) 1]);
	plane(i,:) = N;
    tmpD(M) = -(X(M) * N(1) + Y(M) * N(2) + N(3));
end
% figure(3);imshow(tmpD/0.0087);
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
