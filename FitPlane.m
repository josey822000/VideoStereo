function plane = FitPlane( Map, D)
%CALCOBJMODEL Summary of this function goes here
%   Detailed explanation goes here
warning_state = warning('query', 'all');
warning off MATLAB:divideByZero
warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix
warning off MATLAB:illConditionedMatrix
warning off MATLAB:rankDeficientMatrix

ObjIds = unique(Map);
segNum = numel(ObjIds);
sz = size(Map);
plane = zeros(segNum,3);

[X Y] = meshgrid(1:sz(2), 1:sz(1));
WC = zeros(sz(2)*sz(1), 3);

WC(:,2) = WC(:,3) .* Y(:);
WC(:,1) = WC(:,3) .* X(:);

WC(:,3) = 1 ./ D(:);
rt = 0.1;
dataCnt = 0;
for i=1:segNum
    % Choose a segment
    s = ObjIds(i);
    M = Map == s;
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