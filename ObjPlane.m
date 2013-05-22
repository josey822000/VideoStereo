function  ObjMap = ObjPlane(images, P, disps, R, options,info , Dproposal)
INFO = info.segpln_gen;
% splitTh = 20;
% for i = 1:size(INFO.segments,3)
%     segMap = INFO.segments(:,:,i);
%     segNum = max(max(segMap));
%     segDisp = zeros(segNum,1);
%     if segNum > splitTh
%     for s = 1:segNum  % calc mean disparity
%         segDisp(s) = mean(segMap(find(segMap(:,:) == s)));
%     end
%     
% end



sz = size(R);
if numel(sz) < 3
    sz(3) = 1;
end
%ephoto = @(F) log(2) - log(exp(sum((F-reshape(double(Dproposal), [], sz(3))) .^ 2, 2)*(-1/(options.col_thresh*sz(3))))+1);
ephoto = @(F) sum(F-reshape(double(Dproposal), [], sz(3)), 2);
R = uint8(R);
[X Y] = meshgrid(1:sz(2), 1:sz(1));
WC = ones(sz(1)*sz(2), 3);
WC(:,1) = X(:);
WC(:,2) = Y(:);
X = sz - 2 * options.window;
corr = zeros(X(1), X(2), numel(disps));
filt = fspecial('average', [1 1+2*options.window]);
% For each image...
for a = 1:numel(images)
    % Project the points
    X = WC * P(:,1:3,a)';
    P_ = P(:,4,a);

    % For each disparity...
    for b = 1:numel(disps)
        % Vary image coordinates according to disparity
        d = disps(b) * P_;
        Z = 1 ./ (X(:,3) + d(3));

        % Look up the colours
        Y = squeeze(vgg_interp2(images{a}, (X(:,1) + d(1)) .* Z, (X(:,2) + d(2)) .* Z, 'linear', -1000));

        % Calculate the RSSD
        Y = ephoto(Y);
        Y = conv2(filt, filt', reshape(Y, sz(1:2)), 'valid');
        corr(:,:,b) = corr(:,:,b) + Y;
    end
end
% Normalize
X = ephoto(-1000) * numel(images);
corr = (X(1) - corr) / X(1);
clear X Y Z h1 h2

% Extract highest scoring matches (winner takes all)
[info.corr corr] = max(corr, [], 3);
corr = disps(corr);
corr(info.corr<0.07) = 0;
corr = padarray(corr, options.window([1 1]), 'symmetric'); % Return to original size

% Generate image segmentations
if size(R, 3) == 1
    R = repmat(R, [1 1 3]);
end
segment_params = [1 1.5 10 100];
mults = [1:7 3 5 8 12 24 50 100];		%14 kinds
nMaps = numel(mults);
info.segments = zeros(sz(1), sz(2), nMaps, 'uint32');
for b = 1:nMaps
    sp = segment_params * mults(b);
    if b < 8
        % Segment the image using mean shift
        info.segments(:,:,b) = vgg_segment_ms(R, sp(1), sp(2), sp(3));
    else
        % Segment the image using Felzenszwalb's method
        info.segments(:,:,b) = vgg_segment_gb(R, 0, sp(4), sp(3), 1);
    end
end
clear A

% World coordinates for plane fitting
[X Y] = meshgrid(1:sz(2), 1:sz(1));
WC = zeros(sz(2)*sz(1), 3);
WC(:,3) = 1 ./ corr(:);
WC(:,2) = WC(:,3) .* Y(:);
WC(:,1) = WC(:,3) .* X(:);
clear X Y

% Switch off annoying warnings
warning_state = warning('query', 'all');
warning off MATLAB:divideByZero
warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix
warning off MATLAB:illConditionedMatrix
warning off MATLAB:rankDeficientMatrix

% Generate piecewise-planar disparity maps
D = zeros(sz(1)*sz(2), nMaps);
rt = 0.1; %2 * min(abs(diff(Z(:))));
plane = zeros(max(max(info.segments(:,:,b))),3);
for b = 1:nMaps
    for a = 1:max(max(info.segments(:,:,b)))
        % Choose a segment
        M = info.segments(:,:,b) == a;

        N = WC(M,:);
        N = N(N(:,3)~=0,:);
        if size(N, 1) > 3
            % Ransac to weed out outliers
            M_ = rplane(N, rt);
            N = N(M_,:);
        end

        % Find least squares plane from inliers
        N = N \ repmat(-1, [size(N, 1) 1]);
		plane(a,:) = N;
        [Y X] = ind2sub(sz, find(M));
        D(M,b) = -(X * N(1) + Y * N(2) + N(3));
    end
end

% Reset warnings
warning(warning_state);

D(isnan(D)) = 1e-100;
D = reshape(D, sz(1), sz(2), nMaps);
info.disp = corr;
info.disp(isinf(info.disp)) = 0;
return

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



