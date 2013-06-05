function Proposals = calcModel( savePath, disps, Key, Img ,id)
%CALCOBJMODEL Summary of this function goes here
%   Detailed explanation goes here
warning_state = warning('query', 'all');
warning off MATLAB:divideByZero
warning off MATLAB:singularMatrix
warning off MATLAB:nearlySingularMatrix
warning off MATLAB:illConditionedMatrix
warning off MATLAB:rankDeficientMatrix


sz = size(Key.F);


Proposals = cell(sz(3),1);
[X Y] = meshgrid(1:sz(2), 1:sz(1));
WC = zeros(sz(2)*sz(1), 3);
WC(:,2) = WC(:,3) .* Y(:);
WC(:,1) = WC(:,3) .* X(:);
Img = single(reshape(Img,[],3));
tmp = Img(:,1);
Img(:,1) = Img(:,3);
Img(:,3) = tmp;
for k=1:sz(3)
    D = Key.D(:,:,k);
	WC(:,3) = 1 ./ D(:);
	rt = 0.1;
    ObjIds = unique(Key.segMap(:,:,k));
    segNum = numel(ObjIds);
    Proposals{k}.segMap = Key.segMap(:,:,k);
    Proposals{k}.F = Key.F(:,:,k);
    Proposals{k}.D = Key.D(:,:,k);
    Proposals{k}.parallax = zeros(segNum,ceil((2*numel(disps)+1)/2));
    Proposals{k}.GMM_Name = cell(segNum,1);
    Proposals{k}.ObjPln = FitPlane(Key.D(:,:,k),Key.segMap(:,:,k));
	for i=1:segNum
		% Choose a segment
		s = ObjIds(i);
		M = Key.segMap(:,:,k) == s;
		findM = find(M);
		[Y X] = ind2sub(sz(1:2), findM);
        N = Proposals{k}.ObjPln(i,:);
		OM = -(X * N(1) + Y * N(2) + N(3));
		parallax(i,:) = hist(D(M)-OM(:),-disps(1):2*disps(1)/numel(disps):disps(1))/numel(OM);
		
        GMM_Name{i} = [savePath num2str(k) '_' num2str(i) '.yml'];
		
		Pix = Img(reshape(M,[],1),:);
	   
		TrainGMM(Pix',GMM_Name{i});      % element column major
		
	end
end
% Reset warnings
warning(warning_state);
end
