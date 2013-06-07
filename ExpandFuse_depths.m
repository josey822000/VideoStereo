function [N info energy V] = ExpandFuse_depths(nowModel,newModel, vals, options)


% $Id: ibr_fuse_depths.m,v 1.3 2008/11/17 11:27:35 ojw Exp $
t_start = cputime;
figure(2);imshow(nowModel.D/vals.disps(1));
figure(3);imshow(newModel.D/vals.disps(1));
% figure(4);sc(nowModel.segMap,'rand');
% figure(5);sc(newModel.segMap,'rand');

% Initialize values
Kinf = int32(0);
occl_cost = cast(1000000,class(Kinf));
Kinf = occl_cost + 1; % Smaller value seems to avoid errors in QPBO
info.timings = zeros(3, 1, numel(vals.improve));
info.numbers = zeros(4, 1, numel(vals.improve), 'uint32');
energy = zeros(1, numel(vals.improve));
%vals.I = images(2:end);
num_in = numel(vals.I);
% image sz
sp = size(nowModel.D);
% pixel num
tp = numel(nowModel.D);
planar = size(vals.SEI, 1) == 3;
oobv = cast(25, class(vals.R));
out_unlabel = vals.improve(end) < 0;


seg1 = nowModel.segMap;
seg2 = newModel.segMap;
MDLaux = numel(unique([seg1 seg2]));
% Initialise arrays for the data terms
if vals.visibility
    EI = reshape(repmat(uint32(tp+(0:num_in-1)*2*tp), [2*tp 1]), 1, []);
    EI = [repmat(uint32(1:tp), [1 2*num_in]); EI+repmat(uint32(1:2*tp), [1 num_in])];
    E = repmat(cat(3, [occl_cost 0 0 0], [0 0 occl_cost 0]), [tp 1 1 num_in]);
else
    % Data edges are not needed
    EI = zeros(2, 0, 'uint32');
    E = zeros(4, 0, class(Kinf));
end


% unary
% U = zeros(tp+MDLaux, 2, class(Kinf));
U = zeros(tp, 2, class(Kinf));
TEI = zeros(2, 0, 'uint32');
TE = zeros(4, 0, class(Kinf));
% lambda
lambdaOcc = 25;
lambdaOcoh = 10;
lambdaDcoh = 10;
lambdaColor = 0;
lambdaParallax = 2;
lambdaMDL = 0;
lambdaInf = 1000000;
% For each input image...
sz = size(nowModel.D);

[X Y] = meshgrid(0:sz(2)-1, 0:sz(1)-1);
pts = ones(sz(1)*sz(2), 3);
pts(:,1) = X(:);
pts(:,2) = Y(:);
P = vals.P;
Kf = P.K(:,:,1);
Rf = P.R(:,:,1);
Tf = P.T(:,1);
for a = 1:num_in
    % Calculate the coordinates in the input image
    K2 = P.K(:,:,a+1);
    R2 = P.R(:,:,a+1);
    T2 = P.T(:,a+1);
    tmp_mat = K2 * R2' * Rf;
    tmp_vec = K2 * R2' * ( Tf - T2);
    tmp_vec = repmat(tmp_vec, [1 prod(sz)]);
    epl_pts1 = ( tmp_mat * (Kf\(pts(:,:)')) + repmat(nowModel.D(:)',[3 1]) .* tmp_vec )';
    epl_pts1 = epl_pts1 ./ repmat(epl_pts1(:,3),[1 3]);
    epl_pts2 = ( tmp_mat * (Kf\(pts(:,:)')) + repmat(newModel.D(:)',[3 1]) .* tmp_vec )';
    epl_pts2 = epl_pts2 ./ repmat(epl_pts2(:,3),[1 3]);
    epl_pts1(:,1:2) = epl_pts1(:,1:2)+1;
    epl_pts2(:,1:2) = epl_pts2(:,1:2)+1;
    % Epc
    % D1,D2
    [mpD1 mpD2 Check1 Check2 outB1 outB2] = warp2DepthConsiderNoise(nowModel.D,newModel.D,nowModel,newModel,nowModel.otherView{a},newModel.otherView{a},epl_pts1,epl_pts2);
%     figure(6);
%     imshow(Check1);
%     figure(7);
%     imshow(Check2);

    Check = [reshape(Check1,[],1); reshape(Check2,[],1)] ;
    
    Dp = [reshape(nowModel.D,[],1) ; reshape(newModel.D,[],1)];
    
    Dmp = [reshape(mpD1,[],1) ; reshape(mpD2,[],1)];
    nCheck = ~Check;
    occ = Dp(nCheck)<=Dmp(nCheck);
    nCheck(nCheck>0)=occ;
    occ = nCheck | [outB1(:) ; outB2(:)];
    epl_pts = [epl_pts1 ; epl_pts2];
    clear epl_pts1 epl_pts2
    M = vgg_interp2(vals.I{a}, epl_pts(:,1), epl_pts(:,2), 'linear', oobv);
    M = squeeze(M) - vals.R;
    M = vals.ephoto(M);
    M(M>lambdaOcc-1) = lambdaOcc-1;
    CheckDisRange = [reshape(nowModel.D>=0 & nowModel.D<=vals.d_step,[],1); reshape(newModel.D>=0 & newModel.D<=vals.d_step,[],1)];
    occ = occ & CheckDisRange;
    Check = Check & CheckDisRange;
%     figure(8);imshow(reshape(occ(1:end/2).*lambdaOcc+Check1(:).*M(1:end/2),size(Check1))/max(max(M(1:end/2))));
%     figure(9);imshow(reshape(occ(end/2+1:end).*lambdaOcc+Check2(:).*M(end/2+1:end),size(Check1))/max(max(M(end/2+1:end))));
    Epc = reshape((~Check.*lambdaOcc + Check.*M), tp, 2);
    figure(10);
    imshow(reshape(~(occ(1:end/2)|Check(1:end/2)),size(Check1)));
    figure(11);
    imshow(reshape(~(occ(end/2+1:end)|Check(end/2+1:end)),size(Check1)));
    clear M N Check V1 V2 occ ObjCheck DepCheck CheckDisRange
    % Obj Color term
    Ecol = reshape(lambdaColor*-GMM(vals.R,nowModel.segMap,newModel.segMap,nowModel.GMM_Name,newModel.GMM_Name,nowModel.table,newModel.table),tp,2);
    
    % Obj parallax
    Epar = reshape(lambdaParallax*(1-ParallaxModel(nowModel.D,newModel.D,nowModel,newModel,vals.d_step,vals.ndisps,nowModel.table,newModel.table,P)),tp,2);
    IA = cast((Epc) ,class(Kinf));
    clear Epc Ecol Epar


    % Find interactions
    epl_pts(:,3) = epl_pts(:,3)./[nowModel.D(:); newModel.D(:)];
    [T M] = sortrows(epl_pts);
    N = find_interactions(T, 0.5); % Optimized version
    N = M(N); % Unsort
    N = uint32(N(:,abs(diff(N))~=tp)); % Remove interactions between the same node

    % Add the pixel interactions to the graph
    M = N(1,:) > tp;
    
    %N(2,N(2,:)>tp) = N(2,N(2,:)>tp) + MDLaux;
%     TEI = [TEI [N(1,:)-uint32(tp*M); uint32(tp*2*a-tp)+MDLaux+N(2,:)]];
    TEI = [TEI [N(1,:)-uint32(tp*M); uint32(tp*2*a-tp)+N(2,:)]];
    T = zeros(4, numel(M));
    T(2,~M) = Kinf;
    T(4,M) = Kinf;
    TE = [TE T];

    if vals.visibility
        % not used
        % Set up photoconsistency edges
        E(:,2,1,a) = IA(:,1);
        E(:,4,2,a) = IA(:,2);
        
        if vals.compress_graph
            % Determine the photoconsistency nodes which have no interactions
            M = ones(tp, 2, class(Kinf));
            M(N(2,:)) = 0;

            % Add those photoconsistency terms to the unaries
            U = U + M .* IA;
        end
    else
        % Add the photoconsistency terms to the unaries
        U = U + IA;
    end
    clear T M N
end
clear WC IA

EI_ = EI;
if vals.visibility
    % not used this part cause in our code visibility is false
    E = reshape(permute(E, [2 1 3 4]), 4, []);
    if vals.compress_graph
        % The unary and pairwise energies as they stand are entirely
        % correct, i.e. will give the correct labelling. However, it can be
        % compressed into a smaller but equivalent graph, which will be
        % faster to solve, by removing superfluous nodes and edges.
        [U EI EI_ E N T] = compress_graph(U', EI, E, TEI, TE, tp, num_in);
    else
        %
        U = zeros(2, tp+tp*2*num_in, class(Kinf));
        T = TE;
        N = TEI;
    end
    % Concatenate data and visibility edges
    E = [E T];
    EI_ = [EI_ N];
    clear T N
else
    U = U';
end
TE = TE(2,:) ~= 0;
info.timings(1,:) = cputime - t_start; % Time data term evaluation
% ====================%
% smooth term
% ====================%
% Map = logical([ones(1,size(vals.SEI,2)); zeros(2,size(vals.SEI,2)); ones(1,size(vals.SEI,2))])' ; 
% Map = permute(Map,[3 2 1]);
% Eoc
Oc = [reshape(nowModel.segMap,1,[]) ; reshape(newModel.segMap,1,[])];
Oc = reshape(Oc(:,vals.SEI),4,[]);
Oc = (diff(int32(reshape(Oc([1 3; 1 4; 2 3; 2 4]',:), 2, 4, [])))~=0);
% Oc = lambdaOcoh*Oc + 1.*~Oc;
Oc = lambdaOcoh*Oc;
% Oc = Oc .* Map;
Eoc = reshape(cast(Oc, class(E)), [], size(vals.SEI, 2)); 
% Edc
Dc = [reshape(nowModel.F,1,[]) ; reshape(newModel.F,1,[])];
Dc = reshape(Dc(:,vals.SEI),4,[]);
Dc = diff(int32(reshape(Dc([1 3; 1 4; 2 3; 2 4]',:), 2, 4, [])))~=0;
EqObjNeqDep = ~Oc & Dc;
% EqObjNeqDep = EqObjNeqDep & Map;

calcDeptDis = [reshape(nowModel.D,1,[]) ; reshape(newModel.D,1,[])];
calcDeptDis = reshape(calcDeptDis(:,vals.SEI),4,[]);
calcDeptDis = abs(diff(reshape(calcDeptDis([1 3; 1 4; 2 3; 2 4]',:), 2, 4, [])));
% Edc = (lambdaDcoh/2)*(calcDeptDis<=1 & EqObjNeqDep) + lambdaDcoh*(calcDeptDis>1 & EqObjNeqDep) + 1.*~EqObjNeqDep;
Edc = (lambdaDcoh/2)*(calcDeptDis<=1 & EqObjNeqDep) + lambdaDcoh*(calcDeptDis>1 & EqObjNeqDep);
Edc = reshape(cast(Edc, class(E)), [], size(vals.SEI, 2));
SE = Eoc + Edc;
% SE = ones(4,size(vals.SEI,2),'int32');
clear Eoc Edc Map Oc Dc calcDeptDis 
%

if ~planar
    E = [E SE];
    EI_ = [EI_ vals.SEI];
end

if vals.visibility
    [tmpE tmpEI tmpU] = GenClique(seg1,seg2,tp,num_in,MDLaux,lambdaMDL);
else
    [tmpE tmpEI tmpU] = GenClique(seg1,seg2,tp,0,MDLaux,lambdaMDL);
end
E = [E tmpE];
EI_ = [EI_ tmpEI];
U = [U zeros(2,MDLaux,'int32')]+tmpU;

clear tmpE tmpEI tmpU
info.timings(2,:) = cputime - t_start; % Time smoothness term evaluation
figure(1);
for a = 1:numel(vals.improve)
    t_start = cputime;
    % Fuse the two labellings, using contract and/or improve if desired
%     qpbo_params = int32([tp+MDLaux ((vals.improve(a)==1)+(vals.improve(a)==4)*2) vals.contract(a) vals.contract(a)>0]);
    qpbo_params = int32([tp+MDLaux ((vals.improve(a)==1)+(vals.improve(a)==4)*2) vals.contract(a) vals.contract(a)>0]);
    if vals.improve(a) == 4
        % Add callback function handle
        qpbo_params = {qpbo_params, @(L) (choose_labels(L, U, E, EI, SE, vals.SEI, TE, TEI, num_in, vals.visibility, 2, vals.independent) > 0)};
    end
    try
        if planar
            [M stats] = vgg_qpbo(U, EI_, E, vals.SEI, SE, qpbo_params);
        else
            [M stats] = vgg_qpbo(U, EI_, E, qpbo_params);
        end
    catch
        % Error probably due to probe failure
        stats = [0 0 Inf];
        M = false(tp, 1);
    end
    clear qpbo_params
    info.numbers(2:4,a) = stats;

    if stats(1) && vals.improve(a) >= 2 && vals.improve(a) <= 3
        if nargout > 2 || vals.show_output
            [M info.numbers(3,a) U_ E_ SE_ V] = choose_labels(M, U, E(:,1:size(EI, 2)), EI, SE, vals.SEI, TE, TEI, num_in, vals.visibility, vals.improve(a), vals.independent);
            energy(a) = sum(U_) + sum(E_) + sum(SE_);
        else
            [M info.numbers(3,a)] = choose_labels(M, U, E(:,1:size(EI, 2)), EI, SE, vals.SEI, TE, TEI, num_in, vals.visibility, vals.improve(a), vals.independent);
        end
        N = M > 0;
    elseif nargout > 2 || vals.show_output
        N = M > 0;
        [U_ E_ SE_ V] = calc_vis_energy(N, U, E(:,1:size(EI, 2)), EI, SE, vals.SEI, TE, TEI, num_in, MDLaux);
        
        M = M(1:end-MDLaux);
        energy(a) = sum(U_) + sum(E_) + sum(SE_);
    end
    info.numbers(1,a) = sum(N);
    info.timings(3,a) = cputime - t_start + info.timings(2,a); % Time optimization
end
clear TEI TE E SE EI_

if nargout > 3 || vals.show_output
    % Generate output visibilities
    T = (tp * N(1:end-MDLaux)) + (1:tp)';
    for b = 1:num_in
        V(1:tp,b) = V(T);
        T = T + 2*tp;
    end
    V(tp+1:end,:) = [];
end

if vals.show_output
    % Display the output figures
    N = N(1:end-MDLaux);
    MDLU = U_(end-MDLaux+1:end);
    U_ = U_(1:end-MDLaux);
    U_ = double(U_) + accum(EI(1,:)', E_, [tp 1]);
    U_ = reshape(U_, sp(1), sp(2));
    if vals.visibility
        % Take off the occlusion costs and normalize
        EI = reshape(sum(V(1:tp,:), 2), sp(1), sp(2));
        U_ = U_ - (num_in - EI) * double(occl_cost);
        E_ = EI ~= 0;
        U_(E_) = U_(E_) ./ EI(E_);
    end
    
    set(0, 'CurrentFigure', vals.show_output);
    subplot('Position', [1/3 0.5 1/3 0.5]);
    nowModel.D(N) = newModel.D(N);
    sc(nowModel.D, 'contrast', vals.d_min+[0 vals.d_step]);
    subplot('Position', [2/3 0.5 1/3 0.5]);
    disp(['U:' num2str(sum(sum(U_)))]);    
    disp(['MDLU:' num2str(sum(MDLU))]);
    sc(U_, 'jet');
    subplot('Position', [0 0 1/3 0.5]);
    T = reshape(sc(reshape(M, sp(1), sp(2)), 'prism'), [], 3);
    I = M < 0;
    T(I,:) = 1 - (1 - T(I,:)) * 0.3; % Lighten unlabelled pixels set to 0
    I = M > 1;
    T(I,:) = T(I,:) .* 0.3; % Darken unlabelled pixels set to 1 by optimal splice
    T(M==0,:) = 1;
    T(M==1,:) = 0;
    sc(reshape(T, sp(1), sp(2), 3), [0 1]);
    subplot('Position', [1/3 0 1/3 0.5]);
    sc(reshape(sum(V, 2), sp(1), sp(2)), [0 num_in], 'contrast');
    subplot('Position', [2/3 0 1/3 0.5]);
    U_ = -accum(vals.SEI(2,:)', SE_, [tp 1]);
    disp(['E:' num2str(sum(sum(U_)))]);    
    sc(reshape(U_, sp(1), sp(2)));
    
    figure(12);
    subplot('Position', [1/3 0.5 1/3 0.5]);
    U = U(:,1:end-MDLaux);
    U_ = U((0:tp-1)'*2+1);
    U_ = reshape(U_, sp(1), sp(2));
    sc(U_, 'jet');
    
    subplot('Position', [2/3 0.5 1/3 0.5]);
    U_ = U((0:tp-1)'*2+2);
    U_ = reshape(U_, sp(1), sp(2));
    sc(U_, 'jet');
    drawnow;
end
if out_unlabel
    N = M; % Output unlabelled pixels
end
return

function [M num_regions U_ E_ SE_ V] = choose_labels(M, U, E, EI, SE, SEI, TE, TEI, num_in, visibility, improve, independent)
% Calculate visibilities and regions assuming unlabelled pixels are set
% to 0 then 1.
[U_ E_ SE_ V] = calc_vis_energy(M==1, U, E, EI, SE, SEI, TE, TEI, num_in);
[U2 E2 SE2 V2] = calc_vis_energy(M~=0, U, E, EI, SE, SEI, TE, TEI, num_in);
num_regions = double(-min(M(:)));

if improve == 2
    % We want to do optimal splice.
    sz = [num_regions 1];
    % Merge strongly connected regions that are connected by smoothness
    % cliques
    SEI2 = M(SEI);
    SEI2 = SEI2(:,any(SEI2 < 0));
    SEI2 = SEI2(:,~all(SEI2 >= 0 | ojw_bsxfun(@eq, SEI2, min(SEI2))));
    while independent && ~isempty(SEI2)
        num_regions = num_regions - 1;
        T = SEI2(:,1);
        N = min(T);
        T = T(T < 0 & T ~= N);
        M(M==T(1)) = N;
        SEI2(SEI2==T(1)) = N;
        SEI2 = SEI2(:,~all(SEI2 >= 0 | ojw_bsxfun(@eq, SEI2, min(SEI2))));
    end
    if visibility
        % Merge strongly connected regions that are connected by visibility
        % cliques
        TEI2 = TEI(:,M(TEI(1,:))<0);
        tp = numel(M);
        M = repmat(M, [1+2*num_in 1]);
        M(TEI2(2,:)) = M(TEI2(1,:));
        SEI2 = M(TEI);
        SEI2 = SEI2(:,any(SEI2 < 0));
        SEI2 = SEI2(:,~all(SEI2 >= 0 | ojw_bsxfun(@eq, SEI2, min(SEI2))));
        while independent && ~isempty(SEI2)
            num_regions = num_regions - 1;
            T = SEI2(:,1);
            N = min(T);
            T = T(T < 0 & T ~= N);
            M(M==T(1)) = N;
            SEI2(SEI2==T(1)) = N;
            SEI2 = SEI2(:,~all(SEI2 >= 0 | ojw_bsxfun(@eq, SEI2, min(SEI2))));
        end
        % Go through each region and determine whether a labelling of 1 or 0
        % gives a lower energy, starting with the visibility edges
        EI2 = min(M(EI))';
        T = EI2 < 0;
        E2 = E2(T) - E_(T);
        engy = accum(-EI2(T), E2, sz);
        M = M(1:tp);
    else
        engy = zeros(sz);
    end
    % Go through each region and determine whether a labelling of 1 or 0
    % gives a lower energy
    SEI2 = M(SEI);
    T = any(SEI2 < 0);
    SEI2 = SEI2(:,T);
    SE2 = SE2(T) - SE_(T);
    T = -min(SEI2);
    engy = engy + accum(T(:), SE2, sz);
    T = M < 0;
    U2 = U2(T) - U_(T);
    engy = engy + accum(-M(T), U2, sz);
    update = false;
    for b = 1:sz(1)
        if engy(b) < 0
            M(M==-b) = b + 1;
            update = true;
        end
    end
    if update && nargout > 2
        % Recalculate the visibilities and energies
        [U_ E_ SE_ V] = calc_vis_energy(M > 0, U, E, EI, SE, SEI, TE, TEI, num_in);
    end
else
    % Just choose the label which gives the lowest energy
    if (sum(U_) + sum(E_) + sum(SE_)) >= (sum(U2) + sum(E2) + sum(SE2))
        % Update the labelling
        T = M < 0;
        M(T) = 1 - M(T);
        % Update the energy and visibility
        U_ = U2;
        E_ = E2;
        SE_ = SE2;
        V = V2;
    end
end
return

function [U E SE V] = calc_vis_energy(L, U, E, EI, SE, SEI, TE, TEI, num_in, MDLaux)
% Generate visibility maps
tp = numel(L);
V = true(2*tp, num_in);
V(TEI(2,L(TEI(1,:))~=TE')-(tp-MDLaux)) = false;

% Calculate energies
U = U((0:tp-1)'*2+L+1);
L = [L; V(:)];
E = E((0:size(EI, 2)-1)'*4+L(EI(1,:))*2+L(EI(2,:))+1);
if size(SEI, 1) == 3
    SE = SE((0:size(SEI, 2)-1)'*8+L(SEI(1,:))*4+L(SEI(2,:))*2+L(SEI(3,:))+1);
else
    SE = SE((0:size(SEI, 2)-1)'*4+L(SEI(1,:))*2+L(SEI(2,:))+1);
end
return

function [U EI EI_ E TEI TE] = compress_graph(U, EI, E, TEI, TE, tp, num_in)
% Count the number of interactions per input sample
SE = accum(TEI(2,:)', (1:size(TEI, 2))', [tp+tp*2*num_in 1], @num_first);
SE = SE(tp+1:end);

% Remove single interactions, attaching the photoconsitency edge
% directly to the interacting pixel
M = find(SE > 0);
L = SE(M);
EI(2,M) = TEI(1,L);
M = M(TE(4,L)~=0);
E(:,M) = E([2 1 4 3],M);
TEI(:,L) = [];
TE(:,L) = [];

% Remove the superfluous edges - photoconsistency edges with no
% interactions, that have already been incorporated into the unary term
M = SE ~= 0;
E = E(:,M);
EI = EI(:,M);

% Compress the node indices
M = zeros(tp+2*tp*num_in, 1, 'uint32');
SE = SE < 0;
L = sum(SE);
M([true(tp, 1); SE]) = uint32(1):uint32(L+tp);
EI_ = EI;
EI_(2,:) = M(EI(2,:));
TEI(2,:) = M(TEI(2,:));
U = [U zeros(2, L, class(U))];
return

function B = accum(I, A, sz, varargin)
% Older versions of Matlab can't accumulate integer arrays!
try
    B = accumarray(I, A, sz, varargin{:});
catch
    if isempty(I)
        B = zeros(sz);
    else
        B = accumarray(double(I), double(A), sz, varargin{:});
    end
end
return
function outval = ParallaxModel(D1,D2,ObjModel1,ObjModel2,d_step,ndisps,D1_table,D2_table,P)
    sz = size(D1);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    probMap1 = zeros(sz);
    probMap2 = zeros(sz);
    
    Kf = P.K(:,:,1);
    Rf = P.R(:,:,1);
    Tf = P.T(:,1);
    
    ObjIds = unique(ObjModel1.segMap);
    for i= 1:numel(ObjIds)
        s = ObjIds(i);
        M = ObjModel1.segMap==s;
        N = ObjModel1.ObjPln(D1_table(s),:)';
        Dtmp = -(X(M) * N(1) + Y(M) * N(2) + N(3));
        D = D1(M)-Dtmp;
        parallax = ObjModel1.parallax(D1_table(s),:);
        inBound = (D<d_step & D>-d_step);
        inBound = M(M>0) & inBound;
        M(M>0) = inBound;
        D = D(inBound);
        probMap1(M) = parallax(int32(ceil((D+d_step+d_step/ndisps)/(2*d_step/ndisps))));
    end
    ObjIds = unique(ObjModel2.segMap);
    for i= 1:numel(ObjIds)
        s = ObjIds(i);
        M = ObjModel2.segMap==s;
        N = ObjModel2.ObjPln(D2_table(s),:)';
        Dtmp = -(X(M) * N(1) + Y(M) * N(2) + N(3));
        D = D2(M)-Dtmp;
        parallax = ObjModel2.parallax(D2_table(s),:);
        inBound = (D<d_step & D>-d_step);
        inBound = M(M>0) & inBound;
        M(M>0) = inBound;
        D = D(inBound);
        probMap2(M) = parallax(int32(ceil((D+d_step+d_step/ndisps)/(2*d_step/ndisps))));
    end
%     figure(14); imshow(reshape(probMap1,size(ObjModel1.segMap)));
%     figure(15); imshow(reshape(probMap2,size(ObjModel1.segMap)));
    outval = [reshape(probMap1,[],1) ; reshape(probMap2,[],1)];
return


function [outval] = GMM(R,ObjMap1,ObjMap2,ModelName1,ModelName2,D1_table,D2_table)
    halfR = R(1:end/2,:);
    tmp = halfR(:,1);
    halfR(:,1) = halfR(:,3);
    halfR(:,3) = tmp;
    probMap1 = zeros(size(halfR,1),1);
    probMap2 = zeros(size(halfR,1),1);
    ObjIds = unique(ObjMap1);
    for i= 1:numel(ObjIds)
        s = ObjIds(i);
        samples = halfR(ObjMap1==s,:);
        prob = PredictGMM(samples',ModelName1{D1_table(s)});
        probMap1(ObjMap1==s) = prob;
    end
    ObjIds = unique(ObjMap2);
    for i= 1:numel(ObjIds)
        s = ObjIds(i);
        samples = halfR(ObjMap2==s,:);
        prob = PredictGMM(samples',ModelName2{D2_table(s)});
        probMap2(ObjMap2==s) = prob;
    end
    
    outval = [reshape(probMap1,[],1) ; reshape(probMap2,[],1)];
    middleV = mean(outval(1:end/2));
    outval(outval >= middleV) = middleV;
    outval = outval-middleV;
%     figure(12); imshow(reshape(abs(outval(1:end/2))/6,size(ObjMap2)));
%     figure(13); imshow(reshape(abs(outval(end/2+1:end))/6,size(ObjMap2)));
return
function B = num_first(A)
% Return:        A  if numel(A) == 1
%         -numel(A) otherwise
B = -numel(A);
if B == -1
    B = A;
end
return
function [tmpE tmpEI tmpU] = GenClique(seg1,seg2,tp,num_in,MDLaux,lambdaMDL)
tmpU = zeros(2,tp+tp*2*num_in+MDLaux,'int32');
tmpE = zeros(4,0,'int32');
tmpEI = zeros(2,0,'uint32');
tp2num_in = tp+2*tp*num_in;
% ObjIds = max(unique([seg1 seg2]));
ObjIds(unique([seg1 seg2])) = 1:numel(unique([seg1 seg2]));
ObjIdInSeg1 = unique(seg1);
segNum1 = numel(ObjIdInSeg1);
% ============================== %
% in seg1,        seg1    seg2
%           00  thumbda    0
%           01     0       0
%           10     0    thumbda
%           11     0       0
% ============================== %
%  seg1
for i= 1:segNum1
    s = ObjIdInSeg1(i);
    index = find(seg1==s);
    index(:,2) = tp2num_in+ObjIds(s);
    tmpEI = [tmpEI  index' ];
end
tmpE = [ones(1,tp,'int32')*lambdaMDL; zeros(3,tp,'int32')];



ObjIdInSeg2 = unique(seg2);
segNum2 = numel(ObjIdInSeg2);
% seg2
for i= 1:segNum2
    s = ObjIdInSeg2(i);
    index = find(seg2==s);
    index(:,2) = tp2num_in+ObjIds(s);
    tmpEI = [tmpEI  index' ];
end
tmpE = [tmpE [zeros(2,tp,'int32'); ones(1,tp,'int32')*lambdaMDL; zeros(1,tp,'int32')]];

tmpU(2,tp2num_in+1:tp2num_in+MDLaux)  = lambdaMDL;

return