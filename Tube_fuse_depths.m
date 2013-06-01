function [N info energy V] = Tube_fuse_depths(nowModel, newModel, vals, options)


% $Id: ibr_fuse_depths.m,v 1.3 2008/11/17 11:27:35 ojw Exp $
t_start = cputime;
nowD = [];
newD = [];
segMap = [];
for i=1:size(nowModel,1)
    nowD = [ nowD  nowModel{i}.D];
    newD = [ newD  newModel{i}.D];
    segMap = [segMap [nowModel{i}.segMap; newModel{i}.segMap]];
end
figure(2);imshow([nowD ; newD]/vals.disps(1));
figure(3);sc(segMap,'rand');
clear nowD newD
% figure(3);imshow(D2/vals.disps(1));
% figure(4);sc(D1_Obj.segMap,'rand');
% figure(5);sc(D2_Obj.segMap,'rand');

% lambda
lambdaOcc = 25;
lambdaOcoh = 25;
lambdaDcoh = 25;
lambdaColor = 4;
lambdaParallax = 2;
lambdaMDL = 10000;
lambdaInf = 100000;
num_in= 1;
% Initialize values
Kinf = int32(0);
occl_cost = cast(1000000,class(Kinf));
Kinf = occl_cost + 1; % Smaller value seems to avoid errors in QPBO
info.timings = zeros(3, 1, numel(vals.improve));
info.numbers = zeros(4, 1, numel(vals.improve), 'uint32');
energy = zeros(1, numel(vals.improve));
%vals.I = images(2:end);
% tube sz
tubeSz = size(nowModel,1);
% image sz
sp = size(nowModel{1}.D);
% pixel num
tp = numel(nowModel{1}.D);
planar = size(vals.SEI, 1) == 3;
oobv = cast(-1000, class(vals.R{1}));
out_unlabel = vals.improve(end) < 0;

% Initialise arrays for the data terms
if vals.visibility
    EI = reshape(repmat(uint32(tp*tubeSz), [2*tp*tubeSz 1]), 1, []);
    EI = [repmat(uint32(1:tp*tubeSz), [1 2]); EI+repmat(uint32(1:2*tp*tubeSz), [1 1])];
    E = repmat(cat(3, [occl_cost 0 0 0], [0 0 occl_cost 0]), [tp*tubeSz 1 1 1]);
else
    % Data edges are not needed
    EI = zeros(2, 0, 'uint32');
    E = zeros(4, 0, class(Kinf));
end


% unary
% U = zeros(tp+MDLaux, 2, class(Kinf));
U = zeros(tp*tubeSz, 2, class(Kinf));
TEI = zeros(2, 0, 'uint32');
TE = zeros(4, 0, class(Kinf));
%

P = vals.P;
pts = vals.Pts;
% For each input image...
D1 = [];
D2 = [];
uniqueElement = zeros(2,1);
uniqueObj1 = [];
uniqueObj2 = [];
PcCost = [];
PGmm = [];
Ppar = [];
for w=1:size(newModel,1)
    D1 = [D1 nowModel{w}.D];
    D2 = [D2 newModel{w}.D];
    % ref P
    Kf = P.K(:,:,w);
    Rf = P.R(:,:,w);
    Tf = P.T(:,w);
    % match P
    if w == size(newModel,1)
        m = w-1;
    else
        m = w+1;
    end
    K2 = P.K(:,:,m);
    R2 = P.R(:,:,m);
    T2 = P.T(:,m);
    tmp_mat = K2 * R2' * Rf;
    tmp_vec = K2 * R2' * ( Tf - T2);
    tmp_vec = repmat(tmp_vec, [1 prod(sp)]);
    epl_pts1 = ( tmp_mat * (Kf\(pts{w}')) + repmat(nowModel{w}.D(:)',[3 1]) .* tmp_vec )';
    epl_pts1 = epl_pts1 ./ repmat(epl_pts1(:,3),[1 3]);
    epl_pts2 = ( tmp_mat * (Kf\(pts{w}')) + repmat(newModel{w}.D(:)',[3 1]) .* tmp_vec )';
    epl_pts2 = epl_pts2 ./ repmat(epl_pts2(:,3),[1 3]);
    [mpD1 mpD2 Check1 Check2 outB1 outB2] = warp2DepthConsiderNoise(D1,D2,D1_Obj,D2_Obj,D1_Obj.otherView{a},D2_Obj.otherView{a},epl_pts1,epl_pts2);

    Check = [reshape(Check1,[],1); reshape(Check2,[],1)] ;
    
    Dp = [reshape(D1,[],1) ; reshape(D2,[],1)];
    
    Dmp = [reshape(mpD1,[],1) ; reshape(mpD2,[],1)];
    nCheck = ~Check;
    occ = Dp(nCheck)<=Dmp(nCheck);
    nCheck(nCheck>0)=occ;
    occ = nCheck | [outB1(:) ; outB2(:)];
    epl_pts = [epl_pts1 ; epl_pts2];
    clear epl_pts1 epl_pts2
    M = vgg_interp2(vals.I{m}, epl_pts(:,1), epl_pts(:,2), 'linear', oobv);
    M = squeeze(M) - vals.R{w};
    M = vals.ephoto(M);
%     figure(8);imshow(reshape(M(1:end/2),sp)/max(max(M(1:end/2))));
%     figure(9);imshow(reshape(M(end/2+1:end),sp)/max(max(M(end/2+1:end))));
    PcCost = [PcCost [reshape(M(1:end/2),sp(1),[]) ; reshape(M(end/2+1:end),sp(1),[])]];
    Epc = reshape(M, tp, 2);
    % calc GMM
    Ecol = reshape(lambdaColor*-GMM(vals.R{w},nowModel{w}.segMap,newModel{w}.segMap,nowModel{w}.GMM_Name,newModel{w}.GMM_Name,nowModel{w}.table,newModel{w}.table),tp,2);
    PGmm = [PGmm [reshape(Ecol(1:end/2),sp(1),[]) ; reshape(Ecol(end/2+1:end),sp(1),[])]];
    % calc parallax
    Epar = reshape(lambdaParallax*(1-ParallaxModel(nowModel{w}.D,newModel{w}.D,nowModel{w},newModel{w},vals.d_step,vals.ndisps,nowModel{w}.table,newModel{w}.table)),tp,2);
    Ppar = [Ppar [reshape(Epar(1:end/2),sp(1),[]) ; reshape(Epar(end/2+1:end),sp(1),[])]];
    % data term
    IA = cast((Epc + Ecol+ Epar) ,class(Kinf));
    if vals.visibility
        % not used
        % Set up photoconsistency edges
        E(1+(w-1)*tp:w*tp,2,1,1) = IA(:,1);
        E(1+(w-1)*tp:w*tp,4,2,1) = IA(:,2);
    else
        % Add the photoconsistency terms to the unaries
        U(1+(w-1)*tp:w*tp,:) = U(1+(w-1)*tp:w*tp,:) + IA;
    end
    uniqueElement = uniqueElement + [numel(unique(nowModel{w}.segMap)); numel(unique(newModel{w}.segMap))];
    uniqueObj1 = [uniqueObj1 ; unique(nowModel{w}.segMap)];
    uniqueObj2 = [uniqueObj2 ; unique(newModel{w}.segMap)];
end
figure(4); imshow(PcCost/25);
figure(5); imshow(PGmm/25);
figure(6); imshow(Ppar/2);
uniqueObj1 = unique(uniqueObj1);
uniqueObj2 = unique(uniqueObj2);

EI_ = EI;
if vals.visibility
    % reshape to 4*tp
    E = reshape(permute(E, [2 1 3 4]), 4, []);
    U = zeros(2, tp*tubeSz+tp*tubeSz*2, class(Kinf));
else
    U = U';
end
SE = [];
for w=1:size(newModel,1)
    % calc ocoh
    Oc = [reshape(nowModel{w}.segMap,1,[]) ; reshape(newModel{w}.segMap,1,[])];
    Oc = reshape(Oc(:,vals.SEI),4,[]);
    Oc = (diff(int32(reshape(Oc([1 3; 1 4; 2 3; 2 4]',:), 2, 4, [])))~=0);
    Oc = lambdaOcoh*Oc;
    Eoc = reshape(cast(Oc, class(E)), [], size(vals.SEI, 2)); 
    % calc pcoh
    Dc = [reshape(nowModel{w}.F,1,[]) ; reshape(newModel{w}.F,1,[])];
    Dc = reshape(Dc(:,vals.SEI),4,[]);
    Dc = diff(int32(reshape(Dc([1 3; 1 4; 2 3; 2 4]',:), 2, 4, [])))~=0;
    EqObjNeqDep = ~Oc & Dc;

    calcDeptDis = [reshape(nowModel{w}.D,1,[]) ; reshape(newModel{w}.D,1,[])];
    calcDeptDis = reshape(calcDeptDis(:,vals.SEI),4,[]);
    calcDeptDis = abs(diff(reshape(calcDeptDis([1 3; 1 4; 2 3; 2 4]',:), 2, 4, [])));
    Edc = (lambdaDcoh/2)*(calcDeptDis<=1 & EqObjNeqDep) + lambdaDcoh*(calcDeptDis>1 & EqObjNeqDep);
    Edc = reshape(cast(Edc, class(E)), [], size(vals.SEI, 2));
    SE = [SE Eoc + Edc];
    
    
    EI_ = [EI_ vals.SEI+tp*(w-1)];

end
E = [E SE];
MDLaux = numel(uniqueObj1)+numel(uniqueObj2);
if vals.visibility
    [tmpE tmpEI tmpU] = GenClique(nowModel,newModel,tp,num_in,tubeSz,uniqueObj1,uniqueObj2,lambdaMDL);
else
    [tmpE tmpEI tmpU] = GenClique(nowModel,newModel,tp,0,tubeSz,uniqueObj1,uniqueObj2,lambdaMDL);
end
E = [E tmpE];
EI_ = [EI_ tmpEI];
U = [U zeros(2,numel(uniqueObj1)+numel(uniqueObj2),'int32')]+tmpU;
info.timings(1,:) = cputime - t_start; % Time data term evaluation

clear tmpE tmpEI tmpU
info.timings(2,:) = cputime - t_start; % Time smoothness term evaluation
figure(1);
for a = 1:numel(vals.improve)
    t_start = cputime;
    % Fuse the two labellings, using contract and/or improve if desired
    qpbo_params = int32([tp*tubeSz + MDLaux ((vals.improve(a)==1)+(vals.improve(a)==4)*2) vals.contract(a) vals.contract(a)>0]);
    if vals.improve(a) == 4
        % Add callback function handle
        qpbo_params = {qpbo_params, @(L) (choose_labels(L, U, E, EI, SE, vals.SEI, TE, TEI, num_in, vals.visibility, 2, vals.independent) > 0)};
    end
    try
        % planar == false
        [M stats] = vgg_qpbo(U, EI_, E, qpbo_params);
        
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
        [U_ E_ SE_ V] = calc_vis_energy(N, U, E(:,1:size(EI, 2)), EI, SE, vals.SEI, TE, TEI, num_in,MDLaux);
        M = M(1:end-MDLaux);
        energy(a) = sum(U_) + sum(E_) + sum(SE_);
    end
    info.numbers(1,a) = sum(N);
    info.timings(3,a) = cputime - t_start + info.timings(2,a); % Time optimization
end
clear TEI TE U E SE EI_

if nargout > 3 || vals.show_output
    % Generate output visibilities
    T = (tp* tubeSz* N(1:end-MDLaux)) + (1:tp*tubeSz)';
    for b = 1:num_in
        V(1:tp*tubeSz,b) = V(T);
        T = T + 2*tp*tubeSz;
    end
    V(tp*tubeSz+1:end,:) = [];
end

if vals.show_output
    % Display the output figures
    N = N(1:end-MDLaux);
    MDLU = U_(end-MDLaux+1:end);
    U_ = U_(1:end-MDLaux);
    U_ = double(U_) + accum(EI(1,:)', E_, [tp*tubeSz 1]);
    U_ = reshape(U_, sp(1), sp(2)*tubeSz);
    if vals.visibility
        % Take off the occlusion costs and normalize
        EI = reshape(sum(V(1:tp*tubeSz,:), 2), sp(1), sp(2)*tubeSz);
        U_ = U_ - (num_in - EI) * double(occl_cost);
        E_ = EI ~= 0;
        U_(E_) = U_(E_) ./ EI(E_);
    end
    
    set(0, 'CurrentFigure', vals.show_output);
    subplot('Position', [1/3 0.5 1/3 0.5]);
    D1(N) = D2(N);
    sc(D1, 'contrast', vals.d_min+[0 vals.d_step]);
    subplot('Position', [2/3 0.5 1/3 0.5]);
    disp(['U:' num2str(sum(sum(U_)))]);    
    disp(['MDLU:' num2str(sum(MDLU))]);
    sc(U_, 'jet');
    subplot('Position', [0 0 1/3 0.5]);
    T = reshape(sc(reshape(M, sp(1), sp(2)*tubeSz), 'prism'), [], 3);
    I = M < 0;
    T(I,:) = 1 - (1 - T(I,:)) * 0.3; % Lighten unlabelled pixels set to 0
    I = M > 1;
    T(I,:) = T(I,:) .* 0.3; % Darken unlabelled pixels set to 1 by optimal splice
    T(M==0,:) = 1;
    T(M==1,:) = 0;
    sc(reshape(T, sp(1), sp(2)*tubeSz, 3), [0 1]);
    subplot('Position', [1/3 0 1/3 0.5]);
    sc(reshape(sum(V, 2), sp(1), sp(2)*tubeSz), [0 num_in], 'contrast');
    subplot('Position', [2/3 0 1/3 0.5]);
    U_ = -accum(vals.SEI(2,:)', SE_, [tp*tubeSz 1]);
    disp(['E:' num2str(sum(sum(U_)))]);    
    sc(reshape(U_, sp(1), sp(2)*tubeSz));
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

function [U E SE V] = calc_vis_energy(L, U, E, EI, SE, SEI, TE, TEI, num_in,MDLaux)
% Generate visibility maps
tp = numel(L);
V = true(2*(tp-MDLaux), num_in);
% V(TEI(2,L(TEI(1,:))~=TE')-tp) = false;

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
function outval = ParallaxModel(D1,D2,ObjModel1,ObjModel2,d_step,ndisps,D1_table,D2_table)
    sz = size(D1);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    probMap1 = zeros(sz);
    probMap2 = zeros(sz);

    
    ObjIds = unique(ObjModel1.segMap);
    for i= 1:numel(ObjIds)
        s = ObjIds(i);
        M = ObjModel1.segMap==s;
        N = ObjModel1.plane(D1_table(s),:)';
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
        N = ObjModel2.plane(D2_table(s),:)';
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
    R = R(1:end/2,:);
    probMap1 = zeros(size(R,1),1);
    probMap2 = zeros(size(R,1),1);
    ObjIds = unique(ObjMap1);
    for i= 1:numel(ObjIds)
        s = ObjIds(i);
        samples = R(ObjMap1==s,:);
        prob = PredictGMM(samples',ModelName1{D1_table(s)});
        probMap1(ObjMap1==s) = prob;
    end
    ObjIds = unique(ObjMap2);
    for i= 1:numel(ObjIds)
        s = ObjIds(i);
        samples = R(ObjMap2==s,:);
        prob = PredictGMM(samples',ModelName2{D2_table(s)});
        probMap2(ObjMap2==s) = prob;
    end
    
    outval = [reshape(probMap1,[],1) ; reshape(probMap2,[],1)];
%     outval = outval - max(outval);
%     outval(outval<-20) = -20;
    middleV = mean(outval(1:end/2));
%     disp(['GMM mean obj1:' num2str(middleV)])
    outval(outval >= middleV) = middleV;
    outval = outval-middleV;
    outval(outval < -100) = -100;
    outval = outval/100;
%     figure(13); imshow(reshape(abs(outval(1:end/2))/6,size(ObjMap2)));
%     figure(14); imshow(reshape(abs(outval(end/2+1:end))/6,size(ObjMap2)));
return
function B = num_first(A)
% Return:        A  if numel(A) == 1
%         -numel(A) otherwise
B = -numel(A);
if B == -1
    B = A;
end
return
function [tmpE tmpEI tmpU] = GenClique(nowModel,newModel,tp,num_in,tubeSz,uniqueObj1,uniqueObj2,lambdaMDL)
% if visibility , here have to change index
MDLaux = numel(uniqueObj1) + numel(uniqueObj2);
tmpU = zeros(2,tp*tubeSz+MDLaux,'int32');
tmpE = zeros(4,0,'int32');
tmpEI = zeros(2,0,'uint32');

for w= 1:size(nowModel,1)
    ObjIds = unique(nowModel{w}.segMap);
    segNum1 = numel(ObjIds);
    for i= 1:segNum1
        s = ObjIds(i);
        index = tp*(w-1) + find(nowModel{w}.segMap==s);
        index(:,2) = tp*tubeSz+i;
        tmpEI = [tmpEI  index' ];
    end
    tmpE = [tmpE [ones(1,tp,'int32')*lambdaMDL; zeros(3,tp,'int32')]];
    
    ObjIds = unique(newModel{w}.segMap);
    for i= 1:numel(ObjIds)
        s = ObjIds(i);
        index = tp*(w-1) + find(newModel{w}.segMap==s);
        index(:,2) = tp*tubeSz+numel(uniqueObj1)+i;
        tmpEI = [tmpEI  index' ];
    end
    tmpE = [tmpE [ zeros(3,tp,'int32'); ones(1,tp,'int32')*lambdaMDL]];
    
end
tmpU(2,tp*tubeSz+1:tp*tubeSz+numel(uniqueObj1))  = lambdaMDL;
tmpU(1,tp*tubeSz+numel(uniqueObj1)+1:end) = lambdaMDL;
return