function [N info energy V] = objFuse( D1, D2, vals, D1_Obj, D2_Obj, options )
%OBJFUSE Summary of this function goes here
%   Detailed explanation goes here
t_start = cputime;

    
% Initialize values
Kinf = int32(0);
scale_factor = 1e5 / vals.ephoto(1e6);
occl_cost = cast(scale_factor * vals.occl_val, class(Kinf));
Kinf = occl_cost + 1; % Smaller value seems to avoid errors in QPBO
info.timings = zeros(3, 1, numel(vals.improve));
info.numbers = zeros(4, 1, numel(vals.improve), 'uint32');
energy = zeros(1, numel(vals.improve));
%vals.I = images(2:end);
num_in = numel(vals.I);
% image sz
sp = size(D1);
% pixel num
tp = numel(D1);
planar = size(vals.SEI, 1) == 3;
if planar
   sprintf('err');
   return;
end
oobv = cast(-1000, class(vals.R));
out_unlabel = vals.improve(end) < 0;

% Calculate the homogenous coordinates of our two labellings
X = repmat(1:sp(2), [sp(1) 1]);
Y = repmat((1:sp(1))', [1 sp(2)]); % Faster than meshgrid
%                 X   Y   1   D1
%  WC =  
%                 X   Y   1   D2
WC = [repmat([X(:) Y(:)], [2 1]) ones(2*tp, 1) [D1(:); D2(:)]];
clear X Y

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
U = zeros(tp, 2, class(Kinf));
TEI = zeros(2, 0, 'uint32');
TE = zeros(4, 0, class(Kinf));
%Josey : U(nary) term => Epc,Ecol,Epar
% lambda
lambdaOcc = 25;
lambdaOcoh = 25;
lambdaDcoh = 25;
lambdaColor = 4;
lambdaParallax = 2;
lambdaMdl = 100000;
% 
% ObjMap of Image{a}
ObjCheck = segpln_Obj.otherView{a}.segMap;
DepCheck = segpln_Obj.otherView{a}.F;
% D1,D2
Check = (ObjCheck ~= segpln_Obj.segMap) & (DepCheck ~= segpln_Obj.F);
V1 = [D1_V1;D2_V2];
V2 = [D1_V2;D2_V2];
occ = V1(~Check)<V2(~Check);
M = vgg_interp2(vals.I{a}, T(:,1), T(:,2), 'linear', oobv);
M = squeeze(M) - vals.R;
M = vals.ephoto(M);
M(M<lambdaOcc) = lambdaOcc;
Epc = reshape(cast(scale_factor *(occ.*lambdaOcc + Check.*M) + (~(occ+Check))*realmax, class(Kinf)), tp, 2);
clear M N Check V1 V2 occ ObjCheck DepCheck
% Obj Color term
Ecol = lambdaColor*-log(GMM(vals.R,segpln_Obj));
% Obj parallax
Epar = lambdaParallax*(1-ParallaxModel(D1-Dop,Op));
% Object MDL
%     Emdl = lambdaMDL*(numel(unique(segpln_Obj.SegMap)));
% For each input image...
for a = 1:num_in
    % Calculate the coordinates in the input image
    T = WC * vals.P(:,:,a);
    N = 1 ./ T(:,3);
    T(:,1) = T(:,1) .* N;
    T(:,2) = T(:,2) .* N;
    
    % Calculate photoconsistency -> IA
    M = vgg_interp2(vals.I{a}, T(:,1), T(:,2), 'linear', oobv);
    M = squeeze(M) - vals.R;
    IA = reshape(cast(scale_factor * vals.ephoto(M), class(Kinf)), tp, 2);
    clear M N

    % Find interactions
    T(:,3) = T(:,3) ./ WC(:,4);
    [T M] = sortrows(T);
    N = find_interactions(T, 0.5); % Optimized version
    N = M(N); % Unsort
    N = uint32(N(:,abs(diff(N))~=tp)); % Remove interactions between the same node

    % Add the pixel interactions to the graph
    M = N(1,:) > tp;
    TEI = [TEI [N(1,:)-uint32(tp*M); (tp*2*a-tp)+N(2,:)]];
    T = zeros(4, numel(M));
    T(2,~M) = Kinf;
    T(4,M) = Kinf;
    TE = [TE T];

    if vals.visibility
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
    E = reshape(permute(E, [2 1 3 4]), 4, []);
    if vals.compress_graph
        % The unary and pairwise energies as they stand are entirely
        % correct, i.e. will give the correct labelling. However, it can be
        % compressed into a smaller but equivalent graph, which will be
        % faster to solve, by removing superfluous nodes and edges.
        [U EI EI_ E N T] = compress_graph(U', EI, E, TEI, TE, tp, num_in);
    else
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

% Add surface smoothness constraints
SE = (double([D1(:)'; D2(:)']) - vals.d_min) / vals.d_step;
SE = SE(:,vals.SEI);

% Fronto-parallel prior - Finite differences 1st derivative of
% disparity
SE = reshape(SE, 4, []);
SE = diff(reshape(SE([1 3; 1 4; 2 3; 2 4]',:), 2, 4, []));

% Apply our smoothness weighting
SE = reshape(cast(scale_factor * vals.esmooth(SE(:)), class(E)), [], size(vals.SEI, 2));
if ~planar
    E = [E SE];
    EI_ = [EI_ vals.SEI];
end
info.timings(2,:) = cputime - t_start; % Time smoothness term evaluation

for a = 1:numel(vals.improve)
    t_start = cputime;
    % Fuse the two labellings, using contract and/or improve if desired
    qpbo_params = int32([tp ((vals.improve(a)==1)+(vals.improve(a)==4)*2) vals.contract(a) vals.contract(a)>0]);
    if vals.improve(a) == 4
        % Add callback function handle
        qpbo_params = {qpbo_params, @(L) (choose_labels(L, U, E, EI, SE, vals.SEI, TE, TEI, num_in, vals.visibility, 2, vals.independent) > 0)};
    end
    try
        
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
        [U_ E_ SE_ V] = calc_vis_energy(N, U, E(:,1:size(EI, 2)), EI, SE, vals.SEI, TE, TEI, num_in);
        energy(a) = sum(U_) + sum(E_) + sum(SE_);
    end
    info.numbers(1,a) = sum(N);
    info.timings(3,a) = cputime - t_start + info.timings(2,a); % Time optimization
end
clear TEI TE U E SE EI_

if nargout > 3 || vals.show_output
    % Generate output visibilities
    T = (tp * N) + (1:tp)';
    for b = 1:num_in
        V(1:tp,b) = V(T);
        T = T + 2*tp;
    end
    V(tp+1:end,:) = [];
end

if vals.show_output
    % Display the output figures
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
    D1(N) = D2(N);
    sc(D1, 'contrast', vals.d_min+[0 vals.d_step]);
    subplot('Position', [2/3 0.5 1/3 0.5]);
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
    sc(reshape(U_, sp(1), sp(2)));
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

function [U E SE V] = calc_vis_energy(L, U, E, EI, SE, SEI, TE, TEI, num_in)
% Generate visibility maps
tp = numel(L);
V = true(2*tp, num_in);
V(TEI(2,L(TEI(1,:))~=TE')-tp) = false;

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

function B = num_first(A)
% Return:        A  if numel(A) == 1
%         -numel(A) otherwise
B = -numel(A);
if B == -1
    B = A;
end
return
