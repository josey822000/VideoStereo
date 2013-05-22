    ONum = 20;
    sameTable = false(ONum,ONum);
    
    plane = [1 1 3 2 2 5 5 4 1 1 3 3 3 5 2 1 1 4 2 3];
    plane = plane';
    %% looking for the same object plane
    [T M] = sortrows(plane);
    % find a within aThr
    T = abs(diff(T));
    T = find(T == 0);
    % if the index is continue than this plane is inside a Thr
    diffT = diff(T);
    diffT = find(diffT > 1);
    diffT(end+1) = numel(T);
    start = 0;
    for s = 1:numel(diffT)
        a = M([T(start+1:diffT(s)); T(diffT(s))+1] );
        c = reshape( repmat(a', numel(a), 1), numel(a) * numel(a), 1 );
        d = repmat(a(:), length(a), 1);
        sameTable((c-1)*double(ONum) + d) = 1;
        start = diffT(s);
    end
    