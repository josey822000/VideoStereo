function Proposal = ExpandProposal(Key,plane,oldF)
    oldF = unique(oldF);
    Seg = unique(Key.F);
    Obj = unique(Key.segMap);
    sz = size(Key.D);
    [X Y] = meshgrid(1:sz(2),1:sz(1));
    SegCnt = histc(Key.F(:),1:Seg(end));
    Que = find(SegCnt>500);
    Proposal = cell(numel(Que)+1,1);
    Proposal{1} = Key;
    for s=1:numel(Que)
        segId = Que(s);
        M = Key.F == segId;
        ObjId = Key.segMap(M);
        Proposal{s+1}.segMap = zeros(sz)+ObjId(end);
        Proposal{s+1}.F = zeros(sz)+segId;
        plnId = find(oldF==segId);
        Proposal{s+1}.D = -(X*plane(plnId,1) + Y*plane(plnId,2) + plane(plnId,3));
        Proposal{s+1}.ObjPln = Key.ObjPln(Key.table(ObjId(end)),:);
        Proposal{s+1}.parallax = Key.parallax(Key.table(ObjId(end)),:);
        Proposal{s+1}.GMM_Name = Key.GMM_Name(Key.table(ObjId(end)));
    end
end