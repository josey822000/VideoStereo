function [outD1 outD2 Check1 Check2 outBound1 outBound2] = warp2DepthConsiderDiffLabel( Depth1, Depth2, Obj1, Obj2,otherView1,otherView2,X1,X2,DiffThr)
%WARP2DEPTHCONSIDERNOISE Summary of this function goes here
%   Detailed explanation goes here
    [X Y] = meshgrid(1:size(Depth1,2), 1:size(Depth1,1));
    % For each image...
    outD1 = zeros(size(Depth1,1),size(Depth1,2));
    outD2 = zeros(size(Depth1,1),size(Depth1,2));

    % Project the points


    % Look up the colours
    warpX1 = ceil(X1(:,1));
    warpY1 = ceil(X1(:,2));
    outBound1 = warpX1 <1 | warpX1 > size(Depth1,2) | warpY1 <1 | warpY1 > size(Depth1,1);
    warpX1(warpX1<1 | warpX1 > size(Depth1,2)) = round(X(warpX1<1 | warpX1 > size(Depth1,2)));
    warpY1(warpY1<1 | warpY1 > size(Depth1,1)) = round(Y(warpY1<1 | warpY1 > size(Depth1,1)));

    outD1 = otherView1.Disparity(warpY1+size(Depth1,1)*(warpX1-1));
%     mpSeg1_1 = reshape(otherView1.segMap(warpY1+size(Depth1,1)*(warpX1-1)), size(Depth1));
    mpSeg1_1 = abs(reshape(outD1, size(Depth1))-Depth1);
    warpX1_1 = warpX1-1;
    warpX1_1(warpX1_1<1) = 1;
    mpSeg1_2 = abs(reshape(otherView1.Disparity(warpY1+size(Depth1,1)*(warpX1_1-1)), size(Depth1))-Depth1);
%     CheckObj1 = mpSeg1_1 == Obj1.segMap | mpSeg1_2 == Obj1.segMap ;
    Check1 = mpSeg1_1 < DiffThr | mpSeg1_2 < DiffThr;
%     figure(6);
%     imshow(CheckObj1);
    clear mpSeg1_1 mpSeg1_2 
%     mpF1_1 = reshape(otherView1.F(warpY1+size(Depth1,1)*(warpX1-1)), size(Depth1));
%     mpF1_2 = reshape(otherView1.F(warpY1+size(Depth1,1)*(warpX1_1-1)), size(Depth1));
%     CheckDep1 = mpF1_1 == Obj1.F | mpF1_2 == Obj1.F;
%     Check1 = CheckObj1 & CheckDep1;
%     clear mpF1_1 mpF1_2 CheckObj1 CheckDep1
    %

    % Look up the colours
    warpX2 = ceil(X2(:,1));
    warpY2 = ceil(X2(:,2));
    outBound2 = warpX2 <1 | warpX2 > size(Depth2,2) | warpY2 < 1 | warpY2 > size(Depth2,1);
    warpX2(warpX2<1 | warpX2 > size(Depth2,2)) = round(X(warpX2<1 | warpX2 > size(Depth2,2)));
    warpY2(warpY2<1 | warpY2 > size(Depth2,1)) = round(Y(warpY2<1 | warpY2 > size(Depth2,1)));

    outD2 = otherView2.Disparity(warpY2+size(Depth2,1)*(warpX2-1));

%     mpSeg2_1 = reshape(otherView2.segMap(warpY2+size(Depth2,1)*(warpX2-1)), size(Depth2));
    mpSeg2_1 = abs(reshape(outD2, size(Depth2)) - Depth2);

    warpX2_1 = warpX2-1;
    warpX2_1(warpX2_1<1) = 1;
    mpSeg2_2 = abs(reshape(otherView2.Disparity(warpY2+size(Depth2,1)*(warpX2_1-1)), size(Depth2)) - Depth2);
%     CheckObj2 = mpSeg2_1 == Obj2.segMap | mpSeg2_2 == Obj2.segMap;
    Check2 = mpSeg2_1 < DiffThr | mpSeg2_2 < DiffThr;
    
    clear mpSeg2_1 mpSeg2_2 
%     mpF2_1 = reshape(otherView2.F(warpY2+size(Depth2,1)*(warpX2-1)), size(Depth2));
%     mpF2_2 = reshape(otherView2.F(warpY2+size(Depth2,1)*(warpX2_1-1)), size(Depth2));
% 
%     CheckDep2 = mpF2_1 == Obj2.F | mpF2_2 == Obj2.F ;
%     Check2 = CheckObj2 & CheckDep2;
%     clear mpF2_1 mpF2_2 



end

