function [outD1 outD2 mpSeg1 mpSeg2 mpF1 mpF2 outBound1 outBound2] = warp2Depth( Depth1, Depth2, Obj1, Obj2,P)
%WARP2DEPTH Summary of this function goes here
%   Detailed explanation goes here
    sz = size(Depth1);
    [X Y] = meshgrid(1:sz(2), 1:sz(1));
    WC = ones(sz(1)*sz(2), 4);
    
    
    % For each image...
    outD1 = zeros(size(Depth1,1),size(Depth1,2),size(P,3));
    outD2 = zeros(size(Depth1,1),size(Depth1,2),size(P,3));
    for a = 1:size(P,3)
        % Project the points
        WC(:,1) = X(:);
        WC(:,2) = Y(:);
        WC(:,3) = 1;
        WC(:,4) = Depth1(:);
        X1 = WC * P(:,:,a)';

%         P_ = P(:,4,a);

        
%         D1 = reshape(Depth1, sz(2)*sz(1),[]);
%         D1 = repmat(D1, [1 3]);
%         D1(:,1) = D1(:,1)*P_(1);
%         D1(:,2) = D1(:,2)*P_(2);
%         D1(:,3) = D1(:,3)*P_(3);
%         d1 = D1;
%         Z1 = 1 ./ (X1(:,3) + d1(:,3));

        % Look up the colours
        warpX1 = round(X1(:,1));
        warpY1 = round(X1(:,2));
        outBound1 = warpX1 <1 | warpX1 > size(Depth1,2) | warpY1 < 1 | warpY1 > size(Depth1,1);
        warpX1(warpX1<1 | warpX1 > size(Depth1,2)) = round(X(warpX1<1 | warpX1 > size(Depth1,2)));
        warpY1(warpY1<1 | warpY1 > size(Depth1,1)) = round(Y(warpY1<1 | warpY1 > size(Depth1,1)));
        outD1 = Obj1.Disparity(warpY1+size(Depth1,1)*(warpX1-1));
        mpSeg1 = reshape(Obj1.segMap(warpY1+size(Depth1,1)*(warpX1-1)), size(Depth1));
        mpF1 = reshape(Obj1.F(warpY1+size(Depth1,1)*(warpX1-1)), size(Depth1));
        %
        WC(:,1) = X(:);
        WC(:,2) = Y(:);
        WC(:,3) = 1;
        WC(:,4) = Depth2(:);
        X2 = WC * P(:,:,a)';
%         D2 = reshape(Depth2, sz(2)*sz(1),[]);
%         D2 = repmat(D2, [1 3]);
%         D2(:,1) = D2(:,1)*P_(1);
%         D2(:,2) = D2(:,2)*P_(2);
%         D2(:,3) = D2(:,3)*P_(3);
%         d2 = D2;
%         Z2 = 1 ./ (X2(:,3) + d2(:,3));

        % Look up the colours
        warpX2 = round(X2(:,1));
        warpY2 = round(X2(:,2));
        outBound2 = warpX2 <1 | warpX2 > size(Depth2,2) | warpY2 < 1 | warpY2 > size(Depth2,1);
        warpX2(warpX2<1 | warpX2 > size(Depth2,2)) = round(X(warpX2<1 | warpX2 > size(Depth2,2)));
        warpY2(warpY2<1 | warpY2 > size(Depth2,1)) = round(Y(warpY2<1 | warpY2 > size(Depth2,1)));
        
        outD2 = Obj2.Disparity(warpY2+size(Depth2,1)*(warpX2-1));
%         outD2(outBound2) = -1;
        mpSeg2 = reshape(Obj2.segMap(warpY2+size(Depth2,1)*(warpX2-1)), size(Depth2));
%         mpSeg2(outBound2) = -100;
        mpF2 = reshape(Obj2.F(warpY2+size(Depth2,1)*(warpX2-1)), size(Depth2));
%         mpF2(outBound2) = -100;
        
    end


end

