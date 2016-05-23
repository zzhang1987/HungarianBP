function [KP, KQ, beta, lamQ] = conKnlGphPQDChar(gphs, parKnl)
    gph1 = gphs{1};
    gph2 = gphs{2};
    [n1, m1] = size(gph1.G);
    [n2, m2] = size(gph2.G);
    KP = zeros(n1, n2)';
    % distance
    Dst1 = repmat(gph1.dsts', 1, m2);
    Dst2 = repmat(gph2.dsts, m1, 1);
    Dst = abs(Dst1 - Dst2) / max(max(max(Dst1, Dst2) + eps));

    % angle
    Ang1 = repmat(gph1.angs', 1, m2);
    Ang2 = repmat(gph2.angs, m1, 1);
    Ang = abs(Ang1 - Ang2);
    
    % combine distance and angle
    KQ = exp(-(Dst + Ang) / 2);
    