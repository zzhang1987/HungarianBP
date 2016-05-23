function gph = newGphA(Pt, G1)
% Generate a graph by connecting points.
%
% Remark
%   The edge feature is asymmetric.
%
% Input
%   Pt      -  graph node, d x n
%   parGph  -  parameter for computing the adjacency matrix
%              see gphEg.m for more details
%
% Output
%   gph     -  graph
%     Pt    -  graph node, d x n
%     Eg    -  graph edge, 2 x m
%     vis   -  binary indicator of nodes that have been kept, 1 x n | []
%     G     -  node-edge adjacency (for starting point), n x m
%     H     -  node-edge adjacency (for ending point), n x m
%     PtD   -  edge feature, 2 x m
%     dsts  -  distance, 1 x m
%     angs  -  angle, 1 x m
%
% History
%   create  -  Feng Zhou (zhfe99@gmail.com), 08-11-2011
%   modify  -  Feng Zhou (zhfe99@gmail.com), 03-06-2013

% dimension
n = size(Pt, 2);

% edge
G1 = triu(G1);
idxs = find(G1>0);
[Eg1, Eg2] = ind2sub([n,n], idxs);
Eg = [Eg1',Eg2'; Eg2', Eg1' ];
% incidence for asymmetric edge
[G, H] = gphEg2IncA(Eg, n);

% second-order feature
[PtD, dsts, angs, angAs] = gphEg2Feat(Pt, Eg);

% store
gph.Pt = Pt;
gph.Eg = Eg;
gph.G = G;
gph.H = H;
gph.PtD = PtD;
gph.dsts = dsts;
gph.angs = angs;
gph.angAs = angAs;
