function wsSrc = motorAsgSrc(pFs, nOus)
% Generate CMU Motion source for assignment problem.
%
% Input
%   pFs      -  frame index, 1 x 2
%   nOus     -  #outliers, 1 x 2
%
% Output
%   wsSrc
%     prex   -  prex
%     asgT   -  ground truth assignment
%     Pts    -  graph node set, 1 x mG (cell), 2 x ni
%     ords   -  order, 1 x mG (cell), 1 x ni
%
% History
%   create   -  Zhen Zhang (zhen@zzhang.org), 

% save option
%prex = cellStr(tag, pFs, nIns);
Fname = sprintf('%s_%d.mat', ...
    './Cars_and_Motorbikes_Graph_Matching_Datasets_and_Code/Data_for_Cars_and_Motorbikes/Data_Pairs_Motorbikes/pair', pFs);
prex = sprintf('Motor_pair_%d', pFs);
load(Fname);
nIn = nF1;

Inliers1 = features1(1:nIn, 2:-1:1)';
Inliers2 = features2(1:nIn, 2:-1:1)';

Inliers1 = Inliers1  ;
Inliers2 = Inliers2  ;
% load
[n1,~] = size(features1);
[n2,~] = size(features2);
n1 = n1 - nIn;
n2 = n2 - nIn;
ord1 = randperm(n1);
ord1 = ord1(1:nOus);
Outliers1 = features1(nIn + ord1, 2:-1:1)';
ord2 = randperm(n2);
ord2 = ord2(1:nOus);
Outliers2 = features2(nIn + ord2, 2:-1:1)';

Pts{1} = [Inliers1, Outliers1];
Pts{2} = [Inliers2, Outliers2];

% marker position
%Pts = CMUM.(tag).XTs(pFs);

% ground-truth assignment
XT = [eye(nIn), zeros(nIn, nOus); zeros(nOus, nIn +  nOus)] ;
asgT.alg = 'truth';
asgT.X = XT;
ords = {[1:nIn, ord1 + nIn], [1:nIn, ord2  + nIn]};
Features={features1, features2};
Fs={I1,I2};
% store
wsSrc.prex = prex;
wsSrc.Pts = Pts;
wsSrc.asgT = asgT;
wsSrc.ords = ords;
wsSrc.tag = 'Motor';
wsSrc.pFs = pFs;
wsSrc.nIns = nIn;
wsSrc.Features = Features;
wsSrc.Fs = Fs;
% save


