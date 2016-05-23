clear variables;
global footpath;
%delete(gcp)
%parpool
footpath = cd;
footpath = strcat(footpath,'/fgm');
rng(45678);



rmpath('./HBPMex/');
rmpath(footpath);
rmpath(genpath([footpath '/src']));
rmpath(genpath([footpath '/lib']));


addpath('./HBPMex/');
addpath(footpath);
addpath(genpath([footpath '/src']));
addpath(genpath([footpath '/lib']));
prSet(1);
%% src parameter
tag = 'pas';
MaxInstance = 30;
NofAlgos = 11;
%% src parameter
tag = 'pas';




parKnl = st('alg', 'pas'); % type of affinity: only edge distance
MaxOutliers = 20;
AvgAcc = zeros(MaxOutliers, NofAlgos);
AvgObj = zeros(MaxOutliers, NofAlgos);
StdAcc = zeros(MaxOutliers, NofAlgos);
StdObj = zeros(MaxOutliers, NofAlgos);
AvgTime = zeros(MaxOutliers, NofAlgos);
StdTime = zeros(MaxOutliers, NofAlgos);

for nOut = 0:20
    accs = cell(MaxInstance, 1);
    objs = cell(MaxInstance, 1);
    times = cell(MaxInstance, 1);
    bpoptions.outIter = 1;
    bpoptions.innerIter = 5;
    BaBoptions.MaxIter = 600;
    BaBoptions.bpoptions = bpoptions;
    for kFs=1:MaxInstance
        acc = zeros(1, NofAlgos);
        obj = zeros(1, NofAlgos);
        tm = zeros(1, NofAlgos);
        
        %    nOut = 5 ;%-5 ; % randomly remove 2 nodes
        parKnl = st('alg', 'pas2'); % type of affinity: only edge distance
        %% algorithm parameter
        [pars, algs] = gmPar(2);
        
        %% src
        wsSrc = carAsgSrc(kFs, nOut);
        asgT = wsSrc.asgT;
        
        parG = st('link', 'del'); % Delaunay triangulation for computing the graphs
        parF = st('smp', 'n', 'nBinT', 4, 'nBinR', 3); % not used, ignore it
        wsFeat = motorAsgFeat(wsSrc, parG, parF, 'svL', 1);
        [gphs, XPs, Fs] = stFld(wsFeat, 'gphs', 'XPs', 'Fs');
        
        [KP, KQ] = conKnlGphPQD(gphs, parKnl);
        K = conKnlGphKD(KP, KQ, gphs);
        Ct = ones(size(KP));
        
        %tdd = tic;
        %asgDD = gmDD(K, Ct, asgT,BaBoptions);
        %tm(11) = toc(tdd);
        %acc(11) = asgDD.acc;
        %obj(11) = asgDD.obj;

        tlsm = tic;
        asgLsm = gmLSM(K, Ct, asgT, BaBoptions);
        tm(10) = toc(tlsm);
        acc(10) = asgLsm.acc;
        obj(10) = asgLsm.obj;
        %BP
        tbp = tic;
        asgBP = HungarianBP(K, Ct, asgT,BaBoptions);
        tm(1) = toc(tbp);        
        acc(1) = asgBP.acc;
        obj(1) = asgBP.obj;
        
        
        %% GA
        tGa = tic;
        asgGa = gm(K, Ct, asgT, pars{1}{:});
        tm(2) = toc(tGa);
        acc(2) = asgGa.acc;
        obj(2) = asgGa.obj;
        
        
        %% PM
        tPm = tic;
        asgPm = pm(K, KQ, gphs, asgT);
        tm( 3) = toc(tPm);
        acc(3) = asgPm.acc;
        obj(3) = asgPm.obj;
        
        
        
        %% SM
        tSm = tic;
        asgSm = gm(K, Ct, asgT, pars{3}{:});
        tm( 4) = toc(tSm);
        acc(4) = asgSm.acc;
        obj(4) = asgSm.obj;
        
        
        
        %% SMAC
        tSmac = tic;
        asgSmac = gm(K, Ct, asgT, pars{4}{:});
        tm(5) = toc(tSmac);
        acc(5) = asgSmac.acc;
        obj(5) = asgSmac.obj;
        
        
        %% IPFP-U
        tIpfu = tic;
        asgIpfpU = gm(K, Ct, asgT, pars{5}{:});
        tm(6) = toc(tIpfu);
        acc(6) = asgIpfpU.acc;
        obj(6) = asgIpfpU.obj;
        
        
        
        %% IPFP-S
        tIpfp = tic;
        asgIpfpS = gm(K, Ct, asgT, pars{6}{:});
        tm( 7) = toc(tIpfp);
        acc(7) = asgIpfpS.acc;
        obj(7) = asgIpfpS.obj;
        
        
        %% RRWM
        tRrwm = tic;
        asgRrwm = gm(K, Ct, asgT, pars{7}{:});
        tm( 8) = toc(tRrwm);
        acc(8) = asgRrwm.acc;
        obj(8) = asgRrwm.obj;
        
        
        % FGM-D
        tFgmD = tic;
        asgFgmD = fgmD(KP, KQ, Ct, gphs, asgT, pars{9}{:});
        tm( 9) = toc(tFgmD);
        acc(9) = asgFgmD.acc;
        obj(9) = asgFgmD.obj;
        
        % print information
        times{kFs} = tm;
        accs{kFs} = acc;
        objs{kFs} = obj;
    end
    times = cell2mat(times);
    objs = cell2mat(objs);
    accs = cell2mat(accs);
    maxobjs = max(objs,[], 2);
    normalised_objs = zeros(MaxInstance,NofAlgos);
    for i=1:NofAlgos
        normalised_objs(:,i) = objs(1:1:MaxInstance, i)./maxobjs(1:1:MaxInstance);
    end
    
    for i=1:NofAlgos
        AvgAcc(nOut+1, i) = mean(accs(:,i));
        StdAcc(nOut+1, i) = std(accs(:,i));
        AvgObj(nOut+1, i) = mean(normalised_objs(:, i));
        StdObj(nOut+1, i) = std(normalised_objs(:,i));
        AvgTime(nOut+1, i) = mean(times(:,i));
        StdTime(nOut+1, i) = std(times(:,i));
    end
    save('CarResult.mat');
end
