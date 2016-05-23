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
[pars, algs] = gmPar(2);

parKnl = st('alg', 'char'); % type of affinity: only edge distance
MaxInstance = 45;
NofAlgos = 10;
rng(123456);
AvgAcc =  zeros(4, NofAlgos);
StdAcc = zeros(4, NofAlgos);
AvgObj = zeros(4, NofAlgos);
StdObj = zeros(4, NofAlgos);
AvgTime = zeros(4, NofAlgos);
StdTime = zeros(4, NofAlgos);
for idx=3:4
    kFs = 1;
    times = cell(MaxInstance,1);
    objs = cell(MaxInstance,1);
    accs = cell(MaxInstance,1);
    for idx1 = 1:10
        for idx2 = (idx1+1):10
            idx1base = (idx - 1) * 10;
            data1src = sprintf('data_chrct/%d.mat', idx1 + idx1base);
            data2src = sprintf('data_chrct/%d.mat', idx2 + idx1base);
            d1=load(data1src);
            d2=load(data2src);
            n = size(d1.Pt,1);
            Fs{1} = d1.I;
            Fs{2} = d2.I;
            NN = randperm(n);
            G1 = zeros(n, n);
            Pt1 = zeros(2, n);
            for i1=1:n
                Pt1(:, i1) = [d1.Pt(NN(i1),2);d1.Pt(NN(i1),1)] ;
                for i2=1:n
                    G1(i1, i2) = d1.G(NN(i1), NN(i2));
                end
            end
            Pt2 = d2.Pt';
            g1 = newGphChar(Pt1, G1);
            g2 = newGphChar([Pt2(2,:);Pt2(1,:)], d2.G);
            gphs{1} = g1;
            gphs{2} = g2;
            
            [KP, KQ] = conKnlGphPQDChar(gphs, parKnl);
            K = conKnlGphKD(KP, KQ, gphs);
            Ct = ones(size(KP));
            
            %% undirected graph -> directed graph (for FGM-D)
            gphDs = gphU2Ds(gphs);
            KQD = [KQ, KQ; KQ, KQ];
            
            bpoptions.outIter = 1;
            bpoptions.innerIter = 5;
            BaBoptions.MaxIter = 600;
            BaBoptions.bpoptions = bpoptions;
            %SG
            %asgCutting = gmCutting(K,Ct, asgT, bpoptions);
            asgT.X = zeros(n,n);
            for i1=1:n
                asgT.X(i1, NN(i1)) = 1;
            end
%             tdd = tic;
%             asgDD = gmDD(K, Ct, asgT,BaBoptions);
%             tm(11) = toc(tdd);
%             acc(11) = asgDD.acc;
%             obj(11) = asgDD.obj;
            %DD excluded due to license problems
            acc(11) = 0;
            obj(11) = 0;
            tm(11) = 0;
            
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
            fprintf('BP    : acc %.2f, obj %.2f\n', asgBP.acc, asgBP.obj);
            fprintf('GA    : acc %.2f, obj %.2f\n', asgGa.acc, asgGa.obj);
            fprintf('PM    : acc %.2f, obj %.2f\n', asgPm.acc, asgPm.obj);
            fprintf('SM    : acc %.2f, obj %.2f\n', asgSm.acc, asgSm.obj);
            fprintf('SMAC  : acc %.2f, obj %.2f\n', asgSmac.acc, asgSmac.obj);
            fprintf('IPFP-U: acc %.2f, obj %.2f\n', asgIpfpU.acc, asgIpfpU.obj);
            fprintf('IPFP-S: acc %.2f, obj %.2f\n', asgIpfpS.acc, asgIpfpS.obj);
            fprintf('RRWM  : acc %.2f, obj %.2f\n', asgRrwm.acc, asgRrwm.obj);
            fprintf('FGM-D : acc %.2f, obj %.2f\n', asgFgmD.acc, asgFgmD.obj);
            fprintf('LSM : acc %.2f, obj %.2f\n', asgLsm.acc, asgLsm.obj);
           % fprintf('DD : acc %.2f, obj %.2f\n', asgDD.acc, asgDD.obj);
            times{kFs} = tm;
            accs{kFs} = acc;
            objs{kFs} = obj;
            kFs = kFs + 1;
         %   figure(1);
          %  subplot(2,2,idx);
           % rows = 1; cols = 1;
            %Ax = iniAx(1, rows, cols, [400 * rows, 900 * cols], 'hGap', .1, 'wGap', .1);
            %parCor = st('cor', 'ln', 'mkSiz', 7, 'cls', {'g', 'y', 'b'});
            %shAsgImg(Fs, gphs, asgBP, asgT, parCor, 'ax', Ax{1}, 'ord', 'n');
            %title('result of FGM-D');
        end
    end
    NofAlgos = 11;
    times = cell2mat(times);
    objs = cell2mat(objs);
    accs = cell2mat(accs);
    maxobjs = max(objs,[], 2);
    normalised_objs = zeros(MaxInstance,NofAlgos);
    for i=1:NofAlgos
        normalised_objs(:,i) = objs(1:1:MaxInstance, i)./maxobjs(1:1:MaxInstance);
    end
    
    for i=1:NofAlgos
        AvgAcc(idx, i) = mean(accs(:,i));
        StdAcc(idx, i) = std(accs(:,i));
        AvgObj(idx, i) = mean(normalised_objs(:, i));
        StdObj(idx, i) = std(normalised_objs(:,i));
        AvgTime(idx, i) = mean(times(:,i));
        StdTime(idx, i) = std(times(:,i));
    end
end
selected_algs = [2,4,6,7,5,8,9,10,1];

for i=1:4
    fprintf('\\hline\n');
    fprintf('Character%d(Acc)', i);
    fprintf('&%.4f ', AvgAcc(i,selected_algs)); fprintf('\\\\\n')
    fprintf('Character%d(Obj)', i);
    fprintf('&%.4f ', AvgObj(i,selected_algs)); fprintf('\\\\\n')
    fprintf('Character%d(Time)', i);
    fprintf('&%.4f ', AvgTime(i,selected_algs)); fprintf('\\\\\n')
    fprintf('\\hline\n');
end
