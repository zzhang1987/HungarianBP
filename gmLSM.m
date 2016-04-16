function assign = gmLSM(K, Ct, asgT, options)
    [NofNodes, ~] = size(Ct);
    
    x = ones(NofNodes^2, 1) / NofNodes;
    bestv = -1e20;
    lambda = -1e20;
    Xast = [];
    for i = 1:10000
        Kx = K * x;
        lastv = lambda;
        lambda = x'*Kx;
        
        KxMat = reshape(Kx, [NofNodes, NofNodes]);
        Xmat = reshape(x, [NofNodes, NofNodes]);
        SumX = sum(Xmat); 
        SumXExpand = ones(NofNodes, 1) * SumX;
        SumXExpand = reshape(SumXExpand, [NofNodes^2, 1]);
        x = x.*sqrt(Kx./SumXExpand/lambda);
        
    %    figure(1);imshow(imresize(reshape(x * 255,[NofNodes, NofNodes]),4, 'nearest'), [] );SumXExpand;
     %   fprintf('Iter = %d, relax obj = %12.5f, Obj = %12.5f\n',i,lambda, bestv);
      %  drawnow;
        if(abs(lastv - lambda) < 1e-6)
            break;
        end
    end
    [IntX, cost] = hungarian(-reshape(x, [NofNodes, NofNodes]));
    IntX = IntX(:);
    DecodeObj = IntX' * K * IntX;
    if(DecodeObj > bestv)
        Xast = IntX;
        bestv = DecodeObj;
    end
    assign.X = Xast;
    acc = matchAsg(assign.X , asgT);
    assign.acc = acc;
    assign.obj = bestv;
end