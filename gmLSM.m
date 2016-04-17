function assign = gmLSM(K, Ct, asgT, options)
    %
    %Graph Matching via Local Sparse Model
    %Reference: @paper{AAAI159386,
    % 	author = {Bo Jiang and Jin Tang and Chris Ding and Bin Luo},
    % 	title = {A Local Sparse Model for Matching Problem},
    % 	conference = {AAAI Conference on Artificial Intelligence},
    % 	year = {2015},
    % 	keywords = {feature matching; sparse model; match selection},
    % 	abstract = {Feature matching problem that incorporates pairwise constraints is usually formulated as a quadratic assignment problem (QAP). Since it is NP-hard, relaxation models are required. In this paper, we first formulate the QAP from the match selection point of view; and then propose a local sparse model for matching problem. Our local sparse matching (LSM) method has the following advantages: (1) It is parameter-free; (2) It generates a local sparse solution which is closer to a discrete matrix than most other continuous relaxation methods for the matching problem. (3) The one-to-one matching constraints are better maintained in LSM solution. Promising experimental results show the effectiveness of the Proposed LSM method.},
    % 
    % 	url = {http://www.aaai.org/ocs/index.php/AAAI/AAAI15/paper/view/9386}
    % }
    %
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