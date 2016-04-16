function [x] = BaBMatchingSolver(NofNodes, Factors, Potentials, U1, V1, BaBOptions)
NofStates = NofNodes * ones(1, NofNodes);
GUB = 1e20;
GLB = -1e20;
U = zeros(1, NofNodes);
V = zeros(1, NofNodes);

RootNode.Constraints = [];
RootNode.UB = inf;
RootNode.LB = -inf;
RootNode.Factors = Factors;
RootNode.Potentials = Potentials;
if(~isempty(U1) && ~isempty(V1))
    RootNode.U = U1;
    RootNode.V = V1;
else
    RootNode.U = U;
    RootNode.V = V;
end

%RootNode.Price = 0.1;

Q = [];

Q = [Q,RootNode];
BaBIter = 1;
while(~isempty(Q))
    AllUB = [Q.UB];
    [CUB, idx] = max(AllUB);
    if(CUB <= GLB)
        break;
    end
    CNode = Q(idx);
    Q(idx) = [];
    AllUB(idx) = [];
    nPotentials1 = CNode.Potentials;
    [nc, ~] = size(CNode.Constraints);
    NofAssigns = zeros(1, NofNodes);
    IsFeasible = true;
    
    for i=1:nc
        NodeIdx = CNode.Constraints(i, 1);
        NodeValue = CNode.Constraints(i, 2);
        if(NodeValue >  0)
            if(NofAssigns(NodeIdx) ~= 0 && NofAssigns(NodeIdx)~= NodeValue)
                IsFeasible = false;
                break;
            else
                NofAssigns(NodeIdx) = NodeValue;
                for xi = 1:NofStates
                    if(xi ~= NodeValue)
                        nPotentials1{NodeIdx}(xi) = -100;
                    end
                end
            end
        else
            NodeValue = - NodeValue;
            nPotentials1{NodeIdx}(NodeValue) = -100;
        end
    end
    if(IsFeasible == false)
        continue;
    end
    [xhat, nFactors, nPotentials1, dual, U, V] = HungarianBPMex(NofNodes, NofStates, CNode.Factors, nPotentials1, CNode.U, CNode.V, 1, BaBOptions.bpoptions.innerIter,GLB);
    UB = dual;
    AllUB = [AllUB(:); dual];
    GUB = min(GUB, max(AllUB));
    LB = CluComputeObjMex(xhat, NofStates, Factors, Potentials);
    if(LB >= GLB)
        GLB = LB;
        x = xhat;
    end
    if(UB > GLB)
        NodePotentials = cell2mat(nPotentials1(1:NofNodes));
%        figure(1);imshow(imresize(reshape(exp(NodePotentials),[NofNodes, NofNodes]),4, 'nearest'), [] );
%       drawnow;
        Prices = zeros(1, NofNodes);
        Values = zeros(1, NofNodes);
        for i=1:NofNodes
            [maxv, Values(i)] = max(NodePotentials(i,:));
            offset = abs(NodePotentials(i,:) - maxv);
            Prices(i) = sum(offset < 1e-4);
        end
        [~, NodeToSplit] = max(Prices);
        NodeAssign = Values(NodeToSplit);
        LNode.UB = dual;
        LNode.LB = -1000;
        LNode.Factors = nFactors;
        LNode.Potentials = nPotentials1;
        LNode.U = U;
        LNode.V = V;
        LNode.Constraints = [CNode.Constraints; NodeToSplit, NodeAssign];
        RNode.UB = dual;
        RNode.LB = -1000;
        RNode.Constraints = [CNode.Constraints; NodeToSplit, -NodeAssign];
        RNode.Factors = nFactors;
        RNode.Potentials = nPotentials1;
        RNode.U = U;
        RNode.V = V;
        Q = [Q; LNode; RNode];
    end
    


  %  fprintf('Iter = %d, GUB = %12.7f, GLB= %12.7f, GAP = %12.2f%%\n', BaBIter, GUB, GLB, (GUB - GLB) / GLB * 100);

    if(abs(GUB - GLB)/GLB*100 < 0.5)
        break;
    end
    if(BaBIter >= BaBOptions.MaxIter)
        break;
    end
    BaBIter = BaBIter + 1;
end

end