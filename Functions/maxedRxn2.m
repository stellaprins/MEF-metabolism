function [rxnListLB, rxnListUB] = maxedRxn2(model,FBA)
% takes a model (with lb and ub fields) and an FBA.x list and returns rxns
% that are at the maximum bound

    for i = 1:size(FBA,2)
       diffLB{i} = model.lb - FBA(i).x;
       rxnListLB{i} = model.rxns(diffLB{i} == 0 & FBA(i).x ~= 0);
       diffUB{i} = model.ub - FBA(i).x;
       rxnListUB{i} = model.rxns(diffUB{i} == 0 & FBA(i).x ~= 0);
    end
end

