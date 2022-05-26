function FBA = FBAdummy(model)
% FBAdummy sets up a dummy FBA result for instances where it's required for
% array of FBA results and an error is thrown if there is no solution.


%   INPUT    model              COBRA model

%   OUTPUT   FBAresult          FBA dummy results

    model.lb(:) = -1000;
    model.ub(:) = 1000;
    rxnNumber = size(model.rxns,1);
    FBA = optimizeCbModel(model,'max','zero'); % dummy run to extract fields
        if FBA.stat ~= 1
            disp("Error, FBAdummy requires a model with an FBA solution.");
        end
    FBA.x = zeros([rxnNumber 1]);
    FBA.v = zeros([rxnNumber 1]);
    FBA.full = zeros([rxnNumber 1]);
    
    
end

