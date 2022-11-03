function [rxnTab] = getRxnFlux4Met3(model,FBA,met, nzflag)

    if nargin < 4
        nzflag = 1;
    end

    if met(end) == ']'
        [rxnList, rxnFormulaList] = findRxnsFromMets(model, met);
        index = find(ismember(model.rxns,rxnList));
        fluxes = FBA.x(index);
        lb = model.lb(index);
        ub = model.ub(index);
        rxnTab = table(index, rxnList, fluxes, lb, ub, rxnFormulaList);

    else
        rxnListTot = {};
        rxnFormulaListTot = {};
        fluxesTot = [];
        lbTot = [];
        ubTot = [];
        indexTot = [];

        [defaultCompartmentSymbolList, defaultCompartmentNameList] = getDefaultCompartmentSymbols();

        for i = 1:size(defaultCompartmentSymbolList,2)
            metTemp = strcat(met, '[', defaultCompartmentSymbolList{i}, ']');
            [rxnList, rxnFormulaList] = findRxnsFromMets(model, metTemp);
            index = find(ismember(model.rxns,rxnList));
            fluxes = FBA.x(index);
            lb = model.lb(index);
            ub = model.ub(index);
            rxnListTot = [rxnListTot;rxnList];
            rxnFormulaListTot = [rxnFormulaListTot; rxnFormulaList];
            fluxesTot = [fluxesTot; fluxes];
            lbTot = [lbTot; lb];
            ubTot = [ubTot; ub];
            indexTot = [indexTot; index];
            rxnList = rxnListTot;
            rxnFormulaList = rxnFormulaListTot;
            fluxes = fluxesTot;
            
        end

        indexRet = indexTot;
        rxnList = rxnListTot;
        rxnFormulaList = rxnFormulaListTot;
        fluxes = fluxesTot;
        lb = lbTot;
        ub = ubTot;
        rxnTab = table(indexRet, rxnList, fluxes, lb, ub, rxnFormulaList);
        rxnTab = unique(rxnTab);
        

    end

    if nzflag == 1
        index = rxnTab{:,3} ~= 0;
        rxnTab = rxnTab(index,:);
    end

end

