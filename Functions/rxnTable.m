function rxnTable = rxnTable(modelStore, FBA, index, S)

% rxnTable returns a table of reactions with fluxes and lower and upper
% bounds
% This is helpful in ascertaining how FBA solutions for a set of models are behaving
% The reaction itself is given in the final column for ease of reference

%   INPUT    modelStore         Collection of models (1 x n struct)
%   INPUT    FBA                FBA structures (1 x n struct)  
%   INPUT    index              index of reaction numbers to be used
%   INPUT    S                  Array of model identifiers (1 x n cell)

%   OUTPUT   rxnTable           Table with rows reactions and columns lb f ub sets 
%

% .. Author: - Phil Luthert & Christina Kiel 11/5/22

% convert rxns to strings. Note there is an assumption that all model
% structs are the same size. It would be helpful to add a test

    rxnString = convertCharsToStrings(modelStore(1).rxns);

% set up parameters
    tab = [];
    
% construct table   
    for i = 1:size(modelStore,2)
        tab = horzcat(tab,  modelStore(i).lb(index), ...
                            FBA(i).x(index), ...
                            modelStore(i).ub(index));
    end

% set up parameters
    varString = {};
    
% sort variable names
    for i = 1:size(modelStore,2)
        varTemp = {strcat(S{i}, 'lb'), strcat(S{i}, 'f'), strcat(S{i}, 'ub')}; 
        varString = [varString, varTemp];
    end

% construct table    
    rxnTable = array2table(tab, 'RowNames', rxnString(index), 'VariableNames', varString); 
    formulaTable = table(printRxnFormula(modelStore(1), cellstr(rxnString(index))));
    rxnTable = [rxnTable formulaTable];
end

