function outTab = writeMet(outFileName, modelStore, FBA, S, metab, sheetname)
% writeMet takes xlsx spreadsheet name, model and FBA data and writes
% non-zero fluxes involving metab

index = [];

for i = 1:numel(S)
   [rxnTab] = getRxnFlux4Met3(modelStore(i),FBA(:,i),metab,1);
   index = [index; rxnTab{:,1}];
end

index = unique(index);

outTab = rxnTable(modelStore, FBA, index, S);
writetable(outTab,outFileName, 'Sheet', sheetname, 'WriteRowNames',true);




end

