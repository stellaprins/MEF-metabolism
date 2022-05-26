function [gene, kcat, abundance, abundanceStarKcat] = ATPaseRules(ruleTable, dataGenes, dataExpression, sampleColumn)
% This takes Christina Kiel format of ATPase rules (ruleTable), a gene
% table with gene names (dataGenes) and gene table with expression values
% (dataExpression) and a string defining column within dataExpression to be
% analysed

% In initial version    ATPases = ruleTable
%                       dataRaw = dataGenes
%                       dataAtlasPC = dataExpression
%                       AverageWt = sampleColumn


% In a future version the multiplier vector might be an arrgument as well

geneRules = unique(ruleTable.AVERAGEGROUP_GENERULE);
geneRules = geneRules(2:end);           % lose the leading blank relating to where there are no rules

                                        % presentSUM = ~ismissing(ATPases.SUM);
                                        
multiplier = [2, 1, 1, 16, 1, 8];       % these are defined in the ATPases table rules Group 1 - 6 respectively

% Create a column for all proteins
ruleTable.Protein = strcat(ruleTable.SUM,  ruleTable.AVERAGE);

% Create a column for associated abundances. Need to look up gene in
% dataRaw as gene names not present in dataAtlas etc
for ATPaseCount = 1:size(ruleTable,1)
    ruleTable.Abundance{ATPaseCount} = lookup(dataGenes.NEWSymbol, ruleTable.Protein{ATPaseCount}, dataExpression.(sampleColumn));
end

% Create a column with adjusted abundances according to 'rules' where there
% is no rule
for ATPaseCount = 1:size(ruleTable,1)
    if string(ruleTable.SUM{ATPaseCount}) ~= string('')
        ruleTable.AbundanceMod{ATPaseCount} = ruleTable.Abundance{ATPaseCount};
    end
end    

% Deal with the more complex case where rules are taken into account
for ruleCount = 1:size(geneRules,1)    
    ruleIndex = find(ismember(ruleTable.AVERAGEGROUP_GENERULE, geneRules{ruleCount}));
    groupMean = mean(cell2mat(ruleTable.Abundance(ruleIndex)))*multiplier(ruleCount);
    ruleTable.AbundanceMod(min(ruleIndex)) = num2cell(groupMean);         % insert value for one instance only
    ruleTable.AbundanceMod(ruleIndex(2:end)) = num2cell(0);               % zero remaining values
end

% Calculate product abundance and kcat
for ATPaseCount = 1:size(ruleTable,1)
    ruleTable.kcat(ATPaseCount) = sum([ruleTable.kcat_s_1_(ATPaseCount), ruleTable.kcat_s_1__1(ATPaseCount)], 'omitnan');
    ruleTable.AbundKcat(ATPaseCount)= ruleTable.kcat(ATPaseCount) * ruleTable.AbundanceMod{ATPaseCount};
end

gene = ruleTable.Protein;
kcat = ruleTable.kcat;
abundance = ruleTable.AbundanceMod;
abundanceStarKcat = ruleTable.AbundKcat;

end

