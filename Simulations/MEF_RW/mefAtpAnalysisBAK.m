
% Script that uses preprocessed proteomic expression datasets to inform FBA
% using the COBRA Toolbox

% FUNCTIONS CALLED

%
% INPUT DATASETS. DATAKCAT HAS MISSING VALUES AND MODIFIED ATPase ABUNDANCE
% AND KCAT as well as KCAT*ABUNDANCE VALUES. DATARAW HAS GENENAMES
%

clear;
load('dataKcat');
load('dataRaw');
dataPC = dataKcat;

% lose the kcat column which is only there for reference
dataPC = removevars(dataPC,{'kcat'});

% now select subset of data
% Base data, columns 1-4
% ATPase mod, columns 5-8
% kcat rules included in addition, columns 9-12

dataPC = dataPC(:,9:12); % set column subsetting here
S = dataPC.Properties.VariableNames(1:4); 

%
% CONSTANTS / FLAGS
%

% Squeeze is a parameter for separate runs that determines level of
% compression of expression data
Squeeze = [0.0001 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.25 0.5 -1];

% Multiplier sets the ratio between exchange reaction flux and internal
% reaction fluxes that are set by expression levels
Multiplier = [0.003125 0.00625 0.0125 0.025 0.05 0.1 0.2 0.4 0.8 1.6 3.2];

includeRealFlux = 1;                                                    % set to use measured fluxes to constrain exchange reactions
boundsFileName = 'Opening_reactions2.xlsx';                             % contains exchange reaction bounds
filename_rxns = 'Rxn_Class_Recon3DModel.xlsx';                          % used for classifying results
fieldName = 'ClassMEF';                                                 % used for classifying results


% OVER-ARCHING OUTPUT FILENAME FOR EXPERIMENT
exprFileName = 'MEFcutoff.xlsx'; 

% LOAD MODEL
initCobraToolbox(false);
load('Recon3DModel_301.mat');
if ~exist('model')
    model = Recon3DModel;
end
modelName = 'Recon3DModel';
      
% MODIFY BASIC MODEL
model = Rtotal_patch(model);                    % patch from doi: https://doi.org/10.1101/674150

% now remove 'irrelevant' reactions that create spurious fluxes
model = removeMetabolites(model, {'ggn[c]', 'HC02203[c]', 'HC02205[c]', 'HC02207[c]', 'prostge2[c]'}, true);
model = removeRxns(model, {'r0355', 'DXTRNt','GLDBRAN', 'GLPASE2'});

%
biomassFileName = 'bm_Recon2.txt';              % option to add alternate biomass function
bm = readtable(biomassFileName);
model = bmRead2(model, biomassFileName); 
FBA_dummy = FBAdummy(model);                    % this is a zero'ed FBA solution for use where there is no solution
modelOri=model;

% Subset the main expression dataset for those rows with genes in Rec
gem = readtable('MEFexp_rows.xlsx'); % this includes NaN's and a column 'rows' that contains references to rows in 'data'
gem = gem(1:2248,:);                 % this is required to get rid of NaN row at 2249
geneRows = gem.rows;

% MAIN LOOPS
    
for SqueezeCount = 1:numel(Squeeze)
    
    if Squeeze(SqueezeCount) ~= -1 % saturate expression levels according to 'squeeze'
        dataPCA = table2array(dataPC);
        dataPC_compress = 40*dataPCA(:,:)./(dataPCA(:,:)+Squeeze(SqueezeCount)); % Max set empirically
        % Create a percentage version
        dataPC_compress = 100*bsxfun(@rdivide,dataPC_compress(:,:),sum(dataPC_compress(:,:),1,'omitnan'));

    else   % use linear relationship between expression and magnitude of upper and lower bounds
        dataPC_compress = table2array(dataPC);
        dataPC_compress = 100*bsxfun(@rdivide,dataPC_compress(:,:),sum(dataPC_compress(:,:),1,'omitnan'));
    end

    Exp_tab = dataPC_compress(geneRows,:);                      % subsets Recon genes from total dataset
    Exp_tab = array2table(Exp_tab, 'VariableNames', S);
    %%% S = Exp_tab.Properties.VariableNames;

    % SUBSET AND SUM THE ATPASE ROWS
    ATPases = readtable('EnerSysGO kinetic data.xlsx'); % read in rules
    [gene, kcat, ~, ~] = ATPaseRules(ATPases, dataRaw, dataPC, S{1}); % call simply to get gene list and kcat
    atpaseRows = find(ismember(dataRaw.NEWSymbol,gene));
    atpaseArray = dataPC_compress(atpaseRows,:);
    atpaseSum = sum(atpaseArray,1);
    atpaseTab = array2table(atpaseSum, 'VariableNames', S);

%
%       CONVERT GENE EXPRESSION LEVELS TO RELATIVE REACTION FLUX SIZES
%       ACROSS SELECTED GROUPS GENERATING THE MATRIX fluxMaxTab
    
        % mapExpressionToReactions returns -1 for orphan reactions and no gene
        % expression data generates an NaN that feeds through
        % mapExpressionToReactions

    exprData.gene = model.genes;
    for i = 1:size(Exp_tab,2)
        exprData.value = table2array(Exp_tab(:,i)); 
        [expressionRxns, parsedGPR] = mapExpressionToReactions(model, exprData, 'minSum');
        fluxMaxTab(:,i) = expressionRxns;
    end
  
    Emax = max(fluxMaxTab, [],'all', 'omitnan'); % get maximum expression value entire dataset so all samples scaled in the same way below
    
    for mult = Multiplier

    %
    %   SAMPLE LOOP This is the main FBA loop. A single set of parameters used
    %   for all samples and the output saved as multiple worksheets in a
    %   single Excel workbook
    %

        % SET UP OUTPUT FILENAMES
        outFileName = strcat(exprFileName(1:end-5), '_dmATP_', string(mult), ' ', datestr(now, 'dd-mmm-yyyyHH-MM-SS'),'.xlsx');
        outFileName = convertStringsToChars(outFileName);
        outFileName_run = strcat(outFileName(1:end-5),'.mat');                  % this is used to store workspace
        outFileName_summary = strcat(exprFileName(1:end-5), 'summ.xlsx');       % this is main summary table

        % SET UP VARIABLES
        count = 0;              % this tracks total iterations for output purposes
        rxnListLBTot = [];
        rxnListUBTot =[];

        for sampleNo = 1:size(S,2) % cycle through all the separate experimental samples                               

            count = count+1;                                % increment count


            % SETUP MODEL PARAMETERS
            % 1) GENE EXPRESSION - LOWER AND UPPER BOUNDS
            gexp=fluxMaxTab(:,sampleNo);                % take  single expression col

            gexp(isnan(gexp))=Emax;                     % NaN in expression list. Should never happen in this version
            gexp(gexp==-1)=Emax;                        % orphan reactions
            gexp(gexp==0)=Emax;                         % is this a genuine zero?   
            gexp_scale = gexp/Emax;                     % set up scaling vector


            modelScaled = model;                           % set up scaled model and don't modify model
            modelScaled.lb = modelScaled.lb.*gexp_scale/10;   % increase headroom by factor of 10
            modelScaled.ub = modelScaled.ub.*gexp_scale/10;

            % 2) CLOSE EXTERNAL REACTIONS
            modelClosed = modelClose(modelScaled);


            % 3) OPEN EXTERNAL REACTION LOWER BOUNDS TO ALLOW MEDIA COMPONENTS IN
            EXRs_tab = readtable(boundsFileName);
            modelOpen = modelClosed;
            for k = 1:size(EXRs_tab, 1)
                index = find(ismember(modelOpen.rxns, EXRs_tab{k,1}));
                modelOpen.lb(index) = EXRs_tab{k,2}*mult;
                modelOpen.ub(index) = EXRs_tab{k,3}*mult;
            end  

            % 4) SET OBJECTIVE FUNCTION
            % remember to open bounds if they have already been closed
            modelOpen = changeRxnBounds(modelOpen,'DM_atp_c_', -1000,'l');
            modelOpen = changeRxnBounds(modelOpen,'DM_atp_c_', 1000,'u');
            modelOpen = changeObjective(modelOpen,'DM_atp_c_');

            % 5) OPTION TO SET METABOLITE EXCHANGE CALIBRATED AGAINST DEFAULT DMEM
            % GLUCOSE

            if includeRealFlux == 1
                METs_tab = readtable('MetabExchange.txt');
                METs_tab{:,2:5} = METs_tab{:,2:5}.*(-25/METs_tab{1,2});  % scales to WT glucose uptake of 25       
                METs_tab.Properties.VariableNames = ['rxns' S]; % replace the column names with those in use
                for k = 1:size(METs_tab, 1)
                    index = find(ismember(modelClosed.rxns, METs_tab{k,1}));
                    modelOpen.lb(index) = METs_tab{k,S(sampleNo)}*mult;
                    modelOpen.ub(index) = METs_tab{k,S(sampleNo)}*mult;
                end  
            end 

            % Get fit and store in FBA    
            FBA_temp = optimizeCbModel(modelOpen,'max','zero'); % used to generate list of param outputs
            if FBA_temp.stat == 1 
                FBA(:,count) = optimizeCbModel(modelOpen,'max','zero'); 
            else
                FBA(:,count) = FBA_dummy;
            end

            % Find max'ed reactons
            [rxnListLB, rxnListUB] = maxedRxn2(modelOpen,FBA(1,count));
            rxnListLBTot = [rxnListLBTot rxnListLB];
            rxnListUBTot = [rxnListUBTot rxnListUB];

            % Store FBA results 
            FBA_tab = table(FBA(:,count).x, 'RowNames', modelOpen.rxns);

            % Store model
           modelStore(:,count) = modelOpen;

        end % sample loop

    % NOW OUTPUT TO MULTIPLE WORKBOOK SPREADSHEETS
    % Generate summary table and write to file
    FBA_tab = table(FBA(:,:).x, 'RowNames', modelOpen.rxns);
    writetable(FBA_tab,outFileName, 'Sheet', 'FBAtab', 'WriteRowNames',true); % save flux estimates table

    % Combine results from different cell lines. Note duplicates remain
    rxnListLBTot = vertcat(rxnListLBTot{:});
    rxnListUBTot = vertcat(rxnListUBTot{:});


    % Get indices for upper and lower bounds of max'ed reactions
    indexLB = [];
    indexUB = [];

    for i = 1: numel(rxnListLBTot)
        indexLB = [indexLB, find(strcmp(modelOpen.rxns,rxnListLBTot{i}))];
    end

    for i = 1: numel(rxnListUBTot)
        indexUB = [indexUB, find(strcmp(modelOpen.rxns,rxnListUBTot{i}))];
    end

    index = unique([indexLB, indexUB]);

    % Take copy of modelStore and convert to VMH metabolite names
    modelStoreOut = modelStore;

    % identify EX reactions
    exTabIndex = startsWith(FBA_tab.Properties.RowNames, 'EX_');
    exTabNZ = max(abs(FBA_tab{:,:}),[],2);
    exTabIndexNZ = and(exTabNZ,exTabIndex);
    exTab = FBA_tab(exTabIndexNZ,:);
    exLB = modelOpen.lb(exTabIndexNZ,:);
    exUB = modelOpen.ub(exTabIndexNZ,:);
    exTab = addvars(exTab, exLB, exUB);
    writetable(exTab,outFileName, 'Sheet', 'ExRxn', 'WriteRowNames',true);

    outTab = rxnTable(modelStoreOut, FBA, index, S);
    writetable(outTab,outFileName, 'Sheet', 'MaxRxn', 'WriteRowNames',true);
    writetable(array2table(dataPC_compress), outFileName, 'Sheet', 'compressData', 'WriteRowNames',true);
    writetable(atpaseTab, outFileName, 'Sheet', 'ATPases', 'WriteRowNames',true);


    % Now metabolite tables
    outTab = writeMet(outFileName, modelStoreOut, FBA, S, 'glc_D', 'glc_D_t');
    outTab = writeMet(outFileName, modelStoreOut, FBA, S, 'gln_L', 'gln_L');
    outTab = writeMet(outFileName, modelStoreOut, FBA, S, 'atp', 'atp');
    outTab = writeMet(outFileName, modelStoreOut, FBA, S, 'gtp', 'gtp');
    outTab = writeMet(outFileName, modelStoreOut, FBA, S, 'nad', 'nad');
    outTab = writeMet(outFileName, modelStoreOut, FBA, S, 'nadp', 'nadp');
    outTab = writeMet(outFileName, modelStoreOut, FBA, S, 'pyr', 'pyr');

    end % multiplier loop

end % squeeze loop



