% fillMissingAbundances

% Takes a proteomic dataset, linearises data if necessary, reads GEM gene
% list with cross-referenced row numbers for proteomic dataset, fills in
% missing proteomic data using protein process groups

% Essentially preprocessing for FBA analysis

% 6 September 2021 Reorganisation to give summary tables of results.
% 23 October 2021 Correct problem with mean(temp) = NaN if all values are
% NaNs
% 19 May 2022 simplify around MEF data

clear

% [1] Proteomic data entry 

data = readtable('Wang19300genes.xlsx');
dataRaw = data;

% [1a] Subset data where appropriate 

data = data(:,{'TissueSpecificExpression','Process_1_','Process_2_','Process_3_','AverageWT','AverageG12D','AverageG12V','AverageQ61L'});
sample = data.Properties.VariableNames(5:8);

% [2] Generate summary data for expression levels of proteins linked to
% specific call processes that can be used for filling in missing values in
% the proteomic data

dataOri = data;

P1=[]; % structures to contain summary data about 3 hierarchical levels of process
P2=[];
P3=[];

% pull out individual processes 
P1.list = unique(data.Process_1_);
P2.list = unique(data.Process_2_);
P3.list = unique(data.Process_3_);

% include index for each
for count = 1:numel(P1.list)
P1.idx{count} = find(ismember(data.Process_1_, P1.list{count}));
end

for count = 1:numel(P2.list)
P2.idx{count} = find(ismember(data.Process_2_, P2.list{count}));
end

for count = 1:numel(P3.list)
P3.idx{count} = find(ismember(data.Process_3_, P3.list{count}));
end

% calculate mean and min value for each process for each sample
sampleNo = 1; 
for sampleName = sample
    for count = 1:numel(P1.list)
        temp = table2array(data(cell2mat(P1.idx(count)),sampleName));
        if ~isnan(mean(temp, 'omitnan'))    % test for all temp = NaN in which case risk downstream error
            P1.(strcat(sample{sampleNo},'_Mean'))(count) = mean(temp, 'omitnan');
            P1.(strcat(sample{sampleNo},'_Min'))(count) = min(temp, [], 'omitnan');
        else
            P1.(strcat(sample{sampleNo},'_Mean'))(count) = 0;
            P1.(strcat(sample{sampleNo},'_Min'))(count) = 0;
        end
    end
    sampleNo = sampleNo+1;
end

sampleNo = 1;
for sampleName = sample
    for count = 1:numel(P2.list)
        temp = table2array(data(cell2mat(P2.idx(count)),sampleName));
        if ~isnan(mean(temp, 'omitnan'))    % test for all temp = NaN in which case risk downstream error
            P2.(strcat(sample{sampleNo},'_Mean'))(count) = mean(temp, 'omitnan');
            P2.(strcat(sample{sampleNo},'_Min'))(count) = min(temp, [], 'omitnan');
        else
            P2.(strcat(sample{sampleNo},'_Mean'))(count) = 0;
            P2.(strcat(sample{sampleNo},'_Min'))(count) = 0;
        end
    end
    sampleNo = sampleNo+1;
end

sampleNo = 1;
for sampleName = sample
    for count = 1:numel(P3.list)
        temp = table2array(data(cell2mat(P3.idx(count)),sampleName));
        if ~isnan(mean(temp, 'omitnan'))    % test for all temp = NaN in which case risk downstream error
            P3.(strcat(sample{sampleNo},'_Mean'))(count) = mean(temp, 'omitnan');
            P3.(strcat(sample{sampleNo},'_Min'))(count) = min(temp, [], 'omitnan');
        else
            P3.(strcat(sample{sampleNo},'_Mean'))(count) = 0;
            P3.(strcat(sample{sampleNo},'_Min'))(count) = 0;
        end
    end
    sampleNo = sampleNo+1;
end

% Overall approach now to cycle through values in data table and if NaN
% replace with mean or min for the group. Need nested structure so that if
% one level of process has no value (NaN) we go to the next in order P1, P2
% P3.
% If there is no P3 print to screen no expression value.


data = dataOri;

sampleNo = 1;
for sampleName = sample
    exprList = data(:,sampleName);
    for count = 1:numel(exprList)
        exprLevel = exprList{count,1};
        if isnan(exprLevel)
            testP1 = data.Process_1_{count};                % read off process name
            process1ID = find(ismember(P1.list, testP1));    % find the process ID within the struct
            process1Mean = P1.(strcat(sample{sampleNo},'_Mean'))(process1ID);
            if isnan(process1Mean)
                disp('P1 isnan');                           % place holder for next bit of code
                testP2 = data.Process_2_{count}; 
                process2ID = find(ismember(P2.list, testP2));    % find the process ID within the struct
                process2Mean = P2.(strcat(sample{sampleNo},'_Mean'))(process2ID);
                if isnan(process2Mean)
                    disp('P2 isnan');                           % place holder for next bit of code
                    testP3 = data.Process_3_{count}; 
                    process3ID = find(ismember(P3.list, testP3));    % find the process ID within the struct
                    process3Mean = P3.(strcat(sample{sampleNo},'_Mean'))(process3ID);
                    if isnan(process3Mean)
                        disp('P3 isnan, unassignable expression value'); 
                    else
                        data{count,sampleName} = process3Mean;
                    end
                else
                    data{count,sampleName} = process2Mean;
                end
            else
                data{count,sampleName} = process1Mean;
            end     
        end    
    end
    sampleNo = sampleNo+1;
end

dataMean = data;

data = dataOri;

sampleNo = 1;
for sampleName = sample
    exprList = data(:,sampleName);
    for count = 1:numel(exprList)
        exprLevel = exprList{count,1};
        if isnan(exprLevel)
            testP1 = data.Process_1_{count};                % read off process name
            process1ID = find(ismember(P1.list, testP1));    % find the process ID within the struct
            process1Min = P1.(strcat(sample{sampleNo},'_Min'))(process1ID);
            if isnan(process1Min)
                disp('P1 isnan');                           % place holder for next bit of code
                testP2 = data.Process_2_{count}; 
                process2ID = find(ismember(P2.list, testP2));    % find the process ID within the struct
                process2Min = P2.(strcat(sample{sampleNo},'_Min'))(process2ID);
                if isnan(process2Min)
                    disp('P2 isnan');                           % place holder for next bit of code
                    testP3 = data.Process_3_{count}; 
                    process3ID = find(ismember(P3.list, testP3));    % find the process ID within the struct
                    process3Min = P3.(strcat(sample{sampleNo},'_Min'))(process3ID);
                    if isnan(process3Min)
                        disp('P3 isnan, unassignable expression value'); 
                    else
                        data{count,sampleName} = process3Min;
                    end
                else
                    data{count,sampleName} = process2Min;
                end
            else
                data{count,sampleName} = process1Min;
            end     
        end    
    end
    sampleNo = sampleNo+1;
end

dataMin = data;

data = dataOri;

% And now select min or mean depending on status of Protein Atlas column
% ('TissueSpecificExpression')

sampleNo = 1;
for sampleName = sample
    exprList = data(:,sampleName);
    for count = 1:numel(exprList)
        exprLevel = exprList{count,1};
        tissueSpecificExpression = string(data{count,'TissueSpecificExpression'}); % read off tissue specific expression classification
        if isnan(exprLevel)
            testP1 = data.Process_1_{count};                % read off process name
            process1ID = find(ismember(P1.list, testP1));    % find the process ID within the struct
            process1Min = P1.(strcat(sample{sampleNo},'_Min'))(process1ID);
            process1Mean = P1.(strcat(sample{sampleNo},'_Mean'))(process1ID);
            if isnan(process1Min)
                disp('P1 isnan');                           % place holder for next bit of code
                testP2 = data.Process_2_{count}; 
                process2ID = find(ismember(P2.list, testP2));    % find the process ID within the struct
                process2Min = P2.(strcat(sample{sampleNo},'_Min'))(process2ID);
                process2Mean = P2.(strcat(sample{sampleNo},'_Mean'))(process2ID);
                if isnan(process2Min)
                    disp('P2 isnan');                           % place holder for next bit of code
                    testP3 = data.Process_3_{count}; 
                    process3ID = find(ismember(P3.list, testP3));    % find the process ID within the struct
                    process3Min = P3.(strcat(sample{sampleNo},'_Min'))(process3ID);
                    process3Mean = P3.(strcat(sample{sampleNo},'_Mean'))(process3ID);
                    if isnan(process3Min)
                        disp('P3 isnan, unassignable expression value'); 
                    else
                        switch tissueSpecificExpression
                            case "Detected in all"
                                data{count,sampleName} = process3Mean;
                            case "Detected in many"
                                data{count,sampleName} = process3Min;
                            case "Detected in some"
                                data{count,sampleName} = process3Min;
                            case "Detected in single"
                                data{count,sampleName} = process3Min;
                            case "Not detected"
                                data{count,sampleName} = 0;
                            otherwise
                                data{count,sampleName} = process3Min;
                        end     
                    end
                else
                    switch tissueSpecificExpression
                        case "Detected in all"
                            data{count,sampleName} = process2Mean;
                        case "Detected in many"
                            data{count,sampleName} = process2Min;
                        case "Detected in some"
                            data{count,sampleName} = process2Min;
                        case "Detected in single"
                            data{count,sampleName} = process2Min;
                        case "Not detected"
                            data{count,sampleName} = 0;
                        otherwise
                            data{count,sampleName} = process2Min;
                    end     
                end
            else
                switch tissueSpecificExpression
                    case "Detected in all"
                        data{count,sampleName} = process1Mean;
                    case "Detected in many"
                        data{count,sampleName} = process1Min;
                    case "Detected in some"
                        data{count,sampleName} = process1Min;
                    case "Detected in single"
                        data{count,sampleName} = process1Min;
                    case "Not detected"
                        data{count,sampleName} = 0;
                    otherwise
                        data{count,sampleName} = process1Min;
                end     
            end     
        end    
    end
    sampleNo = sampleNo+1;
end

dataAtlas = data;
data = dataOri;

% Code here as an inset if want to use abundances rules for ATPases and
% include consideration of kcat values

% Create a percentage version of the dataAtlas
dataAtlasPC = 100*bsxfun(@rdivide,dataAtlas{:,sample},sum(dataAtlas{:,sample},1,'omitnan'));
dataAtlasPC = array2table(dataAtlasPC, 'VariableNames', string(sample));

% Creates new columns for kcat and modified abundances (which will only
% change for ATPases.
dataKcat = dataAtlasPC;       % just do for dataAtlasPC at least initially

% read table to calculate median kcat
ATPases = readtable('EnerSysGO kinetic data.xlsx'); % read in rules
[gene, kcat, ~, ~] = ATPaseRules(ATPases, dataRaw, dataAtlasPC, sample{1}); % call simply to get gene list and kcat

% Duplicate kcat columns to later overwrite in modified abundance values (mod)
dataKcat(:,5:8) = dataAtlasPC;
dataKcat(:,9)   = array2table(zeros(size(dataKcat,1),1)); 
dataKcat.Properties.VariableNames(5:8) = strcat(sample,'mod');
dataKcat.Properties.VariableNames(9) = {'kcat'};

% when kcat value is not available assign medium kcat (=1)
dataKcat.kcat(~ismember(dataRaw.NEWSymbol, gene)) = median(kcat); 

% assign kcat values for genes in 'EnerSysGO kinetic data.xlsx'
for g = 1 : length(gene)
    pos = find(ismember(dataRaw.NEWSymbol, gene(g)));
    dataKcat.kcat(pos) = kcat(g);
end

% assign modified abundance values
for s = 1:size(sample, 2)
   [gene, kcat, abundance, abundanceStarKcat] = ATPaseRules(ATPases, dataRaw, dataAtlasPC, sample{s});
    for g = 1 : length(gene);
        pos = find(ismember(dataRaw.NEWSymbol, gene(g)));
        dataKcat(pos,length(sample)+s)= abundance(g);
    end
end

% compute the abundance Kcat products (AsK)
AsK = dataKcat.kcat.*table2array(dataKcat(:,strcat(sample,'mod')));
dataKcat(:,10:13) = array2table(AsK) ;
dataKcat.Properties.VariableNames(10:13) = strcat(sample,'AsK');

save('dataRaw','dataRaw');
save('dataKcat','dataKcat');
