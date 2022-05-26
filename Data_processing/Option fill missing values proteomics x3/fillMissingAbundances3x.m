% Script with following functionality

% Takes a proteomic dataset, linearises data if necessary, reads GEM gene
% list with cross-referenced row numbers for proteomic dataset, fills in
% missing proteomic data using protein process groups

% Essentially preprocessing for FBA analysis

% 6 September 2021 Reorganisation to give summary tables of results.
% 23 October 2021 Correct problem with mean(temp) = NaN if all values are
% NaNs

% 15 November 2021 version to read triplicate datasets

% Plenty of scope to improve execution speed. Currently v slow
% Also generates a lot of data that aren't required for current downstream
% analysis
% Latest review 10 May 2022

%   .. Author: - Phil Luthert 10/5/22

clear

% [1] Proteomic data entry and pre-processing of log transformed expression
% data

data = readtable('Wang19300genes3x.xlsx');
dataRaw = data;

% [1a] Subset data where appropriate to speed execution time

data = data(:,{'TissueSpecificExpression','Process_1_','Process_2_','Process_3_','WT_1_','WT_2_','WT_3_', 'G12D_1_','G12D_2_','G12D_3_','G12V_1_','G12V_2_','G12V_3_','Q61L_1_','Q61L_2_','Q61L_3_'});
sample = data.Properties.VariableNames(5:16);

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
                                data{count,sampleName} = process1Min;
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

% Create a percentage version of the dataAtlas
dataAtlasPC = 100*bsxfun(@rdivide,dataAtlas{:,sample},sum(dataAtlas{:,sample},1,'omitnan'));
dataAtlasPC = array2table(dataAtlasPC, 'VariableNames', string(sample));

% Creates new columns for kcat and modified abundances (which will only
% change for ATPases.
% There will be much faster ways of doing this.
dataKcat = dataAtlasPC;       % just do for dataAtlasPC at least initially

% Add columns to dataKcat. 
dummyCol = zeros(size(dataKcat,1),1); 

sampleID = dataAtlas.Properties.VariableNames(5:16); % to include all abundance data columns
for sampleCount = 1:size(sampleID, 2)
    colName = strcat(sampleID{sampleCount}, 'mod');
    dataKcat = addvars(dataKcat,dummyCol);
    dataKcat.Properties.VariableNames{'dummyCol'} = colName;
end
dataKcat = addvars(dataKcat,dummyCol);
dataKcat.Properties.VariableNames{'dummyCol'} = 'kcat';

% read table to calculate median kcat
ATPases = readtable('EnerSysGO kinetic data.xlsx'); % read in rules
[gene, kcat, ~, ~] = ATPaseRules(ATPases, dataRaw, dataAtlasPC, sampleID{1}); % call simply to get gene list and kcat
kcatMedian = median(kcat);

% go through gene table line by line and populate new table. This is very
% slow and could be faster written
for geneCount = 1:size(dataKcat,1)
  % first complete kcat column
  if ~ismember(dataRaw.NEWSymbol{geneCount}, gene)  
    dataKcat.kcat(geneCount) = kcatMedian;
    % all abundances stay as they were
    for sampleCount = 1:size(sampleID, 2)
        colName = strcat(sampleID{sampleCount}, 'mod');
        dataKcat.(colName)(geneCount) = dataKcat.(sampleID{sampleCount})(geneCount);
    end
  else     
      kcatMod = lookup(gene,dataRaw.NEWSymbol{geneCount}, kcat);
      dataKcat.kcat(geneCount) = kcatMod;
      for sampleCount = 1:size(sampleID, 2)
        [gene, kcat, abundance, abundanceStarKcat] = ATPaseRules(ATPases, dataRaw, dataAtlasPC, sampleID{sampleCount});
        abundanceMod = lookup(gene,dataRaw.NEWSymbol{geneCount}, abundance);
        colName = strcat(sampleID{sampleCount}, 'mod');
        dataKcat.(colName)(geneCount) = cell2mat(abundanceMod);
      end
      kcatMod = lookup(gene,dataRaw.NEWSymbol{geneCount}, kcat);    % kcat needs updating only once
      dataKcat.kcat(geneCount) = kcatMod;
  end
end

% Now compute the abundance Kcat products

for sampleCount = 1:size(sampleID, 2)                   
    colName = strcat(sampleID{sampleCount}, 'AsK');                 % set up column names
    dataKcat = addvars(dataKcat,dummyCol);
    dataKcat.Properties.VariableNames{'dummyCol'} = colName;
    sourceName = strcat(sampleID{sampleCount}, 'mod'); 
    dataKcat.(colName) = dataKcat.kcat.* dataKcat.(sourceName);
end


save('dataRaw','dataRaw');        % required for MEF workflow
save('dataKcat','dataKcat');      % required for MEF workflow

