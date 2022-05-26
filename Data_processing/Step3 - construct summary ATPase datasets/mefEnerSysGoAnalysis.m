clear;

% [1] GENERATE SETS OF COMPRESSED DATASETS WITHOUT FBA
% A modification of the main FBA - script

% Outcome is 'Outarray':
    % Rows = locations
    % Cols = classes of ATPase
    % Z1 = different cell lines
    % Z2 = different levels of compression
    % Z3 = group with base, ATPase or ATPase - Kcat rules
    
% Analysis carried out separately .
% That is 3 datasets will be generated

load('dataRaw');
dataRaw.PROCESSNEW27082021{5428} = '008_Transporter ions'; % this sorts error in original file
dataRaw.PROCESSNEW27082021{5429} = '008_Transporter ions'; % this sorts error in original file
dataTmp = dataRaw;      % use this as has some of the key columns that are lost in prior processing
dataTmp(19301:19303,:)=[];
dataTmp = removevars(dataTmp, {'AverageWT','AverageG12D','AverageG12V','AverageQ61L'}); % this avoids confusion if names are duplicated in dataPC

load('dataKcat');
dataPC = dataKcat; % this is the version of the expression dataset that includes ATPase rules and kcat terms
% lose the kcat column
dataPCcomplete = removevars(dataPC,{'kcat'});
dataPCcomplete(19301:19303,:)=[]; % need to lose dummy rows etc

matrix = readtable('Process_reactions_working2.xlsx'); % this is the matrix with the defining row and column labels to summarise ATPase abundances

% some constants / flags
Squeeze = [0.0001 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.25 0.5 -1];

% OVER-ARCHING OUTPUT FILENAME FOR EXPERIMENT
exprFileName = 'MEFcutoff.xlsx'; 

fileTags = {'base','rules','kcat'}; % string to put into filenames giving condition
samples = {'WT','G12D','G12V','Q61L'};
outArray = zeros(size(matrix,1), size(matrix, 2)-1, size(dataPCcomplete,2)/size(fileTags,2), size(Squeeze, 2), size(fileTags,2));
% MAIN LOOPS

for groupNo = 0:2 % this will take us through normal / ATPase rules / ATPase rules + kcat  
    
    % setup data
    dataPC = dataPCcomplete(:,groupNo*4+1:groupNo*4+4);
    S = dataPC.Properties.VariableNames(1:4); 
        
    for SqueezeCount = 1:numel(Squeeze)
        
        % Compression steps
        if Squeeze(SqueezeCount) ~= -1 % -1 = no compression
            dataPCA = table2array(dataPC);     % convery to array
            dataPC_compress = 40*dataPCA(:,:)./(dataPCA(:,:)+Squeeze(SqueezeCount)); % Saturate expression levels according to 'squeeze'Max set empirically
            dataPC_compress = 100*bsxfun(@rdivide,dataPC_compress(:,:),sum(dataPC_compress(:,:),1,'omitnan')); % Create a percentage version            
        else
            dataPC_compress = table2array(dataPC);
            dataPC_compress = 100*bsxfun(@rdivide,dataPC_compress(:,:),sum(dataPC_compress(:,:),1,'omitnan'));
        end
        dataPCtab = array2table(dataPC_compress, 'VariableNames', S);
        % Set up combined dataTab for current group
        dataTab = addvars(dataTmp, dataPCtab );          % combine data
        dataTab = splitvars(dataTab);                    % split headings

        % Inner loop to generate data 
        count = 0;                                       % required for indexing
        for varName = S
            count = count+1;
            matrixOut = mapExpression2Matrix(dataTab,matrix,varName);
            outArray(:,:,count,SqueezeCount,groupNo+1) = table2array(matrixOut);
        end    % varName
    end
end 

% outArray now available for further analyses

% [2] Revision of initial mefAtpASEAnalysis script to embrace more explicitly GTP as well as ATP

% generate tables 
ATPbool = contains(matrixOut.Properties.VariableNames, 'ATP');
NADbool = contains(matrixOut.Properties.VariableNames, 'NAD');
GTPbool = contains(matrixOut.Properties.VariableNames, 'GTP');
CTPbool = contains(matrixOut.Properties.VariableNames, 'CTP');

outDatProcess = sum(outArray,1, 'omitnan'); % this sums columns to generate list of processes
outDatProcess = squeeze(outDatProcess);

outDatATP = outDatProcess(ATPbool,:,:,:);
outDatNGTPATP = outDatProcess(~ATPbool & ~GTPbool,:,:,:);
outDatNAD = outDatProcess(NADbool,:,:,:);
outDatGTP = outDatProcess(GTPbool,:,:,:);
outDatCTP = outDatProcess(CTPbool,:,:,:);

% CONSTRUCT ATPase and GTPase SUMMARY with nonGTPase and ATPases totalled as final row

% ATPase summary ...
% place base, ATPase rules and kcat datasets side by side.
out1 = zeros(size(outDatATP,1), size(outDatATP,3)*size(outDatATP,4),size(outDatATP,2));
for i = 1:4
    out1(:,:,i) = [squeeze(outDatATP(:,i,:,1)) squeeze(outDatATP(:,i,:,2)) squeeze(outDatATP(:,i,:,3))];
end

% GTPase summary ...
out2 = zeros(size(outDatGTP,1), size(outDatGTP,3)*size(outDatGTP,4),size(outDatGTP,2));
for i = 1:4
    out2(:,:,i) = [squeeze(outDatGTP(:,i,:,1)) squeeze(outDatGTP(:,i,:,2)) squeeze(outDatGTP(:,i,:,3))];
end

% nonATPase totals
NGTPATPtotals = squeeze(sum(outDatNGTPATP,1));
 
out3 = zeros(size(NGTPATPtotals,2)*size(NGTPATPtotals,3),size(NGTPATPtotals,1));
for i = 1:4
    out3(:,i) = [squeeze(NGTPATPtotals(i,:,1)) squeeze(NGTPATPtotals(i,:,2)) squeeze(NGTPATPtotals(i,:,3))];
end

out = vertcat(out1, out2);
% add final row
out((size(out1,1)+size(out2,1)+1),:,:)=out3;

% now 'out' needs to be turned into tables which are written to successive
% sheets in an excel spreadsheet

outRowNames = [matrixOut.Properties.VariableNames(ATPbool) matrixOut.Properties.VariableNames(GTPbool) 'nonGTPATPases'];

outVarNames = {};
for i = 1:3
    for j = 1:size(out,2)/3
       outVarNames = [outVarNames strcat(fileTags{i}, num2str(Squeeze(1,j)))];
    end
end

for i = 1:4
    writetable(array2table(out(:,:,i), 'VariableNames', outVarNames, 'RowNames', outRowNames), 'ATPGTPaseSummary.xlsx', 'sheet', samples{i}, 'WriteRowNames',true);
end 

% Now the same for the nonGTPATPases

% place base, ATPase rules and kcat datasets side by side.
outNGTPATP = zeros(size(outDatNGTPATP,1), size(outDatNGTPATP,3)*size(outDatNGTPATP,4),size(outDatNGTPATP,2));
for i = 1:4
    outNATP(:,:,i) = [squeeze(outDatNGTPATP(:,i,:,1)) squeeze(outDatNGTPATP(:,i,:,2)) squeeze(outDatNGTPATP(:,i,:,3))];
end

outRowNames = matrixOut.Properties.VariableNames(~ATPbool & ~GTPbool);

% outVarNames as above

for i = 1:4
   writetable(array2table(outNATP(:,:,i),  'VariableNames', outVarNames, 'RowNames', outRowNames), 'NGTPATPaseSummary.xlsx', 'sheet', samples{i}, 'WriteRowNames',true);
end


