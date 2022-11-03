clear;

% This is a modification of the main FBA - centric script simply to
% allow generation of complete sets of compressed datasets

% Outcome is array:
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

% To analyse by process type collapse 4D array to 3D summing [1]
% Take out, for instance, tRNA-AA ligases which are column 
% And then create heat map, cell line by compression where color is related
% to magnitude of rRNA-AA ligase abundance (= col 84)

groupNo = 2;

clear 'outdat'
outDat = sum(outArray,1); % this sums columns to generate list of processes
outDat = squeeze(outDat);
outDat = outDat(84,:,:,groupNo); % tRNA-AA ligases
outDat = squeeze(outDat);
image(outDat, 'CDataMapping','scaled');
colorbar

% Now analyse results above in terms of total ATPase abundance to determine
% % ATP flux passing through each subset of ATPases.

ATPbool = contains(matrixOut.Properties.VariableNames, 'ATP');

outDat2 = sum(outArray,1); % this sums columns to generate list of processes
outDat2 = squeeze(outDat2);
outDat2 = outDat2(ATPbool,:,:,:);
outDat2 = sum(outArray,2); 
outDat2 = squeeze(outDat2);
outDat2 = outDat2(1,:,:,groupNo);
outDat2 = squeeze(outDat2);
image(outDat2, 'CDataMapping','scaled');
colorbar

outDat3 = outDat2./outDat;
image(outDat3, 'CDataMapping','scaled');
colorbar

% now start separate analysis.
% generate tables 
ATPbool = contains(matrixOut.Properties.VariableNames, 'ATP');
NADbool = contains(matrixOut.Properties.VariableNames, 'NAD');
GTPbool = contains(matrixOut.Properties.VariableNames, 'GTP');
CTPbool = contains(matrixOut.Properties.VariableNames, 'CTP');

outDat3 = sum(outArray,1, 'omitnan'); % this sums columns to generate list of processes
outDat3 = squeeze(outDat3);

outDatATP = outDat3(ATPbool,:,:,:);
outDatNATP = outDat3(~ATPbool,:,:,:);
outDatNAD = outDat3(NADbool,:,:,:);
outDatGTP = outDat3(GTPbool,:,:,:);
outDatCTP = outDat3(CTPbool,:,:,:);

% CONSTRUCT ATPase SUMMARY with nonATPases totalled as final row

% ATPase summary ...
% place base, ATPase rules and kcat datasets side by side.
out = zeros(size(outDatATP,1), size(outDatATP,3)*size(outDatATP,4),size(outDatATP,2));
for i = 1:4
    out(:,:,i) = [squeeze(outDatATP(:,i,:,1)) squeeze(outDatATP(:,i,:,2)) squeeze(outDatATP(:,i,:,3))];
end

% nonATPase totals
NATPtotals = squeeze(sum(outDatNATP,1));
 
out2 = zeros(size(NATPtotals,2)*size(NATPtotals,3),size(NATPtotals,1));
for i = 1:4
    out2(:,i) = [squeeze(NATPtotals(i,:,1)) squeeze(NATPtotals(i,:,2)) squeeze(NATPtotals(i,:,3))];
end

% add final row
out(41,:,:)=out2;

% now 'out' needs to be turned into tables which are written to successive
% sheets in an excel spreadsheet

outRowNames = [matrixOut.Properties.VariableNames(ATPbool) 'nonATPases'];

outVarNames = {};
for i = 1:3
    for j = 1:size(out,2)/3
       outVarNames = [outVarNames strcat(fileTags{i}, num2str(Squeeze(1,j)))];
    end
end

for i = 1:4
    writetable(array2table(out(:,:,i), 'VariableNames', outVarNames, 'RowNames', outRowNames), 'ATPaseSummary2.xlsx', 'sheet', samples{i}, 'WriteRowNames',true);
end 

% Now the same for the nonATPases

% place base, ATPase rules and kcat datasets side by side.
outNATP = zeros(size(outDatNATP,1), size(outDatNATP,3)*size(outDatNATP,4),size(outDatNATP,2));
for i = 1:4
    outNATP(:,:,i) = [squeeze(outDatNATP(:,i,:,1)) squeeze(outDatNATP(:,i,:,2)) squeeze(outDatNATP(:,i,:,3))];
end

outRowNames = matrixOut.Properties.VariableNames(~ATPbool);

% outVarNames as above

for i = 1:4
    writetable(array2table(outNATP(:,:,i), 'VariableNames', outVarNames, 'RowNames', outRowNames), 'NATPaseSummary2.xlsx', 'sheet', samples{i}, 'WriteRowNames',true);
end

% THE ABOVE IS THE MAIN ANALYSIS. NOW SELECT DATA FOR MEF PAPER ANALYSIS
% Taking kcat rules and compression 0.001 which is at row 3
% outDatATP etc organised cols, lines, compression, rules so...

outDatATPSummary = squeeze(outDatATP(:,:,3,3));
