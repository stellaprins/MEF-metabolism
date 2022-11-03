function modelClosed = modelClose(model)
% modelClose() 
%   Takes a COBRA model and closes it
%   Uses code from COBRA tutorials
%   Also sets .c to zero. Doesn't set OF
%   10/1/2021 modify to close upper bounds as well as lower bounds

modelClosed = model;
    % prepare models for test - these changes are needed for the different
    % recon versions to match the rxn abbr definitions in this script
    %modelClosed.rxns = regexprep(modelClosed.rxns,'\(','\[');
    %modelClosed.rxns = regexprep(modelClosed.rxns,'\)','\]');
    %modelClosed.mets = regexprep(modelClosed.mets,'\(','\[');
    %modelClosed.mets = regexprep(modelClosed.mets,'\)','\]');
    %modelClosed.rxns = regexprep(modelClosed.rxns,'ATPS4mi','ATPS4m');

    %if length(strmatch('EX_glc[e]',modelClosed.rxns))>0
    %    modelClosed.rxns{find(ismember(modelClosed.rxns,'EX_glc[e]'))} = 'EX_glc_D[e]';
    %end
    % add reaction if it does not exist
    %[modelClosed, rxnIDexists] = addReaction(modelClosed,'DM_atp_c_',  'h2o[c] + atp[c]  -> adp[c] + h[c] + pi[c] ');
    %if length(rxnIDexists)>0
    %    modelClosed.rxns{rxnIDexists} = 'DM_atp_c_'; % rename reaction in case that it exists already
    %end



    % close all exchange and sink reactions (lb)
    modelexchanges1 = strmatch('biomass',modelClosed.rxns); 
    modelexchanges4 = strmatch('EX_',modelClosed.rxns);
    modelexchanges2 = strmatch('DM_',modelClosed.rxns);
    modelexchanges3 = strmatch('sink_',modelClosed.rxns);
    % also close biomass reactions
    %BM= (find(~cellfun(@isempty,strfind(lower(modelClosed.mets),'bioma'))));

    selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))'; % redundant in Recon3DModel

    modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
    modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges2))))=0; % close DM
    modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges3))))=0; % close sinks ub
    %modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges4))))=0;
    modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0; % close everything lb
    %modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
    modelClosed.c = zeros(length(modelClosed.rxns),1);
 
    modelClosed;

end

