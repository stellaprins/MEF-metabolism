function modelx = Rtotal_patch(model)
%  Rtotal_patch modifies model to facilitate biomass increase from
%  tissue culture media
% 
%  Modified from code linked to bioarxiv.org: 

% Simultaneous integration of gene expression and nutrient availability for studying metabolism of hepatocellular carcinoma
% Ewelina Weglarz-Tomczak, Thierry D.G.A. Mondeel, Diewertje G.E. Piebes, Hans V. Westerhoff
% doi: https://doi.org/10.1101/674150


%% Patch Recon3D: combine Rtotal species + add synthesis step
% define metabolite IDs that need to change
map_to_Rtotal_c = {'Rtotal2[c]','Rtotal3[c]'};
map_to_Rtotal_e = {'Rtotal2[e]','Rtotal3[e]'};
map_to_Rtotalcoa_c = {'Rtotal2coa[c]','Rtotal3coa[c]'};
change = {'Rtotal3coa[m]','Rtotal3crn[c]','Rtotal3crn[m]'};

% rename some metabolite IDs
for met = change
    indx = findMetIDs(model,met);
    model.mets(indx) = strrep(met,'3','');
end

% Replace stoichiometries to new species
for met = map_to_Rtotal_c
    row = findMetIDs(model,met);
    if row == 0
        break
    end
    newrow = findMetIDs(model,'Rtotal[c]');
    
    active_reactions = find(model.S(row,:));
    stoich = model.S(row, active_reactions);
    model.S(row, active_reactions) = 0;
    model.S(newrow, active_reactions) = stoich;
end

for met = map_to_Rtotal_e
    row = findMetIDs(model,met);
    newrow = findMetIDs(model,'Rtotal[e]');
    
    active_reactions = find(model.S(row,:));
    stoich = model.S(row, active_reactions);
    model.S(row, active_reactions) = 0;
    model.S(newrow, active_reactions) = stoich;
end

for met = map_to_Rtotalcoa_c
    row = findMetIDs(model,met);
    newrow = findMetIDs(model,'Rtotalcoa[c]');
    
    active_reactions = find(model.S(row,:));
    stoich = model.S(row, active_reactions);
    model.S(row, active_reactions) = 0;
    model.S(newrow, active_reactions) = stoich;
end

% Add Rtotal synthesis
modelx = addReaction(model, 'Rtotal_synth', 'metaboliteList',{'stcoa[c]', 'pmtcoa[c]', 'odecoa[c]', 'lneldccoa[c]', 'Rtotalcoa[c]'},...
    'stoichCoeffList',[-1 -1 -1 -1 4], ...
    'reversible',true);

end