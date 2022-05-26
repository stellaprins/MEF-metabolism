function model_out = bmRead2(model, filename)
%   bmRead 
%   takes a file with metabolites in column 1 and coefficients in column 2
%   It returns a model with an extra biomass reaction called biomass_2[c]
%   NOTE: needs checks for presence of metabolites in 'model' and no NaN in
%   column 2
%   NOTE: could add additional argument for new biomass equation name


bm = readtable(filename);

% bm{:,2} = -bm{:,2};

sizebm = size(bm);
metab_string = '';
coeff_string = '';
for i = 1:sizebm(1)
    metab_string = strcat(metab_string, " ", "'", string(bm{i,1}), "' ,");
    coeff_string = strcat(coeff_string, " ", sprintf('%f', bm{i,2}), ", ");
end

metab_string = char(metab_string);
metab_string = metab_string(1:end-1);
coeff_string = char(coeff_string);
coeff_string = coeff_string(1:end-2);

biomass_rxn_string = strcat("addReaction(model, 'biomass_2[c]', 'metaboliteList', {", metab_string, "}, 'stoichCoeffList', [",  coeff_string, " ]);");

model_out = eval(biomass_rxn_string);

end

