function [dataset_proteomics_phenotype] = extract_proteomics_kinases(file_name)
file_proteomics_kinases = file_name;
opts = detectImportOptions(file_proteomics_kinases);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_proteomics_kinases = readtable(file_proteomics_kinases, opts);

ind = find(contains(file_proteomics_kinases.Properties.VariableNames,'HumanDonor M1'));
ind1 = find(ismember(file_proteomics_kinases{:,ind},'True'));
dataset_proteomics_phenotype{1} = file_proteomics_kinases{ind1,1};
ind = find(contains(file_proteomics_kinases.Properties.VariableNames,'HumanDonor M2a expressed'));
ind1 = find(ismember(file_proteomics_kinases{:,ind},'True'));
dataset_proteomics_phenotype{2} = file_proteomics_kinases{ind1,1};
ind = find(contains(file_proteomics_kinases.Properties.VariableNames,'HumanDonor M2c expressed'));
ind1 = find(ismember(file_proteomics_kinases{:,ind},'True'));
dataset_proteomics_phenotype{3} = file_proteomics_kinases{ind1,1};
end
