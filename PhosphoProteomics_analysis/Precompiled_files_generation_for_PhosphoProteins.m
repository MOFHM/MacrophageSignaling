function [] = Precompiled_files_generation_for_PhosphoProteins(root_NetPhorest,file_PhosphoProteomics,Interaction_database,NetPhorest_full_names,loc_prob,p_names,position_protein)
%Filter the data of Phosphoproteomics kinases by the available phosphorylated
%residues from PhosphoSite!
%clear all

opts = detectImportOptions(file_PhosphoProteomics);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_PhosphoProteomics = readtable(file_PhosphoProteomics, opts);

posterior = cellfun(@str2num,file_PhosphoProteomics.Var9);
prior = cellfun(@str2num,file_PhosphoProteomics.Var10);
ind = find(posterior<prior);
file_PhosphoProteomics(ind,:) = [];
posterior = cellfun(@str2num,file_PhosphoProteomics.Var9);
ind = find(posterior<0.035);
file_PhosphoProteomics(ind,:) = [];
ind = find(ismember(file_PhosphoProteomics.Var8,'any_group'));
file_PhosphoProteomics(ind,:) = [];

pat_ini_protein = 'sp|';
pat_fin_protein = '|';
file_PhosphoProteomics{:,1} = cellfun(@(x) extractBetween(x,pat_ini_protein,pat_fin_protein), file_PhosphoProteomics{:,1});
file_proteomics_residuesSorted = cellfun(@(x,y) [x,y],file_PhosphoProteomics{:,3},file_PhosphoProteomics{:,2},'UniformOutput',false);
constructed_proteomics_unique =  cell2table([file_PhosphoProteomics{:,1},file_proteomics_residuesSorted]);
%constructed_proteomics_unique = unique(constructed_proteomics_unique,'rows','stable');

proteins_set_total = cell(1,2);

index_rm = find(cellfun(@str2num,loc_prob)<0.75);
p_names_1 = p_names; p_names_1(index_rm) = [];
position_protein_1= position_protein; position_protein_1(index_rm) = [];

for i = 1:length(p_names_1)
    if(contains(p_names_1(i),"cont",'IgnoreCase',true)==0)
        proteins_set_total{1} = [proteins_set_total{1};extractBetween(p_names_1(i),'|','|')];
        a = strsplit(position_protein_1{i},';');
        proteins_set_total{2} = [proteins_set_total{2};a'];
    end
end
aj3 = sortrows(table(proteins_set_total{1,1},proteins_set_total{1,2}));

[ind] = ismember(constructed_proteomics_unique,aj3);
file_PhosphoProteomics_aux = file_PhosphoProteomics(find(ind),:);
writetable(file_PhosphoProteomics_aux,strcat(root_NetPhorest,'\Results_file_PhosphoProteomics_data_filterd_byPhosphoData.txt'),'Delimiter','\t');
%%   

%Structure where each cell corresponds to one protein containing all the
%kinases that are able to phosphorylate the respective protein further
%splitted by the phosphorylation residues
file_PhosphoProteomics = strcat(root_NetPhorest,'\Results_file_PhosphoProteomics_data_filterd_byPhosphoData.txt');
opts = detectImportOptions(file_PhosphoProteomics);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_PhosphoProteomics = readtable(file_PhosphoProteomics, opts);
file_PhosphoProteomics = sortrows(file_PhosphoProteomics);

n1 = unique(file_PhosphoProteomics.Var1,'stable');
kinase_list_proteomics = cell(length(n1),1);

for j = 1:length(n1)
    ind = find(ismember(file_PhosphoProteomics.Var1,n1(j)));
    kinase_list_proteomics{j,1} = file_PhosphoProteomics(ind,:);
    j
end
clear file_proteomics
for j = 1:length(kinase_list_proteomics)
    n1 = unique(kinase_list_proteomics{j}.Var2);
    as = cell(length(n1),1);
    for jj = 1:length(n1)
        ind = find(ismember(kinase_list_proteomics{j}.Var2,n1(jj)));
        as{jj} = kinase_list_proteomics{j}(ind,:);
    end
    kinase_list_proteomics{j} = as;
end
save(strcat(root_NetPhorest,'file_PhosphoProteomics_splited_by_proteins_byPhosphoPeptideAnalysis'),'kinase_list_proteomics','-v7','-nocompression')

%Dosen't work for nested cell arrays!!
% 
% kk = cell(length(kinase_list_proteomics),1);
% for j = 1:length(kinase_list_proteomics)
%     kk{j} = cellfun(@table2cell,kinase_list_proteomics{j},'UniformOutput',false);
% end
% kinase_list_proteomicsBYTES = hlp_serialize(kk);
% clear kk
% fileID = fopen('file_PhosphoProteomics_splited_by_proteins_byPhosphoPeptideAnalysis.bin','w');
% fwrite(fileID,kinase_list_proteomicsBYTES);
% clear kinase_list_proteomicsBYTES
% fclose(fileID);
%% 


%% 
%Structure where each cell corresponds to one protein containing all the
%kinases that are able to phosphorylate the respective protein further
%splitted by the phosphorylation residues filtered by interaction data
file_PhosphoProteomics = strcat(root_NetPhorest,'\Results_file_PhosphoProteomics_data_filterd_byPhosphoData.txt');
opts = detectImportOptions(file_PhosphoProteomics);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_PhosphoProteomics = readtable(file_PhosphoProteomics, opts);
file_PhosphoProteomics = sortrows(file_PhosphoProteomics);

NetPhorest_full_data = readtable(NetPhorest_full_names,'VariableNamingRule','preserve');
NetPhorest_full_data = sortrows(NetPhorest_full_data);
%Choose the desire database file - with multiple or fewer
%values
file_interactions = Interaction_database;
opts = detectImportOptions(file_interactions);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_interactions = readtable(file_interactions, opts);
ind_interactorA = find(ismember(file_interactions.Properties.VariableNames,'protein1'));
ind_interactorB = find(ismember(file_interactions.Properties.VariableNames,'protein2'));


n1 = table(file_PhosphoProteomics.Var1,file_PhosphoProteomics.Var8);
n1 = unique(n1,'stable');
n_aux = cell(size(n1,1),1);
for j = 1:size(n1,1)
    n_aux{j} = NetPhorest_full_data.Unip_ID(find(ismember(NetPhorest_full_data.NetPhorest,n1{j,2})));
end
n1{:,3} = n_aux;
n_aux = cell(size(n1,1),1);
n_aux1 = cell(size(n1,1),1);
for j = 1:size(n1,1)
    n_aux{j} = repmat(n1{j,1},length(n1{j,3}{:}),1);
    n_aux1{j} = repmat(j,length(n1{j,3}{:}),1);
end
n1{:,4} = n_aux;
n1{:,5} = n_aux1;
ns = cell2table([vertcat(n1{:,4}{:}),vertcat(n1{:,3}{:}),num2cell(vertcat(n1{:,5}{:}))]);

file_interactions_table1 = cell2table([file_interactions{:,ind_interactorA},file_interactions{:,ind_interactorB}]);
file_interactions_table2 = cell2table([file_interactions{:,ind_interactorB},file_interactions{:,ind_interactorA}]);

ind1 = ismember(ns(:,1:2),file_interactions_table1);
ind2 = ismember(ns(:,1:2),file_interactions_table2);
ind = or(ind1,ind2);
interaction_ind=find(ind);
ns = ns(interaction_ind,:);
ind_aux = unique(ns{:,3});
ind = ismember(table(file_PhosphoProteomics.Var1,file_PhosphoProteomics.Var8),n1(ind_aux,1:2));
file_PhosphoProteomics = file_PhosphoProteomics(ind,:);

n1 = unique(file_PhosphoProteomics.Var1,'stable');
kinase_list_proteomics = cell(length(n1),1);

for j = 1:length(n1)
    ind = find(ismember(file_PhosphoProteomics.Var1,n1(j)));
    kinase_list_proteomics{j,1} = file_PhosphoProteomics(ind,:);
    j
end
clear file_proteomics
for j = 1:length(kinase_list_proteomics)
    n1 = unique(kinase_list_proteomics{j}.Var2);
    as = cell(length(n1),1);
    for jj = 1:length(n1)
        ind = find(ismember(kinase_list_proteomics{j}.Var2,n1(jj)));
        as{jj} = kinase_list_proteomics{j}(ind,:);
    end
    kinase_list_proteomics{j} = as;
end

save(strcat(root_NetPhorest,'file_PhosphoProteomics_splited_by_proteins_byPhosphoPeptideAnalysis_fitlered_by_interaction'),'kinase_list_proteomics','-v7','-nocompression')

%Dosen't work for nested cell arrays!!
% z
% kk = cell(length(kinase_list_proteomics),1);
% for j = 1:length(kinase_list_proteomics)
%     kk{j} = cellfun(@table2cell,kinase_list_proteomics{j},'UniformOutput',false);
% end
% kinase_list_proteomicsBYTES = hlp_serialize(kk);
% clear kk
% fileID = fopen('file_PhosphoProteomics_splited_by_proteins_byPhosphoPeptideAnalysis_fitlered_by_interaction.bin','w');
% fwrite(fileID,kinase_list_proteomicsBYTES);
% clear kinase_list_proteomicsBYTES
% fclose(fileID);
%% 


end