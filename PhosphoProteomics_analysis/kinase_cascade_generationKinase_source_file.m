function [] = kinase_cascade_generationKinase_source_file(file_kinases_Phospho,dir_cr,file_kinases_list,file_TF,file_TF_mapping,file_TF_measured,dir_NetPhorest_file,workspace_file)

load(workspace_file)

%Extract all the proteins names
pat_ini_protein = 'sp|';
pat_fin_protein = '|';
pat_ini_gene = 'sp|'+ wildcardPattern + '|';
pat_fin_gene = '_HUMAN' ;


opts = detectImportOptions(file_kinases_Phospho);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_kinases_Phospho1 = readtable(file_kinases_Phospho, opts);

gene_protein_list = {};
for i = 1:length(p_names)
    aa1 = extractBetween(p_names(i),pat_ini_protein,pat_fin_protein);
    aa2 = extractBetween(p_names(i),pat_ini_gene,pat_fin_gene);
    gene_protein_list = [gene_protein_list;[aa1,aa2]];
end
gene_protein_list = unique(cell2table(gene_protein_list),'stable');
gene_protein_list = [gene_protein_list;[file_kinases_Phospho1{:,2},file_kinases_Phospho1{:,5}]];
gene_protein_list = [gene_protein_list;[file_kinases_Phospho1{:,1},file_kinases_Phospho1{:,4}]];
gene_protein_list = unique(gene_protein_list,'stable');

%Load the directories and files
dir1 = dir_cr;

analyzed_file = strcat(dir1,'\Excel_data_Kinase_Phosphorylated_by_Kinase\data__kinase_kinase_interactions_separate_data_group_data.xlsx');
Xlsx_Sheet_analyzed_file_data = sheetnames(analyzed_file);
file_kinase_total = strcat(dir1,'\Excel_data_Kinase_lists_from_data\initial_kinase_list.xlsx');
Xlsx_Sheet_file_kinase_total = sheetnames(file_kinase_total);

file_upregulated_members = strcat(dir1,'\Excel_data_unique_peptide_values\list_proteins_significant_group_comparison.xlsx');
Xlsx_Sheet_file_upregulated_members = sheetnames(file_upregulated_members);

dir_data_Netphorest = dir(dir_NetPhorest_file);

%TF based on each phenotype
interaction_data = cell(length(Xlsx_Sheet_analyzed_file_data),1);
for i = 1:length(Xlsx_Sheet_analyzed_file_data)
    %CHECK HERE!!!!!!!!!!!!!!! instead of 2 it should be i??
    %ind = find(contains(lower({dir_data_Netphorest.name}),lower(erase(Xlsx_Sheet_analyzed_file_data{2},'_'))));
    ind = find(contains(lower(erase({dir_data_Netphorest.name},'_')),lower(erase(Xlsx_Sheet_analyzed_file_data{i},'_'))));
    
    file_kinases_NetPhorest = strcat(dir_data_Netphorest(ind).folder,'\',dir_data_Netphorest(ind).name);
    
    aa1 = strsplit(lower(Xlsx_Sheet_analyzed_file_data{i}),'vs');
    ind_xlsx_sheet = find(ismember(lower(Xlsx_Sheet_file_kinase_total),aa1{1}));
    interaction_data{i} = kinase_cascade_generation(analyzed_file,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,Xlsx_Sheet_file_kinase_total{ind_xlsx_sheet},Xlsx_Sheet_analyzed_file_data{i},file_TF_measured,file_upregulated_members,Xlsx_Sheet_file_upregulated_members{i},file_kinases_list,p_names,p_full_phenotype{ind_xlsx_sheet});
    
    
end

opts = detectImportOptions(file_kinases_Phospho);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_kinases_Phospho = readtable(file_kinases_Phospho, opts);

opts = detectImportOptions(file_TF);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_TF = readtable(file_TF, opts);


opts = detectImportOptions(file_TF_mapping);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_TF_mapping = readtable(file_TF_mapping, opts);

for i = 1:length(interaction_data)
    
    ind_kinase = find(ismember(lower(interaction_data{i}{:,1}),lower(file_TF_mapping{:,2})));
    aux = zeros(size(interaction_data{i},1),1);
    aux(ind_kinase) = 1;
    aux = table(aux);
    aux.Properties.VariableNames = {'Kinase_TF'};
    interaction_data{i} = [interaction_data{i},aux];
    
    ind_substrate = find(ismember(lower(interaction_data{i}{:,2}),lower(file_TF_mapping{:,2})));
    aux = zeros(size(interaction_data{i},1),1);
    aux(ind_substrate) = 1;
    aux = table(aux);
    aux.Properties.VariableNames = {'Substrate_TF'};
    interaction_data{i} = [interaction_data{i},aux];
    
    ind_kinase = cellfun(@(y) find(ismember(lower(gene_protein_list{:,1}),lower(y))),interaction_data{i}{:,1},'UniformOutput',false);
    ind_kinase = cellfun(@(x) x(1),ind_kinase);
    aux = gene_protein_list{ind_kinase,2};
    aux = table(aux);
    aux.Properties.VariableNames = {'Kinase_gene'};
    interaction_data{i} = [interaction_data{i},aux];
    
    ind_substrate = cellfun(@(y) find(ismember(lower(gene_protein_list{:,1}),lower(y))),interaction_data{i}{:,2},'UniformOutput',false);
    ind_substrate = cellfun(@(x) x(1),ind_substrate);
    aux = gene_protein_list{ind_substrate,2};
    aux = table(aux);
    aux.Properties.VariableNames = {'Substrate_gene'};
    interaction_data{i} = [interaction_data{i},aux];
    
    
end

%Give a score of 2 in the upregulated rubric for the proteins that were
%not measured

for i = 1:length(interaction_data)
    ind_cut = [];
    ind = find(interaction_data{i}{:,3}==10);
    ind1 = find(interaction_data{i}{:,3}==210);
    ind2 = find(interaction_data{i}{:,3}==510);
    ind3 = find(interaction_data{i}{:,3}==110);
    ind = union(ind,ind1);
    ind = union(ind,ind2);
    ind = union(ind,ind3);
    ind_k = find(interaction_data{i}{ind,10}==0);
    ind_s = find(interaction_data{i}{ind,11}==0);
    interaction_data{i}{ind(ind_k),5} = 2;
    interaction_data{i}{ind(ind_s),6} = 2;
    interaction_data{i}{ind(ind_s),7} = 2;
end

dir_current = strcat(dir1,'\Kinase_Network_data');
mkdir(dir_current);
save_name = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison.xlsx');
delete(save_name)
for i = 1:length(Xlsx_Sheet_analyzed_file_data)
    writetable(interaction_data{i},save_name,'Sheet',Xlsx_Sheet_analyzed_file_data{i});
end

interaction_data_unique = cell(size(interaction_data));

for i = 1:length(interaction_data)
    [~,ind] = unique(interaction_data{i}(:,1:3),'stable');
    interaction_data_unique{i} = interaction_data{i}(ind,[1:3 5:end]);
end


save_name = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison_unique_value_no_residues.xlsx');
delete(save_name)
for i = 1:length(Xlsx_Sheet_analyzed_file_data)
    writetable(interaction_data_unique{i},save_name,'Sheet',Xlsx_Sheet_analyzed_file_data{i});
end


interaction_data_unique_only_measured = cell(size(interaction_data_unique));

for i = 1:length(interaction_data_unique)
    ind = find(interaction_data_unique{i}{:,3}~=10);
    ind1 = find(interaction_data_unique{i}{:,3}~=210);
    ind2 = find(interaction_data_unique{i}{:,3}~=510);
    ind3 = find(interaction_data_unique{i}{:,3}~=110);
    ind = intersect(ind,ind1);
    ind = intersect(ind,ind2);
    ind = intersect(ind,ind3);
    [~,iind] = unique(interaction_data{i}(ind,1:3),'stable');
    interaction_data_unique_only_measured{i} = interaction_data_unique{i}(ind(iind),:);
end

save_name = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison_unique_value_no_residues_only_measured.xlsx');
delete(save_name)
for i = 1:length(Xlsx_Sheet_analyzed_file_data)
    writetable(interaction_data_unique_only_measured{i},save_name,'Sheet',Xlsx_Sheet_analyzed_file_data{i});
end


%Merged togheter the pathways that have both a value of 2 and a value of 5
%because there are no more residues to be considered (leave the solo
%pathways with 5 alone)
%Do the same for a score of 10 combined with 110 210 or 510
interaction_data_unique_merged_elements = interaction_data_unique;

for i = 1:length(interaction_data_unique_merged_elements)
    ind_cut = [];
    ind = find(ismember(interaction_data_unique_merged_elements{i}{:,3},5));
    aa = table(interaction_data_unique_merged_elements{i}{ind,1},interaction_data_unique_merged_elements{i}{ind,2},repmat(2,length(ind),1));
    for j = 1:size(aa,1)
        ind1 = find(ismember(table(interaction_data_unique_merged_elements{i}{:,1},interaction_data_unique_merged_elements{i}{:,2},interaction_data_unique_merged_elements{i}{:,3}),aa(j,:)));
        if(isempty(ind1)==0)
            ind_cut = [ind_cut,ind(j)];
        end
    end
    interaction_data_unique_merged_elements{i}(ind_cut,:) = [];
end


for i = 1:length(interaction_data_unique_merged_elements)
    ind_cut = [];
    ind = find(ismember(interaction_data_unique_merged_elements{i}{:,3},10));
    aa_110 = table(interaction_data_unique_merged_elements{i}{ind,1},interaction_data_unique_merged_elements{i}{ind,2},repmat(110,length(ind),1));
    aa_210 = table(interaction_data_unique_merged_elements{i}{ind,1},interaction_data_unique_merged_elements{i}{ind,2},repmat(210,length(ind),1));
    aa_510 = table(interaction_data_unique_merged_elements{i}{ind,1},interaction_data_unique_merged_elements{i}{ind,2},repmat(510,length(ind),1));
    for j = 1:length(ind)
        ind1 = find(ismember(table(interaction_data_unique_merged_elements{i}{:,1},interaction_data_unique_merged_elements{i}{:,2},interaction_data_unique_merged_elements{i}{:,3}),aa_110(j,:)));
        ind2 = find(ismember(table(interaction_data_unique_merged_elements{i}{:,1},interaction_data_unique_merged_elements{i}{:,2},interaction_data_unique_merged_elements{i}{:,3}),aa_210(j,:)));
        ind3 = find(ismember(table(interaction_data_unique_merged_elements{i}{:,1},interaction_data_unique_merged_elements{i}{:,2},interaction_data_unique_merged_elements{i}{:,3}),aa_510(j,:)));
        ind1 = union(ind1,ind2);
        ind1 = union(ind1,ind3);
        if(isempty(ind1)==0)
            ind_cut = [ind_cut,ind(j)];
        end
    end
    interaction_data_unique_merged_elements{i}(ind_cut,:) = [];
end

%Remove the self-phosphorylation entries that have a score of 5
%If needed

%  for i = 1:length(interaction_data_unique_merged_elements)
%      ind_cut = [];
%      ind = find(ismember(interaction_data_unique_merged_elements{i}{:,3},5));
%      aa = table(interaction_data_unique_merged_elements{i}{ind,1},interaction_data_unique_merged_elements{i}{ind,2});
%      for j = 1:size(aa,1)
%          ind1 = find(ismember(aa{j,1},aa{j,2}));
%          if(isempty(ind1)==0)
%              ind_cut = [ind_cut,ind(j)];
%          end
%      end
%      interaction_data_unique_merged_elements{i}(ind_cut,:) = [];
%  end



save_name = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison_unique_value_no_residues_merged_elements.xlsx');
delete(save_name)
for i = 1:length(Xlsx_Sheet_analyzed_file_data)
    writetable(interaction_data_unique_merged_elements{i},save_name,'Sheet',Xlsx_Sheet_analyzed_file_data{i});
end





%Prepare file for visualization - starting from procced file

save_name = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison_unique_value_no_residues_merged_elements_Prepared_for_visualization.xlsx');
delete(save_name)

save_name_aux = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison_unique_value_no_residues.xlsx');
delete(save_name_aux)
% 
% file_visualization_name = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison.xlsx');
% 
% Xlsx_Sheet_analyzed_file_data = sheetnames(file_visualization_name);
% for ks = 1:length(Xlsx_Sheet_analyzed_file_data)
%     
%     opts = detectImportOptions(file_visualization_name);
%     opts = setvartype(opts,1:length(opts.VariableNames), 'char');
%     opts.PreserveVariableNames=true;
%     file_visualization = readtable(file_visualization_name, opts,'Sheet',Xlsx_Sheet_analyzed_file_data{ks});
%     
%     col_n = file_visualization.Properties.VariableNames;
%     
%     file_visualization{:,end+1} = {'0'};
%     file_visualization.Properties.VariableNames(end) = {'Kinase_presence(0-not measured;1-measured&upreg; 2-measured)'};
%     file_visualization{find(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Kinase_Measured"))})==1),end} = {'2'};
%     file_visualization{find(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Kinase_Upregulated"))})==1),end} = {'1'};
%     a1 = file_visualization{:,end};
%     a2 = file_visualization{:,13};
%     file_visualization = file_visualization(:,1:14);
%     
%     %Remove all the entries that have a score of 5, 510
%     ind = find(ismember(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Database"))}),5));
%     file_visualization(ind,:) = [];
%     ind = find(ismember(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Database"))}),510));
%     file_visualization(ind,:) = [];
%     
%     %Remove all the entries that have a score of 3 and come from and have
%     %a link to or from an measured and not upregulated
%     %kinase
%     ind = find(ismember(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Database"))}),3));
%     ind1 = ind(find(ismember(cellfun(@str2num,file_visualization{ind,find(ismember(col_n,"Kinase_Upregulated"))}),0)));
%     ind2 = ind(find(ismember(cellfun(@str2num,file_visualization{ind,find(ismember(col_n,"Substrate_Upregulated"))}),0)));
%     ind = union(ind1,ind2);
%     file_visualization(ind,:) = [];
%     
%     %Keep only the PSP links that are on upregulated residues of TF. For
%     %both measured and unemasured kinases.
%     ind = find(ismember(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Database"))}),1));
%     ind = ind(find(ismember(cellfun(@str2num,file_visualization{ind,find(ismember(col_n,"Substrate_Upregulated"))}),0)));
%     file_visualization(ind,:) = [];
%     ind = find(ismember(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Database"))}),110));
%     ind = ind(find(ismember(cellfun(@str2num,file_visualization{ind,find(ismember(col_n,"Substrate_Upregulated"))}),0)));
%     file_visualization(ind,:) = [];
%     
%     %Keep only the auxiliary kinases that are bridges between two other
%     %upregulated elements.
%     ind1 = find(ismember(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Database"))}),10));
%     ind2 = find(ismember(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Database"))}),210));
%     ind3 = find(ismember(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Database"))}),110));
%     ind = union(ind1,ind2);
%     ind = union(ind,ind3);
%     ind1 = find(ismember(cellfun(@str2num,file_visualization{ind,find(ismember(col_n,"Kinase_Upregulated"))}),0));
%     ind2 = find(ismember(cellfun(@str2num,file_visualization{ind,find(ismember(col_n,"Substrate_Upregulated"))}),0));
%     ind3 = union(ind(ind1),ind(ind2));
%     ind(union(ind1,ind2)) = [];
%     
% %     aux = file_visualization(ind,:);
% %     un_aux = {};
% %     for i=1:size(ind,1)
% %         un_aux(i,find(str2num(cell2mat(aux{i,4:5}'))==2)) = aux{i,find(str2num(cell2mat(aux{i,4:5}'))==2)};
% %     end
% %     
% %     un1 = unique(string(vertcat(un_aux{(cellfun(@length,vertcat(un_aux(:,1)))==6),1})));
% %     un2 = unique(string(vertcat(un_aux{(cellfun(@length,vertcat(un_aux(:,2)))==6),2})));
% %     
% %     un = un1(find(ismember(un1,un2)));
% %     
% %     is1 = find(ismember(aux{:,1},un));
% %     is2 = find(ismember(aux{:,2},un));
% %     is = union(is1,is2);
% %     ind(is) = [];
% %     
% %     ind = union(ind,ind3);
%     ind = ind3;
%     file_visualization(ind,:) = [];
%     
%     %Remove entries of kinasese that self-phosphorylates but are not
%     %upregulated
%     
%     file_visualization{:,1} = cellfun(@(x) {string(x)},file_visualization{:,find(ismember(col_n,"Kinase"))});
%     file_visualization{:,2} = cellfun(@(x) {string(x)},file_visualization{:,find(ismember(col_n,"Substrate"))});
%     ind = find(cellfun(@(x,y) ismember(x,y),file_visualization{:,find(ismember(col_n,"Kinase"))},file_visualization{:,find(ismember(col_n,"Substrate"))}));
%     ind1 = ind(find(ismember(cellfun(@str2num,file_visualization{ind,find(ismember(col_n,"Kinase_Upregulated"))}),0)));
%     ind2 = ind(find(ismember(cellfun(@str2num,file_visualization{ind,find(ismember(col_n,"Substrate_Upregulated"))}),0)));
%     ind = intersect(ind1,ind2);
%     file_visualization(ind,:) = [];
%     
%     %Remove auxilary and not upregualted kinases that are dangling.
%     
%     ind1 = find(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Kinase_Upregulated"))})~=1);
%     ind2 = find(cellfun(@str2num,file_visualization{:,find(ismember(col_n,"Substrate_Upregulated"))})~=1);
%     
%     un1 = file_visualization{ind1,find(ismember(col_n,"Kinase"))};
%     un2 = file_visualization{ind2,find(ismember(col_n,"Substrate"))};
%     
%     ind11 = find(ismember(vertcat(un1{:}),vertcat(un2{:})));
%     
%     ind22 = find(ismember(vertcat(un2{:}),vertcat(un1{:})));
%     
%     ind1(ind11) = [];
%     ind2(ind22) = [];
%     
%     ind = union(ind1,ind2);
%     file_visualization(ind,:) = [];
%     
%     
%     file_visualization{:,[3,4,5,6,9,10,11,12]} = cellfun(@(x) {str2num(x)},file_visualization{:,[3,4,5,6,9,10,11,12]});
%     %file_visualization{:,[3,5,6,7,10,11,12,13]} = cellfun(@(x) {str2num(x)},file_visualization{:,[3,5,6,7,10,11,12,13]});
%     Kinase_presence = cell2table([cellfun(@(x) {str2num(x)},a1),a2]);
%     Kinase_presence.Properties.VariableNames = ["Kinase_presence(0-not measured;1-measured&upreg; 2-measured)","Kinase_ID"];
%     
%     writetable(file_visualization,save_name,'Sheet',Xlsx_Sheet_analyzed_file_data{ks});
%     writetable(Kinase_presence,save_name_aux,'Sheet',Xlsx_Sheet_analyzed_file_data{ks});
%     
% end







%Prepare file for visualization - starting from raw file

save_name = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison_unique_value_no_residues_merged_elements_Prepared_for_visualization.xlsx');
delete(save_name)

save_name_aux = strcat(dir_current,'\Kinase_Network_data_aux_info.xlsx');
delete(save_name_aux)

file_visualization_name = strcat(dir_current,'\Kinase_Network_data_phenotype_comaprison.xlsx');

Xlsx_Sheet_analyzed_file_data = sheetnames(file_visualization_name);
for ks = 1:length(Xlsx_Sheet_analyzed_file_data)
    
    opts = detectImportOptions(file_visualization_name);
    opts = setvartype(opts,1:length(opts.VariableNames), 'char');
    opts.PreserveVariableNames=true;
    file_visualization = readtable(file_visualization_name, opts,'Sheet',Xlsx_Sheet_analyzed_file_data{ks});
    
    file_visualization{:,end+1} = {'0'};
    file_visualization.Properties.VariableNames(end) = {'Kinase_presence(0-not measured;1-measured&upreg; 2-measured)'};
    file_visualization{find(cellfun(@str2num,file_visualization{:,10})==1),end} = {'2'};
    file_visualization{find(cellfun(@str2num,file_visualization{:,5})==1),end} = {'1'};
    a1 = file_visualization{:,end};
    a2 = file_visualization{:,14};
    file_visualization = file_visualization(:,1:15);
    
    %Remove all the entries that have a score of 5, 510
    ind = find(ismember(cellfun(@str2num,file_visualization{:,3}),5));
    file_visualization(ind,:) = [];
    ind = find(ismember(cellfun(@str2num,file_visualization{:,3}),510));
    file_visualization(ind,:) = [];
    
    %Remove all the entries that have a score of 3 and come from and have
    %a link to or from an measured and not upregulated
    %kinase
    ind = find(ismember(cellfun(@str2num,file_visualization{:,3}),3));
    ind1 = ind(find(ismember(cellfun(@str2num,file_visualization{ind,5}),0)));
    ind2 = ind(find(ismember(cellfun(@str2num,file_visualization{ind,6}),0)));
    ind = union(ind1,ind2);
    file_visualization(ind,:) = [];
    
    %Keep only the PSP links that are on upregulated residues of TF. For
    %both measured and unemasured kinases. EDIT: ONYL FOR THE UNMEASURED
%     ind = find(ismember(cellfun(@str2num,file_visualization{:,3}),1));
%     ind = ind(find(ismember(cellfun(@str2num,file_visualization{ind,6}),0)));
%     file_visualization(ind,:) = [];
    ind = find(ismember(cellfun(@str2num,file_visualization{:,3}),110));
    ind = ind(find(ismember(cellfun(@str2num,file_visualization{ind,6}),0)));
    file_visualization(ind,:) = [];
    
    %Keep only the auxiliary kinases that are bridges between two other
    %upregulated elements.
    ind1 = find(ismember(cellfun(@str2num,file_visualization{:,3}),10));
    ind2 = find(ismember(cellfun(@str2num,file_visualization{:,3}),210));
    ind3 = find(ismember(cellfun(@str2num,file_visualization{:,3}),110));
    ind = union(ind1,ind2);
    ind = union(ind,ind3);
    ind1 = find(ismember(cellfun(@str2num,file_visualization{ind,5}),0));
    ind2 = find(ismember(cellfun(@str2num,file_visualization{ind,6}),0));
    ind3 = union(ind(ind1),ind(ind2));
    ind(union(ind1,ind2)) = [];
    
%     aux = file_visualization(ind,:);
%     un_aux = {};
%     for i=1:size(ind,1)
%         un_aux(i,find(str2num(cell2mat(aux{i,4:5}'))==2)) = aux{i,find(str2num(cell2mat(aux{i,4:5}'))==2)};
%     end
%     
%     un1 = unique(string(vertcat(un_aux{(cellfun(@length,vertcat(un_aux(:,1)))==6),1})));
%     un2 = unique(string(vertcat(un_aux{(cellfun(@length,vertcat(un_aux(:,2)))==6),2})));
%     
%     un = un1(find(ismember(un1,un2)));
%     
%     is1 = find(ismember(aux{:,1},un));
%     is2 = find(ismember(aux{:,2},un));
%     is = union(is1,is2);
%     ind(is) = [];
%     
%     ind = union(ind,ind3);
    ind = ind3;
    file_visualization(ind,:) = [];
    
    %Remove entries of kinasese that self-phosphorylates but are not
    %upregulated
    
    file_visualization{:,1} = cellfun(@(x) {string(x)},file_visualization{:,1});
    file_visualization{:,2} = cellfun(@(x) {string(x)},file_visualization{:,2});
    ind = find(cellfun(@(x,y) ismember(x,y),file_visualization{:,1},file_visualization{:,2}));
    ind1 = ind(find(ismember(cellfun(@str2num,file_visualization{ind,5}),0)));
    ind2 = ind(find(ismember(cellfun(@str2num,file_visualization{ind,6}),0)));
    ind = intersect(ind1,ind2);
    file_visualization(ind,:) = [];
    
    %Remove auxilary and not upregualted kinases that are dangling. DO it
    %two times!
    
    %Round 1
    ind1 = find(cellfun(@str2num,file_visualization{:,5})~=1);
    ind2 = find(cellfun(@str2num,file_visualization{:,6})~=1);
    un1 = file_visualization{ind1,1};
    un2 = file_visualization{ind2,2}; 
    ind11 = find(ismember(vertcat(un1{:}),vertcat(un2{:})));
    ind22 = find(ismember(vertcat(un2{:}),vertcat(un1{:})));
    ind1(ind11) = [];
    ind2(ind22) = [];
    ind = union(ind1,ind2);
    file_visualization(ind,:) = [];
    
    %Round 2
    ind1 = find(cellfun(@str2num,file_visualization{:,5})~=1);
    ind2 = find(cellfun(@str2num,file_visualization{:,6})~=1);
    un1 = file_visualization{ind1,1};
    un2 = file_visualization{ind2,2}; 
    ind11 = find(ismember(vertcat(un1{:}),vertcat(un2{:})));
    ind22 = find(ismember(vertcat(un2{:}),vertcat(un1{:})));
    ind1(ind11) = [];
    ind2(ind22) = [];
    ind = union(ind1,ind2);
    file_visualization(ind,:) = [];
    
    file_visualization{:,1} = cellfun(@str2mat,file_visualization{:,1} ,'UniformOutput',false);
    file_visualization{:,2} = cellfun(@str2mat,file_visualization{:,2} ,'UniformOutput',false);
    [~,ind] = unique(file_visualization(:,[1,2,5,6,7]),'stable');
    file_visualization = file_visualization(ind,[1:3 5:end]);
    
    file_visualization{:,[3,4,5,6,9,10,11,12]} = cellfun(@(x) {str2num(x)},file_visualization{:,[3,4,5,6,9,10,11,12]});
    Kinase_presence = cell2table([cellfun(@(x) {str2num(x)},a1),a2]);
    Kinase_presence.Properties.VariableNames = ["Kinase_presence_0_not_measured_1_measured_and_upreg_2_measured","Kinase_ID"];
    
    writetable(file_visualization,save_name,'Sheet',Xlsx_Sheet_analyzed_file_data{ks});
    writetable(Kinase_presence,save_name_aux,'Sheet',Xlsx_Sheet_analyzed_file_data{ks});
    
end



end