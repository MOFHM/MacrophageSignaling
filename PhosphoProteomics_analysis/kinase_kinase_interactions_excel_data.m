function [] = kinase_kinase_interactions_excel_data(dir1,dir_current,file_name_kinome,nr_sample_groups,sample_names,list_proteins_significant_non_unique_upreg_complete_kinases,list_proteins_significant_non_unique_downreg_complete_kinase,...
list_proteins_significant_unique_complete_kinases,kinase_original_list)

for i = 1:nr_sample_groups
%Kinase-kinase interaction identification - M1
[dataset_proteomics_phenotype] = extract_proteomics_kinases(file_name_kinome);
aa1 = vertcat(list_proteins_significant_non_unique_upreg_complete_kinases{i,:},list_proteins_significant_non_unique_downreg_complete_kinase{:,i});
%aa1 = cell2table(aa1);
aa1 = unique(aa1,'rows','stable');
aa1 = table2cell(aa1);
aa1 = [aa1(:,1),aa1(:,4)];

if(isempty(list_proteins_significant_unique_complete_kinases{1})==0)
    kinase_list = [[list_proteins_significant_unique_complete_kinases{i}(:,1),list_proteins_significant_unique_complete_kinases{i}(:,4)];aa1];
else
    kinase_list = cell2table(aa1);
end
dataset_phospho = kinase_original_list{i}{:,1}; %If we want only the kinases list
dataset_proteomics = dataset_proteomics_phenotype{i}; %If we want all the kinases measured in the respective phenotype
xlsx_name = strcat('',sample_names{i},'_kinase_kinase_interactions');
[~] = kinase_information_phosphorylation(kinase_list,dataset_phospho,dataset_proteomics,dir1,xlsx_name,strcat(sample_names{i},'VsAll'),1);
end


%Excel with all the informations split into sets
%M1
[dataset_proteomics_phenotype] = extract_proteomics_kinases(file_name_kinome);
for i = 1:nr_sample_groups
aa = {};
if(isempty(list_proteins_significant_unique_complete_kinases{i})==0)
    aa{1} = [list_proteins_significant_unique_complete_kinases{i}(:,1),list_proteins_significant_unique_complete_kinases{i}(:,4)];
end
ind_j = 1;
for j = 1:nr_sample_groups
    aa1 = vertcat(list_proteins_significant_non_unique_upreg_complete_kinases{i,j},list_proteins_significant_non_unique_downreg_complete_kinase{j,i});
    if(isempty(aa1)==0)
        ind_j = ind_j+1;
        aa1 = [aa1(:,1),aa1(:,4)];
        aa{ind_j} = aa1;
    end
end

dataset_phospho = kinase_original_list{i}{:,1}; %If we want only the kinases list
dataset_proteomics = dataset_proteomics_phenotype{i}; %If we want all the proteins measured i nthe respective phenotype
xlsx_name = strcat('',sample_names{i},'_kinase_kinase_interactions_separate_data');
for ii = 1:length(aa)
    if (ii==i)
        Sheet_name = "unique";
    else
        Sheet_name = strcat(sample_names{i},'Vs',sample_names{ii});
    end
    if (ii==1) delete_flag = 1; else delete_flag = 0; end
    [~] = kinase_information_phosphorylation(aa{ii},dataset_phospho,dataset_proteomics,dir1,xlsx_name,Sheet_name,delete_flag);
end
end



%Gropu comparison excel data - filtered by PhosPhoSite!!!!!
delete(strcat(dir_current,'\_kinase_kinase_interactions_separate_data_group_data.xlsx'));
xlsx_name = strcat('_kinase_kinase_interactions_separate_data_group_data');
[dataset_proteomics_phenotype] = extract_proteomics_kinases(file_name_kinome);
for i = 1:size(list_proteins_significant_non_unique_upreg_complete_kinases,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete_kinases,2)
        if(isempty(list_proteins_significant_non_unique_upreg_complete_kinases{i,j})==0)
            aas = [list_proteins_significant_unique_complete_kinases{i};list_proteins_significant_non_unique_upreg_complete_kinases{i,j}];
            aas = aas(:,[1 4]);
            Sheet_name = strcat(sample_names{i},'Vs',sample_names{j});
            delete_flag = 0;
            [~] = kinase_information_phosphorylation(aas,dataset_phospho,dataset_proteomics,dir1,xlsx_name,Sheet_name,delete_flag);
            
        end
    end
end
%The index is inveresed because downregulated means inverse upregulated
xlsx_name = strcat('_kinase_kinase_interactions_separate_data_group_data');
[dataset_proteomics_phenotype] = extract_proteomics_kinases(file_name_kinome);
for i = 1:size(list_proteins_significant_non_unique_downreg_complete_kinase,2)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete_kinase,1)
        if((isempty(list_proteins_significant_non_unique_downreg_complete_kinase{j,i})==0))
            aas = [list_proteins_significant_unique_complete_kinases{i};list_proteins_significant_non_unique_downreg_complete_kinase{j,i}];
            aas = aas(:,[1 4]);
            Sheet_name = strcat(sample_names{i},'Vs',sample_names{j});
            delete_flag = 0;
            [~] = kinase_information_phosphorylation(aas,dataset_phospho,dataset_proteomics,dir1,xlsx_name,Sheet_name,delete_flag);
            
        end
    end
end



end