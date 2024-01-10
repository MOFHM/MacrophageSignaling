function [list] = kinase_information_phosphorylation(kinase_list,dataset_phospho,dataset_proteomics,dir1,xlsx_name,Sheet_name,delete_flag)
%kinase_lsit - contains two columns - the first one contains the acc id of
%the kinases and the second one contains the phosphorylated residues of the
%respective kinases

file_name = ".\Kinase_Substrate_Dataset.txt";
opts = detectImportOptions(file_name);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_kinase_substrate = readtable(file_name, opts);
ind_sub_acc = find(ismember(file_kinase_substrate.Properties.VariableNames,'SUB_ACC_ID'));
ind_sub_rsd = find(ismember(file_kinase_substrate.Properties.VariableNames,'SUB_MOD_RSD'));
ind_kin_acc = find(ismember(file_kinase_substrate.Properties.VariableNames,'KIN_ACC_ID'));
ind_kin_name = find(ismember(file_kinase_substrate.Properties.VariableNames,'KINASE'));
ind_sub_name = find(ismember(file_kinase_substrate.Properties.VariableNames,'SUBSTRATE'));

list = {};
for i = 1:size(kinase_list,1)
    
    ind_name = find(ismember(file_kinase_substrate{:,ind_sub_acc},kinase_list{i,1}));
    ind_residue = find(ismember(file_kinase_substrate{ind_name,ind_sub_rsd},kinase_list{i,2}));
    ind_residue = ind_name(ind_residue);
    
    if(isempty(ind_residue)==0)
        list = [list;[file_kinase_substrate{ind_residue,ind_sub_acc},file_kinase_substrate{ind_residue,ind_sub_name},file_kinase_substrate{ind_residue,ind_sub_rsd},file_kinase_substrate{ind_residue,ind_kin_acc},file_kinase_substrate{ind_residue,ind_kin_name}]];
    end
    
end

presence_dataset_phospho = cell(size(list,1),1);
presence_dataset_proteomics = cell(size(list,1),1);

for i = 1:size(list,1)
    ind = find(ismember(dataset_phospho,list{i,4}));
    if(isempty(ind)==0) presence_dataset_phospho(i) = {'true'}; end
    ind = find(ismember(dataset_proteomics,list{i,4}));
    if(isempty(ind)==0) presence_dataset_proteomics(i) = {'true'}; end
end

list = [list,presence_dataset_phospho,presence_dataset_proteomics];


mkdir(dir1,'Excel_data_Kinase_Phosphorylated_by_Kinase')
dir_current = strcat(dir1,'\Excel_data_Kinase_Phosphorylated_by_Kinase');
recycle on % Send to recycle bin instead of permanently deleting.
if(delete_flag==1) delete(strcat(dir_current,'/data_',xlsx_name,'.xlsx')); end
if(isempty(list)==0)
header = {'Substrate(target_kinase)_ACC','Substrate(target_kinase)_Name','Residue','Kinase_ACC(from_PhosphoSitePlus)','Kinase_Name(from_PhosphoSitePlus)','Presence_in_PhosphoProteomics','Presence_in_Proteomics'};
writecell(header,strcat(dir_current,'/data_',xlsx_name,'.xlsx'),'Sheet',Sheet_name,'Range','A1','AutofitWidth',1);
writecell(list,strcat(dir_current,'/data_',xlsx_name,'.xlsx'),'Sheet',Sheet_name,'Range','A2','AutofitWidth',0);
end

end



