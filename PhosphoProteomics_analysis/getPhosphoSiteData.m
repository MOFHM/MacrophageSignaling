function [list_kinase_proteins,score_kin_sub] = getPhosphoSiteData(p_names,position_protein,PSP_file)

pat_ini_protein = 'sp|';
pat_fin_protein = '|';
proteins_set_total = cell(1,2);
for i = 1:length(p_names)
    if(contains(p_names(i),"cont",'IgnoreCase',true)==0)
    ind1 = size(proteins_set_total{1},1);
    if (isempty(find(contains(p_names,'sp'),1)) == 0)
    proteins_set_total{1} = [proteins_set_total{1};extractBetween(p_names(i),pat_ini_protein,pat_fin_protein)];
    else
         proteins_set_total{1} = [proteins_set_total{1};p_names(i)];
    end
    ind2 = size(proteins_set_total{1},1);
    a = strsplit(position_protein{i},';');
    proteins_set_total{2} = [proteins_set_total{2};a'];
    end
end

file_name = PSP_file;
opts = detectImportOptions(file_name);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_kinase_substrate = readtable(file_name, opts);
ind_sub_acc = find(ismember(file_kinase_substrate.Properties.VariableNames,'SUB_ACC_ID'));
ind_sub_rsd = find(ismember(file_kinase_substrate.Properties.VariableNames,'SUB_MOD_RSD'));
ind_kin_acc = find(ismember(file_kinase_substrate.Properties.VariableNames,'KIN_ACC_ID'));
ind_kin_name = find(ismember(file_kinase_substrate.Properties.VariableNames,'KINASE'));


list_kinase_proteins = {};
for i = 1:size(proteins_set_total{1},1)
    
    ind_name = find(ismember(file_kinase_substrate{:,ind_sub_acc},proteins_set_total{1}{i}));
    ind_residue = find(ismember(file_kinase_substrate{ind_name,ind_sub_rsd},proteins_set_total{2}{i}));
    ind_residue = ind_name(ind_residue);
    
    if(isempty(ind_residue)==0)
    list_kinase_proteins = [list_kinase_proteins;{[file_kinase_substrate{ind_residue,ind_sub_acc},file_kinase_substrate{ind_residue,ind_sub_rsd},file_kinase_substrate{ind_residue,ind_kin_acc},file_kinase_substrate{ind_residue,ind_kin_name}]}];
    end
    
end


if(isempty(list_kinase_proteins)==0)
score_kin_sub = {};       
%Unique after acc
list_kinases_subTABLE = vertcat(list_kinase_proteins{:});
list_kinases_sub = list_kinases_subTABLE(:,3);
s1 = unique(list_kinases_sub,'stable');
%Unique after name
% list_kinases_sub = vertcat(list_kinase_proteins{:});
% list_kinases_sub = list_kinases_sub(:,4);
% s1 = unique(list_kinases_sub,'stable');

        for j = 1:length(s1)
            aa1 = num2cell(length(find(ismember(list_kinases_sub,s1(j)))));
            %Divide the number of identified kinases by the number of
            %analyzed phophopeptides (proteins with a specific
            %phosphorylated site in order to see how active that kinase is
            %in general
            aa2 = num2cell(length(find(ismember(list_kinases_sub,s1(j))))./size(list_kinase_proteins,1));
            aa3 = list_kinases_subTABLE(find(ismember(list_kinases_subTABLE(:,3),s1(j)),1),4);
            score_kin_sub = [score_kin_sub;[s1(j),aa3,aa1,aa2]];
        end
        
        [~,a2] = sort(cell2mat(score_kin_sub(:,3)),'descend');
        score_kin_sub = score_kin_sub(a2,:);
else
    score_kin_sub = {};
end
        
end