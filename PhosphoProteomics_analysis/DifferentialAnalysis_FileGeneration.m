
%Write excel files
mkdir(dir1,'Excel_data_unique_peptide_values')
dir_current = strcat(dir1,'\Excel_data_unique_peptide_values');
recycle on % Send to recycle bin instead of permanently deleting.

for i = 1:size(list_proteins_significant_unique,2)
    delete(strcat(dir_current,'/list_proteins_significant_unique_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_unique,1)
        if(isempty(list_proteins_significant_unique{j,i})==0)
            xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{j},'None','Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
            writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_unique_','p',num2str(i),'.xlsx'),'Range','A1','Sheet',sample_names{j});
            writecell(list_proteins_significant_unique{j,i},strcat(dir_current,'/list_proteins_significant_unique_','p',num2str(i),'.xlsx'),'Range','A2','Sheet',sample_names{j});
       
        end
    end
end

for i = 1:size(list_proteins_significant_non_unique,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique,1)
        for k = 1:size(list_proteins_significant_non_unique,2)
            if(isempty(list_proteins_significant_non_unique{j,k,i})==0)
                xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{j},sample_names{k},'Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
                writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_non_unique_','p',num2str(i),'.xlsx'),'Range','A1','Sheet',strcat(sample_names{j},'vs',sample_names{k}));
                writecell(list_proteins_significant_non_unique{j,k,i},strcat(dir_current,'/list_proteins_significant_non_unique_','p',num2str(i),'.xlsx'),'Range','A2','Sheet',strcat(sample_names{j},'vs',sample_names{k}));
            end
        end
    end
end

for i = 1:size(list_proteins_significant_non_unique_upreg,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_upreg_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique_upreg,1)
        for k = 1:size(list_proteins_significant_non_unique_upreg,2)
            if(isempty(list_proteins_significant_non_unique_upreg{j,k,i})==0)
                xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{j},sample_names{k},'Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
                writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_non_unique_upreg_','p',num2str(i),'.xlsx'),'Range','A1','Sheet',strcat(sample_names{j},'vs',sample_names{k}));
                writecell(list_proteins_significant_non_unique_upreg{j,k,i},strcat(dir_current,'/list_proteins_significant_non_unique_upreg_','p',num2str(i),'.xlsx'),'Range','A2','Sheet',strcat(sample_names{j},'vs',sample_names{k}));
            end
        end
    end
end

for i = 1:size(list_proteins_significant_non_unique_downreg,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_downreg_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique_downreg,1)
        for k = 1:size(list_proteins_significant_non_unique_downreg,2)
            if(isempty(list_proteins_significant_non_unique_downreg{j,k,i})==0)
                xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{j},sample_names{k},'Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
                writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_non_unique_downreg_','p',num2str(i),'.xlsx'),'Range','A1','Sheet',strcat(sample_names{j},'vs',sample_names{k}));
                writecell(list_proteins_significant_non_unique_downreg{j,k,i},strcat(dir_current,'/list_proteins_significant_non_unique_downreg_','p',num2str(i),'.xlsx'),'Range','A2','Sheet',strcat(sample_names{j},'vs',sample_names{k}));
            end
        end
    end
end

delete(strcat(dir_current,'/list_proteins_significant_unique_complete.xlsx'));
for i = 1:length(list_proteins_significant_unique_complete)
    if(isempty(list_proteins_significant_unique_complete{i})==0)
        xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{j},'None','Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
        writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_unique_complete.xlsx'),'Range','A1','Sheet',sample_names{i});
        writecell(list_proteins_significant_unique_complete{i},strcat(dir_current,'/list_proteins_significant_unique_complete.xlsx'),'Range','A2','Sheet',sample_names{i});
    end
end

delete(strcat(dir_current,'/list_proteins_significant_non_unique_downreg_complete.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_downreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete,2)
        if(isempty(list_proteins_significant_non_unique_downreg_complete{i,j})==0)
            xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{i},sample_names{j},'Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
            writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_non_unique_downreg_complete.xlsx'),'Range','A1','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
            writecell(list_proteins_significant_non_unique_downreg_complete{i,j},strcat(dir_current,'/list_proteins_significant_non_unique_downreg_complete.xlsx'),'Range','A2','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end

delete(strcat(dir_current,'/list_proteins_significant_non_unique_upreg_complete.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_upreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete,2)
        if(isempty(list_proteins_significant_non_unique_upreg_complete{i,j})==0)
            xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{i},sample_names{j},'Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
            writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_non_unique_upreg_complete.xlsx'),'Range','A1','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
            writecell(list_proteins_significant_non_unique_upreg_complete{i,j},strcat(dir_current,'/list_proteins_significant_non_unique_upreg_complete.xlsx'),'Range','A2','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end



list_proteins_total = cell(1,nr_sample_groups);
for k = 1:nr_sample_groups
    for i = (k-1)*nr_rep+1:nr_rep*nr_sample_groups:size(int_data,2)
        for j = 1:size(int_data,1)
            vec = int_data(j,i:i+nr_rep-1);
            [vec_sort,initial_ind] = sort(vec,'descend');
            ind = find(vec_sort==0);
            if(length(ind)<=nr_rep-1) list_proteins_total{1,k} = [list_proteins_total{1,k};p_names(j)]; end
        end
    end
end

for i = 1:length(list_proteins_total)
    list_proteins_total{i} =  unique(extractBetween([list_proteins_total{1,i}{:}],'|','|'));
end
list_proteins_total = cellfun(@(x) remove_cont(x),list_proteins_total,'UniformOutput',false);

delete(strcat(dir_current,'/list_proteins_total_each_phenotype.xlsx'));
for i = 1:length(list_proteins_total)
    writecell(list_proteins_total{i},strcat(dir_current,'/list_proteins_total_each_phenotype.xlsx'),'Sheet',strcat(sample_names{i}));
end
for i = 1:length(list_proteins_total)
    delete(strcat(dir_current,'/list_proteins_total_each_phenotype_',sample_names{i},'.txt'));
    writecell(list_proteins_total{i},strcat(dir_current,'/list_proteins_total_each_phenotype_',sample_names{i},'.txt'));
end

%Gropu comparison excel data!!!!!
delete(strcat(dir_current,'/list_proteins_significant_group_comparison.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_upreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete,2)
        if(isempty(list_proteins_significant_non_unique_upreg_complete{i,j})==0)
            aas = [list_proteins_significant_unique_complete{i};list_proteins_significant_non_unique_upreg_complete{i,j}];
            xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{i},sample_names{j},'Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
            writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_group_comparison.xlsx'),'Range','A1','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
            writecell(aas,strcat(dir_current,'/list_proteins_significant_group_comparison.xlsx'),'Range','A2','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end

%The index is inveresed because downregulated means inverse upregulated
for i = 1:size(list_proteins_significant_non_unique_downreg_complete,2)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete,1)
        if(isempty(list_proteins_significant_non_unique_downreg_complete{j,i})==0)
            aas = [list_proteins_significant_unique_complete{i};list_proteins_significant_non_unique_downreg_complete{j,i}];
            xlx_header = {'Protein','Gene','Peptide','Residue.Both','FDR','FC','baseMean','Pvalue',sample_names{i},sample_names{j},'Localization.Probability','STDerr','DF',"ProteinCompensation","Number_Of_Phosphorylations"};
             writecell(xlx_header,strcat(dir_current,'/list_proteins_significant_group_comparison.xlsx'),'Range','A1','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
            writecell(aas,strcat(dir_current,'/list_proteins_significant_group_comparison.xlsx'),'Range','A2','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end



%The identified kinases for all the conditions
mkdir(dir1,'Excel_data_unique_peptide_values_identified_kinases')
dir_current = strcat(dir1,'\Excel_data_unique_peptide_values_identified_kinases');
recycle on % Send to recycle bin instead of permanently deleting.
list_proteins_significant_unique_identified_kinases = cell(size(list_proteins_significant_unique));
for j =  1:size(list_proteins_significant_unique,2)
    for i = 1:size(list_proteins_significant_unique,1)
        if(isempty(list_proteins_significant_unique{i,j})==0)
            list_proteins_significant_unique_identified_kinases{i,j} = get_list_kinases_origianl_data(list_proteins_significant_unique{i,j}(:,1));
            if(isempty(list_proteins_significant_unique_identified_kinases{i,j})==0)
                list_proteins_significant_unique_identified_kinases{i,j} = sortrows(list_proteins_significant_unique_identified_kinases{i,j},1);  
                aa1 = {};
                [~,n1] = unique(list_proteins_significant_unique_identified_kinases{i,j}(:,1),'stable');
                for k = 1:length(n1)
                    ind = find(ismember(list_proteins_significant_unique{i,j}(:,1),list_proteins_significant_unique_identified_kinases{i,j}{n1(k),1}));
                    aa1 = [aa1;list_proteins_significant_unique{i,j}(ind,4:end-1)];
                end
            list_proteins_significant_unique_identified_kinases{i,j} = [list_proteins_significant_unique_identified_kinases{i,j},aa1];
            list_proteins_significant_unique_identified_kinases{i,j}.Properties.VariableNames(end-2) = "Phosphorylated Residue";
            list_proteins_significant_unique_identified_kinases{i,j}.Properties.VariableNames(end-1) = "Q-Value";
            list_proteins_significant_unique_identified_kinases{i,j}.Properties.VariableNames(end) = "FC";
            end
        end
    end
end

for i = 1:size(list_proteins_significant_unique_identified_kinases,2)
    delete(strcat(dir_current,'/list_proteins_significant_unique_identified_kinases_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_unique_identified_kinases,1)
        if(isempty(list_proteins_significant_unique_identified_kinases{j,i})==0)
            writetable(list_proteins_significant_unique_identified_kinases{j,i},strcat(dir_current,'/list_proteins_significant_unique_identified_kinases_','p',num2str(i),'.xlsx'),'Sheet',sample_names{j});
        end
    end
end

list_proteins_significant_non_unique_identified_kinases = cell(size(list_proteins_significant_non_unique));
for i = 1:size(list_proteins_significant_non_unique,3)
    for j = 1:size(list_proteins_significant_non_unique,1)
        for k = 1:size(list_proteins_significant_non_unique,2)
            if(isempty(list_proteins_significant_non_unique{j,k,i})==0)
            list_proteins_significant_non_unique_identified_kinases{j,k,i} = get_list_kinases_origianl_data(list_proteins_significant_non_unique{j,k,i}(:,1));
              if(isempty(list_proteins_significant_non_unique_identified_kinases{j,k,i})==0)
                list_proteins_significant_non_unique_identified_kinases{j,k,i} = sortrows(list_proteins_significant_non_unique_identified_kinases{j,k,i},1);  
                aa1 = {};
                [~,n1] = unique(list_proteins_significant_non_unique_identified_kinases{j,k,i}(:,1),'stable');
                for jj = 1:length(n1)
                    ind = find(ismember(list_proteins_significant_non_unique{j,k,i}(:,1),list_proteins_significant_non_unique_identified_kinases{j,k,i}{n1(jj),1}));
                    aa1 = [aa1;list_proteins_significant_non_unique{j,k,i}(ind,4:end-1)];
                end
            list_proteins_significant_non_unique_identified_kinases{j,k,i} = [list_proteins_significant_non_unique_identified_kinases{j,k,i},aa1];
            list_proteins_significant_non_unique_identified_kinases{j,k,i}.Properties.VariableNames(end-2) = "Phosphorylated Residue";
            list_proteins_significant_non_unique_identified_kinases{j,k,i}.Properties.VariableNames(end-1) = "Q-Value";
            list_proteins_significant_non_unique_identified_kinases{j,k,i}.Properties.VariableNames(end) = "FC";
              end 
            end
        end
    end
end

for i = 1:size(list_proteins_significant_non_unique_identified_kinases,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_identified_kinases_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique_identified_kinases,1)
        for k = 1:size(list_proteins_significant_non_unique_identified_kinases,2)
            if(isempty(list_proteins_significant_non_unique_identified_kinases{j,k,i})==0)
                writetable(list_proteins_significant_non_unique_identified_kinases{j,k,i},strcat(dir_current,'/list_proteins_significant_non_unique_identified_kinases_','p',num2str(i),'.xlsx'),'Sheet',strcat(sample_names{j},'vs',sample_names{k}));
            end
        end
    end
end


list_proteins_significant_non_unique_upreg_identified_kinases = cell(size(list_proteins_significant_non_unique_upreg));
for i = 1:size(list_proteins_significant_non_unique_upreg,3)
    for j = 1:size(list_proteins_significant_non_unique_upreg,1)
        for k = 1:size(list_proteins_significant_non_unique_upreg,2)
            if(isempty(list_proteins_significant_non_unique_upreg{j,k,i})==0)
                list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i} = get_list_kinases_origianl_data(list_proteins_significant_non_unique_upreg{j,k,i}(:,1));
                if(isempty(list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i})==0)
                    list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i} = sortrows(list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i},1);
                    aa1 = {};
                    [~,n1] = unique(list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i}(:,1),'stable');
                    for jj = 1:length(n1)
                        ind = find(ismember(list_proteins_significant_non_unique_upreg{j,k,i}(:,1),list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i}{n1(jj),1}));
                        aa1 = [aa1;list_proteins_significant_non_unique_upreg{j,k,i}(ind,4:end-1)];
                    end
                    list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i} = [list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i},aa1];
                    list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i}.Properties.VariableNames(end-2) = "Phosphorylated Residue";
                    list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i}.Properties.VariableNames(end-1) = "Q-Value";
                    list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i}.Properties.VariableNames(end) = "FC";
                end
            end
        end
    end
end

for i = 1:size(list_proteins_significant_non_unique_upreg_identified_kinases,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_upreg_identified_kinases_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique_upreg_identified_kinases,1)
        for k = 1:size(list_proteins_significant_non_unique_upreg_identified_kinases,2)
            if(isempty(list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i})==0)
                writetable(list_proteins_significant_non_unique_upreg_identified_kinases{j,k,i},strcat(dir_current,'/list_proteins_significant_non_unique_upreg_identified_kinases_','p',num2str(i),'.xlsx'),'Sheet',strcat(sample_names{j},'vs',sample_names{k}));
            end
        end
    end
end

list_proteins_significant_non_unique_downreg_identified_kinases = cell(size(list_proteins_significant_non_unique_downreg));
for i = 1:size(list_proteins_significant_non_unique_downreg,3)
    for j = 1:size(list_proteins_significant_non_unique_downreg,1)
        for k = 1:size(list_proteins_significant_non_unique_downreg,2)
            if(isempty(list_proteins_significant_non_unique_downreg{j,k,i})==0)
                list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i} = get_list_kinases_origianl_data(list_proteins_significant_non_unique_downreg{j,k,i}(:,1));
                if(isempty(list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i})==0)
                    list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i} = sortrows(list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i},1);
                    aa1 = {};
                    [~,n1] = unique(list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i}(:,1),'stable');
                    for jj = 1:length(n1)
                        ind = find(ismember(list_proteins_significant_non_unique_downreg{j,k,i}(:,1),list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i}{n1(jj),1}));
                        aa1 = [aa1;list_proteins_significant_non_unique_downreg{j,k,i}(ind,4:end-1)];
                    end
                    list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i} = [list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i},aa1];
                    list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i}.Properties.VariableNames(end-2) = "Phosphorylated Residue";
                    list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i}.Properties.VariableNames(end-1) = "Q-Value";
                    list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i}.Properties.VariableNames(end) = "FC";
                end
            end
        end
    end
end


for i = 1:size(list_proteins_significant_non_unique_downreg_identified_kinases,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_downreg_identified_kinases_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique_downreg_identified_kinases,1)
        for k = 1:size(list_proteins_significant_non_unique_downreg_identified_kinases,2)
            if(isempty(list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i})==0)
                writetable(list_proteins_significant_non_unique_downreg_identified_kinases{j,k,i},strcat(dir_current,'/list_proteins_significant_non_unique_downreg_identified_kinases_','p',num2str(i),'.xlsx'),'Sheet',strcat(sample_names{j},'vs',sample_names{k}));
            end
        end
    end
end

list_proteins_significant_unique_complete_kinases = cell(size(list_proteins_significant_unique_complete));
for i = 1:length(list_proteins_significant_unique_complete)
    if(isempty(list_proteins_significant_unique_complete{i})==0)
        list_proteins_significant_unique_complete_kinases{i} = get_list_kinases_origianl_data(list_proteins_significant_unique_complete{i}(:,1));
        if(isempty(list_proteins_significant_unique_complete_kinases{i})==0)
            list_proteins_significant_unique_complete_kinases{i} = sortrows(list_proteins_significant_unique_complete_kinases{i},1);
            aa1 = {};
            [~,n1] = unique(list_proteins_significant_unique_complete_kinases{i}(:,1),'stable');
            for jj = 1:length(n1)
                ind = find(ismember(list_proteins_significant_unique_complete{i}(:,1),list_proteins_significant_unique_complete_kinases{i}{n1(jj),1}));
                aa1 = [aa1;list_proteins_significant_unique_complete{i}(ind,4:end-1)];
            end
            list_proteins_significant_unique_complete_kinases{i} = [list_proteins_significant_unique_complete_kinases{i},aa1];
            list_proteins_significant_unique_complete_kinases{i}.Properties.VariableNames(end-2) = "Phosphorylated Residue";
            list_proteins_significant_unique_complete_kinases{i}.Properties.VariableNames(end-1) = "Q-Value";
            list_proteins_significant_unique_complete_kinases{i}.Properties.VariableNames(end) = "FC";
        end
    end
end


delete(strcat(dir_current,'/list_proteins_significant_unique_complete_identified_kinases.xlsx'));
for i = 1:length(list_proteins_significant_unique_complete_kinases)
    if(isempty(list_proteins_significant_unique_complete_kinases{i})==0)
        writetable(list_proteins_significant_unique_complete_kinases{i},strcat(dir_current,'/list_proteins_significant_unique_complete_identified_kinases.xlsx'),'Sheet',sample_names{i});
    end
end

list_proteins_significant_non_unique_downreg_complete_kinase = cell(size(list_proteins_significant_non_unique_downreg_complete));
for i = 1:size(list_proteins_significant_non_unique_downreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete,2)
        if(isempty(list_proteins_significant_non_unique_downreg_complete{i,j})==0)
            list_proteins_significant_non_unique_downreg_complete_kinase{i,j} = get_list_kinases_origianl_data(list_proteins_significant_non_unique_downreg_complete{i,j}(:,1));
            if(isempty(list_proteins_significant_non_unique_downreg_complete_kinase{i,j})==0)
            list_proteins_significant_non_unique_downreg_complete_kinase{i,j} = sortrows(list_proteins_significant_non_unique_downreg_complete_kinase{i,j},1);
            aa1 = {};
            [~,n1] = unique(list_proteins_significant_non_unique_downreg_complete_kinase{i,j}(:,1),'stable');
            for jj = 1:length(n1)
                ind = find(ismember(list_proteins_significant_non_unique_downreg_complete{i,j}(:,1),list_proteins_significant_non_unique_downreg_complete_kinase{i,j}{n1(jj),1}));
                aa1 = [aa1;list_proteins_significant_non_unique_downreg_complete{i,j}(ind,4:end-1)];
            end
            list_proteins_significant_non_unique_downreg_complete_kinase{i,j} = [list_proteins_significant_non_unique_downreg_complete_kinase{i,j},aa1];
            list_proteins_significant_non_unique_downreg_complete_kinase{i,j}.Properties.VariableNames(end-2) = "Phosphorylated Residue";
            list_proteins_significant_non_unique_downreg_complete_kinase{i,j}.Properties.VariableNames(end-1) = "Q-Value";
            list_proteins_significant_non_unique_downreg_complete_kinase{i,j}.Properties.VariableNames(end) = "FC";
        end
        end
    end
end

delete(strcat(dir_current,'/list_proteins_significant_non_unique_downreg_complete_identified_kinases.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_downreg_complete_kinase,1)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete_kinase,2)
        if(isempty(list_proteins_significant_non_unique_downreg_complete_kinase{i,j})==0)
            writetable(list_proteins_significant_non_unique_downreg_complete_kinase{i,j},strcat(dir_current,'/list_proteins_significant_non_unique_downreg_complete_identified_kinases.xlsx'),'Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end

list_proteins_significant_non_unique_upreg_complete_kinases = cell(size(list_proteins_significant_non_unique_upreg_complete));
for i = 1:size(list_proteins_significant_non_unique_upreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete,2)
        if(isempty(list_proteins_significant_non_unique_upreg_complete{i,j})==0)
            list_proteins_significant_non_unique_upreg_complete_kinases{i,j} = get_list_kinases_origianl_data(list_proteins_significant_non_unique_upreg_complete{i,j}(:,1));
            if(isempty(list_proteins_significant_non_unique_upreg_complete_kinases{i,j})==0)
            list_proteins_significant_non_unique_upreg_complete_kinases{i,j} = sortrows(list_proteins_significant_non_unique_upreg_complete_kinases{i,j},1);
            aa1 = {};
            [~,n1] = unique(list_proteins_significant_non_unique_upreg_complete_kinases{i,j}(:,1),'stable');
            for jj = 1:length(n1)
                ind = find(ismember(list_proteins_significant_non_unique_upreg_complete{i,j}(:,1),list_proteins_significant_non_unique_upreg_complete_kinases{i,j}{n1(jj),1}));
                aa1 = [aa1;list_proteins_significant_non_unique_upreg_complete{i,j}(ind,4:end-1)];
            end
            list_proteins_significant_non_unique_upreg_complete_kinases{i,j} = [list_proteins_significant_non_unique_upreg_complete_kinases{i,j},aa1];
            list_proteins_significant_non_unique_upreg_complete_kinases{i,j}.Properties.VariableNames(end-2) = "Phosphorylated Residue";
            list_proteins_significant_non_unique_upreg_complete_kinases{i,j}.Properties.VariableNames(end-1) = "Q-Value";
            list_proteins_significant_non_unique_upreg_complete_kinases{i,j}.Properties.VariableNames(end) = "FC";
            end  
        end
    end
end

delete(strcat(dir_current,'/list_proteins_significant_non_unique_upreg_complete_identified_kinases.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_upreg_complete_kinases,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete_kinases,2)
        if(isempty(list_proteins_significant_non_unique_upreg_complete_kinases{i,j})==0)
            writetable(list_proteins_significant_non_unique_upreg_complete_kinases{i,j},strcat(dir_current,'/list_proteins_significant_non_unique_upreg_complete_identified_kinases.xlsx'),'Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end



%Gropu comparison excel data!!!!!
delete(strcat(dir_current,'/list_proteins_significant_group_comparison_kinases.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_upreg_complete_kinases,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete_kinases,2)
        if(isempty(list_proteins_significant_non_unique_upreg_complete_kinases{i,j})==0)
            aas = [list_proteins_significant_unique_complete_kinases{i};list_proteins_significant_non_unique_upreg_complete_kinases{i,j}];
            writetable(aas,strcat(dir_current,'/list_proteins_significant_group_comparison_kinases.xlsx'),'Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end

%The index is inveresed because downregulated means inverse upregulated
for i = 1:size(list_proteins_significant_non_unique_downreg_complete_kinase,2)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete_kinase,1)
        if(isempty(list_proteins_significant_non_unique_downreg_complete_kinase{j,i})==0)
            aas = [list_proteins_significant_unique_complete_kinases{i};list_proteins_significant_non_unique_downreg_complete_kinase{j,i}];
            writetable(aas,strcat(dir_current,'/list_proteins_significant_group_comparison_kinases.xlsx'),'Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end






