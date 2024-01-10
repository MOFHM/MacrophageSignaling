function [list_proteins_significant_non_unique,list_proteins_significant_non_unique_upreg,list_proteins_significant_non_unique_downreg,list_proteins_significant_unique,list_proteins_significant_unique_complete,list_proteins_significant_non_unique_upreg_complete,list_proteins_significant_non_unique_downreg_complete] = list_generation_proteinsKSEA1(nr_sample_groups,nr_analyzed_phos,data_unique,nr_rep,p_names_unique,p_names_p,position_protein_p,position_protein_unique,fc,p_values_total,baseMean,baseMean_unique,significant_peptides_index,significant_peptides_index_upreg,significant_peptides_index_downreg,Consider_unique_sepparately,p_values_noncorrected,measured_values,measured_values_unique,loc_prob_p,loc_prob_unique,peptide_sequence_p,peptide_sequence_unique,stder_total,df_pooled_total,prot_comp_flag)
%Function that generated the matrix that contains the list of significantly
%upreugalted proteins per conditions

pat_ini_protein = 'sp|';
pat_fin_protein = '|';
pat_ini_gene = 'sp|'+ wildcardPattern + '|';
pat_fin_gene = '_HUMAN' ;

%The list of proteins that are significant between the comparison groups
list_proteins_significant_non_unique = cell(nr_sample_groups,nr_sample_groups,nr_analyzed_phos);
for i = 1:nr_analyzed_phos
    iji = 0;
    for j = 1:nr_sample_groups
        for k = j+1:nr_sample_groups
            if(isempty(significant_peptides_index{j,k,i})==0)
                iji = iji+1;
                for kk = 1:length(significant_peptides_index{j,k,i})
                    if(contains(p_names_p{i}{significant_peptides_index{j,k,i}(kk)},"cont",'IgnoreCase',true)==0)
                        el1 = extractBetween(p_names_p{i}{significant_peptides_index{j,k,i}(kk)},pat_ini_protein,pat_fin_protein);
                        el2 = extractBetween(p_names_p{i}{significant_peptides_index{j,k,i}(kk)},pat_ini_gene,pat_fin_gene);
                        ind1 = size(el2,1);
                        el3 = repmat(peptide_sequence_p{i}{significant_peptides_index{j,k,i}(kk)},length(el1),1);
                        el4 = strsplit(position_protein_p{i}{significant_peptides_index{j,k,i}(kk)},';')';
                        el5 = repmat(p_values_total{i,iji}(significant_peptides_index{j,k,i}(kk)),ind1,1);
                        el6 = repmat({fc{i,iji}(significant_peptides_index{j,k,i}(kk))},ind1,1);
                        el7 = repmat({baseMean{i,iji}(significant_peptides_index{j,k,i}(kk))},ind1,1);
                        el8 = repmat(p_values_noncorrected{i,iji}(significant_peptides_index{j,k,i}(kk)),ind1,1);
                        ams = strsplit(cell2mat(measured_values{i,iji}(significant_peptides_index{j,k,i}(kk))),'||');
                        el9 = repmat(ams{1},ind1,1);
                        el10 = repmat(ams{2},ind1,1);
                        el11 = repmat(loc_prob_p{i}{significant_peptides_index{j,k,i}(kk)},ind1,1);
                        el12 = repmat(stder_total{i,iji}(significant_peptides_index{j,k,i}(kk)),ind1,1);
                        el13 = repmat(df_pooled_total{i,iji}(significant_peptides_index{j,k,i}(kk)),ind1,1);
                        el14 = repmat(prot_comp_flag{i,iji}(significant_peptides_index{j,k,i}(kk)),ind1,1);
                        list_proteins_significant_non_unique{j,k,i} = [list_proteins_significant_non_unique{j,k,i};[el1,el2,el3,el4,el5,el6,el7,el8,el9,el10,el11,el12,el13,el14]];
                    end
                end
            end
        end
    end
end

list_proteins_significant_non_unique_upreg = cell(nr_sample_groups,nr_sample_groups,nr_analyzed_phos);
for i = 1:nr_analyzed_phos
    iji = 0;
    for j = 1:nr_sample_groups
        for k = j+1:nr_sample_groups
            if(isempty(significant_peptides_index_upreg{j,k,i})==0)
                iji = iji+1;
                for kk = 1:length(significant_peptides_index_upreg{j,k,i})
                    if(contains(p_names_p{i}{significant_peptides_index_upreg{j,k,i}(kk)},"cont",'IgnoreCase',true)==0)
                        el1 = extractBetween(p_names_p{i}{significant_peptides_index_upreg{j,k,i}(kk)},pat_ini_protein,pat_fin_protein);
                        el2 = extractBetween(p_names_p{i}{significant_peptides_index_upreg{j,k,i}(kk)},pat_ini_gene,pat_fin_gene);
                        ind1 = size(el2,1);
                        el3 = repmat(peptide_sequence_p{i}{significant_peptides_index_upreg{j,k,i}(kk)},length(el1),1);
                        el4 = strsplit(position_protein_p{i}{significant_peptides_index_upreg{j,k,i}(kk)},';')';
                        el5 = repmat(p_values_total{i,iji}(significant_peptides_index_upreg{j,k,i}(kk)),ind1,1);
                        el6 = repmat({fc{i,iji}(significant_peptides_index_upreg{j,k,i}(kk))},ind1,1);
                        el7 = repmat({baseMean{i,iji}(significant_peptides_index_upreg{j,k,i}(kk))},ind1,1);
                        el8 = repmat(p_values_noncorrected{i,iji}(significant_peptides_index_upreg{j,k,i}(kk)),ind1,1);
                        ams = strsplit(cell2mat(measured_values{i,iji}(significant_peptides_index_upreg{j,k,i}(kk))),'||');
                        el9 = repmat(ams{1},ind1,1);
                        el10 = repmat(ams{2},ind1,1);
                        el11 = repmat(loc_prob_p{i}{significant_peptides_index_upreg{j,k,i}(kk)},ind1,1);
                        el12 = repmat(stder_total{i,iji}(significant_peptides_index_upreg{j,k,i}(kk)),ind1,1);
                        el13 = repmat(df_pooled_total{i,iji}(significant_peptides_index_upreg{j,k,i}(kk)),ind1,1);
                        el14 = repmat(prot_comp_flag{i,iji}(significant_peptides_index_upreg{j,k,i}(kk)),ind1,1);
                        el15 = repmat(num2str(i),ind1,1);
                        list_proteins_significant_non_unique_upreg{j,k,i} = [list_proteins_significant_non_unique_upreg{j,k,i};[el1,el2,el3,el4,el5,el6,el7,el8,el9,el10,el11,el12,el13,el14,el15]];
                    end
                end
            end
        end
    end
end

list_proteins_significant_non_unique_downreg = cell(nr_sample_groups,nr_sample_groups,nr_analyzed_phos);
for i = 1:nr_analyzed_phos
    iji = 0;
    for j = 1:nr_sample_groups
        for k = j+1:nr_sample_groups
            if(isempty(significant_peptides_index_downreg{j,k,i})==0)
                iji = iji+1;
                for kk = 1:length(significant_peptides_index_downreg{j,k,i})
                    if(contains(p_names_p{i}{significant_peptides_index_downreg{j,k,i}(kk)},"cont",'IgnoreCase',true)==0)
                        el1 = extractBetween(p_names_p{i}{significant_peptides_index_downreg{j,k,i}(kk)},pat_ini_protein,pat_fin_protein);
                        el2 = extractBetween(p_names_p{i}{significant_peptides_index_downreg{j,k,i}(kk)},pat_ini_gene,pat_fin_gene);
                        ind1 = size(el2,1);
                        el3 = repmat(peptide_sequence_p{i}{significant_peptides_index_downreg{j,k,i}(kk)},length(el1),1);
                        el4 = strsplit(position_protein_p{i}{significant_peptides_index_downreg{j,k,i}(kk)},';')';
                        el5 = repmat(p_values_total{i,iji}(significant_peptides_index_downreg{j,k,i}(kk)),ind1,1);
                        el6 = repmat({fc{i,iji}(significant_peptides_index_downreg{j,k,i}(kk))},ind1,1);
                        el7 = repmat({baseMean{i,iji}(significant_peptides_index_downreg{j,k,i}(kk))},ind1,1);
                        el8 = repmat(p_values_noncorrected{i,iji}(significant_peptides_index_downreg{j,k,i}(kk)),ind1,1);
                        ams = strsplit(cell2mat(measured_values{i,iji}(significant_peptides_index_downreg{j,k,i}(kk))),'||');
                        el9 = repmat(ams{2},ind1,1);
                        el10 = repmat(ams{1},ind1,1);
                        el11 = repmat(loc_prob_p{i}{significant_peptides_index_downreg{j,k,i}(kk)},ind1,1);
                        el12 = repmat(stder_total{i,iji}(significant_peptides_index_downreg{j,k,i}(kk)),ind1,1);
                        el13 = repmat(df_pooled_total{i,iji}(significant_peptides_index_downreg{j,k,i}(kk)),ind1,1);
                        el14 = repmat(prot_comp_flag{i,iji}(significant_peptides_index_downreg{j,k,i}(kk)),ind1,1);
                        el15 = repmat(num2str(i),ind1,1);
                        list_proteins_significant_non_unique_downreg{j,k,i} = [list_proteins_significant_non_unique_downreg{j,k,i};[el1,el2,el3,el4,el5,el6,el7,el8,el9,el10,el11,el12,el13,el14,el15]];
                    end
                end
            end
        end
    end
end


%The list of proteins that are unique for only one specific condition
if(Consider_unique_sepparately==1)
    list_proteins_significant_unique = cell(nr_sample_groups,nr_analyzed_phos);
    list_positions_significant_unique = cell(nr_sample_groups,nr_analyzed_phos);
    baseMean_unique1 = cell(nr_sample_groups,nr_analyzed_phos);
    measured_values_unique = cell(nr_sample_groups,nr_analyzed_phos);
    peptide_sequence_unique1 = cell(nr_sample_groups,nr_analyzed_phos);
    loc_prob_unique1 = cell(nr_sample_groups,nr_analyzed_phos);
    for i = 1:nr_analyzed_phos
        for j = 1:nr_rep:size(data_unique{i},2)
            ind_j = (j-1)/nr_rep+1;
            for k = 1:size(data_unique{i},1)
                vec = data_unique{i}(k,j:j+nr_rep-1);
                ind = find(vec==0);
                if(length(ind)<nr_rep)
                    list_proteins_significant_unique{ind_j,i} = [list_proteins_significant_unique{ind_j,i};p_names_unique{i}(k)];
                    list_positions_significant_unique{ind_j,i} = [list_positions_significant_unique{ind_j,i};position_protein_unique{i}(k)];
                    baseMean_unique1{ind_j,i} = [baseMean_unique1{ind_j,i};baseMean_unique{i}(k)];
                    measured_values_unique{ind_j,i} = [measured_values_unique{ind_j,i};measured_values_unique{i}(k)];
                    peptide_sequence_unique1{ind_j,i} = [peptide_sequence_unique{ind_j,i};peptide_sequence_unique{i}(k)];
                    loc_prob_unique1{ind_j,i} = [loc_prob_unique{ind_j,i};loc_prob_unique{i}(k)];
                end
            end
        end
    end
    
    ll = cell(size(list_proteins_significant_unique));
    for i = 1:size(list_proteins_significant_unique,1)
        for j = 1:size(list_proteins_significant_unique,2)
            if(isempty(list_proteins_significant_unique{i,j})==0)
                for kk = 1:length(list_proteins_significant_unique{i,j})
                    if(contains(list_proteins_significant_unique{i,j}{kk},"cont",'IgnoreCase',true)==0)
                        el1 = extractBetween(list_proteins_significant_unique{i,j}{kk},pat_ini_protein,pat_fin_protein);
                        el2 = extractBetween(list_proteins_significant_unique{i,j}{kk},pat_ini_gene,pat_fin_gene);
                        ind1 = size(el2,1);
                        el3 = repmat(peptide_sequence_unique1{i,j}{kk},length(el1),1);
                        el4 = strsplit(list_positions_significant_unique{i,j}{kk},';')';
                        el5 = repmat({'NULL'},length(el1),1);
                        el6 = repmat({-1},ind1,1);
                        el7 = repmat({baseMean_unique1{i,j}(kk)},length(el1),1);
                        el8 = repmat({'NULL'},length(el1),1);
                        el9 = repmat({measured_values_unique{i,j}(kk)},length(el1),1);
                        el10 = repmat({'NULL'},length(el1),1);
                        el11 = repmat(loc_prob_unique1{i,j}{kk},ind1,1);
                        el12 = repmat({'NULL'},length(el1),1);
                        el13 = repmat({'NULL'},length(el1),1);
                        el14 = repmat({'NULL'},length(el1),1);
                        el15 = repmat(num2str(j),ind1,1);
                        ll{i,j} = [ll{i,j};[el1,el2,el3,el4,el5,el6,el7,el8,el9,el10,el11,el12,el13,el14,el15]];
                    end
                end
            end
        end
    end
    list_proteins_significant_unique = ll;
    list_proteins_significant_unique_complete = cell(1,size(list_proteins_significant_unique,1));
    iji = 1;
    pval_aux = [];
    for i = 1:size(list_proteins_significant_unique,1)
        for j = 1:size(list_proteins_significant_unique,2)
            if(isempty(list_proteins_significant_non_unique_upreg{i,j,k})==0)
                pval_aux = [pval_aux;p_values_noncorrected{k,iji}];
            end
            list_proteins_significant_unique_complete{i} = [list_proteins_significant_unique_complete{i};list_proteins_significant_unique{i,j}];
            %Give a value of zero (do not compesate) for the residues that
            %comes from a protein that have at least one residue for which
            %the compensation is not neccessary
            aas = list_proteins_significant_unique_complete{i};
            ar = find(cellfun(@str2num,aas(:,14))==1);
            ar1 = ar(find(ismember(aas(ar,1),aas(find(cellfun(@str2num,aas(:,14))==0),1))));
            ar2 = aas(:,14); ar2(ar1) = {'0'};
            list_proteins_significant_unique_complete{i} = [aas,ar2];
        end
        if(isempty(list_proteins_significant_unique_complete{i})==0)
            pval_aux(:,2) = num2cell(mafdr([pval_aux{:,1}],'BHFDR','TRUE'));
            list_proteins_significant_unique_complete{i,j}(:,5) = pval_aux(find(ismember([pval_aux{:,1}],[list_proteins_significant_unique_complete{i,j}{:,8}])),2);
            list_proteins_significant_unique_complete{i}(find([list_proteins_significant_unique_complete{i}{:,5}]>0.05),:)=[];
            pval_aux = [];
            iji = iji + 1;
        end
    end
    
else list_proteins_significant_unique_complete = cell(1,size(list_proteins_significant_non_unique_upreg,1)); list_proteins_significant_unique = cell(nr_sample_groups,nr_analyzed_phos);
end


% for i = 1:size(list_proteins_significant_unique_complete,2)
%  if(isempty(list_proteins_significant_unique_complete{i})==0)
%         [~,aa,~] = unique(list_proteins_significant_unique_complete{i}(:,1),'stable');
% list_proteins_significant_unique_complete{i} = list_proteins_significant_unique_complete{i}(aa,:);
%  end
% end

list_proteins_significant_non_unique_upreg_complete = cell(size(list_proteins_significant_non_unique_upreg,1),size(list_proteins_significant_non_unique_upreg,2));
iji = 1;
pval_aux = [];
for i = 1:size(list_proteins_significant_non_unique_upreg,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg,2)
        for k = 1:size(list_proteins_significant_non_unique_upreg,3)
            if(isempty(list_proteins_significant_non_unique_upreg{i,j,k})==0)
                pval_aux = [pval_aux;p_values_noncorrected{k,iji}];
            end
            list_proteins_significant_non_unique_upreg_complete{i,j} = [list_proteins_significant_non_unique_upreg_complete{i,j};list_proteins_significant_non_unique_upreg{i,j,k}];
        end
        if(isempty(list_proteins_significant_non_unique_upreg_complete{i,j})==0)
            %Give a value of zero (do not compesate) for the residues that
            %comes from a protein that have at least one residue for which
            %the compensation is not neccessary
            aas = list_proteins_significant_non_unique_upreg_complete{i,j};
            ar = find(cellfun(@str2num,aas(:,14))==1);
            ar1 = ar(find(ismember(aas(ar,1),aas(find(cellfun(@str2num,aas(:,14))==0),1))));
            ar2 = aas(:,14); ar2(ar1) = {'0'};
            list_proteins_significant_non_unique_upreg_complete{i,j} = [aas,ar2];
            pval_aux(:,2) = num2cell(mafdr([pval_aux{:,1}],'BHFDR','TRUE'));
            list_proteins_significant_non_unique_upreg_complete{i,j}(:,5) = pval_aux(find(ismember([pval_aux{:,1}],[list_proteins_significant_non_unique_upreg_complete{i,j}{:,8}])),2);
            list_proteins_significant_non_unique_upreg_complete{i,j}(find([list_proteins_significant_non_unique_upreg_complete{i,j}{:,5}]>0.05),:)=[];
            pval_aux = [];
            iji = iji + 1;
        end
    end
end
% for i = 1:size(list_proteins_significant_non_unique_upreg_complete,1)
%     for j = 1:size(list_proteins_significant_non_unique_upreg_complete,2)
%  if(isempty(list_proteins_significant_non_unique_upreg_complete{i,j})==0)
%         [~,aa,~] = unique(list_proteins_significant_non_unique_upreg_complete{i,j}(:,1),'stable');
% list_proteins_significant_non_unique_upreg_complete{i,j} = list_proteins_significant_non_unique_upreg_complete{i,j}(aa,:);
%  end
%     end
% end

list_proteins_significant_non_unique_downreg_complete = cell(size(list_proteins_significant_non_unique_downreg,1),size(list_proteins_significant_non_unique_downreg,2));
iji = 1;
pval_aux = [];
for i = 1:size(list_proteins_significant_non_unique_downreg,1)
    for j = 1:size(list_proteins_significant_non_unique_downreg,2)
        for k = 1:size(list_proteins_significant_non_unique_downreg,3)
            if(isempty(list_proteins_significant_non_unique_upreg{i,j,k})==0)
                pval_aux = [pval_aux;p_values_noncorrected{k,iji}];
            end
            list_proteins_significant_non_unique_downreg_complete{i,j} = [list_proteins_significant_non_unique_downreg_complete{i,j};list_proteins_significant_non_unique_downreg{i,j,k}];
        end
        if(isempty(list_proteins_significant_non_unique_downreg_complete{i,j})==0)
            %Give a value of zero (do not compesate) for the residues that
            %comes from a protein that have at least one residue for which
            %the compensation is not neccessary
            aas = list_proteins_significant_non_unique_downreg_complete{i,j};
            ar = find(cellfun(@str2num,aas(:,14))==1);
            ar1 = ar(find(ismember(aas(ar,1),aas(find(cellfun(@str2num,aas(:,14))==0),1))));
            ar2 = aas(:,14); ar2(ar1) = {'0'};
            list_proteins_significant_non_unique_downreg_complete{i,j} = [aas,ar2];
            pval_aux(:,2) = num2cell(mafdr([pval_aux{:,1}],'BHFDR','TRUE'));
            list_proteins_significant_non_unique_downreg_complete{i,j}(:,5) = pval_aux(find(ismember([pval_aux{:,1}],[list_proteins_significant_non_unique_downreg_complete{i,j}{:,8}])),2);
            list_proteins_significant_non_unique_downreg_complete{i,j}(find([list_proteins_significant_non_unique_downreg_complete{i,j}{:,5}]>0.05),:)=[];
            pval_aux = [];
            iji = iji + 1;
        end
    end
end
% for i = 1:size(list_proteins_significant_non_unique_downreg_complete,1)
%     for j = 1:size(list_proteins_significant_non_unique_downreg_complete,2)
%         if(isempty(list_proteins_significant_non_unique_downreg_complete{i,j})==0)
%         [~,aa,~] = unique(list_proteins_significant_non_unique_downreg_complete{i,j}(:,1),'stable');
% list_proteins_significant_non_unique_downreg_complete{i,j} = list_proteins_significant_non_unique_downreg_complete{i,j}(aa,:);
%         end
%     end
% end

for i = 1:size(list_proteins_significant_non_unique_upreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete,2)
        
        if(isempty(list_proteins_significant_non_unique_upreg_complete{i,j})==0)
            a1 = cell2table(list_proteins_significant_non_unique_upreg_complete{i,j});
            [~,aa,~] = unique(a1(:,[1,4]),'stable');
            a1 = a1(aa,:);
            list_proteins_significant_non_unique_upreg_complete{i,j} = table2cell(a1);
        end
        
        if(isempty(list_proteins_significant_non_unique_downreg_complete{i,j})==0)
            a2 = cell2table(list_proteins_significant_non_unique_downreg_complete{i,j});
            [~,aa,~] = unique(a2(:,[1,4]),'stable');
            a2 = a2(aa,:);
            list_proteins_significant_non_unique_downreg_complete{i,j} = table2cell(a2);
        end
        
        if(union(isempty(list_proteins_significant_non_unique_upreg_complete{i,j}), isempty(list_proteins_significant_non_unique_downreg_complete{i,j}))==0)
            a1 = cell2table(list_proteins_significant_non_unique_upreg_complete{i,j});
            a2 = cell2table(list_proteins_significant_non_unique_downreg_complete{i,j});
            ind1 = find(ismember(a1(:,[1,4]),a2(:,[1,4])));
            ind2 = find(ismember(a2(:,[1,4]),a1(:,[1,4])));
            a1(ind1,:) = [];
            a2(ind2,:) = [];
            list_proteins_significant_non_unique_upreg_complete{i,j} = table2cell(a1);
            list_proteins_significant_non_unique_downreg_complete{i,j} = table2cell(a2);
        end
        
    end
end


end