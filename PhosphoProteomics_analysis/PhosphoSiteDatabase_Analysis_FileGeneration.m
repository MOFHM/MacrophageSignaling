mkdir(dir1,'Excel_data_PhosphoSite_database_analysis')
dir_current = strcat(dir1,'\Excel_data_PhosphoSite_database_analysis');
recycle on % Send to recycle bin instead of permanently deleting.


delete(strcat(dir_current,'/list_proteins_significant_unique_complete.xlsx'));
for i = 1:size(list_proteins_significant_unique_complete,2)
    if(isempty(list_proteins_significant_unique_complete{i})==0)
        %p_names_comaprison_external_list = {};
        %position_protein_external = {};
        p_names_comparison_external_list = list_proteins_significant_unique_complete{i}(:,1);
        position_protein_external = list_proteins_significant_unique_complete{i}(:,4);
        if(PSP_Loc_prob==1)
            index_to_remove1 = find(cellfun(@str2num,list_proteins_significant_unique_complete{i}(:,11))<0.75);
        else index_to_remove1 = []; end
        if(PSP_protein_compensation==1)
            index_to_remove2 = find(cellfun(@str2num,list_proteins_significant_unique_complete{i}(:,14))==1);
        else index_to_remove2 = []; end
        index_to_remove = union(index_to_remove1,index_to_remove2);
        index_to_remove3 = find(cell2mat(list_proteins_significant_unique_complete{i}(:,5))>PSP_p_val);
        index_to_remove = union(index_to_remove,index_to_remove3);
        position_protein_external(index_to_remove) = [];
        p_names_comparison_external_list(index_to_remove) = [];
        index_to_remove = find(cell2mat(list_proteins_significant_unique_complete{i}(:,5))>PSP_p_val);
        [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
        [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
        score_kin_subExternal = [score_kin_subExternal,statistics];
        header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
        writecell(header,strcat(dir_current,'/list_kinases_significant_unique_complete.xlsx'),'Range','A1','Sheet',sample_names{i},'AutofitWidth',1);
        writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_unique_complete.xlsx'),'Range','A2','Sheet',sample_names{i},'AutofitWidth',0);
    end
end

delete(strcat(dir_current,'/list_proteins_significant_non_unique_upreg_complete.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_upreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete,2)
        if(isempty(list_proteins_significant_non_unique_upreg_complete{i,j})==0)
            p_names_comparison_external_list = list_proteins_significant_non_unique_upreg_complete{i,j}(:,1);
            position_protein_external = list_proteins_significant_non_unique_upreg_complete{i,j}(:,4);
            if(PSP_Loc_prob==1)
                index_to_remove1 = find(cellfun(@str2num,list_proteins_significant_non_unique_upreg_complete{i,j}(:,11))<0.75);
            else index_to_remove1 = []; end
            if(PSP_protein_compensation==1)
                index_to_remove2 = find(cellfun(@str2num,list_proteins_significant_non_unique_upreg_complete{i,j}(:,14))==1);
            else index_to_remove2 = []; end
            index_to_remove = union(index_to_remove1,index_to_remove2);
            index_to_remove3 = find(cell2mat(list_proteins_significant_non_unique_upreg_complete{i,j}(:,5))>PSP_p_val);
            index_to_remove = union(index_to_remove,index_to_remove3);
            position_protein_external(index_to_remove) = [];
            p_names_comparison_external_list(index_to_remove) = [];
            [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
            [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
            score_kin_subExternal = [score_kin_subExternal,statistics];
            header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
            writecell(header,strcat(dir_current,'/list_kinases_significant_non_unique_upreg_complete.xlsx'),'Range','A1','Sheet',strcat(sample_names{i},'vs',sample_names{j}),'AutofitWidth',1);
            writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_non_unique_upreg_complete.xlsx'),'Range','A2','Sheet',strcat(sample_names{i},'vs',sample_names{j}),'AutofitWidth',0);
        end
    end
end

delete(strcat(dir_current,'/list_proteins_significant_non_unique_downreg_complete.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_downreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete,2)
        if(isempty(list_proteins_significant_non_unique_downreg_complete{i,j})==0)
            p_names_comparison_external_list = list_proteins_significant_non_unique_downreg_complete{i,j}(:,1);
            position_protein_external = list_proteins_significant_non_unique_downreg_complete{i,j}(:,4);
            if(PSP_Loc_prob==1)
                index_to_remove1 = find(cellfun(@str2num,list_proteins_significant_non_unique_downreg_complete{i,j}(:,11))<0.75);
            else index_to_remove1 = []; end
            if(PSP_protein_compensation==1)
                index_to_remove2 = find(cellfun(@str2num,list_proteins_significant_non_unique_downreg_complete{i,j}(:,14))==1);
            else index_to_remove2 =[]; end
            index_to_remove = union(index_to_remove1,index_to_remove2);
            index_to_remove3 = find(cell2mat(list_proteins_significant_non_unique_downreg_complete{i,j}(:,5))>PSP_p_val);
            index_to_remove = union(index_to_remove,index_to_remove3);
            position_protein_external(index_to_remove) = [];
            p_names_comparison_external_list(index_to_remove) = [];
            [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
            [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
            score_kin_subExternal = [score_kin_subExternal,statistics];
            header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
            writecell(header,strcat(dir_current,'/list_kinases_significant_non_unique_downreg_complete.xlsx'),'Range','A1','Sheet',strcat(sample_names{i},'vs',sample_names{j}),'AutofitWidth',1);
            writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_non_unique_downreg_complete.xlsx'),'Range','A2','Sheet',strcat(sample_names{i},'vs',sample_names{j}),'AutofitWidth',0);
        end
    end
end


for i = 1:size(list_proteins_significant_non_unique_downreg,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_downreg_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique_downreg,1)
        for k = 1:size(list_proteins_significant_non_unique_downreg,2)
            if(isempty(list_proteins_significant_non_unique_downreg{j,k,i})==0)
                p_names_comparison_external_list = list_proteins_significant_non_unique_downreg{j,k,i}(:,1);
                position_protein_external = list_proteins_significant_non_unique_downreg{j,k,i}(:,4);
                if(PSP_Loc_prob==1)
                    index_to_remove1 = find(cellfun(@str2num,list_proteins_significant_non_unique_downreg{j,k,i}(:,11))<0.75);
                else index_to_remove1 = []; end
                if(PSP_protein_compensation==1)
                    index_to_remove2 = find(cellfun(@str2num,list_proteins_significant_non_unique_downreg{j,k,i}(:,14))==1);
                else index_to_remove2 = []; end
                index_to_remove = union(index_to_remove1,index_to_remove2);
                index_to_remove3 = find(cellfun(@str2num,list_proteins_significant_non_unique_downreg{j,k,i}(:,5))>PSP_p_val);
                index_to_remove = union(index_to_remove,index_to_remove3);
                position_protein_external(index_to_remove) = [];
                p_names_comparison_external_list(index_to_remove) = [];
                [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
                [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
                score_kin_subExternal = [score_kin_subExternal,statistics];
                header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
                writecell(header,strcat(dir_current,'/list_kinases_significant_non_unique_downreg_','p',num2str(i),'.xlsx'),'Range','A1','Sheet',strcat(sample_names{j},'vs',sample_names{k}),'AutofitWidth',1);
                writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_non_unique_downreg_','p',num2str(i),'.xlsx'),'Range','A2','Sheet',strcat(sample_names{j},'vs',sample_names{k}),'AutofitWidth',0);
            end
        end
    end
end

for i = 1:size(list_proteins_significant_non_unique_upreg,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_upreg_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique_upreg,1)
        for k = 1:size(list_proteins_significant_non_unique_upreg,2)
            if(isempty(list_proteins_significant_non_unique_upreg{j,k,i})==0)
                p_names_comparison_external_list = list_proteins_significant_non_unique_upreg{j,k,i}(:,1);
                position_protein_external = list_proteins_significant_non_unique_upreg{j,k,i}(:,4);
                if(PSP_Loc_prob==1)
                    index_to_remove1 = find(cellfun(@str2num,list_proteins_significant_non_unique_upreg{j,k,i}(:,11))<0.75);
                else index_to_remove1 = []; end
                if(PSP_protein_compensation==1)
                    index_to_remove2 = find(cellfun(@str2num,list_proteins_significant_non_unique_upreg{j,k,i}(:,14))==1);
                else index_to_remove2 = []; end
                index_to_remove = union(index_to_remove1,index_to_remove2);
                index_to_remove3 = find(cellfun(@str2num,list_proteins_significant_non_unique_upreg{j,k,i}(:,5))>PSP_p_val);
                index_to_remove = union(index_to_remove,index_to_remove3);
                position_protein_external(index_to_remove) = [];
                p_names_comparison_external_list(index_to_remove) = [];
                [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
                [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
                score_kin_subExternal = [score_kin_subExternal,statistics];
                header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
                writecell(header,strcat(dir_current,'/list_kinases_significant_non_unique_upreg_','p',num2str(i),'.xlsx'),'Range','A1','Sheet',strcat(sample_names{j},'vs',sample_names{k}),'AutofitWidth',1);
                writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_non_unique_upreg_','p',num2str(i),'.xlsx'),'Range','A2','Sheet',strcat(sample_names{j},'vs',sample_names{k}),'AutofitWidth',0);
            end
        end
    end
end

for i = 1:size(list_proteins_significant_non_unique,3)
    delete(strcat(dir_current,'/list_proteins_significant_non_unique_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_non_unique,1)
        for k = 1:size(list_proteins_significant_non_unique,2)
            if(isempty(list_proteins_significant_non_unique{j,k,i})==0)
                p_names_comparison_external_list = list_proteins_significant_non_unique{j,k,i}(:,1);
                position_protein_external = list_proteins_significant_non_unique{j,k,i}(:,4);
                if(PSP_Loc_prob==1)
                    index_to_remove1 = find(cellfun(@str2num,list_proteins_significant_non_unique{j,k,i}(:,11))<0.75);
                else index_to_remove1 = []; end
                if(PSP_protein_compensation==1)
                    index_to_remove2 = find(cellfun(@str2num,list_proteins_significant_non_unique{j,k,i}(:,14))==1);
                else index_to_remove2 = []; end
                index_to_remove = union(index_to_remove1,index_to_remove2);
                index_to_remove3 = find(cellfun(@str2num,list_proteins_significant_non_unique{j,k,i}(:,5))>PSP_p_val);
                index_to_remove = union(index_to_remove,index_to_remove3);
                position_protein_external(index_to_remove) = [];
                p_names_comparison_external_list(index_to_remove) = [];
                [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
                [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
                score_kin_subExternal = [score_kin_subExternal,statistics];
                header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
                writecell(header,strcat(dir_current,'/list_kinases_significant_non_unique_','p',num2str(i),'.xlsx'),'Range','A1','Sheet',strcat(sample_names{j},'vs',sample_names{k}),'AutofitWidth',1);
                writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_non_unique_','p',num2str(i),'.xlsx'),'Range','A2','Sheet',strcat(sample_names{j},'vs',sample_names{k}),'AutofitWidth',0);
            end
        end
    end
end

for i = 1:size(list_proteins_significant_unique,2)
    delete(strcat(dir_current,'/list_proteins_significant_unique_','p',num2str(i),'.xlsx'));
    for j = 1:size(list_proteins_significant_unique,1)
        if(isempty(list_proteins_significant_unique{j,i})==0)
            p_names_comparison_external_list = list_proteins_significant_unique{j,i}(:,1);
            position_protein_external = list_proteins_significant_unique{j,i}(:,4);
            if(PSP_Loc_prob==1)
                index_to_remove1 = find(cellfun(@str2num,list_proteins_significant_unique{j,i}(:,11))<0.75);
            else index_to_remove1 = []; end
            if(PSP_protein_compensation==1)
                index_to_remove2 = find(cellfun(@str2num,list_proteins_significant_unique{j,i}(:,14))==1);
            else index_to_remove2 = []; end
            index_to_remove = union(index_to_remove1,index_to_remove2);
            index_to_remove3 = find(cellfun(@str2num,list_proteins_significant_unique{j,i}(:,5))>PSP_p_val);
            index_to_remove = union(index_to_remove,index_to_remove3);
            position_protein_external(index_to_remove) = [];
            p_names_comparison_external_list(index_to_remove) = [];
            [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
            [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
            score_kin_subExternal = [score_kin_subExternal,statistics];
            header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
            writecell(header,strcat(dir_current,'/list_kinases_significant_unique_','p',num2str(i),'.xlsx'),'Range','A1','Sheet',sample_names{j},'AutofitWidth',1);
            writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_unique_','p',num2str(i),'.xlsx'),'Range','A2','Sheet',sample_names{j},'AutofitWidth',0);
        end
    end
end


%Gropu comparison excel data!!!!!
delete(strcat(dir_current,'/list_proteins_significant_group_comparison.xlsx'));
for i = 1:size(list_proteins_significant_non_unique_upreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete,2)
        if(isempty(list_proteins_significant_non_unique_upreg_complete{i,j})==0)
            amx = [list_proteins_significant_unique_complete{i};list_proteins_significant_non_unique_upreg_complete{i,j}];
            p_names_comparison_external_list = amx(:,1);
            position_protein_external = amx(:,4);
            if(PSP_Loc_prob==1)
                index_to_remove1 = find(cellfun(@str2num,amx(:,11))<0.75);
            else index_to_remove1 = []; end
            if(PSP_protein_compensation==1)
                index_to_remove2 = find(cellfun(@str2num,amx(:,14))==1);
            else index_to_remove2 = []; end
            index_to_remove = union(index_to_remove1,index_to_remove2);
            index_to_remove3 = find(vertcat(amx{:,5})>PSP_p_val);
            index_to_remove = union(index_to_remove,index_to_remove3);
            position_protein_external(index_to_remove) = [];
            p_names_comparison_external_list(index_to_remove) = [];         
            [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
            [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
            score_kin_subExternal = [score_kin_subExternal,statistics];
            header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
            writecell(header,strcat(dir_current,'/list_kinases_significant_group_comparison.xlsx'),'Range','A1','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
            writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_group_comparison.xlsx'),'Range','A2','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end

%The index is inveresed because downregulated means inverse upregulated
for i = 1:size(list_proteins_significant_non_unique_downreg_complete,2)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete,1)
        if(isempty(list_proteins_significant_non_unique_downreg_complete{j,i})==0)
            amx = [list_proteins_significant_unique_complete{i};list_proteins_significant_non_unique_downreg_complete{j,i}];
            p_names_comparison_external_list = amx(:,1);
            position_protein_external = amx(:,4);
            if(PSP_Loc_prob==1)
                index_to_remove1 = find(cellfun(@str2num,amx(:,11))<0.75);
            else index_to_remove1 = []; end
            if(PSP_protein_compensation==1)
                index_to_remove2 = find(cellfun(@str2num,amx(:,14))==1);
            else index_to_remove2 = []; end
            index_to_remove = union(index_to_remove1,index_to_remove2);
            index_to_remove3 = find(vertcat(amx{:,5})>PSP_p_val);
            index_to_remove = union(index_to_remove,index_to_remove3);
            position_protein_external(index_to_remove) = [];
            p_names_comparison_external_list(index_to_remove) = [];
            [list_kinase_proteinsExternal,score_kin_subExternal] = getPhosphoSiteData(p_names_comparison_external_list,position_protein_external,PSP_file);
            [statistics] = pValue_PhosphoSitePlus_data(score_kin_subBackground,score_kin_subExternal);
            score_kin_subExternal = [score_kin_subExternal,statistics];
            header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Chi_test_pValue_Storey_corrected','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_Storey_corrected'};
            writecell(header,strcat(dir_current,'/list_kinases_significant_group_comparison.xlsx'),'Range','A1','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
            writecell(score_kin_subExternal,strcat(dir_current,'/list_kinases_significant_group_comparison.xlsx'),'Range','A2','Sheet',strcat(sample_names{i},'vs',sample_names{j}));
        end
    end
end