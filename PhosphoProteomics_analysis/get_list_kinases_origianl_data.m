function [kinase_original_list] = get_list_kinases_origianl_data(p_names)

%Check if human_kinome has unique rows!! Apply this function for the initial
%proteins split by phenotype. 

human_kinome = ".\Human_kinome.txt";
opts = detectImportOptions(human_kinome);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
human_kinome = readtable(human_kinome, opts);

pat_ini_protein = 'sp|';
pat_fin_protein = '|';



aa = {};
for i = 1:length(p_names)
    if(contains(p_names(i),"cont",'IgnoreCase',true)==0)
        if (isempty(find(contains(p_names,'sp'),1)) == 0)
            aa = [aa;extractBetween(p_names(i),pat_ini_protein,pat_fin_protein)];
        else
            aa = [aa;p_names(i)];
        end
    end
end
p_names = aa;

human_kinome = sortrows(human_kinome,1);
kinase_original_list = cellfun(@(x) human_kinome(find(ismember(human_kinome{:,1},x)),:),p_names,'UniformOutput',false);
kinase_original_list = kinase_original_list(~cellfun(@isempty,kinase_original_list));
kinase_original_list = cellfun(@(x) x(:,[1,2,6]),kinase_original_list,'UniformOutput',false);
kinase_original_list = vertcat(kinase_original_list{:});

end


