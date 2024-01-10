function [] = RNAseq_DGE_files(working_dir,folder_prefix,phenotype_names,phenotype_comparison,DGE_fileName,DRIMSeq_StageR_fileName)

dir1 = dir(working_dir);



isub = find(~[dir1(:).isdir]); %# returns logical vector
dir1(isub) = [];

isub = find(contains({dir1.name},folder_prefix)); %# returns logical vector
dir1 = dir1(isub);

%Set Phenotypes and number of experiments
phenotype_comparison_aux = cellfun(@(x) strsplit(x,'vs'),phenotype_comparison,'UniformOutput',false);
exp_name = {};

%Number of experiments
nr_experiments = length(dir1);

%The background and results matrix
background = cell(1,nr_experiments);
vec_common_DESeq2 = cell(1,nr_experiments);
vec_common_StageR = cell(1,nr_experiments);


for id = 1:size(dir1)
    
    %Keep the name of the experiment
    dir_project = dir(strcat(dir1(id).folder,'\',dir1(id).name));
    exp_name{id} = dir1(id).name;
    
    %Read DESeq2 DGE data - read each xlsx sheet and put it in a cell
    DESeq2_data_name = strcat(dir_project(1).folder,'\',DGE_fileName);
    xls_sheets = sheetnames(DESeq2_data_name);
    DESeq2_data = {};
    for j = 1:length(xls_sheets)
        %Read each sheet and replace empty entries with NaN
        opts = detectImportOptions(DESeq2_data_name,"Sheet",xls_sheets(j));
        opts = setvartype(opts,1:length(opts.VariableNames), 'char');
        opts.PreserveVariableNames=true;
        opts.Sheet = xls_sheets(j);
        aux = readtable(DESeq2_data_name, opts);
        a1 = aux{:,:};
        a1(cellfun('isempty',a1)) = {'NaN'};
        aux{:,:} = a1;
        
        %Extract the phenotype comaprison group that is mentioned in the
        %sheet name
        a1 = strsplit(xls_sheets(j),'_');
        a1 = a1(2);
        a1 = strsplit(a1,'vs');
        ind_phen = find(cell2mat(cellfun(@(x) isequal(x,a1), phenotype_comparison_aux, 'UniformOutput', false)));
        %Extract only the columns of interest
        ind_pval = find(ismember(aux.Properties.VariableNames,"padj"));
        DESeq2_data{ind_phen} = aux(:,[1:5,ind_pval,6]);
        %Construct the background of the measured genes in all experiments
        background{ind_phen} = vertcat(background{ind_phen},aux(:,[1,3,4]));
        %pval - extract only adjusted p-values below 0.05
        p_val = cellfun(@(x)str2num(x) ,aux{:,ind_pval},'UniformOutput',false);
        ind_c = find(cellfun(@isempty,p_val));
        p_val(ind_c) = num2cell(1);
        ind_c = find(vertcat(p_val{:})<=0.05);
        DESeq2_data{ind_phen} = DESeq2_data{ind_phen}(ind_c,:);
        %fold change - extract only the fold changes above 2
        fc = cellfun(@(x)str2num(x) ,DESeq2_data{ind_phen}{:,end},'UniformOutput',false);
        ind_c = find(cellfun(@isempty,fc));
        fc(ind_c) = num2cell(0);
        ind_c = find(vertcat(fc{:})>2);
        DESeq2_data{ind_phen} = DESeq2_data{ind_phen}(ind_c,:);
    end
    
    %Read StageR DTU data - read each xlsx sheet and put it in a cell
    StageR_data_name = strcat(dir_project(1).folder,'\',DRIMSeq_StageR_fileName);
    xls_sheets = sheetnames(StageR_data_name);
    StageR_data = {};
    for j = 1:length(xls_sheets)
        %Read the data and replace each empty entry with NaN
        opts = detectImportOptions(StageR_data_name,"Sheet",xls_sheets(j));
        opts = setvartype(opts,1:length(opts.VariableNames), 'char');
        opts.PreserveVariableNames=true;
        opts.Sheet = xls_sheets(j);
        aux = readtable(StageR_data_name, opts);
        a1 = aux{:,:};
        a1(cellfun('isempty',a1)) = {'NaN'};
        aux{:,:} = a1;
        
        %Select the columns of interest
        a1 = strsplit(xls_sheets(j),'_');
        a1 = a1(2);
        a1 = strsplit(a1,'vs');
        ind_phen = find(cell2mat(cellfun(@(x) isequal(x,a1), phenotype_comparison_aux, 'UniformOutput', false)));
        ind_pval = find(ismember(aux.Properties.VariableNames,"padj_gene"));
        StageR_data{ind_phen} = aux(:,[1,3:6]);
    end
    
    dir_project = dir(strcat(dir1(id).folder,'\',dir1(id).name));
    %Write the fitlered DESeq2 data
    am = strsplit(DGE_fileName,'.');
    xls_name = strcat(dir_project(1).folder,'\',am{1},'_filtered.xlsx');
    delete(xls_name);
    for j = 1:size(DESeq2_data,2)
        if(~isempty(DESeq2_data{j}))
            aux = DESeq2_data{j};
            aux{:,1} = cellfun(@(x) strsplit(x,'.'),aux{:,1},'UniformOutput',false);
            s1 = vertcat(aux{:,1}{:});
            aux{:,1} = s1(:,1);
            writetable(aux,xls_name,'Range','A1','Sheet',phenotype_comparison{j},'WriteMode','overwritesheet');
        end
    end
    %Put the different experiments into the same matrix for each phenotypic
    %comparison
    for j = 1:size(DESeq2_data,2)
        if(~(isempty(DESeq2_data{j})))
            vec_common_DESeq2{j} = vertcat(vec_common_DESeq2{j},unique(DESeq2_data{j}{:,1}));
            vec_common_StageR{j} = vertcat(vec_common_StageR{j},unique(StageR_data{j}{:,1}));
        end
    end
    
    
end


%Write the background for each comaprison group
am = strsplit(DGE_fileName,'.');
am = am{1};
am = strsplit(am,'\\');
am = am{length(am)};
xls_name = strcat(dir1(1).folder,'\',am,'_background.xlsx');
delete(xls_name);
for j = 1:size(background,2)
    %aux_names = background{j}.Properties.VariableNames;
    %background{j} = table2cell(background{j});
    ajt = cellfun(@(x) strsplit(x,'.'),background{j}{:,1},'UniformOutput',false);
    ajt = vertcat(ajt{:});
    background{j}{:,1} = ajt(:,1);
    background{j} = unique(background{j},'rows');
    writetable(background{j},xls_name,'Range','A1','Sheet',phenotype_comparison{j},'WriteMode','overwritesheet');
end
%Extract the genes that are upregulated in at least half of the experiments
%(all the genes for each phenotypic comaprisons are in the same matrix so
%we count the number of appearances for each gene in a specific phenotypic
%comparison
un = cell(1,nr_experiments);
for j = 1:size(vec_common_DESeq2,2)
    
    %Count the number of hits for each unique gene
    aux = vec_common_DESeq2{j};
    un{j} = unique(aux);
    un_number_hits = cellfun(@(x) length(find(ismember(aux,x))),un{j});
    
    %Select the required "half of experiments" threshold
    if(j~=1 && j~=4)
        ind = find(un_number_hits>=1);
    else
        ind = find(un_number_hits>=3);
    end
    %Remove the "." from the names
    un{j} = un{j}(ind);
    un{j} = cellfun(@(x) strsplit(x,'.'),un{j},'UniformOutput',false);
    un{j} = vertcat(un{j}{:});
    un{j} = un{j}(:,1);
    is = find(ismember(background{j}{:,1},un{j}));
    un{j} = cell2table(un{j});
    un{j}(:,[2,3]) = background{j}(is,[2,3]);
    un{j}.Properties.VariableNames = background{j}.Properties.VariableNames;
    
end

%Write the DESeq2 data
am = strsplit(DGE_fileName,'.');
am = am{1};
am = strsplit(am,'\\');
am = am{length(am)};
xls_name = strcat(dir1(1).folder,'\',am,'_significat_OverlapExp.xlsx');
delete(xls_name);
for j = 1:size(un,2)
    writetable(un{j},xls_name,'Range','A1','Sheet',phenotype_comparison{j},'WriteMode','overwritesheet');
end

deseq2_upregulated = un;

%Extract the DTU that are upregulated in at least half of the experiments
%(all the genes for each phenotypic comaprisons are in the same matrix so
%we count the number of appearances for each gene in a specific phenotypic
%comparison
un = cell(1,nr_experiments);
for j = 1:size(vec_common_StageR,2)
    
    %Count the number of hits for each unique DTU
    aux = vec_common_StageR{j};
    un{j} = unique(aux);
    un_number_hits = cellfun(@(x) length(find(ismember(aux,x))),un{j});
    
    %Select the required "half of experiments" threshold
    if(j~=1 && j~=4)
        ind = find(un_number_hits>=1);
    else
        ind = find(un_number_hits>=3);
        
    end
    %Remove the "." from the names
    un{j} = un{j}(ind);
    un{j} = cellfun(@(x) strsplit(x,'.'),un{j},'UniformOutput',false);
    un{j} = vertcat(un{j}{:});
    un{j} = un{j}(:,1);
    is = find(ismember(background{j}{:,1},un{j}));
    un{j} = cell2table(un{j});
    un{j}(:,[2,3]) = background{j}(is,[2,3]);
    un{j}.Properties.VariableNames = background{j}.Properties.VariableNames;
end

%Write the DTU data
am = strsplit(DRIMSeq_StageR_fileName,'.');
am = am{1};
am = strsplit(am,'\\');
am = am{length(am)};
xls_name = strcat(dir1(1).folder,'\',am,'_OverlapExp.xlsx');
delete(xls_name);
for j = 1:size(un,2)
    writetable(un{j},xls_name,'Range','A1','Sheet',phenotype_comparison{j},'WriteMode','overwritesheet');
end



end