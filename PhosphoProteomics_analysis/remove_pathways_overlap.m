function [status] = remove_pathways_overlap(name_dir,source,source_databases,effective_size,members_input_overlap,pvalue,padj_values,removal_method,maximum_pathway_dim,min_number_members)

if  (isempty(name_dir) || ~isfolder(name_dir))
    error('name_dir was not correctly provided.');
end

if  (isempty(source) || class(source)~="char")
    error('source was not correctly provided.');
end

if  (isempty(source_databases))
    error('source_databases was not correctly provided.');
end

if  (isempty(effective_size) || class(effective_size)~="char")
    error('effective_size was not correctly provided.');
end

if  (isempty(members_input_overlap) || class(members_input_overlap)~="char")
    error('members_input_overlap was not correctly provided.');
end

if  (isempty(pvalue) || class(pvalue)~="char")
    error('pvalue was not correctly provided.');
end

if  (isempty(padj_values) || class(padj_values)~="char")
    error('padj_values was not correctly provided.');
end

if  (isempty(removal_method) || class(removal_method)~="double")
    error('removal_method was not correctly provided.');
end

if  (isempty(maximum_pathway_dim) || class(maximum_pathway_dim)~="double")
    error('maximum_pathway_dim was not correctly provided.');
end

if  (isempty(min_number_members) || class(min_number_members)~="double")
    error('min_number_members was not correctly provided.');
end

%Removal method - 1 - remove pathways with less than two different memebers; 2 - remove the pathwhays that overlap more than 80% having the reference the second pathway.

dir1 = dir(name_dir);

ind = [];
for i = 1:length(dir1)
    if(length(dir1(i).name)>4)
        if (dir1(i).name(end-4:end)==".xlsx") ind = [ind,i]; end
    end
end

dir1 = dir1(ind);

for i = 1:length(dir1)
    
    file_name = strcat(dir1(i).folder,'\',dir1(i).name);
    xls_sheets = sheetnames(file_name);
    
    data = {};
    
    for ii = 1:length(xls_sheets)
        
        opts = detectImportOptions(file_name,"Sheet",xls_sheets(ii));
        opts = setvartype(opts,1:length(opts.VariableNames), 'char');
        opts.PreserveVariableNames=true;
        data{ii} = readtable(file_name, opts);
    end
    
    for ii = 1:length(data)
        
        xs = data{ii}{:,find(ismember(data{ii}.Properties.VariableNames, members_input_overlap))};
        xs1 = cellfun(@(x) length(strsplit(x,';')),xs);
        ind1 = find(xs1>=min_number_members);
        data{ii} = data{ii}(ind1,:);
        
        ind1 = find(ismember(lower(data{ii}{:,find(ismember(lower(data{ii}.Properties.VariableNames),lower(source)))}),lower(source_databases)));
        data{ii} = data{ii}(ind1,:);
        
        ind = find(ismember(data{ii}.Properties.VariableNames,effective_size));
        if(~isempty(ind))
            ind1 = cellfun(@(x) str2num(x), (data{ii}{:,ind}))<maximum_pathway_dim;
            data{ii} = data{ii}(ind1,:);
        end
        
        ind = find(ismember(data{ii}.Properties.VariableNames, padj_values));
        ind1 = cellfun(@(x) str2num(x), (data{ii}{:,ind}))>0;
        data1 = data{ii}(ind1,:);
        
        ind1 = cellfun(@(x) str2num(x), (data{ii}{:,ind}))<0;
        data2 = data{ii}(ind1,:);
        
        ind = find(ismember(data{ii}.Properties.VariableNames,members_input_overlap));
        
        aux = table2cell(data1(:,ind));
        aux = cellfun(@(x) strsplit(x,';'),aux,'UniformOutput',false);
        
        ind = [];
        for j = 1:length(aux)-1
            for k = j+1:length(aux)
                ind1 = length(find(ismember(aux{k},aux{j})));
                %if (ind1>length(aux{k})-2) ind = [ind ;[j,k]]; end
                %If we have a pathway that contains 3 terms and another
                %pathway with a smaller q_value that contains the same
                %terms plus some more in addition - we would like to keep
                %both in order to have a good picture overview of what is
                %happening inside the cell. If the first pathway (with the
                %lowest q_value) has 10 elements and there are pathways
                %below it that have a 5 elements all common to the previous
                %pathway we want to discard it. This is equivalent of
                %always using as the length referece value the second term
                %aux(k)
                if(length(aux{k})>length(aux{j}))
                    ref = max(length(aux{k}),length(aux{j}));
                else
                    ref = min(length(aux{k}),length(aux{j}));
                end
                if (removal_method==2)
                    if (ind1>=floor(0.8*ref))
                        ind = [ind ;[j,k]];
                    end
                end
                if (removal_method==1)
                    if (ind1>=(min(length(aux{k}),length(aux{j}))-2))
                        ind = [ind ;[j,k]];
                    end
                end
            end
        end
        
        if(~isempty(ind)) data1(ind(:,2),:) = []; end
        
        ind = find(ismember(data{ii}.Properties.VariableNames,members_input_overlap));
        aux = table2cell(data2(:,ind));
        aux = cellfun(@(x) strsplit(x,';'),aux,'UniformOutput',false);
        
        ind = [];
        for j = 1:length(aux)-1
            for k = j+1:length(aux)
                ind1 = length(find(ismember(aux{k},aux{j})));
                %if (ind1>length(aux{k})-2) ind = [ind ;[j,k]]; end
                if(length(aux{k})>length(aux{j}))
                    ref = max(length(aux{k}),length(aux{j}));
                else
                    ref = min(length(aux{k}),length(aux{j}));
                end
                if (removal_method==2)
                    if (ind1>=floor(0.8*ref))
                        ind = [ind ;[j,k]];
                    end
                end
                if (removal_method==1)
                    if (ind1>=(min(length(aux{k}),length(aux{j}))-2))
                        ind = [ind ;[j,k]];
                    end
                end
            end
        end
        
        if(~isempty(ind)) data2(ind(:,2),:) = []; end
        
        %         if(isempty(data1)==0)
        %         aux = num2cell(mafdr(cellfun(@(x) str2num(x), data1{:,find(data{ii}.Properties.VariableNames == pvalue)}),'BHFDR','TRUE'));
        %         aux = cellfun(@num2str,aux,'UniformOutput',false);
        %         data1{:,find(data{ii}.Properties.VariableNames == padj_values)} = aux;
        %         end
        
        %         if(isempty(data2)==0)
        %         aux = num2cell(mafdr(cellfun(@(x) str2num(x), data2{:,find(data{ii}.Properties.VariableNames == pvalue)}),'BHFDR','TRUE'));
        %         aux = cellfun(@num2str,aux,'UniformOutput',false);
        %         data2{:,find(data{ii}.Properties.VariableNames == padj_values)} = aux;
        %         end
        data{ii} = [data1;data2];

    end

    mkdir(strcat(dir1(i).folder,'\Filtered_UsedOne'))
    name_new = strcat(dir1(i).folder,'\Filtered_UsedOne\',dir1(i).name);
    delete(name_new);
    
    for ii = 1:length(data)
        writetable(data{ii},name_new,'Sheet',xls_sheets{ii},'AutofitWidth',0);
    end
    
end
status=0;
end