function [interaction_table] = kinase_cascade_generation(analyzed_file,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,Xlsx_Sheet_file_kinase_total,Xlsx_Sheet_analyzed_file_data,file_TF_measured,file_upregulated_members,Xlsx_Sheet_file_upregulated_members,file_kinases_list,p_names,p_full_phenotype)


pat_ini_protein = 'sp|';
pat_fin_protein = '|';
pat_ini_gene = 'sp|'+ wildcardPattern + '|';
pat_fin_gene = '_HUMAN' ;


opts = detectImportOptions(file_kinases_Phospho);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_kinases_Phospho = readtable(file_kinases_Phospho, opts);
% ind_sub_acc = find(ismember(file_kinases_Phospho.Properties.VariableNames,'SUB_ACC_ID'));
% ind_sub_rsd = find(ismember(file_kinases_Phospho.Properties.VariableNames,'SUB_MOD_RSD'));
% ind_kin_acc = find(ismember(file_kinases_Phospho.Properties.VariableNames,'KIN_ACC_ID'));
% ind_kin_name = find(ismember(file_kinases_Phospho.Properties.VariableNames,'KINASE'));
amsx = {};
if(size(file_kinases_Phospho,2)<10)
    amsx(:,3) = file_kinases_Phospho{:,1};
    amsx(:,7) = file_kinases_Phospho{:,2};
    amsx(:,10) = file_kinases_Phospho{:,3};
end
amsx = cell2table(amsx);
amsx.Properties.VariableNames([3,7,10]) = {'KIN_ACC_ID','SUB_ACC_ID','SUB_MOD_RSD'};
file_kinases_Phospho = amsx;

opts = detectImportOptions(file_TF);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_TF = readtable(file_TF, opts);
% ind_is_TF = find(ismember(lower(file_TF.Properties.VariableNames),lower('Is TF?')));
% ind = find(ismember(lower(file_TF{:,ind_is_TF}),lower('No')));
% file_TF(ind,:) = [];


opts = detectImportOptions(file_TF_mapping);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_TF_mapping = readtable(file_TF_mapping, opts);

opts = detectImportOptions(file_TF_measured);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_TF_measured = readtable(file_TF_measured, opts);

opts = detectImportOptions(file_kinases_NetPhorest);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_kinases_NetPhorest = readtable(file_kinases_NetPhorest, opts);

opts = detectImportOptions(analyzed_file);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
analyzed_file = readtable(analyzed_file, opts,'Sheet',Xlsx_Sheet_analyzed_file_data);
ind1 = find(ismember(analyzed_file{:,6},'true'));
ind2 = find(ismember(analyzed_file{:,7},'true'));
ind = union(ind1,ind2);
analyzed_file = analyzed_file(ind,:);

opts = detectImportOptions(file_kinase_total);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_kinase_total = readtable(file_kinase_total, opts,'Sheet',Xlsx_Sheet_file_kinase_total);

opts = detectImportOptions(file_upregulated_members);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_upregulated_members = readtable(file_upregulated_members, opts,'Sheet',Xlsx_Sheet_file_upregulated_members);


opts = detectImportOptions(file_kinases_list);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
file_kinases_list = readtable(file_kinases_list, opts);


interaction_table = table('Size',[1 4],'VariableTypes',["cell","cell","double","cell"],'VariableNames',["Kinase","Substrate","Database","Residues"]);

seed_list_substrate = table('Size',[1 2],'VariableTypes',["cell","cell"],'VariableNames',["Var1","Var2"]);
seed_list_kin = table('Size',[1 1],'VariableTypes',"cell",'VariableNames',"Var1");


%1 - TF, 2 - PSP, 3 - NetPhorest, 10 - PSP identified but not found in our
%dataset
for i = 1:size(analyzed_file,1)
    tic
    kin_substrate = analyzed_file{i,1};
    kin_kinase = analyzed_file{i,4};
    residue_kin_substrate = analyzed_file{i,3};
    %The table is first kinase and then substrate
    interaction_table = [interaction_table;[kin_kinase,kin_substrate,2,residue_kin_substrate]];
    seed_list_substrate = [seed_list_substrate;[kin_substrate,residue_kin_substrate]];
    seed_list_kin = [seed_list_kin;kin_kinase];
    if (i==1) interaction_table(1,:) = []; seed_list_substrate(1,:) = []; seed_list_kin(1,:) = []; end
    [list1,seed_list_kin,seed_list_substrate] = kinase_cascade_generationSubstrate(kin_substrate,seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,0,kin_kinase,residue_kin_substrate);
    [list2,seed_list_substrate,seed_list_kin] = kinase_cascade_generationKinase(kin_substrate,seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,0,kin_kinase,residue_kin_substrate);
    [list3,seed_list_kin,seed_list_substrate] = kinase_cascade_generationSubstrate(kin_kinase,seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,0,kin_substrate,residue_kin_substrate);
    [list4,seed_list_substrate,seed_list_kin] = kinase_cascade_generationKinase(kin_kinase,seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,0,kin_substrate,residue_kin_substrate);
    interaction_table = [interaction_table;[list1;list2;list3;list4]];
    toc
    
end


interaction_table = unique(interaction_table,'stable');

%Upregulated only by protein without taking into consideration the residues
ind_kinase_upregulated = find(ismember(lower(interaction_table{:,1}),lower(file_upregulated_members{:,1})));
aux = zeros(size(interaction_table,1),1);
aux(ind_kinase_upregulated) = 1;
aux = table(aux);
aux.Properties.VariableNames = {'Kinase_Upregulated'};
interaction_table = [interaction_table,aux];

ind_substrate_upregulated = find(ismember(lower(table(interaction_table{:,2},interaction_table{:,4})),lower(table(file_upregulated_members{:,1},file_upregulated_members{:,4}))));
aux = zeros(size(interaction_table,1),1);
aux(ind_substrate_upregulated) = 1;
aux = table(aux);
aux.Properties.VariableNames = {'Substrate_Upregulated'};
interaction_table = [interaction_table,aux];

ind_substrate_upregulated_n = find(ismember(lower(interaction_table{:,2}),lower(file_upregulated_members{:,1})));
aux = zeros(size(interaction_table,1),1);
aux(ind_substrate_upregulated_n) = 1;
aux = table(aux);
aux.Properties.VariableNames = {'Substrate_Upregulated_Not_considering_residues'};
interaction_table = [interaction_table,aux];

aux = cellfun(@(x) cell2mat(file_kinases_list{find(ismember(file_kinases_list{:,1},x),1),6}),interaction_table{:,1},'UniformOutput',false);
aux = table(aux);
aux.Properties.VariableNames = {'Kinase_Family_Name'};
interaction_table = [interaction_table,aux];

aux = cellfun(@(x) cell2mat(file_kinases_list{find(ismember(file_kinases_list{:,1},x),1),6}),interaction_table{:,2},'UniformOutput',false);
aux = table(aux);
aux.Properties.VariableNames = {'Substrate_Family_Name'};
interaction_table = [interaction_table,aux];


% gene_protein_list = {};
% for i = 1:length(p_names)
%     if(contains(p_names(i),"cont",'IgnoreCase',true)==0)
%         aa1 = extractBetween(p_names(i),pat_ini_protein,pat_fin_protein);
%         aa2 = extractBetween(p_names(i),pat_ini_gene,pat_fin_gene);
%         gene_protein_list = [gene_protein_list;[aa1,aa2]];
%     end
% end
% gene_protein_list = unique(cell2table(gene_protein_list),'stable');

gene_protein_list = cell2table(p_full_phenotype(:,1));

ind_kinase_measured = cell2table(num2cell(double(ismember(lower(interaction_table{:,1}),lower(gene_protein_list{:,1})))));
ind_kinase_measured.Properties.VariableNames = {'Kinase_Measured'};
interaction_table = [interaction_table,ind_kinase_measured];

ind_sub_measured = cell2table(num2cell(double(ismember(lower(interaction_table{:,2}),lower(gene_protein_list{:,1})))));
ind_sub_measured.Properties.VariableNames = {'Substrate_Measured'};
interaction_table = [interaction_table,ind_sub_measured];

%Cut exactly the source of the dangling node
% aa1 = interaction_table(find([interaction_table{:,3}]==101),:);
% for i = 1:size(aa1,1)
%     aa2 = aa1(i,1:2);
%     aa2 = [aa2,aa1(i,4)];
%     ind = find(ismember(interaction_table(:,[1 2 4]),aa2));
%     aa22 = aa1{i,[2 1]};
%     aa22 = [aa22,aa1(i,4)];
%     ind1 = find(ismember(cell2table(interaction_table{:,[1 2 4]}),cell2table(aa22{:,:})));
%     ind2 = union(ind,ind1);
%     iind = find(ismember(interaction_table{ind2,3},210));
%     iind1 = find(ismember(interaction_table{ind2,3},510));
%     iind = union(iind,iind1);
%     if(isempty(iind))
%         interaction_table(ind2,:)=[];
%     end
%
% end

%Cut anything related to the dangling node
aa1 = interaction_table(find([interaction_table{:,3}]==101),:);
for i = 1:size(aa1,1)
    aa2 = aa1{i,1:2};
    %Find the elements that are not measured
    ind1 = find(ismember(aa2{1},file_kinase_total{:,1}));
    ind11 = find(ismember(aa2{1},file_TF_measured{:,2}));
    ind1 = union(ind1,ind11);
    ind2 = find(ismember(aa2{2},file_kinase_total{:,1}));
    ind22 = find(ismember(aa2{2},file_TF_measured{:,2}));
    ind2 = union(ind2,ind22);
    
    if(isempty(ind1)) aa2 = aa2{1};
    elseif(isempty(ind2)) aa2 = aa2{2};
    end
    
    ind1 = find(ismember(interaction_table{:,1},aa2));
    ind2 = find(ismember(interaction_table{:,2},aa2));
    ind = union(ind1,ind2);
    
    iind = find(ismember(interaction_table{ind,3},210));
    iind1 = find(ismember(interaction_table{ind,3},510));
    iind = union(iind,iind1);
    
    if(isempty(iind))
        interaction_table(ind,:)=[];
    end
    
end

ind = find([interaction_table{:,3}]==101);
interaction_table(ind,:) = [];


%Remove a node that has a score of 10 and is connected to only one other
%node (it is a solo identified transcription factor)
aa1 = interaction_table(find([interaction_table{:,3}]==10),:);
for i = 1:size(aa1,1)
    aa2 = aa1{i,1:2};
    %Find the elements that are not measured
    ind1 = find(ismember(aa2{1},file_kinase_total{:,1}));
    ind11 = find(ismember(aa2{1},file_TF_measured{:,2}));
    ind1 = union(ind1,ind11);
    ind2 = find(ismember(aa2{2},file_kinase_total{:,1}));
    ind22 = find(ismember(aa2{2},file_TF_measured{:,2}));
    ind2 = union(ind2,ind22);
    
    if(isempty(ind1)) aa2 = aa2{1};
    elseif(isempty(ind2)) aa2 = aa2{2};
    end
    
    ind1 = find(ismember(interaction_table{:,1},aa2));
    ind2 = find(ismember(interaction_table{:,2},aa2));
    ind = union(ind1,ind2);
    
    table_TF_dangling = unique(interaction_table(ind,1:2));
    if(size(table_TF_dangling,1)<=1)
        interaction_table(ind,:) = [];
    end
    
end