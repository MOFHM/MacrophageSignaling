function [] = NetPhorest_Kinase_extraction2_byPhosphoPeptides(naming_convention_comparison,LocProb_threshold,FDR_threshold,compensate_for_protein_abundances,dir_original_fasta_files,date_dir,dir_original_fasta_files_linux,NetPhorest_script,root_NetPhorest,file_PhosphoProteomics,Interaction_database,NetPhorest_full_names,workspace_file)

load(workspace_file)

list_proteins_total_per_phenotype = cell(1,nr_sample_groups);
for k = 1:nr_sample_groups
    for i = (k-1)*nr_rep+1:nr_rep*nr_sample_groups:size(int_data,2)
        for j = 1:size(int_data,1)
            vec = int_data(j,i:i+nr_rep-1);
            [vec_sort,initial_ind] = sort(vec,'descend');
            ind = find(vec_sort==0);
            if(length(ind)<nr_rep) list_proteins_total_per_phenotype{1,k} = [list_proteins_total_per_phenotype{1,k};[p_names(j),position_protein(j),peptide_sequence(j)]]; end
        end
    end
    aa1 = unique(cell2table(list_proteins_total_per_phenotype{1,k}));
    aa1 = table2cell(aa1);
    list_proteins_total_per_phenotype{1,k} = aa1;
end



file_upregulated_members = strcat(dir1,'\Excel_data_unique_peptide_values\list_proteins_significant_group_comparison.xlsx');
Xlsx_Sheet_file_upregulated_members = sheetnames(file_upregulated_members);
up_reg_P = cell(1,length(Xlsx_Sheet_file_upregulated_members));

%Filter the measured residues!
for i = 1:length(Xlsx_Sheet_file_upregulated_members)
    opts = detectImportOptions(file_upregulated_members);
    opts = setvartype(opts,1:length(opts.VariableNames), 'char');
    opts.PreserveVariableNames=true;
    aux = readtable(file_upregulated_members, opts,'Sheet',Xlsx_Sheet_file_upregulated_members(i));
    ind_s = intersect(find(cellfun(@str2num,(aux{:,11}))>=LocProb_threshold),find(cellfun(@str2num,(aux{:,5}))<=FDR_threshold));
    if(compensate_for_protein_abundances==1)
        ind_s = intersect(ind_s,find(cellfun(@str2num,aux{:,16})==0));
    end
    up_reg_P{i} = aux{ind_s,[1,4,3]};
    
end
Xlsx_Sheet_file_upregulated_members = strrep(Xlsx_Sheet_file_upregulated_members,'_','');



dir11 = dir(dir_original_fasta_files);

dir_linux_save = strcat('/home/totu/NetPhorest1/',date_dir);

st11 = system(strcat('wsl mkdir /home/totu/NetPhorest1/'));
st11 = system(strcat('wsl mkdir /home/totu/NetPhorest1/',date_dir));
st11 = system(strcat("wsl mkdir ",dir_linux_save,"/NetPhorest_results"));

dir_linux_results = strcat(dir_linux_save,"/NetPhorest_results");

st1 = system(strcat("wsl cp ", dir_original_fasta_files_linux, "/Original_fasta_files ",dir_linux_save, " ","-r"));
st2 = system('wsl chown -R $USER ~');
st3 = system(strcat("wsl cc -O3 -o netphorest"," ", NetPhorest_script,"/netphorest.c", " -lm"));
for i = 3:size(dir11,1)
    a = strsplit(dir11(i).name,'.');
    st4 = system(strcat("wsl bash -c  ",'"',"cat ",dir_linux_save,"/Original_fasta_files/",dir11(i).name," | ",NetPhorest_script,"/netphorest > ",dir_linux_results,"/Results_file_",a{1},".txt",'"'));
end
st5 = system(strcat("wsl \cp ",dir_linux_results," ",dir_original_fasta_files_linux,"/NetPhorest_results/  -r"));

Precompiled_files_generation_for_PhosphoProteins(root_NetPhorest,file_PhosphoProteomics,Interaction_database,NetPhorest_full_names,loc_prob,p_names,position_protein)

dir_results = dir(strcat(dir1,'\NetPhorest\NetPhorest_results'));

for i = 3:size(dir_results,1)
    
    if(isequal(dir_results(i).name(end-3:end),'.txt'))
        
        if(~isempty(find(contains(dir_results(i).name,"Phospho"))))
            continue
        end
        
        NetPhorest_full_data = readtable(NetPhorest_full_names,'VariableNamingRule','preserve');
        NetPhorest_full_data = sortrows(NetPhorest_full_data);
        
        
        for index_interaction = 1:2
            
            
            %Read the analyzed kinase data
            file = readtable(strcat(dir_results(i).folder,'\',dir_results(i).name));
            file.Var1 = extractBetween(file.Var1,'|','|');
            Protein_names = file.Var1;
            
            %Delete rows that do not fit the threshold
            posterior = file.Var9;
            prior = file.Var10;
            ind = find(posterior<prior);
            file(ind,:) = [];
            posterior = file.Var9;
            ind = find(posterior<0.035);
            file(ind,:) = [];
            ind = find(ismember(file.Var8,'any_group'));
            file(ind,:) = [];
            
            if(index_interaction == 1)
                %Load the precompiled list of kinases
                
                clear kinase_list_proteomics
                %                 fileID = fopen('file_proteomics_splited_by_proteins_byPhosphoPeptideAnalysis_fitlered_by_interaction.bin');
                %                 %fileID = fopen('file_PhosphoProteomics_splited_by_proteins_byPhosphoPeptideAnalysis_fitlered_by_interaction.bin');
                %                 kk = fread(fileID);
                %                 kinase_list_proteomics = hlp_deserialize(kk);
                %                 clear kk
                %                 fclose(fileID);
                %                 for j = 1:length(kinase_list_proteomics)
                %                     kinase_list_proteomics{j} = cellfun(@cell2table,kinase_list_proteomics{j},'UniformOutput',false);
                %                      for ijk = 1:length(kinase_list_proteomics{j})
                %                         kinase_list_proteomics{j}{ijk}{:,9:10} = cell2mat(kinase_list_proteomics{j}{ijk}{:,9:10});
                %                      end
                %                 end
                load(strcat(root_NetPhorest,'file_PhosphoProteomics_splited_by_proteins_byPhosphoPeptideAnalysis_fitlered_by_interaction'))
                
                
                %Choose the desire database file - with multiple or fewer
                %values
                file_interactions = Interaction_database;
                opts = detectImportOptions(file_interactions);
                opts = setvartype(opts,1:length(opts.VariableNames), 'char');
                opts.PreserveVariableNames=true;
                file_interactions = readtable(file_interactions, opts);
                ind_interactorA = find(ismember(file_interactions.Properties.VariableNames,'protein1'));
                ind_interactorB = find(ismember(file_interactions.Properties.VariableNames,'protein2'));
                
                n1 = table(file.Var1,file.Var8);
                n1 = unique(n1,'stable');
                n_aux = cell(size(n1,1),1);
                for j = 1:size(n1,1)
                    n_aux{j} = NetPhorest_full_data.Unip_ID(find(ismember(NetPhorest_full_data.NetPhorest,n1{j,2})));
                end
                n1{:,3} = n_aux;
                n_aux = cell(size(n1,1),1);
                n_aux1 = cell(size(n1,1),1);
                for j = 1:size(n1,1)
                    n_aux{j} = repmat(n1{j,1},length(n1{j,3}{:}),1);
                    n_aux1{j} = repmat(j,length(n1{j,3}{:}),1);
                end
                n1{:,4} = n_aux;
                n1{:,5} = n_aux1;
                ns = cell2table([vertcat(n1{:,4}{:}),vertcat(n1{:,3}{:}),num2cell(vertcat(n1{:,5}{:}))]);
                
                file_interactions_table1 = cell2table([file_interactions{:,ind_interactorA},file_interactions{:,ind_interactorB}]);
                file_interactions_table2 = cell2table([file_interactions{:,ind_interactorB},file_interactions{:,ind_interactorA}]);
                
                ind1 = ismember(ns(:,1:2),file_interactions_table1);
                ind2 = ismember(ns(:,1:2),file_interactions_table2);
                ind = or(ind1,ind2);
                interaction_ind=find(ind);
                ns = ns(interaction_ind,:);
                ind_aux = unique(ns{:,3});
                ind = ismember(table(file.Var1,file.Var8),n1(ind_aux,1:2));
                file = file(ind,:);
            else
                clear kinase_list_proteomics
                %                 fileID = fopen('file_proteomics_splited_by_proteins_byPhosphoPeptideAnalysis.bin');
                %                 %fileID = fopen('file_PhosphoProteomics_splited_by_proteins_byPhosphoPeptideAnalysis.bin');
                %                 kk = fread(fileID);
                %                 kinase_list_proteomics = hlp_deserialize(kk);
                %                 clear kk
                %                 fclose(fileID);
                %                 for j = 1:length(kinase_list_proteomics)
                %                     kinase_list_proteomics{j} = cellfun(@cell2table,kinase_list_proteomics{j},'UniformOutput',false);
                %                     for ijk = 1:length(kinase_list_proteomics{j})
                %                         kinase_list_proteomics{j}{ijk}{:,9:10} = cell2mat(kinase_list_proteomics{j}{ijk}{:,9:10});
                %                     end
                %                 end
                load(strcat(root_NetPhorest,'file_PhosphoProteomics_splited_by_proteins_byPhosphoPeptideAnalysis'))
            end
            
            
            kinase_target = file.Var1;
            kinase_location_number = file.Var2;
            kinase_location_residue = file.Var3;
            kinase_peptide = file.Var4;
            
            n1 = unique(kinase_target,'stable');
            
            
            %Remove entries that are not found in the experimental data
            kinase_list = {};
            %Choose the protein_set_total
            %              dir_name_evaluated = dir_results(i).name;
            %             if(length(find(contains(dir_name_evaluated,"_")))>0)
            %                 dir_name_evaluated = erase(dir_name_evaluated,"_");
            %             end
            %                                 index = cellfun(@(x) strfind(dir_name_evaluated,x),naming_convention_comparison,'UniformOutput',false);
            %                                 ind = cellfun(@isempty,index);
            %                                 index{ind} = 10000;
            %                                 index = cell2mat(index);
            %                                 index = find(index==min(index));
            %                                 proteins_set_total = protein_set_total_function_generation(list_proteins_total_per_phenotype{index}(:,1),list_proteins_total_per_phenotype{index}(:,2),list_proteins_total_per_phenotype{index}(:,3));
            dir_name_evaluated = dir_results(i).name;
            if(length(find(contains(dir_name_evaluated,"_")))>0)
                dir_name_evaluated = erase(dir_name_evaluated,"_");
            end
            index =  find(~cellfun(@isempty,(cellfun(@(x) find(contains(lower(dir_name_evaluated),lower(x))),Xlsx_Sheet_file_upregulated_members,'UniformOutput',false))));
            proteins_set_total = {up_reg_P{index}(:,1),cellfun(@(x) x(1),up_reg_P{index}(:,2)),cellfun(@(x) x(2:end),up_reg_P{index}(:,2),'UniformOutput',false),up_reg_P{index}(:,3)};
            
            
            for j = 1:length(n1)
                ind1 = find(contains(proteins_set_total{1},n1(j)));
                %kinase_list{j,1}{1} = n1(j);
                ind = [];
                ind_total = [];
                i_first = find(ismember(kinase_target,n1(j)),1,'first');
                i_last = find(ismember(kinase_target,n1(j)),1,'last');
                %If the search is done on the peptide sequence as well the nubmer
                %of hits increases but the actual residue of interest is no longer
                %alligned with the results
                %  ind_peptide1 = cellfun(@(x) find(contains(upper(proteins_set_total{4}(ind1)),upper(x))),kinase_peptide(i_first:i_last),'UniformOutput',false);
                %  ind_peptide2 = cellfun(@(x) find(contains(upper(x),upper(proteins_set_total{4}(ind1)))),kinase_peptide(i_first:i_last),'UniformOutput',false);
                ind_residue = cellfun(@ (x) find(ismember(proteins_set_total{2}(ind1),x)), kinase_location_residue(i_first:i_last), 'UniformOutput', false);
                ind_number = cellfun(@ (x) find(ismember(proteins_set_total{3}(ind1),num2str(x))), num2cell(kinase_location_number(i_first:i_last)), 'UniformOutput', false);
                ind_total = cellfun(@(x,y) intersect(x,y,'stable'), ind_residue,ind_number,'UniformOutput',false);
                % ind_total = cellfun(@(x,y) union(x,y,'stable'), ind_total,ind_peptide1,'UniformOutput',false);
                % ind_total = cellfun(@(x,y) union(x,y,'stable'), ind_total,ind_peptide2,'UniformOutput',false);
                ind = find(cellfun(@(x) isempty(x)==0,ind_total)==1);
                ind = ind+i_first-1;
                kinase_list{j,1} = file(ind,:);
                %         for k = find(ismember(kinase_target,n1(j)),1,'first'):find(ismember(kinase_target,n1(j)),1,'last')
                %
                %            ind_residue = find(ismember(proteins_set_total{2}(ind1),kinase_location_residue(k)));
                %            ind_number =  find(ismember(proteins_set_total{3}(ind1),num2str(kinase_location_number(k))));
                %            if (isempty(intersect(ind_residue,ind_number,'stable'))==0)
                %            ind_total = [ind_total,intersect(ind_residue,ind_number,'stable')];
                %            ind = [ind,k];
                %            end
                %         end
                %                kinase_list{j}{2} = proteins_set_total{2}(ind1(ind_total));
                %                kinase_list{j}{3} = proteins_set_total{3}(ind1(ind_total));
                %                kinase_list{j}{4} = proteins_set_total{4}(ind1(ind_total));
                %                kinase_list{j}{5} = file(ind,:);
                %                kinase_list{j,1} = file(ind,:);
            end
            
            for j = 1:length(kinase_list)
                n1 = unique(kinase_list{j}.Var2);
                as = {};
                for jj = 1:length(n1)
                    ind = find(ismember(kinase_list{j}.Var2,n1(jj)));
                    as{jj} = kinase_list{j}(ind,:);
                end
                kinase_list{j} = as;
            end
            
            
            
            %Sort by the posterior probability and then select the frist three
            %kinase families.
            kinase_list_sorted = cell(length(kinase_list),1);
            for jj = 1:length(kinase_list)
                [~,a2] = cellfun(@(x) sort(x.Var9,'descend'),kinase_list{jj},'UniformOutput',false);
                kinase_list_sorted{jj} = cellfun(@(x,y) x(y,:),kinase_list{jj},a2,'UniformOutput',false);
                var = cellfun(@(x) unique(x.Var8,'stable'),kinase_list_sorted{jj},'UniformOutput',false);
                
                for j = 1:length(var)
                    if(isempty(var{j})==0&&length(var{j})>=3)
                        ind1 = find(ismember(kinase_list_sorted{jj}{j}.Var8,var{j}(1)),1,'first');
                        ind2 = find(ismember(kinase_list_sorted{jj}{j}.Var8,var{j}(2)),1,'first');
                        ind3 = find(ismember(kinase_list_sorted{jj}{j}.Var8,var{j}(3)),1,'first');
                        
                        kinase_list_sorted{jj}{j}(1,:) = kinase_list_sorted{jj}{j}(ind1,:);
                        kinase_list_sorted{jj}{j}(2,:) = kinase_list_sorted{jj}{j}(ind2,:);
                        kinase_list_sorted{jj}{j}(3,:) = kinase_list_sorted{jj}{j}(ind3,:);
                        kinase_list_sorted{jj}{j}(4:end,:) = [];
                    end
                    
                    %For extracting until the first element of the third group is
                    %extracted
                    %kinase_list_sorted{j} = kinase_list_sorted{j}(1:find(ismember(kinase_list_sorted{j}.Var8,var{j}(3)),1,'first'),:);
                end
            end
            
            %Sort by the posterior probability and then select the frist three
            %kinase families for the full list of background proteins
            
            for jj = 1:length(kinase_list_proteomics)
                kinase_list_proteomics{jj} = cellfun(@(x) sortrows(x,9,'descend'),kinase_list_proteomics{jj},'UniformOutput',false);
                
                for j = 1:length(kinase_list_proteomics{jj})
                    var =  unique(kinase_list_proteomics{jj}{j}.Var8,'stable');
                    if(isempty(var)==0&&length(var)>=3)
                        ind1 = find(ismember(kinase_list_proteomics{jj}{j}.Var8,var(1)),1,'first');
                        ind2 = find(ismember(kinase_list_proteomics{jj}{j}.Var8,var(2)),1,'first');
                        ind3 = find(ismember(kinase_list_proteomics{jj}{j}.Var8,var(3)),1,'first');
                        
                        kinase_list_proteomics{jj}{j}(1,:) = kinase_list_proteomics{jj}{j}(ind1,:);
                        kinase_list_proteomics{jj}{j}(2,:) = kinase_list_proteomics{jj}{j}(ind2,:);
                        kinase_list_proteomics{jj}{j}(3,:) = kinase_list_proteomics{jj}{j}(ind3,:);
                        kinase_list_proteomics{jj}{j}(4:end,:) = [];
                    end
                    
                    %For extracting until the first element of the third group is
                    %extracted
                    %kinase_list_sorted{j} = kinase_list_sorted{j}(1:find(ismember(kinase_list_sorted{j}.Var8,var{j}(3)),1,'first'),:);
                end
            end
            
            
            %Compute the three associated protein kinases metrics
            score_kinase = {};
            for jj = 1:length(kinase_list_sorted)
                for j = 1:length(kinase_list_sorted{jj})
                    if(isempty(kinase_list_sorted{jj}{j}.Var8)==0)
                        asa1 = arrayfun(@num2str,kinase_list_sorted{jj}{j}.Var2,'UniformOutput',false);
                        asa2 = cellfun(@(x,y) strcat(x,y),kinase_list_sorted{jj}{j}.Var3,asa1,'UniformOutput',false);
                        score_kinase = [score_kinase;[kinase_list_sorted{jj}{j}.Var8,num2cell(length(kinase_list_sorted{jj}{j}.Var8):-1:1)',kinase_list_sorted{jj}{j}.Var1,asa2]];
                    end
                end
            end
            
            s1 = unique(score_kinase(:,1),'stable');
            score_kinase_final = {};
            for j = 1:length(s1)
                iind = find(ismember(score_kinase(:,1),s1(j)));
                aa1 = num2cell(length(iind));
                aa2 = num2cell(length(iind)./(length(score_kinase(:,1))/3));
                aa3 = num2cell(sum(cell2mat(score_kinase(iind,2)))./(length(score_kinase(:,1))/3));
                aa4 = score_kinase(iind,3);
                aa5 = score_kinase(iind,4);
                aa6 = unique(table(aa4,aa5),'stable');
                aa4 = aa6{:,1};
                aa5 = aa6{:,2};
                aa4 = cellfun(@(x) strcat(x,'||'),aa4,'UniformOutput',false);
                aa4 = [aa4{:}];
                aa4 = aa4(1:end-2);
                aa5 = cellfun(@(x) strcat(x,'||'),aa5,'UniformOutput',false);
                aa5 = [aa5{:}];
                aa5 = aa5(1:end-2);
                score_kinase_final = [score_kinase_final;[s1(j),aa1,aa2,aa3,aa4,aa5]];
            end
            
            [~,a2] = sort(cell2mat(score_kinase_final(:,3)),'descend');
            score_kinase_final = score_kinase_final(a2,:);
            
            %Compute the three associated protein kinases metrics for the full
            %list of background proteins
            score_kinase_proteomics = {};
            for jj = 1:length(kinase_list_proteomics)
                for j = 1:length(kinase_list_proteomics{jj})
                    if(isempty(kinase_list_proteomics{jj}{j}.Var8)==0)
                        score_kinase_proteomics = [score_kinase_proteomics;[kinase_list_proteomics{jj}{j}.Var8,num2cell(length(kinase_list_proteomics{jj}{j}.Var8):-1:1)']];
                    end
                end
            end
            
            s1 = unique(score_kinase_proteomics(:,1),'stable');
            score_kinase_finalProteomics = {};
            for j = 1:length(s1)
                aa1 = num2cell(length(find(ismember(score_kinase_proteomics(:,1),s1(j)))));
                aa2 = num2cell(length(find(ismember(score_kinase_proteomics(:,1),s1(j))))./(length(score_kinase_proteomics(:,1))/3));
                aa3 = num2cell(sum(cell2mat(score_kinase_proteomics(find(ismember(score_kinase_proteomics(:,1),s1(j))),2)))./(length(score_kinase_proteomics(:,1))/3));
                score_kinase_finalProteomics = [score_kinase_finalProteomics;[s1(j),aa1,aa2,aa3]];
            end
            
            [~,a2] = sort(cell2mat(score_kinase_finalProteomics(:,3)),'descend');
            score_kinase_finalProteomics = score_kinase_finalProteomics(a2,:);
            
            
            %Add the data related to the family and kinase names for the
            %identified kinases.
            for j = 1:size(score_kinase_final,1)
                ind = find(ismember(NetPhorest_full_data.NetPhorest,score_kinase_final{j,1}));
                score_kinase_final{j,7} = unique(NetPhorest_full_data.Family(ind),'stable');
                score_kinase_final{j,8} = unique(NetPhorest_full_data.Kin_name(ind),'stable');
                score_kinase_final{j,9} = unique(NetPhorest_full_data.Unip_ID(ind),'stable');
            end
            
            %Write the data
            header = {'NetPhorest_kinase','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Sum_of_rank/Number_of_PhosphoPeptides','Substrate_members','Substrate_members_residues','Family_name','Associated_kinases','Associated_kinases_uniport'};
            
            if(index_interaction == 1)
                dir_results_netphorest = strcat(dir_results(1).folder,'\','Excel_results','\','With_interactions');
                mkdir(dir_results_netphorest)
                xlx_name = strcat(dir_results_netphorest,'\score_kinase_final_byPhosphoPeptides',dir_results(i).name,'_interactions','.xlsx');
            else
                dir_results_netphorest = strcat(dir_results(1).folder,'\','Excel_results','\','No_interactions');
                mkdir(dir_results_netphorest)
                xlx_name = strcat(dir_results_netphorest,'\score_kinase_final_byPhosphoPeptides',dir_results(i).name,'.xlsx');
            end
            delete(xlx_name);
            writecell(header,xlx_name,'Range','A1:I1','AutofitWidth',1);
            writecell(score_kinase_final(:,1:end-3),xlx_name,'Range','A2','AutofitWidth',0);
            aa1 = score_kinase_final(:,end-2);
            aa1 = cellfun(@(x) strcat(x,'||'),aa1,'UniformOutput',false);
            aa1 = cellfun(@(x) [x{:}],aa1,'UniformOutput',false);
            aa1 = cellfun(@(x) x(1:end-2),aa1,'UniformOutput',false);
            writecell(aa1,xlx_name,'Range','G2','AutofitWidth',1);
            aa1 = score_kinase_final(:,end-1);
            aa1 = cellfun(@(x) strcat(x,'||'),aa1,'UniformOutput',false);
            aa1 = cellfun(@(x) [x{:}],aa1,'UniformOutput',false);
            aa1 = cellfun(@(x) x(1:end-2),aa1,'UniformOutput',false);
            writecell(aa1,xlx_name,'Range','H2','AutofitWidth',1);
            aa1 = score_kinase_final(:,end);
            aa1 = cellfun(@(x) strcat(x,'||'),aa1,'UniformOutput',false);
            aa1 = cellfun(@(x) [x{:}],aa1,'UniformOutput',false);
            aa1 = cellfun(@(x) x(1:end-2),aa1,'UniformOutput',false);
            writecell(aa1,xlx_name,'Range','I2','AutofitWidth',1);
            
            statistics = {};
            for j = 1:size(score_kinase_final,1)
                if(isempty(find(ismember(score_kinase_finalProteomics(:,1),score_kinase_final(j,1))))==0)
                    chi_kinase = score_kinase_final{find(ismember(score_kinase_final(:,1),score_kinase_final(j,1))),2};
                    chi_background = score_kinase_finalProteomics{find(ismember(score_kinase_finalProteomics(:,1),score_kinase_final(j,1))),2};
%                     al1 = cellfun(@(x) strsplit(x,"||"),score_kinase_final(:,5),'UniformOutput',false);
%                     al2 = cellfun(@(x) strsplit(x,"||"),score_kinase_final(:,6),'UniformOutput',false);
%                     al3 = cellfun(@(x,y) strcat(x,y),al1,al2,'UniformOutput',false);
%                     al3 = [al3{:}];
                    chi_kinase_total = round(length([kinase_list{:}]));
                    %chi_kinase_total = round(length(score_kinase(:,1))/3);
                    chi_background_total = round(length(vertcat(kinase_list_proteomics{:})));
                    %chi_background_total = round(length(score_kinase_proteomics(:,1))/3);
                    
                    %Pooled estimate of proportions:
                    p0 = (chi_kinase+chi_background) / (chi_kinase_total+chi_background_total);
                    
                    %Expected counts under the hypothesis of H0:
                    chi_kinase_H0 = chi_kinase_total * p0;
                    chi_background_H0 = chi_background_total * p0;
                    
                    observed = [chi_kinase, chi_kinase_total-chi_kinase, chi_background, chi_background_total-chi_background];
                    expected = [chi_kinase_H0, chi_kinase_total-chi_kinase_H0, chi_background_H0, chi_background_total-chi_background_H0];
                    
                    chiScore = sum((observed-expected).^2 ./ expected);
                    p_value_chi = 1-chi2cdf(chiScore,1);
                    statistics{j,1} = chiScore;
                    statistics{j,2} = p_value_chi;
                    
                    as1 = [repmat('a',chi_kinase_total,1); repmat('b',chi_background_total,1)];
                    as2 = [repmat(1,chi_kinase,1); repmat(2,chi_kinase_total-chi_kinase,1); repmat(1,chi_background,1); repmat(2,chi_background_total-chi_background,1)];
                    
                    [tbl,chi2,pval_tab,labels] = crosstab(as1,as2);
                    [h,p_fisher,stats] = fishertest(tbl);
                    
                    statistics{j,3} = stats.OddsRatio;
                    statistics{j,4} = p_fisher;
                    
                end
                
            end
            
            %Correct for the multiple testing
            statistics_aux = {};
            statistics_aux(:,1) = statistics(:,1);
            statistics_aux(:,2) = statistics(:,2);
            ind = find(cellfun(@isempty,statistics(:,2)));
            s1 = num2cell(mafdr([statistics{:,2}]','BHFDR','true'));
            %                 try
            %                     s1 = num2cell(mafdr([statistics{:,2}]'));
            %                 catch
            %                     s1 = num2cell(mafdr([statistics{:,2}]','BHFDR','true'));
            %                 end
            for k1 = 1:length(ind)
                s1 = [s1(1:ind(k1)-1);'NaN';s1(ind(k1):end)];
            end
            statistics_aux = [statistics_aux,s1];
            statistics_aux(:,4) = statistics(:,3);
            statistics_aux(:,5) = statistics(:,4);
            ind = find(cellfun(@isempty,statistics(:,4)));
            s1 = num2cell(mafdr([statistics{:,4}]','BHFDR','true'));
            %                 try
            %                     s1 = num2cell(mafdr([statistics{:,4}]'));
            %                 catch
            %                     s1 = num2cell(mafdr([statistics{:,4}]','BHFDR','true'));
            %                 end
            for k1 = 1:length(ind)
                s1 = [s1(1:ind(k1)-1);'NaN';s1(ind(k1):end)];
            end
            statistics_aux = [statistics_aux,s1];
            
            header_statistics = {'Chi_test_value','Chi_test_pValue','Chi_test_pValue_BH_adjusted','Fisher_test_OddsRatio','Fisher_test_pValue','Fisher_test_pValue_BH_adjusted'};
            writecell(header_statistics,xlx_name,'Range','J1:O1','AutofitWidth',1);
            writecell(statistics_aux,xlx_name,'Range','J2','AutofitWidth',1);
            
        end
        
    end
end



end