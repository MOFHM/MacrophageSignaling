function [list,seed_list_kin,seed_list_substrate] = kinase_cascade_generationSubstrate(sub,seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,index_stop,upstream_kinase,upstream_residue)


list = {};

% ind_TF = [];
% if(index_stop==0)
% ind_TF = find(ismember(lower(file_TF_mapping{:,2}),lower(sub)));
% end
%
% if(isempty(ind_TF)==0)
%     list = [list;[sub,sub,num2cell(1),{'aaa'}]];
% else
if(index_stop==0)
    ind_PSP = find(ismember(lower(file_kinases_Phospho{:,7}),lower(sub)));
    
    if(isempty(ind_PSP)==0)
        %Extract all the kinases from PSP
        a_uniport_residue = table(file_kinases_Phospho{ind_PSP,3},file_kinases_Phospho{ind_PSP,10});
        a2 = ones(size(a_uniport_residue,1),1).*10;
        %Consider as hit only when the kinase was measured as well as the
        %specific phosphorylated site!
        ind_uniport_residue = find(ismember(lower(a_uniport_residue{:,1}),lower(file_kinase_total{:,1})));
        a2(ind_uniport_residue) = 5;
        ind_aux = find(ismember(lower(table(repmat(sub,size(a_uniport_residue{:,2},1),1),a_uniport_residue{:,2})),lower(table(file_kinase_total{:,1},file_kinase_total{:,4}))));
        ind_uniport_residue = intersect(ind_uniport_residue,ind_aux);
        a2(ind_uniport_residue) = 2;
        ind_TF = find(ismember(lower(a_uniport_residue{:,1}),lower(file_TF_measured{:,2})));
        a2(ind_TF) = 1;
        ll = [a_uniport_residue{:,1},repmat(sub,size(a_uniport_residue,1),1),num2cell(a2),a_uniport_residue{:,2}];
        list = [list;ll];
        %Check the correctness
        %         iind_TF = find(ismember(lower(file_kinases_Phospho{ind_PSP,3}),lower(file_TF_mapping{:,2})));
        %         if(isempty(iind_TF)==0)
        %             list = [list;[file_kinases_Phospho{ind_PSP(iind_TF),3},repmat(sub,length(iind_TF),1),repmat(num2cell(1),length(iind_TF),1),repmat({'aaa'},length(iind_TF),1)]];
        %         end
        %Check the correctenss
        
        for i = 1:size(a_uniport_residue,1)
            if(isempty(find(ismember(lower(seed_list_kin{:,1}),lower(a_uniport_residue{i,1})))))
                if(a2(i)==10) indexx = 1; else indexx = 0; end
                seed_list_kin = [seed_list_kin;a_uniport_residue{i,1}];
                seed_list_substrate = [seed_list_substrate;a_uniport_residue{i,:}];
                [ll,sk,ss] = kinase_cascade_generationSubstrate(a_uniport_residue{i,1},seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,indexx,sub,a_uniport_residue{i,2});
                list = [list;ll];
                list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                seed_list_kin = [seed_list_kin;sk];
                seed_list_kin = unique(seed_list_kin,'stable');
                seed_list_substrate = [seed_list_substrate;ss];
                seed_list_substrate = unique(seed_list_substrate,'stable');
                [ll,ss,sk] = kinase_cascade_generationKinase(a_uniport_residue{i,1},seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,indexx,sub,a_uniport_residue{i,2});
                list = [list;ll];
                list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                seed_list_substrate = [seed_list_substrate;ss];
                seed_list_substrate = unique(seed_list_substrate,'stable');
                seed_list_kin = [seed_list_kin;sk];
                seed_list_kin = unique(seed_list_kin,'stable');
            end
        end
    end
end



if(index_stop==1)
ind_TF = find(ismember(lower(file_TF_mapping{:,2}),lower(sub)));
    ind_PSP = find(ismember(lower(file_kinases_Phospho{:,7}),lower(sub)));

    if(isempty(ind_PSP)&&isempty(ind_TF)) list = [list;[upstream_kinase,sub,num2cell(101),upstream_residue]]; end
        
    if(isempty(ind_PSP)==0)
        
        a_uniport_residue = table(file_kinases_Phospho{ind_PSP,3},file_kinases_Phospho{ind_PSP,10});
        %Consider as hit only when the kinase was measured as well as the
        %specific phosphorylated site!
        ind = find(ismember(a_uniport_residue{:,1},upstream_kinase));
        a_uniport_residue(ind,:) = [];
        ind_TF1 = find(ismember(lower(a_uniport_residue{:,1}),lower(file_TF_measured{:,2})));
        if(isempty(ind_TF1)==0)
            ll = [a_uniport_residue{ind_TF1,1},repmat(sub,size(a_uniport_residue(ind_TF1,1),1),1),repmat(num2cell(110),size(a_uniport_residue(ind_TF1,1),1),1),a_uniport_residue{ind_TF1,2}];
            list = [list;ll];
        end
        ind_uniport_residue = find(ismember(lower(a_uniport_residue{:,1}),lower(file_kinase_total{:,1})));
        a_uniport_residue = a_uniport_residue(ind_uniport_residue,:);
        %Check to do not identify the upstream kinase and introduce the
        %same link with two different scores - one 10 and one 210
        a2 = ones(size(a_uniport_residue,1),1).*210;
        
        if(isempty(ind_uniport_residue)&&isempty(ind_TF)&&isempty(ind_TF1)) list = [list;[upstream_kinase,sub,num2cell(101),upstream_residue]]; end
        
        ll = [a_uniport_residue{:,1},repmat(sub,size(a_uniport_residue,1),1),num2cell(a2),a_uniport_residue{:,2}];
        list = [list;ll];
        %Check the correctness
        %         iind_TF = find(ismember(lower(file_kinases_Phospho{ind_PSP,3}),lower(file_TF_mapping{:,2})));
        %         if(isempty(iind_TF)==0)
        %             list = [list;[file_kinases_Phospho{ind_PSP(iind_TF),3},repmat(sub,length(iind_TF),1),repmat(num2cell(1),length(iind_TF),1),repmat({'aaa'},length(iind_TF),1)]];
        %         end
        %Check the correctenss
        
        for i = 1:size(a_uniport_residue,1)
            if(isempty(find(ismember(lower(seed_list_kin{:,1}),lower(a_uniport_residue{i,1})))))
                indexx = 0;
                seed_list_kin = [seed_list_kin;a_uniport_residue{i,1}];
                seed_list_substrate = [seed_list_substrate;a_uniport_residue{i,:}];
                [ll,sk,ss] = kinase_cascade_generationSubstrate(a_uniport_residue{i,1},seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,indexx,sub,a_uniport_residue{i,2});
                list = [list;ll];
                list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                seed_list_kin = [seed_list_kin;sk];
                seed_list_kin = unique(seed_list_kin,'stable');
                seed_list_substrate = [seed_list_substrate;ss];
                seed_list_substrate = unique(seed_list_substrate,'stable');
                [ll,ss,sk] = kinase_cascade_generationKinase(a_uniport_residue{i,1},seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,indexx,sub,a_uniport_residue{i,2});
                list = [list;ll];
                list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                seed_list_substrate = [seed_list_substrate;ss];
                seed_list_substrate = unique(seed_list_substrate,'stable');
                seed_list_kin = [seed_list_kin;sk];
                seed_list_kin = unique(seed_list_kin,'stable');
            end
        end
    end
end










if(index_stop==0)
    %Filter the file kinases_NetPhorest by the posterior probability
    aa1 = file_kinases_NetPhorest{:,5};
    aa1 = cellfun(@(x) strsplit(x,'||'),aa1,'UniformOutput',false);
    aa2 = file_kinases_NetPhorest{:,6};
    aa2 = cellfun(@(x) strsplit(x,'||'),aa2,'UniformOutput',false);
    NetPhorest_uniport_kinase = table(aa1,aa2);
    ind_NET_substrate = cellfun(@(x) find(ismember(lower(x),lower(sub))),NetPhorest_uniport_kinase{:,1},'UniformOutput',false);
    %Select only the kinase families that phosphorylated the specific site
    %on sub
    ind_NET = ind_NET_substrate;
    %     ind_NET_residue = cellfun(@(x) find(ismember(lower(x),lower(sub_residue))),NetPhorest_uniport_kinase{:,2},'UniformOutput',false);
    %     ind_NET = cellfun(@(x,y) intersect(x,y),ind_NET_substrate,ind_NET_residue,'UniformOutput',false);
    iind_NET = cellfun(@isempty,ind_NET);
    iind_NET = find(iind_NET==0);
    NET_substrate_residues = cellfun(@(x,y) {aa2{x}{[y]}},num2cell(iind_NET),ind_NET(iind_NET),'UniformOutput',false);
    NetPhorest_uniport_kinase_associated_substrates = file_kinases_NetPhorest{iind_NET,9};
    NetPhorest_uniport_kinase_associated_substrates = cellfun(@(x) strsplit(x,'||'),NetPhorest_uniport_kinase_associated_substrates,'UniformOutput',false);
    NET_substrate_residues = cellfun(@(x,y) repmat(x,length(y),1),NET_substrate_residues,NetPhorest_uniport_kinase_associated_substrates,'UniformOutput',false);
    col_dim = max(cellfun(@(x) size(x,2),NET_substrate_residues));
    for ix = 1:size(NET_substrate_residues)
        if(size(NET_substrate_residues{ix},2)<col_dim)
            NET_substrate_residues{ix}(:,size(NET_substrate_residues{ix},2)+1:col_dim) = repmat(NET_substrate_residues{ix}(:,1),1,col_dim-size(NET_substrate_residues{ix},2));
        end
    end
    NET_substrate_residues = vertcat(NET_substrate_residues{:});
    NetPhorest_uniport_kinase_associated_substrates = [NetPhorest_uniport_kinase_associated_substrates{:}]';
    if (isempty(NetPhorest_uniport_kinase_associated_substrates)==0)
        
        %     NetPhorest_uniport_kinase_associated_substrates = unique(NetPhorest_uniport_kinase_associated_substrates,'stable');
        iind  = find(ismember(lower(NetPhorest_uniport_kinase_associated_substrates),lower(file_kinase_total{:,1})));
        NetPhorest_uniport_kinase_associated_substrates =  NetPhorest_uniport_kinase_associated_substrates(iind);
        NET_substrate_residues = NET_substrate_residues(iind,:);
        if(isempty(NetPhorest_uniport_kinase_associated_substrates)==0)
            aa1 = NetPhorest_uniport_kinase_associated_substrates;
            aa2 = NET_substrate_residues;
            ll = {};
            for ixj = 1:size(aa2,2)
                ll = [ll;[aa1,repmat(sub,size(aa1)),num2cell(repmat(3,size(aa1))),aa2(:,ixj)]];
            end
            list = [list;ll];
            for i = 1:length(NetPhorest_uniport_kinase_associated_substrates)
                if(isempty(find(ismember(lower(seed_list_kin{:,1}),lower(NetPhorest_uniport_kinase_associated_substrates{i})))))
                    seed_list_kin = [seed_list_kin;NetPhorest_uniport_kinase_associated_substrates(i)];
                    seed_list_substrate = [seed_list_substrate;[NetPhorest_uniport_kinase_associated_substrates(i),NET_substrate_residues(i)]];
                    [ll,sk,ss] = kinase_cascade_generationSubstrate(NetPhorest_uniport_kinase_associated_substrates(i),seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,0,sub,'hope_not_appear');
                    list = [list;ll];
                    list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                    seed_list_kin = [seed_list_kin;sk];
                    seed_list_kin = unique(seed_list_kin,'stable');
                    seed_list_substrate = [seed_list_substrate;ss];
                    seed_list_substrate = unique(seed_list_substrate,'stable');
                    [ll,ss,sk] = kinase_cascade_generationKinase(NetPhorest_uniport_kinase_associated_substrates(i),seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,0,sub,'hope_not_appear');
                    list = [list;ll];
                    list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                    seed_list_substrate = [seed_list_substrate;ss];
                    seed_list_substrate = unique(seed_list_substrate,'stable');
                    seed_list_kin = [seed_list_kin;sk];
                    seed_list_kin = unique(seed_list_kin,'stable');
                end
            end
        end
    end
end
%end
