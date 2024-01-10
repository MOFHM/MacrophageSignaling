function [list,seed_list_substrate,seed_list_kin] = kinase_cascade_generationKinase(kin,seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,index_stop,upstream_kinase,upstream_residue)

%Put an initial condition of index_stop == 1 for trating differently the
%data - if a2(i) == 10 then select only the data that will have a2(i) = 2
%or 5

list = {};

if(index_stop==0)
    ind_PSP = find(ismember(lower(file_kinases_Phospho{:,3}),lower(kin)));
    
    if(isempty(ind_PSP)==0)
        %Extract all the substrates from PSP and then select only the ones
        %that are kinases or TF
        a_uniport_residue = table(file_kinases_Phospho{ind_PSP,7},file_kinases_Phospho{ind_PSP,10});
        ind1 = (find(ismember(a_uniport_residue{:,1},file_kinases_list{:,1})));
        ind2 = (find(ismember(a_uniport_residue{:,1},file_TF_mapping{:,2})));
        ind = union(ind1,ind2);
        a_uniport_residue = a_uniport_residue(ind,:);
        a2 = ones(size(a_uniport_residue,1),1).*10;
        %Consider as hit only when the kinase was measured as well as the
        %specific phosphorylated site!
        %Verify that the substrate is in the dataset (without considering
        %the residue)
        ind_uniport_residue = find(ismember(lower(a_uniport_residue(:,1)),lower(table(file_kinase_total{:,1}))));
        a2(ind_uniport_residue) = 5;
        %
        ind_uniport_residue = find(ismember(lower(a_uniport_residue),lower(table(file_kinase_total{:,1},file_kinase_total{:,4}))));
        ind_aux = find(ismember(lower(repmat(kin,size(a_uniport_residue,1),1)),lower(file_kinase_total{:,1})));
        ind_uniport_residue = intersect(ind_uniport_residue,ind_aux);
        a2(ind_uniport_residue) = 2;
        ind_TF = find(ismember(lower(a_uniport_residue{:,1}),lower(file_TF_measured{:,2})));
        a2(ind_TF) = 1;
        %If it comes from a 10 value dataset it should remain
        %with that value eventhough the other kinase/substrate
        %was measured
        ll = [repmat(kin,size(a_uniport_residue,1),1),a_uniport_residue{:,1},num2cell(a2),a_uniport_residue{:,2}];
        list = [list;ll];
        %Check the correctness
        %         iind_TF = find(ismember(lower(file_kinases_Phospho{ind_PSP,7}),lower(file_TF_mapping{:,2})));
        %         if(isempty(iind_TF)==0)
        %             list = [list;[repmat(kin,length(iind_TF),1),file_kinases_Phospho{ind_PSP(iind_TF),7},repmat(num2cell(1),length(iind_TF),1),repmat({'aaa'},length(iind_TF),1)]];
        %         end
        %Check the correctenss
        
        for i = 1:size(a_uniport_residue,1)
            if(isempty(find(ismember(lower(seed_list_substrate),lower(a_uniport_residue(i,:))))))
                if(a2(i)==10) indexx = 1; else indexx = 0; end
                seed_list_substrate = [seed_list_substrate;a_uniport_residue{i,:}];
                seed_list_kin = [seed_list_kin;a_uniport_residue{i,1}];
                [ll,ss,sk] = kinase_cascade_generationKinase(a_uniport_residue{i,1},seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,indexx,kin,a_uniport_residue{i,2});
                list = [list;ll];
                list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                seed_list_substrate = [seed_list_substrate;ss];
                seed_list_substrate = unique(seed_list_substrate,'stable');
                seed_list_kin = [seed_list_kin;sk];
                seed_list_kin = unique(seed_list_kin,'stable');
                [ll,sk,ss] = kinase_cascade_generationSubstrate(a_uniport_residue{i,1},seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,indexx,kin,a_uniport_residue{i,2});
                list = [list;ll];
                list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                seed_list_kin = [seed_list_kin;sk];
                seed_list_kin = unique(seed_list_kin,'stable');
                seed_list_substrate = [seed_list_substrate;ss];
                seed_list_substrate = unique(seed_list_substrate,'stable');
            end
        end
    end
end


if(index_stop==1)
ind_TF = find(ismember(lower(file_TF_mapping{:,2}),lower(kin)));
    ind_PSP = find(ismember(lower(file_kinases_Phospho{:,3}),lower(kin)));
    
    if(isempty(ind_PSP)&&isempty(ind_TF)) list = [list;[kin,upstream_kinase,num2cell(101),upstream_residue]]; end
    
    if(isempty(ind_PSP)==0)
        
        a_uniport_residue = table(file_kinases_Phospho{ind_PSP,7},file_kinases_Phospho{ind_PSP,10});
        %Consider as hit only when the kinase was measured as well as the
        %specific phosphorylated site!
        %Verify that the substrate is in the dataset (without considering
        %the residue)
        ind = find(ismember(a_uniport_residue{:,1},upstream_kinase));
        a_uniport_residue(ind,:) = [];
          ind_TF1 = find(ismember(lower(a_uniport_residue{:,1}),lower(file_TF_measured{:,2})));
        if(isempty(ind_TF1)==0)
            ll = [repmat(kin,size(a_uniport_residue(ind_TF1,1),1),1),a_uniport_residue{ind_TF1,1},repmat(num2cell(110),size(a_uniport_residue(ind_TF1,1),1),1),a_uniport_residue{ind_TF1,2}];
            list = [list;ll];
        end
        ind_uniport_residue = find(ismember(lower(a_uniport_residue(:,1)),lower(table(file_kinase_total{:,1}))));
        a_uniport_residue = a_uniport_residue(ind_uniport_residue,:);
        a2 = ones(size(a_uniport_residue,1),1).*510;
        if(isempty(ind_uniport_residue)&&isempty(ind_TF)&&isempty(ind_TF1)) list = [list;[kin,upstream_kinase,num2cell(101),upstream_residue]]; end
        %
        ind_uniport_residue = find(ismember(lower(a_uniport_residue),lower(table(file_kinase_total{:,1},file_kinase_total{:,4}))));
        a2(ind_uniport_residue) = 210;
        %If it comes from a 10 value dataset it should remain
        %with that value eventhough the other kinase/substrate
        %was measured
        ll = [repmat(kin,size(a_uniport_residue,1),1),a_uniport_residue{:,1},num2cell(a2),a_uniport_residue{:,2}];
        list = [list;ll];
        
        %Check the correctness
        %         iind_TF = find(ismember(lower(file_kinases_Phospho{ind_PSP,7}),lower(file_TF_mapping{:,2})));
        %         if(isempty(iind_TF)==0)
        %             list = [list;[repmat(kin,length(iind_TF),1),file_kinases_Phospho{ind_PSP(iind_TF),7},repmat(num2cell(1),length(iind_TF),1),repmat({'aaa'},length(iind_TF),1)]];
        %         end
        %Check the correctenss
        
        for i = 1:size(a_uniport_residue,1)
            if(isempty(find(ismember(lower(seed_list_substrate),lower(a_uniport_residue(i,:))))))
                indexx = 0;
                seed_list_substrate = [seed_list_substrate;a_uniport_residue{i,:}];
                seed_list_kin = [seed_list_kin;a_uniport_residue{i,1}];
                [ll,ss,sk] = kinase_cascade_generationKinase(a_uniport_residue{i,1},seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,indexx,kin,a_uniport_residue{i,2});
                list = [list;ll];
                list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                seed_list_substrate = [seed_list_substrate;ss];
                seed_list_substrate = unique(seed_list_substrate,'stable');
                seed_list_kin = [seed_list_kin;sk];
                seed_list_kin = unique(seed_list_kin,'stable');
                [ll,sk,ss] = kinase_cascade_generationSubstrate(a_uniport_residue{i,1},seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,indexx,kin,a_uniport_residue{i,2});
                list = [list;ll];
                list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                seed_list_kin = [seed_list_kin;sk];
                seed_list_kin = unique(seed_list_kin,'stable');
                seed_list_substrate = [seed_list_substrate;ss];
                seed_list_substrate = unique(seed_list_substrate,'stable');
            end
        end
    end
end














if(index_stop==0)
    NetPhorest_uniport_kinase = file_kinases_NetPhorest{:,9};
    NetPhorest_uniport_kinase = cellfun(@(x) strsplit(x,'||'),NetPhorest_uniport_kinase,'UniformOutput',false);
    ind_NET = cellfun(@(x) find(ismember(lower(x),lower(kin))),NetPhorest_uniport_kinase,'UniformOutput',false);
    iind = cellfun(@isempty,ind_NET);
    iind = find(iind==0);
    aa1 = file_kinases_NetPhorest{iind,5};
    aa1 = cellfun(@(x) strsplit(x,'||'),aa1,'UniformOutput',false);
    aa2 = file_kinases_NetPhorest{iind,6};
    aa2 = cellfun(@(x) strsplit(x,'||'),aa2,'UniformOutput',false);
    NetPhorest_uniport_kinase_associated_substrates = table(aa1,aa2);
    NetPhorest_uniport_kinase_associated_substrates = NetPhorest_uniport_kinase_associated_substrates{:,:};
    aaj = {};
    for ixj = 1:size(NetPhorest_uniport_kinase_associated_substrates,1)
        aaj = [aaj;vertcat(NetPhorest_uniport_kinase_associated_substrates{ixj,:})'];
    end
    if(isempty(aaj)==0)
        NetPhorest_uniport_kinase_associated_substrates = table(aaj(:,1),aaj(:,2));
        iind  = find(ismember(lower(NetPhorest_uniport_kinase_associated_substrates),lower(table(file_kinase_total{:,1},file_kinase_total{:,4}))));
        %NetPhorest_uniport_kinase_associated_substrates = table2cell(NetPhorest_uniport_kinase_associated_substrates);
        NetPhorest_uniport_kinase_associated_substrates =  NetPhorest_uniport_kinase_associated_substrates(iind,:);
        
        
        if(isempty(NetPhorest_uniport_kinase_associated_substrates)==0)
            len = size(NetPhorest_uniport_kinase_associated_substrates,1);
            ll = [repmat(kin,len,1),NetPhorest_uniport_kinase_associated_substrates{:,1},num2cell(repmat(3,len,1)),NetPhorest_uniport_kinase_associated_substrates{:,2}];
            list = [list;ll];
            for i = 1:len
                if(isempty(find(ismember(lower(seed_list_substrate),lower(NetPhorest_uniport_kinase_associated_substrates(i,:))))))
                    seed_list_substrate = [seed_list_substrate;NetPhorest_uniport_kinase_associated_substrates{i,:}];
                    seed_list_kin = [seed_list_kin;NetPhorest_uniport_kinase_associated_substrates{i,1}];
                    [ll,ss,sk] = kinase_cascade_generationKinase(NetPhorest_uniport_kinase_associated_substrates{i,1},seed_list_substrate,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_kin,0,kin,'hope_not_appear');
                    list = [list;ll];
                    list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                    seed_list_substrate = [seed_list_substrate;ss];
                    seed_list_substrate = unique(seed_list_substrate,'stable');
                    seed_list_kin = [seed_list_kin;sk];
                    seed_list_kin = unique(seed_list_kin,'stable');
                    [ll,sk,ss] = kinase_cascade_generationSubstrate(NetPhorest_uniport_kinase_associated_substrates{i,1},seed_list_kin,file_kinases_Phospho,file_TF,file_TF_mapping,file_kinases_NetPhorest,file_kinase_total,file_TF_measured,file_kinases_list,seed_list_substrate,0,kin,'hope_not_appear');
                    list = [list;ll];
                    list = cell2table(list); list = unique(list,'stable'); list = table2cell(list);
                    seed_list_kin = [seed_list_kin;sk];
                    seed_list_kin = unique(seed_list_kin,'stable');
                    seed_list_substrate = [seed_list_substrate;ss];
                    seed_list_substrate = unique(seed_list_substrate,'stable');
                end
            end
        end
    end
end
