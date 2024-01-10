function [fc,pp,ind] = protein_logfc_compensation(protein_logfc_compensation_file_name,sample_names,p_names,fc,stder,df_pooled,pp,up_down_val,j,ind_i,ind_k,ind_fc)

Protein_lfc_comp = protein_logfc_compensation_file_name;
opts = detectImportOptions(Protein_lfc_comp);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
Protein_lfc_comp = readtable(Protein_lfc_comp, opts);
v1 = strrep(sample_names{ind_i},'_','');
v2 = strrep(sample_names{ind_k+ind_i-1},'_','');
str1 = strcat(v1,'vs',v2);
ind1 = find(ismember(Protein_lfc_comp.Properties.VariableNames,'UniProt ID'));
ind2 = find(contains(Protein_lfc_comp.Properties.VariableNames,strcat(v1,'_R')));
ind3 = find(contains(Protein_lfc_comp.Properties.VariableNames,strcat(v2,'_R')));
ind_aux = find(contains(Protein_lfc_comp.Properties.VariableNames,str1));
ind4 = ind_aux(find(contains(Protein_lfc_comp.Properties.VariableNames(ind_aux),'log2FC')));
ind5 = ind_aux(find(contains(Protein_lfc_comp.Properties.VariableNames(ind_aux),'FDR')));
%index for the unique hits
ind6 = find(contains(Protein_lfc_comp.Properties.VariableNames,strcat(v1,'_Group')));
ind7 = find(contains(Protein_lfc_comp.Properties.VariableNames,strcat(v2,'_Group')));
Protein_lfc_comp = Protein_lfc_comp(:,[ind1,ind4,ind5,ind2,ind3,ind6,ind7]);
Protein_lfc_comp{:,2:end} = num2cell(cellfun(@(x) str2num(x), Protein_lfc_comp{:,2:end}));

if(up_down_val==1)
    ind1 = find(cell2mat(Protein_lfc_comp{:,2})>0);
%     ind2 = find(cell2mat(Protein_lfc_comp{:,3})<0.05);
     ind3 = find(cell2mat(Protein_lfc_comp{:,end})>0);
end
if(up_down_val==-1)
    ind1 = find(cell2mat(Protein_lfc_comp{:,2})<0);
 %    ind2 = find(cell2mat(Protein_lfc_comp{:,3})<0.05);
     ind3 = find(cell2mat(Protein_lfc_comp{:,end-1})>0);
end
% ind = intersect(ind1,ind2);
 %Protein_lfc_comp = Protein_lfc_comp(intersect(ind,ind3),:);
Protein_lfc_comp = Protein_lfc_comp(intersect(ind1,ind3),:);
vec_ref = [];
for irk = 1:size(Protein_lfc_comp,1)
    vec_comp1 = cell2mat(Protein_lfc_comp{irk,find(contains(Protein_lfc_comp.Properties.VariableNames,strcat(v1,'_R')))});
    vec_comp2 = cell2mat(Protein_lfc_comp{irk,find(contains(Protein_lfc_comp.Properties.VariableNames,strcat(v2,'_R')))});
    std1 = std(vec_comp1,'omitnan');
    std2 = std(vec_comp2,'omitnan');
    n1 = length(isnan(vec_comp1));
    n2 = length(isnan(vec_comp2));
    sp = sqrt((std1^2*(n1-1)+std2^2*(n2-1))/(n1+n2-2));
    %std_ref = sqrt((std1/sqrt(n1-1))^2 + (std2/sqrt(n2-1))^2);
    std_ref = sp*sqrt(1/n1+1/n2);
    %df_ref = WelchSatterwhaite(std1,std2,n1,n2);
    df_ref = n1+n2-2;
    lgfc_ref = mean(vec_comp1,'omitnan')-mean(vec_comp2,'omitnan');
    if(std_ref == 0)  vec_ref = [vec_ref;[0,0,0]];
    else
        vec_ref = [vec_ref;[lgfc_ref,std_ref,df_ref]];
    end
end
vec_ref = [Protein_lfc_comp(:,1),cell2table(num2cell(vec_ref))];
vec_ref.Properties.VariableNames(2:end) = [{'lgfc'},{'std'},{'df'}];

vec_u = [table(p_names),cell2table(num2cell([fc{j,ind_fc}',stder',df_pooled']))];
vec_u.Properties.VariableNames = vec_ref.Properties.VariableNames;
if(up_down_val==1)
%ind1 = intersect(find(vec_u{:,2}>1),find(pp<0.05));
ind1 = find(vec_u{:,2}>0);
end
if(up_down_val==-1)
%ind1 = intersect(find(vec_u{:,2}<-1),find(pp<0.05));
ind1 = find(vec_u{:,2}<0);
end
ind2 = find(contains(vec_u{:,1},vec_ref{:,1}));
ind = intersect(ind1,ind2);
vec_u1 = vec_u(ind,:);
pp_correction_val = [];
amx = [];
for irk = 1:size(vec_u1,1)
    name_p = extractBetween(vec_u1{irk,1},'|','|');
    ref = vec_ref(find(ismember(vec_ref{:,1},name_p)),:);
    log2fc = abs(vec_u1{irk,2}) - abs(ref{1,2});
    amx(irk) = log2fc;
    s2 = (vec_u1{irk,3}*sqrt(2))^2;
    s2_ref = (ref{1,3}*sqrt(2))^2;
    n_s2 = vec_u1{irk,4};
    n_s2_ref = ref{1,4};
    if(s2==0 || s2_ref==0) pp_correction_val = [pp_correction_val,1];
    else
    %Method 1 for equal variance assumption between groups
    %- decent assumption
    sp = sqrt((s2*(n_s2-1)+s2_ref*(n_s2_ref-1))/(n_s2+n_s2_ref-2));
    stderr = sp*sqrt(1/n_s2+1/n_s2_ref);
    df = n_s2+n_s2_ref-2;
    %Method 2 for unequal variance assumption - more
    %drastic
   % stderr = sqrt((vec_u1{irk,3}/sqrt(n_s2-1))^2 + (ref{1,3}/sqrt(n_s2_ref-1))^2);
%    df = WelchSatterwhaite(vec_u1{irk,3},ref{1,3},vec_u1{irk,4},ref{1,4});
%     Method 3 - from MSStatPTM but becauese we use t-test insetead of lm it is too drastic
%     stderr = sqrt(s2+s2_ref);
%     numer = (s2 + s2_ref)^2;
%     denom = (s2^2/vec_u1{irk,4} + s2_ref^2/ref{1,4});
%     df = numer/denom;
    tval = log2fc/stderr;
    % If we want to consider the cases when the proteomics is higher than
    % the phosphoproteomics and remove these by default
    if(tval<0) pp_correction_val = [pp_correction_val,1];
    else
      %pp_correction_val = [pp_correction_val,(1-tcdf(abs(tval),df))];  if
      %we would like a double sided test we would use 2*(...). However, our
      %test is unilateral in this case.
    pp_correction_val = [pp_correction_val,(1-tcdf(abs(tval),df))];
    end
    end
end
pp_correction_val = mafdr(pp_correction_val,'BHFDR','true');
ind = ind(find(pp_correction_val>0.05));
is1 = find(ismember(vec_u{ind,1},vec_u1{find(pp_correction_val<0.05),1}));
ind(is1) = [];
fc{j,ind_fc}(ind) = 0;
pp(ind) = 1;

end