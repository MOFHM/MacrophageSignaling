function [p_names_prot,p_names_LFQmean] = normalize_abundance_v1(protein_file,nr_rep,nr_sample_groups,p_LFQ_columns,header_name_string)

%Read the proteomcis file
opts = detectImportOptions(protein_file);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
protein = readtable(protein_file, opts);

%Remove contaminants
ind = find(contains(lower(protein.Properties.VariableNames),lower('Potential contaminant')));
if(isempty(ind)) error('The column "Potential contaminant" does not exist in the proteomics file used of scaling. Introduce it, or contact the authors'); end
a=protein{:,ind};
mark = a(find(~cellfun(@isempty,a),1));
ind = [];
for i = 1:length(a)
    if (isequal(a{i},mark{:})) ind = [ind;i]; end
end
protein_cont = protein;
protein_cont(ind,:)=[];

%Remove reverse
ind = find(contains(lower(protein.Properties.VariableNames),lower('Reverse')));
if(isempty(ind)) error('The column "Reverse" does not exist in the proteomics file used of scaling. Introduce it, or contact the authors'); end
a=protein_cont{:,ind};
mark = a(find(~cellfun(@isempty,a),1));
ind = [];
for i = 1:length(a)
    if (isequal(a{i},mark{:})) ind = [ind;i]; end
end
protein_cont_rev = protein_cont;
protein_cont_rev(ind,:)=[];

%Remove identifide by site
ind = find(contains(lower(protein.Properties.VariableNames),lower('Only identified by site')));
if(isempty(ind)) error('The column "Only identified by site" does not exist in the proteomics file used of scaling. Introduce it, or contact the authors'); end
a=protein_cont_rev{:,ind};
mark = a(find(~cellfun(@isempty,a),1));
ind = [];
for i = 1:length(a)
    if (isequal(a{i},mark{:})) ind = [ind;i]; end
end
 protein_cont_rev_site = protein_cont_rev;
protein_cont_rev_site(ind,:)=[];

%Extract the proteins and the measured intensities
col_names = protein.Properties.VariableNames;
ind = find(ismember(lower(protein.Properties.VariableNames),lower('Protein IDs')));
if(isempty(ind)) error('The column "Protein IDs" does not exist in the proteomics file used of scaling. Introduce it, or contact the authors'); end
p_names_aux = protein_cont_rev_site{:,ind};
p_LFQ = protein_cont_rev_site{:,p_LFQ_columns};
p_LFQ = cellfun(@str2num,p_LFQ);
p_names_prot = {};
p_names_LFQ = [];

%Rearange the header
header = col_names(p_LFQ_columns);
for i = 1:length(header)
    ind = strfind(header{i},header_name_string);
    header{i} = header{i}(ind:end);
end

[~,ind] = sort(header);
header = header(ind);
col_names = col_names(ind);
p_LFQ = p_LFQ(:,ind);

%Extract all the proteins and copy the corresponding intensity levels (in
%one cell there are initially multiple proteins)
for i = 1:length(p_names_aux)
size_ini = length(p_names_prot);
names = split(p_names_aux(i),';');
p_names_prot(size_ini+1:size_ini+length(names),1) = names;
p_names_LFQ(size_ini+1:size_ini+length(names),:) = repmat(p_LFQ(i,:),length(names),1);
end

%Cut all the proteins that were not measured in at least two replicates of
%a phenotype - when all three phenotypes are missing replace them with the mean
%of all the measurements, when only one phenotype is measured copy its values to all the other phenotypes,
%when two or more phenotypes are measured but not all, compute the mean of
%the mean intensitites of the measured phenotypes and atribute it to all of
%them (both measured and not measured). Only for the proteins that are
%measured in all phenotypes compute different means values - this is done
%such that to alter the phosphoproteomics results only when it is correct
%to do so.
p_names_LFQfin = p_names_LFQ;
ind_cut = [];
for i = 1:size(p_names_LFQ,1)
    ind = [];
for j = 1:nr_rep:size(p_names_LFQ,2)
vec = p_names_LFQfin(i,j:j+nr_rep-1);
ind((j-1)/nr_rep+1) = length(find(vec==0));
if (ind((j-1)/nr_rep+1)>0.5*nr_rep) p_names_LFQfin(i,j:j+nr_rep-1) = 0; end
end
nn = find(ind<=0.5*nr_rep);
if (isempty(nn)) ind_cut = [ind_cut,i]; end
end
p_names_LFQfin(ind_cut,:) = 0;
%p_names_prot(ind_cut,:) = [];

for i = 1:size(p_names_LFQfin,1)
    vec = {};
    ind = [];
for j = 1:nr_rep:size(p_names_LFQfin,2)
     vec{(j-1)/nr_rep+1} = p_names_LFQfin(i,j:j+nr_rep-1);
     ind((j-1)/nr_rep+1) = length(find(vec{(j-1)/nr_rep+1}==0));
end
nr_zeros = find(ind==nr_rep);
nr_ref = find(ind==min(ind));
nr_mean = find(ind~=nr_rep);

if(length(nr_zeros)==nr_sample_groups) p_names_LFQfin(i,:) = mean(p_names_LFQ(find(p_names_LFQ~=0))); end
if(length(nr_zeros)==nr_sample_groups-1)
    for j = 1:length(nr_zeros)
        ind = (nr_zeros(j)-1)*nr_rep+1;
        ind_ref = (nr_ref-1)*nr_rep+1;
         p_names_LFQfin(i,ind:ind+nr_rep-1) = p_names_LFQfin(i,ind_ref:ind_ref+nr_rep-1);
    end
end
if((length(nr_zeros)<nr_sample_groups-1)&&(length(nr_zeros)>0))
    mean_val = 0;
    for j = 1:length(nr_mean)
        ind = (nr_mean(j)-1)*nr_rep+1;
        pp = p_names_LFQfin(i,ind:ind+nr_rep-1);
        mean_val = mean_val+mean(pp(find(pp~=0)));
    end
    mean_val = mean_val/j;
    p_names_LFQfin(i,:) = mean_val;
end
end
            
%Reduce all the values such that after division to do not have values
%problems with the further analysis
p_names_LFQfin = p_names_LFQfin/(2^13);

%The normalization should be done on a replicate level!! However, if there
%is no actual link between the phosphoprotein and protein replicate it
%would be wrong to use a normalization based on replicate.

p_names_LFQfin(find(p_names_LFQfin==0)) = NaN;
p_names_LFQmean = [];
for i = 1:size(p_names_LFQfin,1)
    aux = [];
    for j = 1:nr_rep:size(p_names_LFQfin,2)
        aux = [aux,mean(p_names_LFQfin(i,j:j+nr_rep-1),'omitnan')];
    end
        p_names_LFQmean(i,:) = aux;
end

    
end
