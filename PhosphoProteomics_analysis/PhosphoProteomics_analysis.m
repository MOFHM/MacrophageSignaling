%__________________________________________________________________________
% Phosphoproteomics analyisis master file
% Last update: 26/07/2023
%__________________________________________________________________________
%The main scrip for running a PhosphoProteomics analysis
% - reads the MaxQuant files, compensate for the protein abundances using two different methods,
% extract only the most relevant measurements, transform and normalizes the
% data, create the quality control plots, do the differential analysis and
% make the upstream kinase analysis and network file generation
%__________________________________________________________________________
%
%

%
%__________________________________________________________________________
%


function [dir1,Perc_fals_poz,best_tr,test_res,boxcox_param,dist_names_group,parameter_mean_values_group,dist_performance_ind_group] = PhosphoProteomics_analysis(file_name_Phospho,protein_file,protein_logfc_compensation_file_name,dir_results,...
    Analysis_considering_NrPhosSites,Protein_abundance_normalization_moment,Intensity_column_of_interest,p_LFQ_columns,Protein_logfc_compensation_flag,statistical_test_used,Min_number_of_unique_replicates,...
    choose_Protein_or_Proteins,filter_by_evidence,filter_by_evidenceThreshold,normalization_method,Consider_unique_sepparately,...
    test_p_value,test_fc_value,nr_rep,nr_sample_groups,loc_prob_threshold,PEP_threshold,header_name_string,sample_names,imputation_method,imputation_remove,...
    imputation_DreamAI,distirbution_imputation_penalty,distirbution_imputation_normal,Meanshift_imputation_penalty,STDshift_imputation_penalty,...
    Meanshift_imputation_normal,STDshift_imputation_normal,PSP_file,PSP_Loc_prob,PSP_protein_compensation,PSP_p_val,Potential_contaminant_string,Rev_string,Loc_prob_string,PEP_string,...
    MS_MS_ID_string,Evidence_IDs_string,Positions_within_proteins_String,Amino_acid_String,Sequence_window_String,data_transformation,file_name_kinome,R_exe_location,R_library_location)


%Random number seed
rng(123)
%Colors used for the plotting of data
Colors = linspecer(12);
%Read the PhosphoProteomics data from MAXQuant by considering all the
%columns as Char arrays

if (~isfile(file_name_Phospho))
    error('The phosphoproteomics file does not exist')
end
if (~isfile(protein_file))
    error('The first proteomics file used for the phosphoresidues abundance compensation does not exist')
end
if (~isfile(protein_logfc_compensation_file_name))
    error('The second phosphoproteomics file does not exist')
end
if (~isfolder(dir_results))
    error('The results directory does not exist')
end
if (~isfile(PSP_file))
    error('The Curated database file does not exist')
end
if (~isfile(R_exe_location))
    error('The R executable does not exist')
end
if (~isfolder(R_library_location))
    error('The R library folder does not exist')
end


opts = detectImportOptions(file_name_Phospho);
opts = setvartype(opts,1:length(opts.VariableNames), 'char');
opts.PreserveVariableNames=true;
phospho = readtable(file_name_Phospho, opts);

%The saving directory created in conformity with the current date
mkdir(dir_results,datestr(today('datetime')))
dir1 = strcat(dir_results,'\',datestr(today('datetime')));
%The parameters that describes the experiment:
%nr_rep = number of sample replicates
%nr_sample_groups = number of analyzed phenotypes
%missing_values = the maximum number of imputed value per phenotype
%nr_analyzed_phos = maximum number of considered phosphorylated residues
%loc_prob_threshold = the threshold used to filter the MaxQuant output
%based on the localization probability
%PEP_threshold - the threshold used to filter the MaxQuant output
%based on the posterior probability
%int_data_columns - the columns associated with the measure intensitites
%for each phenotype per replciate and per the nubmer of phosphorylated
%residues
%header_name_string = a letter or string that defines the phenotype
%nomenclature
%nnn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%General parameters

if  (isempty(Analysis_considering_NrPhosSites) || ~isscalar(Analysis_considering_NrPhosSites) || (Analysis_considering_NrPhosSites<0 || Analysis_considering_NrPhosSites>2))
    warning('Analysis_considering_NrPhosSites was not correctly provided. Using the default value of 1, the analysis being done on the individual phosphoresidues intensities');
    Analysis_considering_NrPhosSites = 1;
end

if  (isempty(Intensity_column_of_interest))
    error('Please provide the  proper intensities columns')
end

if (Analysis_considering_NrPhosSites==1)
    int_data_columns = Intensity_column_of_interest;
    nr_analyzed_phos = 3;
    if(nr_analyzed_phos*nr_rep*length(sample_names) ~= length(Intensity_column_of_interest))
        error('The number of intensity columns does not match the number of replicates, samples and selected type of analysis(merged intensities or separated)')
    end
end
if (Analysis_considering_NrPhosSites==0)
    int_data_columns = Intensity_column_of_interest;
    nr_analyzed_phos = 1;
    if(nr_analyzed_phos*nr_rep*length(sample_names) ~= length(Intensity_column_of_interest))
        error('The number of intensity columns does not match the number of replicates, samples and selected type of analysis(merged intensities or separated)')
    end
end
if (Analysis_considering_NrPhosSites==2)
    ind = find(contains(phospho.Properties.VariableNames,"_p___1"));
    int_data_columns = ind;
    nr_analyzed_phos = 1;
    if(nr_analyzed_phos*nr_rep*length(sample_names) ~= length(Intensity_column_of_interest) || (ind~=Intensity_column_of_interest))
        error('The number of intensity columns does not match the number of replicates, samples and selected type of analysis(merged intensities or separated)')
    end
end

%Function parameters

if  (isempty(Protein_abundance_normalization_moment) || ~isscalar(Protein_abundance_normalization_moment) || (Protein_abundance_normalization_moment<0 || Protein_abundance_normalization_moment>1))
    warning('Protein_abundance_normalization_moment was not provided. Using the default value of -1, the phosphoresidues intensities will not be scaled to any exernal proteomics data');
    Protein_abundance_normalization_moment = -1;
end

if  (isempty(Protein_logfc_compensation_flag) || ~isscalar(Protein_logfc_compensation_flag) || (Protein_logfc_compensation_flag<1 || Protein_logfc_compensation_flag>2))
    warning('Protein_logfc_compensation_flag was not provided. Using the default value of 2, the phosphoresidues intensities will be compensated after the DEA analysis');
    Protein_logfc_compensation_flag = 2;
end

if  (isempty(statistical_test_used) || ~ismember(statistical_test_used,["moderate_t_test","t_test"]))
    warning('The DEA statistical test was not correctly provided. Using the default moderated t-test');
    statistical_test_used = "moderate_t_test";
end

if  (isempty(Min_number_of_unique_replicates) || (Min_number_of_unique_replicates<=0 || Min_number_of_unique_replicates>=1))
    warning('The Min_number_of_unique_replicates was not correctly provided. Using the default value of 0.5');
    Min_number_of_unique_replicates = 0.5;
end

if  (isempty(choose_Protein_or_Proteins) || ~ismember(choose_Protein_or_Proteins,["Protein","Leading proteins","Proteins"]))
    warning('The choose_Protein_or_Proteins was not correctly provided. Using the default column "Protein", hence not considering any homologs');
    choose_Protein_or_Proteins = "Protein";
end

if  (isempty(filter_by_evidence) ||  (filter_by_evidence<0 || filter_by_evidence>1))
    warning('The filter_by_evidence was not correctly provided. Using the default value of "0", not filtering the data based on the number of evidences');
    filter_by_evidence = 0;
end

if  (isempty(filter_by_evidenceThreshold) ||  (filter_by_evidenceThreshold<0))
    warning('The filter_by_evidenceThreshold was not correctly provided. Using the default value of "2", and use it only if filter_by_evidence was set to 1');
    filter_by_evidenceThreshold=2;
end

if  (isempty(normalization_method) ||  (normalization_method<0 || normalization_method>3))
    warning('The normalization_method was not correctly provided. Using the default value of "1", performing quantile normalization');
    normalization_method = 1;
end

if  (isempty(Consider_unique_sepparately) ||  (Consider_unique_sepparately<0 || Consider_unique_sepparately>1))
    warning('The Consider_unique_sepparately was not correctly provided. Using the default value of "0", the unique cases being imputed and not treated separately');
    Consider_unique_sepparately = 0;
end

if  (isempty(test_p_value) ||  (test_p_value<0 || test_p_value>1))
    warning('The test_p_value was not correctly provided. Using the default DEA statistical confidence of 95%');
    test_p_value = 0.05;
end

if  (isempty(test_fc_value) ||  (test_fc_value<0))
    warning('The test_fc_value was not correctly provided. Using the default DEA log2FoldChange threshold of 1');
    test_fc_value = 1;
end


if  (isempty(nr_rep) ||  (nr_rep<0))
    error('The nr_rep was not correctly provided. Please provide the number of replicates');
end

if  (isempty(nr_sample_groups) ||  (nr_sample_groups<0))
    error('The nr_sample_groups was not correctly provided. Please provide the number of sample groups');
end

if  (isempty(loc_prob_threshold) ||   (test_p_value<0 || test_p_value>1))
    warning('The loc_prob_threshold was not correctly provided. Using the default localization probablity of 0');
    loc_prob_threshold = 0;
end

if  (isempty(PEP_threshold) ||   (PEP_threshold<0 || PEP_threshold>1))
    warning('The PEP_threshold was not correctly provided. Using the default posterior error probability threshold of 1');
    PEP_threshold = 1;
end

if  (isempty(PEP_threshold) ||   (PEP_threshold<0 || PEP_threshold>1))
    warning('The PEP_threshold was not correctly provided. Using the default posterior error probability threshold of 1');
    PEP_threshold = 1;
end

if  (isempty(header_name_string) ||   isempty(find(contains(phospho.Properties.VariableNames,header_name_string))))
    error('The header_name_string was not correctly provided. Please provide a valid identifier for all the samples names codining in the phosphotable header');
end

if  (isempty(sample_names))
    error('Please provide the sample names');
end

if  (isempty(PEP_threshold) ||   (PEP_threshold<0 || PEP_threshold>1))
    warning('The PEP_threshold was not correctly provided. Using the default posterior error probability threshold of 1');
    PEP_threshold = 1;
end


% Imputation analysis characteristics:

if  (isempty(imputation_method) ||   (imputation_method<0 || imputation_method>1))
    warning('The imputation_method was not correctly provided. Using the default imputation method based on DreamAI and random draw from a normal distribution');
    imputation_method = 1;
end

if  (isempty(imputation_remove) ||   (min(imputation_remove) <-1 || max(imputation_remove)>nr_rep))
    warning('The imputation_remove was not correctly provided. Do not remove any comparison conditions');
    imputation_remove = [-1];
end

if  (isempty(imputation_DreamAI) ||   (min(imputation_DreamAI) <-1 || max(imputation_DreamAI)>nr_rep))
    warning('The imputation_DreamAI was not correctly provided. Do not impute based on DreamAI method any condition');
    imputation_DreamAI = [-1];
end

if  (isempty(distirbution_imputation_penalty) ||   (min(distirbution_imputation_penalty) <-1 || max(distirbution_imputation_penalty)>nr_rep))
    warning('The distirbution_imputation_penalty was not correctly provided. Using the default imputation penalty for the cases when there are nr_rep-1 and nr_rep missing values');
    distirbution_imputation_penalty = [3,4];
end

if  (isempty(distirbution_imputation_normal) ||   (min(distirbution_imputation_normal) <-1 || max(distirbution_imputation_normal)>nr_rep))
    warning('The distirbution_imputation_normal was not correctly provided. Using the default distirbution_imputation_normal for the cases that have from 0 to nr_rep-2 missing values');
    distirbution_imputation_normal = [0,1,2];
end

if  (isempty(Meanshift_imputation_penalty) ||   (Meanshift_imputation_penalty <0))
    warning('The Meanshift_imputation_penalty was not correctly provided. Using the default mean shift for the penalty cases value of 1.8');
    Meanshift_imputation_penalty = 1.8;
end

if  (isempty(STDshift_imputation_penalty) ||   (STDshift_imputation_penalty <0))
    warning('The STDshift_imputation_penalty was not correctly provided. Using the default std for the penalty distributtion cases of 0.3');
    STDshift_imputation_penalty = 0.3;
end

if  (isempty(Meanshift_imputation_normal) ||   (Meanshift_imputation_normal <0))
    warning('The Meanshift_imputation_normal was not correctly provided. Using the default std for normaly imputed cases value of 0.5');
    Meanshift_imputation_normal = 0.5;
end

if  (isempty(STDshift_imputation_normal) ||   (STDshift_imputation_normal <0))
    warning('The STDshift_imputation_normal was not correctly provided. Using the default std for the normaly imputed distributtion cases of 0.3');
    STDshift_imputation_normal = 0.3;
end


%Curated Data analysis characteristics:

if  (isempty(PSP_Loc_prob) ||  (PSP_Loc_prob<0 || PSP_Loc_prob>1))
    warning('The PSP_Loc_prob was not correctly provided. Using the default value of "1", hence filter the data for the upstream kinase assessment based on the localization probability');
    PSP_Loc_prob = 1;
end

if  (isempty(PSP_protein_compensation) ||  (PSP_protein_compensation<0 || PSP_protein_compensation>1))
    warning('The PSP_protein_compensation was not correctly provided. Using the default value of "1", hence compensate the phosphoresidues abundances with respect to the proteomics abundances during the upstream kinases assessment');
    PSP_protein_compensation = 1;
end


if  (isempty(PSP_p_val) ||  (PSP_p_val<0 || PSP_p_val>1))
    warning('The PSP_p_val was not correctly provided. Using the default value of "0.01" representing a confidence level of 95% for the upstream kinase assessmen algorithm');
    PSP_p_val = 0.01;
end

if  (isempty(p_LFQ_columns))
    warning('The p_LFQ_columns was not correctly provided. Using the default value of 106:117');
    p_LFQ_columns = 106:117;
end

if  (ischar(R_exe_location)==0)
    warning('R_exe_location was string. Converting to char');
    R_exe_location = convertStringsToChars(R_exe_location);
end

if  (ischar( R_library_location)==0)
    warning(' R_library_location was string. Converting to char');
     R_library_location = convertStringsToChars( R_library_location);
end

missing_values_selfImputation = 3;

%make_PSP_protein_compensation_on_Protein_level = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Remove contaminants - first find the corresponding name
ind = find(contains(lower(phospho.Properties.VariableNames),lower(Potential_contaminant_string)));
if(isempty(ind)) error('The column "Potential contaminant" does not exist. Introduce it, or contact the authors'); end
a=phospho{:,ind};
%Mark for contaminants/reverse etc
mark = a(find(~cellfun(@isempty,a),1));
ind = [];
for i = 1:length(a)
    if (isequal(a{i},mark{:})) ind = [ind;i]; end
end
phospho_cont = phospho;
phospho_cont(ind,:)=[];

%Remove reverse - first find the corresponding name
ind = find(contains(lower(phospho.Properties.VariableNames),lower(Rev_string)));
if(isempty(ind)) error('The column "Reverse" does not exist. Introduce it, or contact the authors'); end
a=phospho_cont{:,ind};
%Mark for contaminants/reverse etc
mark = a(find(~cellfun(@isempty,a),1));
ind = [];
for i = 1:length(a)
    if (isequal(a{i},mark{:})) ind = [ind;i]; end
end
phospho_cont_rev = phospho_cont;
phospho_cont_rev(ind,:)=[];

%Filter by the localization probability - first find the corresponding name
ind = find(ismember(lower(phospho.Properties.VariableNames),lower(Loc_prob_string)));
if(isempty(ind)) error('The column "Localization prob" does not exist. Introduce it, or contact the authors'); end
a=phospho_cont_rev{:,ind};
ind = [];
%Choose the treshold as the defined value
%nnn
for i = 2:length(a)
    if (str2num(a{i})<loc_prob_threshold) ind = [ind;i]; end
end
phospho_cont_rev_loc = phospho_cont_rev;
phospho_cont_rev_loc(ind,:)=[];

%Filter by the PEP (posterior error probability) - first find the corresponding name
ind = find(ismember(lower(phospho.Properties.VariableNames),lower(PEP_string)));
if(isempty(ind)) error('The column "PEP" does not exist. Introduce it, or contact the authors'); end
a=phospho_cont_rev_loc{:,ind};
ind = [];
for i = 1:length(a)
    if (str2num(a{i})>PEP_threshold) ind = [ind;i]; end
end
phospho_cont_rev_loc_PEP = phospho_cont_rev_loc;
phospho_cont_rev_loc_PEP(ind,:)=[];
%Fitlere out possible remaining contaminants
ind = find(ismember(lower(phospho.Properties.VariableNames),lower('Leading proteins')));
a=phospho_cont_rev_loc_PEP{:,ind};
ind = find(contains(a,"Cont"));
if(isempty(ind)==0)
    phospho_cont_rev_loc_PEP(ind,:)=[];
end
%Percentage of fals pozitives - due to the contaminants and too light
%restrictions
Perc_fals_poz = abs((length(phospho_cont_rev_loc_PEP{:,1})-length(phospho{:,1}))/length(phospho{:,1}));

%%%%%nn
%Remove the Phosphoresidues that come from a number of spectra less that
%the chosen threshold
if(filter_by_evidence==1)
    evi_id = find(contains(lower(phospho.Properties.VariableNames),lower(Evidence_IDs_string)));
    if(isempty(evi_id)) error('The column "Evidence IDs" does not exist. Introduce it, or contact the authors'); end
    aux = phospho_cont_rev_loc_PEP{:,evi_id};
    aux = cellfun(@(x) split(x,";"),aux,'UniformOutput',false);
    ind = cellfun(@(x) size(x,1),aux);
    ind = find(ind<filter_by_evidenceThreshold);
    phospho_cont_rev_loc_PEP(ind,:) = [];
    msms_id = find(contains(lower(phospho.Properties.VariableNames),lower(MS_MS_ID_string)));
    if(isempty(msms_id)) error('The column "MS/MS IDs" does not exist. Introduce it, or contact the authors'); end
    aux = phospho_cont_rev_loc_PEP{:,msms_id};
    aux = cellfun(@(x) split(x,";"),aux,'UniformOutput',false);
    ind = cellfun(@(x) size(x,1),aux);
    ind = find(ind<filter_by_evidenceThreshold);
    phospho_cont_rev_loc_PEP(ind,:) = [];
end
%%%%%



%%%%nnn
%Names of the proteins - it is correct to take the leading protein or the
%first one of them
ind = find(ismember(lower(phospho.Properties.VariableNames),lower('Proteins')));
if(isempty(ind)) error('The column "Proteins" does not exist. Introduce it, or contact the authors. Error:378'); end
p_names = phospho_cont_rev_loc_PEP{:,ind};
%Position of phosphorylated residue
ind_residue_number = find(contains(lower(phospho.Properties.VariableNames),lower(Positions_within_proteins_String)));
if(isempty(ind_residue_number)) error('The column "Positions within proteins" does not exist. Introduce it, or contact the authors'); end
ind_residue_name = find(contains(lower(phospho.Properties.VariableNames),lower(Amino_acid_String)));
if(isempty(ind_residue_name)) error('The column "Amino acid" does not exist. Introduce it, or contact the authors'); end
position_protein = cellfun(@strcat,phospho_cont_rev_loc_PEP{:,ind_residue_name},strrep(cellstr(phospho_cont_rev_loc_PEP{:,ind_residue_number}),' ',''),'UniformOutput',false);
position_protein = cellfun(@(x) insertAfter(x,';',x(1)),position_protein,'UniformOutput',false);
ind_sequence = find(contains(lower(phospho.Properties.VariableNames),lower(Sequence_window_String)));
if(isempty(ind_sequence)) error('The column "Sequence window" does not exist. Introduce it, or contact the authors'); end
peptide_sequence = phospho_cont_rev_loc_PEP{:,ind_sequence};

ind_locProb = find(ismember(lower(phospho.Properties.VariableNames),lower(Loc_prob_string)));
loc_prob = phospho_cont_rev_loc_PEP{:,ind_locProb};


if(choose_Protein_or_Proteins=="Protein")
    ind_aux = find(ismember(lower(phospho.Properties.VariableNames),lower('Protein')));
    if(isempty(ind_sequence)) error('The column "Protein" does not exist. Introduce it, or contact the authors'); end
end
if(choose_Protein_or_Proteins=="Leading proteins")
    ind_aux = find(ismember(lower(phospho.Properties.VariableNames),lower('Leading proteins')));
    if(isempty(ind_sequence)) error('The column "Leading proteins" does not exist. Introduce it, or contact the authors'); end
end
if(exist('ind_aux','var')==1)
    p_names_aux = phospho_cont_rev_loc_PEP{:,ind_aux};
    a1 = cellfun(@(x) strsplit(x,';'),p_names,'UniformOutput',false);
    a2 = cellfun(@(x) strsplit(x,';'),p_names_aux,'UniformOutput',false);
    a3 = cellfun(@(x) strsplit(x,';'),position_protein,'UniformOutput',false);
    for i = 1:length(a1)
        ind_pos = find(ismember(a1{i},a2{i}));
        a3{i} = a3{i}(ind_pos);
    end
    as = cellfun(@(x) horzcat(x{:}),a3,'UniformOutput',false);
    as1 = cellfun(@(x) [x(1),insertBefore(x(2:end),x(1),';')],as,'UniformOutput',false);
    position_protein = as1;
    p_names = p_names_aux;
    a1 = cellfun(@(x) strsplit(x,';'),peptide_sequence,'UniformOutput',false);
    peptide_sequence = cellfun(@(x) x(1),a1);
end
%%%%%%%%%%%%%%%%


%Analyzed the data of interest
col_names = phospho.Properties.VariableNames;
int_data = phospho_cont_rev_loc_PEP{:,int_data_columns};
int_data = cellfun(@str2num,int_data);


%Sort and arrange the header such that to have all the replicated for all
%the phenotypes sorted by the number of phsphorylated values (ex. for the
%default values the first 12 columns are associated with one phosphorylated
%site, the next 12 with the following and so on.
header = col_names(int_data_columns);
for i = 1:length(header)
    ind = strfind(header{i},header_name_string);
    header{i} = header{i}(ind:end);
end


[~,ind] = sort(header);
header = header(ind);
col_names = col_names(ind);
int_data = int_data(:,ind);
if(nr_analyzed_phos==3)
    %Reorder the data such that to have _1 for all the groups followed by _2
    header1 = cell(size(header));
    col_names1 = cell(size(col_names));
    int_data1 = zeros(size(int_data));
    for i = 1:3:length(header)
        header1{(i-1)/3+1} = header{i};
        col_names1{(i-1)/3+1} = col_names{i};
        int_data1(:,(i-1)/3+1) = int_data(:,i);
    end
    for i = 2:3:length(header)
        header1{(i-2)/3+1+length(header)/3} = header{i};
        col_names1{(i-2)/3+1+length(header)/3} = col_names{i};
        int_data1(:,(i-2)/3+1+length(header)/3) = int_data(:,i);
    end
    for i = 3:3:length(header)
        header1{(i-3)/3+1+2*length(header)/3} = header{i};
        col_names1{(i-3)/3+1+2*length(header)/3} = col_names{i};
        int_data1(:,(i-3)/3+1+2*length(header)/3) = int_data(:,i);
    end
    
    int_data = int_data1;
    header = header1;
    col_names = col_names1;
end


%Protein abundance normalization. Scale the phosphopeptides based on the
%measrued abudnacnes of the proteins from another dataset or experiment.
if(Protein_abundance_normalization_moment==0)
    [p_names_prot,p_names_LFQmean] = normalize_abundance_v1(protein_file,nr_rep,nr_sample_groups,p_LFQ_columns,header_name_string);
    p_names_LFQmean = p_names_LFQmean/floor(log2(min(p_names_LFQmean(:))));
    for i = 1:size(int_data,1)
        names_protein_phospho = split(p_names(i),';');
        for j = 1:nr_rep:size(int_data,2)
            div_val = [];
            ff = fix(((j-1)/nr_rep+1)/nr_rep);
            ind_j = mod((j-1)/nr_rep+1,nr_rep)+ff;
            ind_j(find(ind_j==nr_rep))=1;
            for k = 1:length(names_protein_phospho)
                aux = find(strcmp(p_names_prot,names_protein_phospho(k))==1,1);
                if(isempty(aux)==0) div_val = [div_val,p_names_LFQmean(aux,ind_j)]; end
            end
            if(isempty(div_val)) div_val = median(p_names_LFQmean(:,:)); end
            div_val_mean = mean(div_val);
            int_data(i,j:j+nr_rep-1) = int_data(i,j:j+nr_rep-1)./div_val_mean;
        end
    end
end

%Cut the rows that do not have more than the required number of replicates in at least one condition
%considering all the phosphates
ind_cut = cell(1,nr_analyzed_phos);
for i = 1:size(int_data,1)
    for jj = 1:nr_analyzed_phos
        ind = [];
        ss = (jj-1)*nr_sample_groups*nr_rep+1;
        for j = ss:nr_rep:jj*nr_sample_groups*nr_rep
            vec = int_data(i,j:j+nr_rep-1);
            ind((j-ss)/nr_rep+1) = length(find(vec==0));
        end
        nn = find(ind<=(1-Min_number_of_unique_replicates)*nr_rep);
        if (isempty(nn)) ind_cut{1,jj} = [ind_cut{1,jj},i]; end
    end
end
if(size(ind_cut,2)>1)
    iind = intersect(ind_cut{1},ind_cut{2});
    for i = 3:length(ind_cut)
        iind = intersect(iind,ind_cut{i});
    end
else iind = ind_cut{1,1};
end
int_data(iind,:) = [];
p_names(iind,:) = [];
position_protein(iind,:) = [];
peptide_sequence(iind,:) = [];
loc_prob(iind,:) = [];

if (Analysis_considering_NrPhosSites==1)
    %Replace the semi-rows that do not have more than two replicates in one condition
    %considering one phosphate at a time (replace the measured values with
    %zeros)
    for i = 1:size(int_data,1)
        for jj = 1:nr_analyzed_phos
            ind = [];
            ss = (jj-1)*nr_analyzed_phos*nr_rep+1;
            for j = ss:nr_rep:jj*nr_analyzed_phos*nr_rep
                vec = int_data(i,j:j+nr_rep-1);
                ind((j-ss)/nr_rep+1) = length(find(vec==0));
            end
            nn = find(ind<=(1-Min_number_of_unique_replicates)*nr_rep);
            if (isempty(nn)) int_data(i,ss:jj*nr_analyzed_phos*nr_rep) = 0; end
        end
    end
end
%List the initial identified proteins before further analysis
initial_identified_unique_proteins = unique(extractBetween([p_names{:}],'|','|'));



%Best transformation finding
best_tr = {};
test_res = {};
boxcox_param = {};

%For each individual sample that was measured
for i = 1:size(int_data,2)
    if(length(find(int_data(:,i)>0))>10)
        [best_tr{i},test_res{i},boxcox_param{i}] = find_best_transform(int_data(:,i));
    end
end
%Input data transfromation. Log2 transform very good
try
    int_data = feval(data_transformation,int_data);
catch
    error('The "data_transfromation" string that was provided is invalig. You can use "log2" as an input data transfroamtion')
end

%Inf and NaN values handling
[aa,bb] = find(int_data==-Inf);
for i = 1:length(aa)
    int_data(aa(i),bb(i)) = 0;
end


%Data imputation
[dist_names_group,parameter_mean_values_group,dist_performance_ind_group] = find_dist(int_data,nr_rep);

%Impute data based on the selected distribtuion
% Normal
% birnbaussander
%nnn
distribution = 'normal';
position_of_mean_in_param = 1;
position_of_std_in_param = 2;
param = {};

for i = 1:length(parameter_mean_values_group)
    nn = dist_names_group{i};
    ind = find(strcmp(nn,distribution));
    param{i} = parameter_mean_values_group{i}(ind,:);
end
%Attention! Change Normal to other type of desired distribution
%nnn
distribution = 'Normal';
%Penalty factor is the proportion of the mean difference that is
%substracted from the imputed values when there a three such imputation per
%a condition - between 0.3 and 1 should be reasonable
penalty_factor = 0.5;
[data_fin,nr_missing_val] = impute_data2(int_data,nr_rep,nr_sample_groups,penalty_factor,distribution,param,missing_values_selfImputation,position_of_mean_in_param,position_of_std_in_param);

%Merged the data
[merged_data] = merge_data(data_fin,missing_values_selfImputation,nr_missing_val,nr_rep);


%Extract the lines that have a peptide completely missing in other
%conditions - take care to the difference between multiple phosphorilation
%sites - extract real unique peptides - it was chosen to consider different
%number of phosphorylation sites as different conditions - thus it is
%possible that the same peptide to be unique in two different phenotypes
%because it is unique with regard to the number of phosphorylated sites
used_data = merged_data{1,3};

if (Consider_unique_sepparately==1)
    ind_extracted_rows = cell(nr_analyzed_phos,1);
    ind_extracted_rows_p123 = cell(nr_analyzed_phos,1);
    for i = 1:size(used_data,1)
        for j = 1:nr_rep*nr_sample_groups:size(used_data,2)
            ind = [];
            for k = j:nr_rep:j+nr_rep*nr_sample_groups-1
                vec = used_data(i,k:k+nr_rep-1);
                ind = [ind,sum(vec==0)];
            end
            ind1 = find(ind==nr_rep);
            %>=1 for considering unique values the ones that are found missing in at
            %least one condition or >1 that are found missing in at least two
            %conditions and are measured in only one
            if (length(ind1)>1&&length(ind1)<3) ind_extracted_rows{(j-1)/(nr_rep*nr_sample_groups)+1} = [ind_extracted_rows{(j-1)/(nr_rep*nr_sample_groups)+1},i]; ind_extracted_rows_p123{(j-1)/(nr_rep*nr_sample_groups)+1} = [ind_extracted_rows_p123{(j-1)/(nr_rep*nr_sample_groups)+1},k];
                % else if (length(ind1)<2) break;
                %     end
            end
        end
    end
end

%Split the data into the corresponding groups
final_data_p = cell(1,nr_analyzed_phos);
p_names_p = cell(1,nr_analyzed_phos);
header_p = cell(1,nr_analyzed_phos);
position_protein_p = cell(1,nr_analyzed_phos);
final_data_p_initial = cell(1,nr_analyzed_phos);
peptide_sequence_p = cell(1,nr_analyzed_phos);
loc_prob_p = cell(1,nr_analyzed_phos);
for i = 1:nr_rep*nr_sample_groups:size(used_data,2)
    inn = (i-1)/(nr_rep*nr_sample_groups)+1;
    final_data_p{inn} = used_data(:,i:i+nr_rep*nr_sample_groups-1);
    final_data_p_initial{inn} = int_data(:,i:i+nr_rep*nr_sample_groups-1);
    p_names_p{inn} = p_names;
    peptide_sequence_p{inn} = peptide_sequence;
    loc_prob_p{inn} = loc_prob;
    header_p{inn} =  header(i:i+nr_rep*nr_sample_groups-1);
    position_protein_p{inn} = position_protein;
    if(Consider_unique_sepparately==1)
        p_names_p{inn}(ind_extracted_rows{inn}) = [];
        final_data_p{inn}(ind_extracted_rows{inn},:) =  [];
        final_data_p_initial{inn}(ind_extracted_rows{inn},:) =  [];
        position_protein_p{inn}(ind_extracted_rows{inn}) = [];
        peptide_sequence_p{inn}(ind_extracted_rows{inn}) = [];
        loc_prob_p{inn}(ind_extracted_rows{inn}) = [];
    end
end

data_unique = {};
data_unique_initial = {};
p_names_unique = {};
position_protein_unique = {};
baseMean_unique = {};
peptide_sequence_unique = {};
measured_values_unique = {};
loc_prob_unique = {};

if (Consider_unique_sepparately==1)
    %%%%%Here is modified to indentify the unique elements
    for i = 1:nr_rep*nr_sample_groups:size(used_data,2)
        inn = (i-1)/(nr_rep*nr_sample_groups)+1;
        if(Min_number_of_unique_replicates>0.5)
            index_bun = cellfun(@(x) find(length(find(int_data(x,i:i+nr_rep*nr_sample_groups-1)>0))>=nr_rep*Min_number_of_unique_replicates),num2cell(ind_extracted_rows{inn}),'UniformOutput',false);
            index_bun = find(cellfun(@(x) ~isempty(x),index_bun));
            index_bun = ind_extracted_rows{inn}(index_bun);
        else  index_bun = ind_extracted_rows{inn};
        end
        data_unique{inn} = used_data(index_bun,i:i+nr_rep*nr_sample_groups-1);
        data_unique_initial{inn} = int_data(index_bun,i:i+nr_rep*nr_sample_groups-1);
        p_names_unique{inn} = p_names(index_bun);
        position_protein_unique{inn} = position_protein(index_bun);
        baseMean_unique{inn} = mean(nr_sample_groups.*data_unique{inn},2);
        peptide_sequence_unique{inn} = peptide_sequence(index_bun);
        loc_prob_unique{inn} = loc_prob(index_bun);
    end
    measured_values_unique = cell(size(data_unique_initial,2),1);
    for j = 1:size(data_unique_initial,2)
        ds = data_unique_initial{j};
        measured_values_unique{j} = cellfun(@(x) length(find(x~=0)),num2cell(ds,2));
    end
    
end

%Remove the rows with zeros
for i = 1:size(final_data_p,2)
    ind = [];
    for j = 1:size(final_data_p{i},1)
        vec = final_data_p{i}(j,:);
        if (sum(vec==0)==nr_rep*nr_sample_groups) ind = [ind,j]; end
    end
    final_data_p{i}(ind,:) = [];
    p_names_p{i}(ind) = [];
    position_protein_p{i}(ind) = [];
    final_data_p_initial{i}(ind,:) = [];
    peptide_sequence_p{i}(ind,:) = [];
    loc_prob_p{i}(ind,:) = [];
end

%Chose the normalization method
if(normalization_method==0)
    test_name = 'median';
else if (normalization_method==1)
        test_name = 'qnorm';
    else if (normalization_method==2)
            test_name = 'zscore';
        else if (normalization_method==3)
                test_name = 'mean';
            end
        end
    end
end

%Data normalization
final_data_p_transform = {};
for i = 1:size(final_data_p,2)
    aa1 = final_data_p{i}(:,:);
    ind = find(aa1==0);
    aa1(ind) = NaN;
    
    if(isequal(test_name,'qnorm'))
        final_data_p_transform{i} = quantilenorm(aa1);
    end
    if(isequal(test_name,'zscore'))
        final_data_p_transform{i} = mean(std(aa1,1,1,'omitnan'),'omitnan').*normalize(aa1,1);
    end
    if(isequal(test_name,'median'))
        final_data_p_transform{i} = aa1 - repmat(median(aa1,1,'omitnan'),size(aa1,1),1);
    end
    if(isequal(test_name,'mean'))
        final_data_p_transform{i} = aa1 - repmat(mean(aa1,1,'omitnan'),size(aa1,1),1);
    end
    
end

measured_values = cell(size(final_data_p_initial,2),nchoosek(nr_sample_groups,2));
for j = 1:size(final_data_p_initial,2)
    ind_ms = 0;
    for i = 1:nr_rep:size(final_data_p_initial{j},2)-nr_rep
        for k = i+nr_rep:nr_rep:size(final_data_p_initial{j},2)
            close all
            ind_i = (i-1)/nr_rep+1;
            ind_k = (k-i)/nr_rep+1;
            ind_ms = ind_ms+1;
            ds = final_data_p_initial{j}(:,[i:i+nr_rep-1 k:k+nr_rep-1]);
            v1 = cellfun(@(x) length(find(x~=0)),num2cell(ds(:,1:nr_rep),2));
            v2 = cellfun(@(x) length(find(x~=0)),num2cell(ds(:,nr_rep+1:2*nr_rep),2));
            v = cellfun(@(x,y) strcat(x,'||',y),string(v1),string(v2),'UniformOutput',false);
            measured_values{j,ind_ms} = v;
        end
    end
end


nr_miss_total_p = cell(nr_rep+1,2,size(int_data,2)/nr_rep);
nr_miss_total_p(:,1,:) = {0};

for kk = 1:nr_rep+1
    for i = 1:nr_rep:size(int_data,2)
        for j = 1:size(int_data,1)
            vec = int_data(j,i:i+nr_rep-1);
            [vec_sort,initial_ind] = sort(vec,'descend');
            ind = find(vec_sort==0);
            if(isempty(ind)) ind = []; end
            if(length(ind) == kk-1)
                nr_miss_total_p{kk,1,(i-1)/nr_rep+1} = nr_miss_total_p{kk,1,(i-1)/nr_rep+1}+1;
                nr_miss_total_p{kk,2,(i-1)/nr_rep+1} = [nr_miss_total_p{kk,2,(i-1)/nr_rep+1},j];
            end
        end
    end
end

%Plot the data - extract significant peptides
mkdir(dir1,'Images_generated')
dir_current = strcat(dir1,'\Images_generated');

%Create all the plots used for the quality control
mkdir(dir_current,'Quality_control')
dir_current1 = strcat(dir_current,'\Quality_control');
plot_generation_validation(final_data_p,final_data_p_transform,int_data,header,header_p,p_names_p,position_protein_p,nr_rep,nr_miss_total_p,Colors,sample_names,nr_sample_groups,dir_current1)

%Extract the differentialy expresed peptides
mkdir(dir_current,'Differentialy_expressed_peptides')
dir_current2 = strcat(dir_current,'\Differentialy_expressed_peptides');

Imputation_script;
[list_proteins_significant_non_unique,list_proteins_significant_non_unique_upreg,list_proteins_significant_non_unique_downreg,list_proteins_significant_unique,list_proteins_significant_unique_complete,list_proteins_significant_non_unique_upreg_complete,list_proteins_significant_non_unique_downreg_complete] = list_generation_proteins(nr_sample_groups,nr_analyzed_phos,data_unique,nr_rep,p_names_unique,p_names_p,position_protein_p,position_protein_unique,fc,p_values_total,baseMean,baseMean_unique,significant_peptides_index,significant_peptides_index_upreg,significant_peptides_index_downreg,Consider_unique_sepparately,p_values_noncorrected,measured_values,measured_values_unique,loc_prob_p,loc_prob_unique,peptide_sequence_p,peptide_sequence_unique,stder_total,df_pooled_total,prot_comp_flag);
DifferentialAnalysis_FileGeneration;

%Coefficient of variation
%method == 1 for std and 0 for std/mean
[CV_normal] = CV_comp(cellfun(@(x)2.^x,final_data_p,'UniformOutput',false),nr_rep,Colors,'non-transformed',0,sample_names,dir_current2);
[CV_transform] = CV_comp(final_data_p_transform,nr_rep,Colors,test_name,1,sample_names,dir_current2);

data_full = cell(1,3);
fac = nr_rep*nr_sample_groups;
for ii = 1:fac:size(int_data,2)
    dist_indices = [];
    dist_parameters = {};
    for j = 1:size(int_data,1)
        vec = int_data(j,ii:ii+fac-1);
        if (isempty(find(vec==0)))
            data_full{(ii-1)/fac+1} = [data_full{(ii-1)/fac+1};vec];
        end
    end
end

[CV_normal_full] = CV_comp(cellfun(@(x)2.^x,data_full,'UniformOutput',false),nr_rep,Colors,' unique non-transformed',0,sample_names,dir_current2);



%Extract full data - clustergram on fully measured data - do the
%differntial analysis on this data as well
final_data_p_transformed_full_measured = {};
for j = 1:nr_analyzed_phos-1
    g1_ind = [nr_miss_total_p{1,2,(j-1)*nr_sample_groups+1},nr_miss_total_p{2,2,(j-1)*nr_sample_groups+1}];
    g2_ind = [nr_miss_total_p{1,2,(j-1)*nr_sample_groups+2},nr_miss_total_p{2,2,(j-1)*nr_sample_groups+2}];
    g3_ind = [nr_miss_total_p{1,2,(j-1)*nr_sample_groups+3},nr_miss_total_p{2,2,(j-1)*nr_sample_groups+3}];
    int = intersect(g1_ind,g2_ind);
    int = intersect(int,g3_ind);
    gr_comp = int_data(int,:);
    gr_comp = int_data(int,(j-1)*12+1:j*12);
    gr_comp(find(gr_comp==0)) = NaN;
    for k = 1:4:size(gr_comp,2)
        gr_comp(:,k:k+3) = knnimpute(gr_comp(:,k:k+3));
    end
    p_comp{j} = p_names(int);
    pos_comp{j} = position_protein(int);
    if(isequal(test_name,'qnorm'))
        final_data_p_transformed_full_measured{j} = quantilenorm(gr_comp);
    end
    if(isequal(test_name,'zscore'))
        final_data_p_transformed_full_measured{j} = mean(std(gr_comp,1,1,'omitnan'),'omitnan').*normalize(gr_comp,1);
    end
    if(isequal(test_name,'median'))
        final_data_p_transformed_full_measured{j} = gr_comp - repmat(median(gr_comp,1,'omitnan'),size(gr_comp,1),1);
    end
    if(isequal(test_name,'mean'))
        final_data_p_transformed_full_measured{j} = gr_comp - repmat(mean(gr_comp,1,'omitnan'),size(gr_comp,1),1);
    end
end

mkdir(dir1,'Images_generated_full_measured')
dir_current = strcat(dir1,'\Images_generated_full_measured');
%Change qnorm with zscore!!!!!!!!!

[significant_peptides_index_full_measured,significant_peptides_index_upreg_full_measured,significant_peptides_index_downreg_full_measured,fc_full_measured] =  plot_generation2_fullMeasure(final_data_p,final_data_p_transformed_full_measured,int_data,header,header_p,p_comp,pos_comp,nr_rep,nr_miss_total_p,Colors,sample_names,nr_sample_groups,dir_current);


%[list_kinase_proteinsBackground,score_kin_subBackground] = getPhosphoSiteData(p_names_PROTEOMICS,position_protein_PROTEOMICS);
if (PSP_Loc_prob==1)
    index_to_remove = find(cellfun(@str2num,loc_prob)<0.75);
    p_names_background = p_names; p_names_background(index_to_remove) = [];
    position_protein_background = position_protein; position_protein_background(index_to_remove) = [];
else p_names_background = p_names; position_protein_background = position_protein;
end
[list_kinase_proteinsBackground,score_kin_subBackground] = getPhosphoSiteData(p_names_background,position_protein_background,PSP_file);
PhosphoSiteDatabase_Analysis_FileGeneration;


ind_full_phenotype = cell(1,3);
%p_full_phenotype = cell(1,3);
fac = nr_rep*nr_sample_groups;
for i = 1:nr_sample_groups
    for ii = (i-1)*nr_rep+1:fac:size(int_data,2)
        for j = 1:size(int_data,1)
            vec = int_data(j,ii:ii+nr_rep-1);
            if (length(find(vec==0))<=3)
                %There will be many replicate values because some are measured for
                %multiple phosphates change p_nmeas with j and then extract unique
                %elements of them in order to get an easy measuremet of the number
                %of identified phosphopeptides
                ind_full_phenotype{i} = [ind_full_phenotype{i};j];
                %p_full_phenotype{i} = [p_full_phenotype{i};p_names(j)];
            end
        end
    end
end
ind_full_phenotype = cellfun(@unique,ind_full_phenotype,'UniformOutput',false);%This may be commented when p_full_phenotype{i} = [p_full_phenotype{i};p_names(j)];
p_full_phenotype = cellfun(@(x) [p_names(x),position_protein(x)],ind_full_phenotype,'UniformOutput',false);%This may be commented when p_full_phenotype{i} = [p_full_phenotype{i};p_names(j)];
%p_full_phenotype = cellfun(@unique,p_full_phenotype,'UniformOutput',false); Uncommet when p_full_phenotype{i} = [p_full_phenotype{i};p_names(j)];

pat_ini_protein = 'sp|';
pat_fin_protein = '|';
aa1 = cell(size(p_full_phenotype));
aa2 = cell(size(p_full_phenotype));
for i = 1:length(p_full_phenotype)
    for ii = 1:size(p_full_phenotype{i},1)
        if(contains(p_full_phenotype{i}(ii),"cont",'IgnoreCase',true)==0)
            if (isempty(find(contains(p_full_phenotype{i},'sp'),1)) == 0)
                aa1{i} = [aa1{i};extractBetween(p_full_phenotype{i}{ii,1},pat_ini_protein,pat_fin_protein)];
            else
                aa1{i} = [aa1{i};p_full_phenotype{i}(ii,1)];
            end
            a = strsplit(p_full_phenotype{i}{ii,2},';');
            aa2{i} = [aa2{i};a'];
        end
    end
end
p_full_phenotype = cellfun(@(x,y) [x,y],aa1,aa2,'UniformOutput',false);

kinase_original_list = cellfun(@(x) get_list_kinases_origianl_data(x(:,1)),p_full_phenotype,'UniformOutput',false);
kinase_original_list = cellfun(@(x) sortrows(x,1),kinase_original_list,'UniformOutput',false);
aa1 = cell(size(p_full_phenotype));
for i = 1:length(kinase_original_list)
    [~,n1] = unique(kinase_original_list{i}{:,1},'stable');
    for j = 1:length(n1)
        ind = find(ismember(p_full_phenotype{i}(:,1),kinase_original_list{i}{n1(j),1}));
        aa1{i} = [aa1{i};p_full_phenotype{i}(ind,2)];
    end
end
kinase_original_list = cellfun(@(x,y) [x,y],kinase_original_list,aa1,'UniformOutput',false);
for i = 1:length(kinase_original_list)
    kinase_original_list{i}.Properties.VariableNames(end) = "Phosphorylated Residue";
end



mkdir(dir1,'Excel_data_Kinase_lists_from_data')
dir_current = strcat(dir1,'\Excel_data_Kinase_lists_from_data');
recycle on % Send to recycle bin instead of permanently deleting.
delete(strcat(dir_current,'/initial_kinase_list.xlsx'));
for i = 1:length(kinase_original_list)
    %header = {'Kinase_ACC','Kinase_Name','Number_of_appearences','Number_of_appearences/Numer_of_PhosphoPeptides','Chi_test_value','Chi_test_pValue','Fisher_test_OddsRatio','Fisher_test_pValue'};
    %writecell(header,strcat(dir_current,'/initial_kinase_list.xlsx'),'Range','A1','Sheet',sample_names{i},'AutofitWidth',1);
    writetable(kinase_original_list{i},strcat(dir_current,'/initial_kinase_list.xlsx'),'Sheet',sample_names{i},'AutofitWidth',1);
end


%It is done on three comparison groups!!!
file_name_kinome = ".\Human_kinome.txt";
kinase_kinase_interactions_excel_data(dir1,dir_current,file_name_kinome,nr_sample_groups,sample_names,list_proteins_significant_non_unique_upreg_complete_kinases,list_proteins_significant_non_unique_downreg_complete_kinase,...
    list_proteins_significant_unique_complete_kinases,kinase_original_list)



%Comment the code from up because it should have been only with kinase
%lists and not the default lists.

save(strcat(dir1,'\workspace.mat'))

mkdir(strcat(dir1,'\Group_comparison_pathways'))

writecell(unique(extractBetween([p_names{:}],'|','|')),strcat(dir1,'\Group_comparison_pathways\Total_Phospho.txt'),"Delimiter","\t");


mkdir(strcat(dir1,'\Kinase_prediction_PSP'))
for i = 1:size(list_proteins_significant_non_unique_upreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_upreg_complete,2)
        if(~isempty(list_proteins_significant_non_unique_upreg_complete{i,j}))
            amx = list_proteins_significant_non_unique_upreg_complete{i,j};
            amx = amx(intersect(intersect(find(cellfun(@str2num,amx(:,11))>0.75),find(vertcat(amx{:,5})<0.01)),find(cellfun(@str2num,amx(:,16))==0)),3);
            writecell(amx,strcat(dir1,"\Kinase_prediction_PSP\",sample_names{i},"vs",sample_names{j},"_001pval.txt"),"Delimiter","\t")
        end
    end
end

for i = 1:size(list_proteins_significant_non_unique_downreg_complete,1)
    for j = 1:size(list_proteins_significant_non_unique_downreg_complete,2)
        if(~isempty(list_proteins_significant_non_unique_downreg_complete{i,j}))
            amx = list_proteins_significant_non_unique_downreg_complete{i,j};
            amx = amx(intersect(intersect(find(cellfun(@str2num,amx(:,11))>0.75),find(vertcat(amx{:,5})<0.01)),find(cellfun(@str2num,amx(:,16))==0)),3);
            writecell(amx,strcat(dir1,"\Kinase_prediction_PSP\",sample_names{j},"vs",sample_names{i},"_001pval.txt"),"Delimiter","\t")
        end
    end
end

% a4 = list_proteins_significant_non_unique_downreg_complete{1,2}
% a5 = list_proteins_significant_non_unique_upreg_complete{1,2}
% a6 = [a4;a5]
% a6 = a6(intersect(intersect(find(cellfun(@str2num,a6(:,11))>0.75),find(vertcat(a6{:,5})<0.05)),find(cellfun(@str2num,a4(:,16))==0)),[3,6])
% writecell(a6,".\Kinase_prediction_PSP\M1ComparedtoM2a_005pval.txt","Delimiter","\t")
% 


end