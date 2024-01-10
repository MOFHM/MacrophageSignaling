###This is the main file to run the phosphoproteomics, transcriptomics and multi-omics network integration pipelines.
###One may need to setup (download and copy) specific files in the right folder locations. These are marked between two lines with "!!!"




###General parameters

tmp_bin_folder <- c("./Location_of_TmpFolder")
if(file.exists(tmp_bin_folder)==FALSE){dir.create(tmp_bin_folder)}



###PhosphoProteomics analysis
library(R.matlab)
library(stringr)

#Set the directroy where all the script are found
setwd(".")
root_dir <- getwd()

#Set the directory where all the results will be save

dir_results_R <- paste0(getwd(),"/Results")
if(file.exists(dir_results_R)==FALSE){dir.create(dir_results_R)}

Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)
setOption(matlab, "readResult/interval", 1200)
setOption(matlab, "readResult/maxTries", 1500 * (60 / 10))


# !!! Download and save in the correct folders the files for the following variables !!!
###All the following filepaths must be with "\\" instead of "/" because this is for matlab.
#The MaxQuant phosphoproteomics file Phospho__STY_Sites.txt
file_name_Phospho = 'Path to Phospho__STY_Sites.txt'
#The MaxQuant proteomics file proteinGroups.txt
protein_file = "Path to proteinGroups.txt"
#The proteomics DEA file
protein_logfc_compensation_file_name = "Path to proteomics DEA file"

# /END !!! Download and save in the correct folders the files for the following variables !!!

#Matlab parameters - You need to have the R.matlab package, the DreamAI package and the limma package (as well as R version 4.1.2) for matlab code to work.
R_exe_location  = "C:/Programs/R/R-4.1.2/bin/R.exe" #Provide the R.exe location!
R_library_location = "C:/Programs/R/R-4.1.2/library" #Provide the R library location!

Intensity_column_of_interest = 112:147 #Provide the column numbers of the intensities columns coresponding for each analyzed sample from the PhosphoProteomics files
p_LFQ_columns = 106:117 #Provide the column numbers of the intensities (LFQ) columns coresponding for each analyzed sample from the Proteomics files

Analysis_considering_NrPhosSites = 1 #Consider the intensities of each individual phosposites - 1, consider the aggregated intensities - 0, consider only the the intensities of one phosphorylation site - 2
Protein_logfc_compensation_flag = 2 #1 - compensate from the biggining for DEA of proteomics; 2 - compensate later for DEA specifically for upstream kinases enrichment
statistical_test_used = "moderate_t_test" #Choose "t_test" or "moderate_t_test" fot the DEA test
Min_number_of_unique_replicates = 0.5 #What is the maximum ratio of missing values for each analysed conditions
choose_Protein_or_Proteins = "Protein" #Chose "Protein" to keep only the first homolog associated with the respective phosphopeptide, "Proteins" to keep all the homologs or "Leading proteins" only to keep the leadining assocaited proteins homologs.
filter_by_evidence = 0 #Should the Phospho talbe be filtered to include only peptides with at least a specific number of evidences? 0 - no ,1 - yes
filter_by_evidenceThreshold=2 #The minimum number of evidences required for each identified phosphopeptide.
normalization_method = 3 #Choose the data normalization method: 0 - Median, 1 - quantile, 2 - zscore, 3 - mean
Consider_unique_sepparately = 0  #0 - the unique values (no measurements in one phenotype) will be imputed; 1 - the unique values will be treated separately.
test_p_value = 0.05 #The FDR threshold for the DEA analysis
test_fc_value = 1 #The log2FC threshold value for the DEA analysis
nr_rep = 4 #Number of sample replicates
nr_sample_groups = 3 #Number of analyzed phenotypes
loc_prob_threshold = 0 #The threshold used to filter the MaxQuant output based on the localization probability
PEP_threshold = 1 #The threshold used to filter the MaxQuant output based on the posterior probability
header_name_string = "M" #Samples information - starting letter for the header and the phenotypes names
sample_names = c('M_1','M_2_a','M_2_c') #The names of the samples from the table header

Protein_abundance_normalization_moment = -1 #0 for normalizing with raw values, 1 for normalizing with imputed values.
#make_PSP_protein_compensation_on_Protein_level = 1

#Imputation analysis characteristics:
#For each imputation method specify for how many missing replciates it
#should apply. When the imputation should not be used then put -1.
#Moreover, specify when to remove from the imputation and analysis scheme
#(for DreamAI 3vs1, 4vs1 and 2vs1, as well as vs0 cannot be computed -
#DreamAI does not impute for less that 2 measurements per group.
#Set as well the parameters for the random draw imputation method - the mean shift and the std.

imputation_method = 1 #Imputation method to use 0 - property method, 1 - DreamAI method and random draw from a distribution
imputation_remove = c(-1)
imputation_DreamAI = c(-1)
distirbution_imputation_penalty = c(3,4)
distirbution_imputation_normal = c(0,1,2)
Meanshift_imputation_penalty = 1.8
STDshift_imputation_penalty = 0.3
Meanshift_imputation_normal = 0.5
STDshift_imputation_normal = 0.3


PSP_Loc_prob =  1 #0 - do the Upstream kinases identification based on curated data based on the data as it is; 1 - filter for loc_prob below 0.75
PSP_protein_compensation = 1 #0 - do the Upstream kinases identification based on curated data based on the data as it is; 1 - remove the phosphopeptides that do not have a high enought log2FC when comapred to the proteomics.

PSP_file = ".\\PTM_Merged_databases.txt" #Database of all kinase-substrate interactions used for the upstream kinase activity inference based on currated data.
PSP_p_val = 0.01 #FDR threshold used to filter the DEA results before performing the upstream kinase enrichment algorithm
Potential_contaminant_string = "Potential contaminant" #The column name that contains the information related to the measured contaminants (if a protein is considered a contaminant or not)
Rev_string = "Reverse" #The column name that points to the flag which identifies peptides identified in reverse
Loc_prob_string = "Localization prob" #The column name that contains the localization probablity values
PEP_string = "PEP" #The column name that contains the PEP values
MS_MS_ID_string = "MS/MS IDs" #The column name that contains the MS/MS Ids
Evidence_IDs_string = "Evidence IDs" #The column name that contains the Evidence Ids
Positions_within_proteins_String = "Positions within proteins" #The column name that contains the position of the identified phosphorylated residues
Amino_acid_String = "Amino acid" #The column name that contains the the symbol of the identified phosphorylated residues
Sequence_window_String = "Sequence window" #The column name that contains the identified sequence window
data_transformation = "log2" #What kind of data transfromation to use on the Intensities values to ensure the normality of the measurementes
file_name_kinome = ".\\Human_kinome.txt" #Database that contains all the human kinases

#The script directory to be run for the initial Phosphoproteomics analysis, including the DEA, upstream kinase analysis based on curated datasets and quality control
ScriptsDir <- paste0(getwd(),"\\PhosphoProteomics_analysis")
ScriptsDir <- gsub("/","\\\\",ScriptsDir)

dir_results = gsub("/","\\\\",dir_results_R)
dir_results <- paste0(dir_results)

setVariable(matlab,ScriptsDir=ScriptsDir)
evaluate(matlab,"cd(ScriptsDir)")

setVariable(matlab, file_name_Phospho=file_name_Phospho,protein_file=protein_file,protein_logfc_compensation_file_name=protein_logfc_compensation_file_name,
            dir_results=dir_results,Analysis_considering_NrPhosSites=Analysis_considering_NrPhosSites,Protein_abundance_normalization_moment=Protein_abundance_normalization_moment,
            Intensity_column_of_interest=Intensity_column_of_interest,p_LFQ_columns=p_LFQ_columns,Protein_logfc_compensation_flag=Protein_logfc_compensation_flag,
            statistical_test_used=statistical_test_used,Min_number_of_unique_replicates=Min_number_of_unique_replicates,choose_Protein_or_Proteins=choose_Protein_or_Proteins,
            filter_by_evidence=filter_by_evidence,filter_by_evidenceThreshold=filter_by_evidenceThreshold,normalization_method=normalization_method,Consider_unique_sepparately=Consider_unique_sepparately,
            test_p_value=test_p_value,test_fc_value=test_fc_value,nr_rep=nr_rep,nr_sample_groups=nr_sample_groups,loc_prob_threshold=loc_prob_threshold,PEP_threshold=PEP_threshold,
            header_name_string=header_name_string,sample_names=sample_names,imputation_method=imputation_method,imputation_remove=imputation_remove,
            imputation_DreamAI=imputation_DreamAI,distirbution_imputation_penalty=distirbution_imputation_penalty,distirbution_imputation_normal=distirbution_imputation_normal,
            Meanshift_imputation_penalty=Meanshift_imputation_penalty,STDshift_imputation_penalty=STDshift_imputation_penalty,Meanshift_imputation_normal=Meanshift_imputation_normal,
            STDshift_imputation_normal=STDshift_imputation_normal,PSP_file=PSP_file,PSP_Loc_prob=PSP_Loc_prob,PSP_protein_compensation=PSP_protein_compensation,PSP_p_val=PSP_p_val,
            Potential_contaminant_string=Potential_contaminant_string,Rev_string=Rev_string,Loc_prob_string=Loc_prob_string,PEP_string=PEP_string,
            MS_MS_ID_string=MS_MS_ID_string,Evidence_IDs_string=Evidence_IDs_string,Positions_within_proteins_String=Positions_within_proteins_String,Amino_acid_String=Amino_acid_String,
            Sequence_window_String=Sequence_window_String,data_transformation=data_transformation,file_name_kinome=file_name_kinome,R_exe_location=R_exe_location,R_library_location=R_library_location)
evaluate(matlab,"[dir1,Perc_fals_poz,best_tr,test_res,boxcox_param,dist_names_group,parameter_mean_values_group,dist_performance_ind_group] = PhosphoProteomics_analysis(file_name_Phospho,protein_file,protein_logfc_compensation_file_name,dir_results,Analysis_considering_NrPhosSites,Protein_abundance_normalization_moment,Intensity_column_of_interest,p_LFQ_columns,Protein_logfc_compensation_flag,statistical_test_used,Min_number_of_unique_replicates,choose_Protein_or_Proteins,filter_by_evidence,filter_by_evidenceThreshold,normalization_method,Consider_unique_sepparately,test_p_value,test_fc_value,nr_rep,nr_sample_groups,loc_prob_threshold,PEP_threshold,header_name_string,sample_names,imputation_method,imputation_remove,imputation_DreamAI,distirbution_imputation_penalty,distirbution_imputation_normal,Meanshift_imputation_penalty,STDshift_imputation_penalty,Meanshift_imputation_normal,STDshift_imputation_normal,PSP_file,PSP_Loc_prob,PSP_protein_compensation,PSP_p_val,Potential_contaminant_string,Rev_string,Loc_prob_string,PEP_string,MS_MS_ID_string,Evidence_IDs_string,Positions_within_proteins_String,Amino_acid_String,Sequence_window_String,data_transformation,file_name_kinome,R_exe_location,R_library_location)")

aux <- getVariable(matlab,c("dir1","Perc_fals_poz","best_tr","test_res","boxcox_param","dist_names_group","parameter_mean_values_group","dist_performance_ind_group"))

dir_results_R <- aux$dir1
dir_results_R <- gsub("\\\\","/",dir_results_R)

###Pathway ORA parameter
dir_save <- paste0(dir_results_R,'/Group_comparison_pathways')

setwd(paste0(root_dir,"/PhosphoProteomics_analysis"))

filter_by_pvalue <- 0 #Do you want to filter the DEA tables for a particular FDR value, 0-No, 1- Yes (The tables are already filtered based on the above selected FDR value
filter_by_pvalue_FDR_threshold <- 0.05 #FDR threshold to filter the DEA tables for ORA analysis
filter_by_proteinAbundance <- 0 #Do you want to filter the DEA tables to exclude the PhosphoProteins whose phosphorylation levels is not higher than the one expected from proteomics, 0-No, 1- Yes
file_CPDBPathways <- paste0(getwd(),"/CPDB_pathways_genes.tab") #Set the path to the CPDB pathways data
Choose_Pathways_databases <- c("KEGG","Reactome") #Select which Databases to use for the ORA enrichment
DEA_file <- paste0(dir_results_R,"/Excel_data_unique_peptide_values/list_proteins_significant_group_comparison.xlsx")
universe_file <- c("Total_Phospho.txt")

source('./Pathways_phospho_CPDBDatabaseUsed.R')

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Pathways_Phospho(dir_save,filter_by_pvalue,filter_by_pvalue_FDR_threshold,filter_by_proteinAbundance,file_CPDBPathways,
                          Choose_Pathways_databases,DEA_file,universe_file)
  )
} 

###Filter the pathways ORA results to remove the pathways that overlap too much
setwd(paste0(root_dir,"/PhosphoProteomics_analysis"))

close(matlab,port=9999)
Matlab$startServer(port=9998)
matlab <- Matlab(host="localhost", port=9998)
isOpen <- open(matlab,port=9998)
setOption(matlab, "readResult/interval", 300)
setOption(matlab, "readResult/maxTries", 600 * (60 / 10))
setVariable(matlab,ScriptsDir=ScriptsDir)
evaluate(matlab,"cd(ScriptsDir)")

name_dir <- gsub("/","\\\\",dir_save)
#General parameters that describe the pathways tables
source = "source"
source_databases = c('reactome','kegg')
effective_size = "effective_size"
members_input_overlap = "members_input_overlap"
pvalue = "pvalue"
padj_values = "padj_values"

removal_method = 1 #1 - remove pathways with less than two different memebers; 2 - remove the pathwhays that overlap more than 80% having the reference the second pathway
maximum_pathway_dim = 300 #Maximum number of elements that a pathway should have in theory
min_number_members = 5 #Minimum number of elements that were found to be assoicated with a pathway

setVariable(matlab, name_dir=name_dir,source=source,source_databases=source_databases,
            effective_size=effective_size,members_input_overlap=members_input_overlap,
            pvalue=pvalue,padj_values=padj_values,removal_method=removal_method,
            maximum_pathway_dim=maximum_pathway_dim,min_number_members=min_number_members)

evaluate(matlab, "[status] = remove_pathways_overlap(name_dir,source,source_databases,effective_size,members_input_overlap,pvalue,padj_values,removal_method,maximum_pathway_dim,min_number_members)")
close(matlab,port=9998)


### Run the NetPhorest analysis for predicting the most active kinases. For this, WSL 2 must be installed and working, as well as NetPhorest installed!
setwd(paste0(root_dir,"/PhosphoProteomics_analysis"))


DEA_file <- paste0(dir_results_R,"/Excel_data_unique_peptide_values/list_proteins_significant_group_comparison.xlsx")
universe_file <- paste0(dir_results_R,'/Group_comparison_pathways/Total_Phospho.txt')

dir.create(paste0(dir_results_R,"/NetPhorest"))
dir.create(paste0(dir_results_R,"/NetPhorest/Original_fasta_files"))
save_dir <- paste0(dir_results_R,"/NetPhorest/Original_fasta_files")

#Extract the fastafiles of the proteins for each contrast of interest
source('FastaFile_extraction.R')


status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- FastaFile_extraction(DEA_file,universe_file,save_dir)
  )
} 

naming_convention_comparison = c('M1','M2a','M2c'); #The phneotype names
LocProb_threshold <- 0.75 #The minimum localization probability used for filtering the DEA tables
FDR_threshold <- 0.05 #The maximum FDR used for filtering the DEA tables
compensate_for_protein_abundances <- 1 #Should the phosphopeptides be compensated for the corresponding protein abundance
dir_original_fasta_files <- paste0(dir_results_R,"/NetPhorest/Original_fasta_files")
dir_original_fasta_files <- gsub("/", "\\", dir_original_fasta_files, fixed=TRUE)
dir_original_fasta_files_linux <- paste0("/mnt/",dir_results_R,"/NetPhorest")
dir_original_fasta_files_linux <- str_replace_all(dir_original_fasta_files_linux,':','')
dir_original_fasta_files_linux <- str_replace_all(dir_original_fasta_files_linux,'C','c')
date_dir <- strsplit(dir_results_R,'/')
date_dir <-date_dir[[1]][lengths(date_dir)]
NetPhorest_script <- "Netphorest scripts path under the Linux wsl2 - it should be something similar to /mnt/userName/NetPhorest/scripts"

workspace_file <- paste0(dir_results_R,"/workspace")
workspace_file <- gsub("/", "\\", workspace_file, fixed=TRUE)

file_PhosphoProteomics <- paste0(dir_results_R,"/NetPhorest/NetPhorest_results/Results_file_Phospho_background.txt")
file_PhosphoProteomics <- gsub("/", "\\", file_PhosphoProteomics, fixed=TRUE)
root_NetPhorest <- paste0(dir_results_R,"/NetPhorest/NetPhorest_results")
Interaction_database <- c(".\\Interaction_Lookup_tableUniprot.txt")
NetPhorest_full_names <- c(".\\NetPhorest_full_kinase.txt")


Matlab$startServer(port=9995)
matlab <- Matlab(host="localhost", port=9995)
isOpen <- open(matlab,port=9995)
setOption(matlab, "readResult/interval", 3000)
setOption(matlab, "readResult/maxTries", 6000 * (60 / 10))
setVariable(matlab,ScriptsDir=ScriptsDir)
evaluate(matlab,"cd(ScriptsDir)")


setVariable(matlab,naming_convention_comparison=naming_convention_comparison,LocProb_threshold=LocProb_threshold,FDR_threshold=FDR_threshold,
            compensate_for_protein_abundances=compensate_for_protein_abundances,dir_original_fasta_files=dir_original_fasta_files,
            dir_original_fasta_files_linux=dir_original_fasta_files_linux,date_dir=date_dir,NetPhorest_script=NetPhorest_script,
            file_PhosphoProteomics=file_PhosphoProteomics,root_NetPhorest=root_NetPhorest,Interaction_database=Interaction_database,
            NetPhorest_full_names=NetPhorest_full_names,workspace_file=workspace_file)

evaluate(matlab,"NetPhorest_Kinase_extraction2_byPhosphoPeptides(naming_convention_comparison,LocProb_threshold,FDR_threshold,compensate_for_protein_abundances,dir_original_fasta_files,date_dir,dir_original_fasta_files_linux,NetPhorest_script,root_NetPhorest,file_PhosphoProteomics,Interaction_database,NetPhorest_full_names,workspace_file)")

close(matlab,port=9995)

files_NetPh <- list.files(paste0(root_NetPhorest,"/Excel_results/No_interactions"),".xlsx",full.names = TRUE)
xlsx_name_NP <- paste0(root_NetPhorest,"/Excel_results/No_interactions/NetPhorest_results.xlsx")

if (file.exists(xlsx_name_NP)) {
  #Delete file if it exists
  file.remove(xlsx_name_NP)
}

for (i in files_NetPh){
  aux <- openxlsx::read.xlsx(i,colNames = TRUE)
  xlsx::write.xlsx(aux,append=TRUE,xlsx_name_NP,sheetName=gsub(".*Results_file_(.+).xlsx.*", "\\1", i),row.names = FALSE)
}


files_NetPh <- list.files(paste0(root_NetPhorest,"/Excel_results/With_interactions"),".xlsx",full.names = TRUE)
xlsx_name_NP <- paste0(root_NetPhorest,"/Excel_results/With_interactions/NetPhorest_results.xlsx")

if (file.exists(xlsx_name_NP)) {
  #Delete file if it exists
  file.remove(xlsx_name_NP)
}

for (i in files_NetPh){
  aux <- openxlsx::read.xlsx(i,colNames = TRUE)
  xlsx::write.xlsx(aux,append=TRUE,xlsx_name_NP,sheetName=gsub(".*Results_file_(.+).xlsx.*", "\\1", i),row.names = FALSE)
}

###Create the kinase-kinase interaction networks
Matlab$startServer(port=9993)
matlab <- Matlab(host="localhost", port=9993)
isOpen <- open(matlab,port=9993)
setOption(matlab, "readResult/interval", 300)
setOption(matlab, "readResult/maxTries", 600 * (60 / 10))
setVariable(matlab,ScriptsDir=ScriptsDir)
evaluate(matlab,"cd(ScriptsDir)")


file_kinases_Phospho = ".\\PTM_Merged_databases.txt"
dir_cr = gsub("/", "\\", dir_results_R, fixed=TRUE)
file_kinases_list = ".\\Human_kinome.txt"
file_TF = ".\\_TF.txt"
file_TF_mapping = ".\\TF_entrez_uniport_list_mapping_filtered.txt"
file_TF_measured = ".\\file_TF_measured.txt" #This should be the list of measured TF

dir_NetPhorest_file = paste0(dir_results_R,"/NetPhorest/NetPhorest_results/Excel_results/No_interactions")
dir_NetPhorest_file = gsub("/", "\\", dir_NetPhorest_file, fixed=TRUE)
workspace_file <- paste0(dir_results_R,"/workspace")
workspace_file <- gsub("/", "\\", workspace_file, fixed=TRUE)

setVariable(matlab,file_kinases_Phospho=file_kinases_Phospho,dir_cr=dir_cr,
            file_kinases_list=file_kinases_list,file_TF=file_TF,file_TF_mapping=file_TF_mapping,
            file_TF_measured=file_TF_measured,dir_NetPhorest_file=dir_NetPhorest_file,
            workspace_file=workspace_file)
#Create the kinase-kinase interaction tables using the two-stage reccursive algorithms
evaluate(matlab,"kinase_cascade_generationKinase_source_file(file_kinases_Phospho,dir_cr,file_kinases_list,file_TF,file_TF_mapping,file_TF_measured,dir_NetPhorest_file,workspace_file)")

close(matlab,port=9993)

#Plot the networks - You need to have Cytoscape installed!
setwd(paste0(root_dir,"/PhosphoProteomics_analysis"))

working_dir <- paste0(dir_results_R,"/Kinase_Network_data")
Cytoscape_dir <- "Path to \\Cytoscape.exe"
Cytoscape_style_file <- paste0(root_dir,"/Kinase_KInase_Networks_Style3_BLUE.xml")
network_directed <- TRUE

source("Kinase_Kinase_networks.R")

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Kinase_Kinase_networks(working_dir,Cytoscape_dir,Cytoscape_style_file,network_directed)
  )
} 

Sys.sleep(10)


#dir.create(paste0(root_dir,"/Network_analysis"))
file.copy(from=paste0(dir_results_R,"/Excel_data_unique_peptide_values/list_proteins_significant_group_comparison.xlsx"),
          to="./Network_analysis/list_proteins_significant_group_comparison_PhosphoProteomics.xlsx")



##############

#Transcriptomics analysis
library(R.matlab)

#Perform the DGE analysis for each of the considered experiments
setwd(paste0(root_dir,"/Transcriptomics_analysis"))


# !!! Download and save in the correct folders the files for the following variables !!!

working_dir <- paste0(root_dir,"/Transcriptomics_analysis/Results") # You need to download the SRA files and analyze them. Contact the authors for the processed files and prepare a secure file transfer method. A template is already present
phenotype_names <- c("M1","M2a","M2c") # The sample phenotypes names
experiments_folders <- list.dirs(working_dir,recursive = FALSE)
sample_info_fileName <- "Samples_info.csv" #The contrast file
Transcripts_fileName <- "Short_salmon_file_ONLY_Proteins.txt" #The filtered RNAseq file to include only protein-coding isoforms

# /END !!! Download and save in the correct folders the files for the following variables !!!


source("RNAseq_differential_isoform_DRIMSeq_DGE.R")
RNAseq_DEA_DTU(working_dir,phenotype_names,experiments_folders,sample_info_fileName,Transcripts_fileName)

#Merge the results from all the different experiments together
setwd(paste0(root_dir,"/Transcriptomics_analysis"))
ScriptsDir <- paste0(getwd())
ScriptsDir <- gsub("/","\\\\",ScriptsDir)
Matlab$startServer(port=9991)
matlab <- Matlab(host="localhost", port=9991)
isOpen <- open(matlab,port=9991)
setOption(matlab, "readResult/interval", 300)
setOption(matlab, "readResult/maxTries", 600 * (60 / 10))
setVariable(matlab,ScriptsDir=ScriptsDir)
evaluate(matlab,"cd(ScriptsDir)")


working_dir <- ".\\Results"
folder_prefix <- "PRJNA"
phenotype_names <- c('M1','M2a','M2c')
phenotype_comparison <- c('M1vsM2a','M1vsM2c','M2avsM2c','M2avsM1','M2cvsM1','M2cvsM2a')
DGE_fileName <- "DGE_DESeq2\\FoldChange_DESeq2_Gene.xlsx"
DRIMSeq_StageR_fileName <- "DTU_DRIMSeq\\FoldChange_DRIMSeq_StageR_Uniform.xlsx"


setVariable(matlab,working_dir=working_dir,folder_prefix=folder_prefix,phenotype_names=phenotype_names,
            phenotype_comparison=phenotype_comparison,DGE_fileName=DGE_fileName,
            DRIMSeq_StageR_fileName=DRIMSeq_StageR_fileName)

evaluate(matlab,"RNAseq_DGE_files(working_dir,folder_prefix,phenotype_names,phenotype_comparison,DGE_fileName,DRIMSeq_StageR_fileName)")

close(matlab,port=9991)
#####


###Network analysis - The interaction databases must be downloaded a priori and the trasncriptomics nad phosphoproteomics folder structure as was created above

#Set the input files locations!

setwd(paste0(root_dir,"/Network_analysis"))

dir_results_R <- paste0(root_dir,"/Network_analysis")

# !!! Extract and save in the correct folders the files for the following variables !!!

working_dir <- getwd()
Results_dir <- paste0(getwd(),"/Results")
BioGRID_data_file <- "./Databases_data/BIOGRID-MV-Physical-4.4.218.mitab.txt"
STRING_data_file <- "./Databases_data/9606.protein.links.full.v11.5.txt"
IntAct_data_file <- "./Databases_data/IntAct/IntAct_24022022_07Confidence.txt"

file_DEA_names <- c("list_proteins_significant_group_comparison_PhosphoProteomics.xlsx",
                    "list_proteins_significant_group_comparison_Proteomics.xlsx",
                    "FoldChange_DESeq2_Gene_26_01_2023_significat_OverlapExp - Modified.xlsx",
                    "FoldChange_DRIMSeq_StageR_Uniform_26_1_2023_OverlapExp - Modified.xlsx")

# !!! /END Download and save in the correct folders the files for the following variables !!!

phenotype_names <- c('M1','M2a','M2c')
phenotype_comparison <- c('M1vsM2a','M1vsM2c','M2avsM2c','M2avsM1','M2cvsM1','M2cvsM2a')
splicing_file_name <- "DTU"

#Network construction and centrality computation
source('Network_analysis_BioGRID_STRING_IntAct.R')

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Network_analysis(working_dir,Results_dir,BioGRID_data_file,
                               STRING_data_file,IntAct_data_file,file_DEA_names,
                               phenotype_names,phenotype_comparison,splicing_file_name)
  )
}

setwd(paste0(root_dir,"/Network_analysis"))

#MONET decomposition part - You need to have MONET installed !!!
Results_dir <- paste0(getwd(),"/Results")
edge_file_path <- "edge_files_STRINGBioGRIDIntAct/Symbol"
remote_server_working_dir <- "/home/Username/tmp" #Change username with your valid user
monet_path <- "/home/Username/.monet/monet" #Provide the monet path
results_remote_location <- "/home/Username/tmp/MONET"
Monet_method_string <- "--method=M1 --container=docker --avgk=10 --linksdir=undirected"


full_remote_server_working_dir <- paste(server_client,remote_server_working_dir,sep=":")

if(dir.exists(paste0(Results_dir,"/",edge_file_path))==1){
  
  dir_MONET <- paste0(Results_dir,"/",edge_file_path,"/*")
  dir_MONET <- gsub("/", "\\", dir_MONET, fixed=TRUE)
  
  dir.create(paste0(Results_dir,"/MONET_analysis"))
  save_dir <- paste0(Results_dir,"/MONET_analysis")
  
    library(RCurl)

    
    system(paste("wsl mkdir ~/tmp"))
    dir_MONET <- gsub("\\", "/", dir_MONET, fixed=TRUE)
    dir_MONET_wsl <- paste0("/mnt/",dir_MONET)
    dir_MONET_wsl <- str_remove(dir_MONET_wsl,":")
    dir_MONET_wsl <- str_replace(dir_MONET_wsl,"C","c")
    system(paste("wsl cp -r",dir_MONET_wsl,remote_server_working_dir))
    system(paste("wsl mkdir",results_remote_location))
    
    ff <- list.files(paste0(Results_dir,"/",edge_file_path))
    for (i in 1:length(ff)){
      command_monet <- paste0(monet_path," --input=",remote_server_working_dir,"/",ff[i]," ",Monet_method_string," --output=",results_remote_location)
      system(paste("wsl bash -c",command_monet))
    }
    save_dir_wsl <- paste0("/mnt/",save_dir)
    save_dir_wsl <- str_remove(save_dir_wsl,":")
    save_dir_wsl <- str_replace(save_dir_wsl,"C","c")
    system(paste0("wsl /cp -r ", results_remote_location," ", save_dir_wsl))
  
  
}else{
  stop("Edge file directory does not exist")
}


#Pathway extraction part
setwd(paste0(root_dir,"/Network_analysis"))

Results_dir <- paste0(getwd(),"/Results")
working_dir <- paste0(Results_dir,"/MONET_analysis")
CPDB_databases <- c("KEGG","Reactome")
MONET_background_file <- paste0(Results_dir,"Background_total.xlsx")
phenotype_names <- c('M1','M2a','M2c')
phenotype_comparison <- c('M1vsM2a','M1vsM2c','M2avsM2c','M2avsM1','M2cvsM1','M2cvsM2a')
CPDB_database_file <- "./CPDB_pathways_genes.tab"

source("MONET_pathways_extraction_CPDB.R")


status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- MONET_Pathways_extraction(working_dir,CPDB_database_file,CPDB_databases,
                                        MONET_background_file,phenotype_names,phenotype_comparison)
  )
}

#Image construction part
setwd(paste0(root_dir,"/Network_analysis"))
Results_dir <- paste0(getwd(),"/Results")

working_dir <- paste0(Results_dir,"/MONET_analysis")
phenotype_names <- c('M1','M2a','M2c')
phenotype_comparison <- c('M1vsM2a','M1vsM2c','M2avsM2c','M2avsM1','M2cvsM1','M2cvsM2a')
files_edges_path <- paste0(Results_dir,"edge_files_STRINGBioGRIDIntAct/Symbol")


file_extension <- "Total"

source("Extract_edge_files_for_clusters.R")


status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Extract_edges_tables(working_dir,phenotype_names,phenotype_comparison,files_edges_path)
  )
}

setwd(paste0(root_dir,"/Network_analysis"))
Results_dir <- paste0(getwd(),"/Results")

working_dir <- paste0(Results_dir,"/MONET_analysis")
centralities_file <- paste0(Results_dir,"/STRINGBioGRIDIntAct_centralities_values_CINNA_Total.xlsx")
phenotype_names <- c("M1","M2a","M2c")

source("Extract_TopCentralitiesInModuels.R")

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Extract_CentralitiesValues(working_dir,centralities_file,phenotype_names)
  )
}



setwd(paste0(root_dir,"/Network_analysis"))
Results_dir <- paste0(getwd(),"/Results")

working_dir <- paste0(Results_dir,"/MONET_analysis")
edges_dir <- paste0(Results_dir,"\edge_files_STRINGBioGRIDIntAct\Symbol")
TF_Database <- c("./_TF.txt")

source("Circos_plots.R")

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Cricos_plots(working_dir,edges_dir,TF_Database)
  )
}



setwd(paste0(root_dir,"/Network_analysis"))
Results_dir <- paste0(getwd(),"/Results")
pathways_dir <- paste0(Results_dir,"/MONET_analysis","/Results/CPDB_Pathways")

source("Monet_cluster_pathways_image_joint.R")

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- MONET_cluster_pathways_joint_image(pathways_dir)
  )
}


######