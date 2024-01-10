Integrative_Network_analysis <- function(working_dir,BioGRID_data_file,STRING_data_file,IntAct_data_file,file_DEA_names,
                                 phenotype_names,phenotype_comparison,splicing_file_name){
  ###Network analysis
  
  #Set the input files locations!
  
  setwd(working_dir)
  Results_dir <- paste0(getwd(),"/Results")
  
  # 
  # dir_results_R <- "C:/work/Files/Algorithms/Published_scripts/Network_analysis"
  # 
  # working_dir <- getwd()
  # Results_dir <- paste0(getwd(),"/Results")
  # BioGRID_data_file <- "C:/work/Files/General_data_precompiled_files_databases_files/Databases_data/BIOGRID-MV-Physical-4.4.218.mitab.txt"
  # STRING_data_file <- "C:/work/Files/General_data_precompiled_files_databases_files/Databases_data/9606.protein.links.full.v11.5.txt"
  # IntAct_data_file <- "C:/work/Files/General_data_precompiled_files_databases_files/Databases_data/IntAct/IntAct_24022022_07Confidence.txt"
  # 
  # file_DEA_names <- c("list_proteins_significant_group_comparison_PhosphoProteomics.xlsx",
  #                     "list_proteins_significant_group_comparison_Proteomics.xlsx",
  #                     "FoldChange_DESeq2_Gene_26_01_2023_significat_OverlapExp - Modified.xlsx",
  #                     "FoldChange_DRIMSeq_StageR_Uniform_26_1_2023_OverlapExp - Modified.xlsx")
  # 
  # phenotype_names <- c('M1','M2a','M2c')
  # phenotype_comparison <- c('M1vsM2a','M1vsM2c','M2avsM2c','M2avsM1','M2cvsM1','M2cvsM2a')
  # splicing_file_name <- "DTU
  
  
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
  
}


MONET_analysis <- function(working_dir,edge_file_path,server_connection,server_client,server_psw,
                           remote_server_working_dir,plink_path,monet_path,results_remote_location,
                            Monet_method_string,tmp_bin_folder,Windows_Linux){
  
  
  setwd(working_dir)
  Results_dir <- paste0(getwd(),"/Results")
  
  #
  # tmp_bin_folder <- c("C:/tmp_bin_folder")
  #
  #
  # edge_file_path <- "edge_files_STRINGBioGRIDIntAct/Symbol"
  # server_connection <- 1
  # server_client <- "totu@sgl07282"
  #server_psw <- c("C:/work/Files/General_data_precompiled_files_databases_files/SSH_psw.txt")
  server_psw <- gsub("/", "\\", server_psw, fixed=TRUE)
  #remote_server_working_dir <- "/home/totu/tmp"
  #plink_path <- "C:/uitls/Putty/plink.exe"
  plink_path <- gsub("/", "\\", plink_path, fixed=TRUE)
  #monet_path <- "/home/totu/.monet/monet"
  #results_remote_location <- "/home/totu/tmp/MONET"
  #Monet_method_string <- "--method=M1 --container=docker --avgk=10 --linksdir=undirected"
  
  
  full_remote_server_working_dir <- paste(server_client,remote_server_working_dir,sep=":")
  
  if(dir.exists(paste0(Results_dir,"/",edge_file_path))){
    
    dir_MONET <- paste0(Results_dir,"/",edge_file_path,"/*")
    dir_MONET <- gsub("/", "\\", dir_MONET, fixed=TRUE)
    
    dir.create(paste0(Results_dir,"/MONET_analysis"))
    save_dir <- paste0(Results_dir,"/MONET_analysis")
    
    if(Windows_Linux==0){
    if(server_connection){
      library(RCurl)
      
      system(paste("pscp -pwfile", server_psw,"-r", dir_MONET,full_remote_server_working_dir))
      #system(paste0("sshpass -f ", server_psw, " ", "scp -r ",dir_MONET," ", server_client,":",full_remote_server_working_dir))
      system(paste(plink_path,"-ssh",server_client,"-pwfile",server_psw,"mkdir",results_remote_location))
      ff <- list.files(paste0(Results_dir,"/",edge_file_path))
      for (i in 1:length(ff)){
        command_monet <- paste0("echo ",readLines(server_psw)," | sudo -S -i; ",monet_path," --input=",remote_server_working_dir,"/",ff[i]," ",Monet_method_string," --output=",results_remote_location)
        tmp_folder <- paste0(tmp_bin_folder,"/Monet.txt")
        file.remove(tmp_folder)
        write(command_monet,tmp_folder)
        system(paste(plink_path," -ssh ",server_client," -pwfile ",server_psw,"-t -m ", gsub("/", "\\", tmp_folder, fixed=TRUE)))
        file.remove(tmp_folder)
      }
      system(paste0("pscp -pwfile ", server_psw," -r ",server_client,":", results_remote_location," ", save_dir))
      
    }else{
      
      system(paste("wsl mkdir ~/tmp"))
      dir_MONET <- gsub("\\", "/", dir_MONET, fixed=TRUE)
      dir_MONET_wsl <- paste0("/mnt/",dir_MONET)
      dir_MONET_wsl <- str_remove(dir_MONET_wsl,":")
      dir_MONET_wsl <- str_replace(dir_MONET_wsl,"C","c")
      system(paste("wsl cp -r",dir_MONET_wsl,remote_server_working_dir))
      system(paste("wsl mkdir",results_remote_location))
      
      ff <- list.files(paste0(Results_dir,"/",edge_file_path))
      for (i in 1:length(ff)){
        command_monet <- paste0("echo ",readLines(server_psw)," | sudo -S -i; ",monet_path," --input=",remote_server_working_dir,"/",ff[i]," ",Monet_method_string," --output=",command_location)
        system(paste("wsl bash -c",command_monet))
      }
      save_dir_wsl <- paste0("/mnt/",save_dir)
      save_dir_wsl <- str_remove(save_dir_wsl,":")
      save_dir_wsl <- str_replace(save_dir_wsl,"C","c")
      system(paste0("wsl /cp -r ", results_remote_location," ", save_dir_wsl))
    }
    }else{
      system(paste("wsl mkdir ~/tmp"))
      dir_MONET <- gsub("\\", "/", dir_MONET, fixed=TRUE)
      dir_MONET_wsl <- paste0("/mnt/",dir_MONET)
      dir_MONET_wsl <- str_remove(dir_MONET_wsl,":")
      dir_MONET_wsl <- str_replace(dir_MONET_wsl,"C","c")
      system(paste("wsl cp -r",dir_MONET_wsl,remote_server_working_dir))
      system(paste("wsl mkdir",results_remote_location))
      
      ff <- list.files(paste0(Results_dir,"/",edge_file_path))
      for (i in 1:length(ff)){
        command_monet <- paste0("echo ",readLines(server_psw)," | sudo -S -i; ",monet_path," --input=",remote_server_working_dir,"/",ff[i]," ",Monet_method_string," --output=",command_location)
        system(paste("wsl bash -c",command_monet))
      }
      save_dir_wsl <- paste0("/mnt/",save_dir)
      save_dir_wsl <- str_remove(save_dir_wsl,":")
      save_dir_wsl <- str_replace(save_dir_wsl,"C","c")
      system(paste0("wsl /cp -r ", results_remote_location," ", save_dir_wsl))
      
    }
    
    
  }else{
    stop("Edge file directory does not exist")
  }
  
  
}




MONET_pathways <- function(working_dir,CPDB_databases,MONET_background_file=NULL,phenotype_names,
                           phenotype_comparison,CPDB_database_file){
  
  
  
  setwd(working_dir)
  
   Results_dir <- paste0(getwd(),"/Results")
   working_dir <- paste0(Results_dir,"/MONET_analysis/MONET")
  # CPDB_databases <- c("KEGG","Reactome")
  if(is.null(MONET_background_file)){
  MONET_background_file <- paste0(Results_dir,"/Background_total.xlsx")
  }
  # phenotype_names <- c('M1','M2a','M2c')
  # phenotype_comparison <- c('M1vsM2a','M1vsM2c','M2avsM2c','M2avsM1','M2cvsM1','M2cvsM2a')
  # CPDB_database_file <- "C:/work/Files/General_data_precompiled_files_databases_files/Databases_data/CPDB_pathways_genes.tab"
  # 
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

  
}



Circos_and_auxiliary <- function(working_dir,phenotype_names,phenotype_comparison,files_edges_path=NULL,
                                 centralities_file=NULL,TF_Database,file_extension=NULL){

  wkd <- working_dir
  
setwd(wkd)
Results_dir <- paste0(getwd(),"/Results")

working_dir <- paste0(Results_dir,"/MONET_analysis/MONET")
# 
# phenotype_names <- c('M1','M2a','M2c')
# phenotype_comparison <- c('M1vsM2a','M1vsM2c','M2avsM2c','M2avsM1','M2cvsM1','M2cvsM2a')

if(is.null(files_edges_path)){
files_edges_path <- paste0(Results_dir,"/edge_files_STRINGBioGRIDIntAct/Symbol")
}

if(is.null(file_extension)){file_extension = "Total"}

source("Extract_edge_files_for_clusters.R")



status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Extract_edges_tables(working_dir,phenotype_names,phenotype_comparison,files_edges_path,file_extension)
  )
}

if(is.null(centralities_file)){
centralities_file <- paste0(Results_dir,"/STRINGBioGRIDIntAct_centralities_values_CINNA_Total.xlsx")
}
 
setwd(wkd)


source("Extract_TopCentralitiesInModuels.R")

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Extract_CentralitiesValues(working_dir,centralities_file,phenotype_names,file_extension)
  )
}

setwd(wkd)

# TF_Database <- c("C:/work/Files/General_data_precompiled_files_databases_files/Databases_data/_TF.txt")

source("Circos_plots.R")

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- Cricos_plots(working_dir,files_edges_path,TF_Database,file_extension)
  )
}

setwd(wkd)


pathways_dir <- paste0(working_dir,"/CPDB_Pathways")
print(pathways_dir)

source("Monet_cluster_pathways_image_joint.R")

status <- NULL
attempt <- 1
while( is.null(status) && attempt <= 5 ) {
  attempt <- attempt + 1
  try(
    status <- MONET_cluster_pathways_joint_image(pathways_dir,file_extension)
  )
}

}
