FastaFile_extraction <- function(DEA_file,universe_file,save_dir){
  
  library(Rcpi)
  library(openxlsx)
  # library(biomaRt)
  # 
  # mart <- useMart('ENSEMBL_MART_ENSEMBL')
  # mart <- useDataset('hsapiens_gene_ensembl', mart)
  
  
  options(timeout = max(1000, getOption("timeout")))
  set.seed(21)
  
  name <- DEA_file
  xls_sheet <- readxl::excel_sheets(path = name)
  
  for (k in 1:length(xls_sheet)){
    
    proteins <- openxlsx::read.xlsx(name,sheet = xls_sheet[k])
    proteins <- unique(proteins$Protein)

    
    fastas <- c()
    count <- 1
    for (ij in proteins){
      
      attempt <- 1
      aux <- c()
      while(length(grep('sp',aux))==0 && attempt <= 5 ) {
        attempt <- attempt + 1
        Sys.sleep(0.2)
        count <- count+1
        try(
          aux <- getFASTAFromUniProt(ij,parallel = 1)
        )
      } 
      if(length(grep('sp',aux))==1){
      fastas <- c(fastas,aux)
      }
      if(count > 50){
        Sys.sleep(5)
        count <- 1
      }
    }
    
    file.remove(paste0(save_dir,"/",xls_sheet[k],'.txt'))
    write(fastas,paste0(save_dir,"/",xls_sheet[k],'.txt'))
    
  }
  

  
  universe_phospho = read.table(universe_file)
  universe_phospho <- unique(universe_phospho[,1])
  
  fastas <- c()
  for (ij in universe_phospho){
    
    attempt <- 1
    aux <- c()
    while(length(grep('sp',aux))==0 && attempt <= 5 ) {
      attempt <- attempt + 1
      Sys.sleep(0.2)
      count <- count+1
      try(
        aux <- getFASTAFromUniProt(ij,parallel = 1)
      )
    } 
    
    if(length(grep('sp',aux))==1){
    fastas <- c(fastas,aux)
    }
    if(count > 50){
      Sys.sleep(5)
      count <- 1
    }
    
  }
  
  file.remove(paste0(save_dir,"/Phospho_background.txt"))
  write(fastas,paste0(save_dir,"/Phospho_background.txt"))
  return(1)
}
