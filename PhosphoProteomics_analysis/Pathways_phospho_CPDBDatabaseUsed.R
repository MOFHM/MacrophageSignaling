Pathways_Phospho <- function(dir_save,filter_by_pvalue,filter_by_pvalue_FDR_threshold,filter_by_proteinAbundance,file_CPDBPathways,
							Choose_Pathways_databases,DEA_file,universe_file){

library(readxl)
library(openxlsx)
library(xlsx)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(rbioapi)
library(KEGGREST)
library(msigdbr)
library(magrittr)
library(AnnotationDbi)
library(stringr)

options(timeout = max(1000, getOption("timeout")))
# setOptions <- options
# setOptions(clusterProfler.download.method = "wininet")
set.seed(21)

setwd(dir_save)


cpdb_db <- read.table(file_CPDBPathways,fill=TRUE,sep="\t",header=TRUE)
cpdb_db <- cpdb_db[which(cpdb_db$source %in% Choose_Pathways_databases),]
cpdb_db_aux <- data.frame(matrix(nrow=1, ncol=4))
colnames(cpdb_db_aux) <- colnames(cpdb_db)


for(i in 1:nrow(cpdb_db)){
  aux <- unlist(str_split(cpdb_db[i,4],','))
  cpdb_db_aux <- rbind(cpdb_db_aux,data.frame(pathway=rep(cpdb_db[i,1],length(aux)),
                                              external_id=rep(cpdb_db[i,2],length(aux)),
                                              source=rep(cpdb_db[i,3],length(aux)),
                                              entrez_gene_ids = aux))
}

cpdb_db_KeggReactome <- cpdb_db_aux[2:nrow(cpdb_db_aux),]

universe_folder <- getwd()


name = DEA_file

xls_sheet <- excel_sheets(path = name)

for (k in 1:length(xls_sheet)){
  proteins <- read.xlsx(name,sheetName = xls_sheet[k])
  
  if(filter_by_pvalue==1){
    
    proteins <- proteins[which(proteins$FDR<filter_by_pvalue_FDR_threshold),]
    
  }
  
  if(filter_by_proteinAbundance==1){
    
    proteins <- proteins[which(proteins$ProteinCompensation_proteinLevel==0),]
    
  }
  
  
  proteins <- as.data.frame(proteins$Protein)#for all the measured phosphopeptides change saved name as well!
  # proteins <- as.data.frame(proteins$Protein[proteins$FC==-1])#only for unique ones change saved name as well!
  
  
  proteins_uniprot <- AnnotationDbi::select(org.Hs.eg.db, 
                                            keys = as.character(proteins[,]),
                                            columns = c("UNIPROT", "SYMBOL"),
                                            keytype = "UNIPROT")
  
  
  name_universe = paste(universe_folder,'/',universe_file,sep="")
  universe_phospho = read.table(name_universe)
  
  ind <- which(proteins_uniprot$UNIPROT %in% intersect(proteins_uniprot$UNIPROT, universe_phospho$V1))
  proteins_uniprot <- proteins_uniprot[ind,]
  
  
  
  proteins_ENTREZ <- AnnotationDbi::select(org.Hs.eg.db, 
                                           keys = as.character(proteins_uniprot$UNIPROT),
                                           columns = c("ENTREZID", "UNIPROT","SYMBOL"),
                                           keytype = "UNIPROT")
  
  universe_phospho_ENTREZ <- AnnotationDbi::select(org.Hs.eg.db, 
                                                   keys = as.character(universe_phospho$V1),
                                                   columns = c("ENTREZID", "UNIPROT","SYMBOL"),
                                                   keytype = "UNIPROT")
  
  #universe_phospho_ENTREZ <- universe_phospho_ENTREZ[-which(universe_phospho_ENTREZ$ENTREZID%in%proteins_ENTREZ$ENTREZID),]
  
  # proteins_ENTREZ <- proteins_ENTREZ[!is.na(proteins_ENTREZ$ENTREZID),]
  # universe_phospho_ENTREZ <- universe_phospho_ENTREZ[!is.na(universe_phospho_ENTREZ$ENTREZID),]
  KeggReactome_ora_results <- c()
  KeggReactome_ora_results <- enricher(
    gene = proteins_ENTREZ$ENTREZID, # A vector of your genes of interest
    pvalueCutoff = 1, # Can choose a FDR cutoff,
    qvalueCutoff = 1,
    pAdjustMethod = "BH", # Method to be used for multiple testing correction
    universe = universe_phospho_ENTREZ$ENTREZID, # A vector containing your background set genes
    # The pathway information should be a data frame with a term name or
    # identifier and the gene identifiers
    TERM2GENE = dplyr::select(
      cpdb_db_KeggReactome,
      pathway,
      entrez_gene_ids
    ),
    minGSSize = round(1*length(proteins_ENTREZ$ENTREZID)/100),maxGSSize = 5000
  )
  
  if(is.null(KeggReactome_ora_results)==FALSE){
    ind <- which(KeggReactome_ora_results@result$p.adjust<1)
    KeggReactome_ora_results <- data.frame(Pathway = KeggReactome_ora_results@result$ID[ind],
                                           Source = cpdb_db_KeggReactome$source[match(KeggReactome_ora_results@result$ID[ind], cpdb_db_KeggReactome$pathway)],
                                           External_id = cpdb_db_KeggReactome$external_id[match(KeggReactome_ora_results@result$ID[ind], cpdb_db_KeggReactome$pathway)],
                                           pvalue = KeggReactome_ora_results@result$pvalue[ind],
                                           padj_values = KeggReactome_ora_results@result$p.adjust[ind],
                                           members_input_overlap = KeggReactome_ora_results@result$geneID[ind],
                                           geneRatio = KeggReactome_ora_results@result$GeneRatio[ind],
                                           bcgRatio = KeggReactome_ora_results@result$BgRatio[ind],
                                           Number_of_components = sapply(KeggReactome_ora_results@result$ID,FUN = function(x){length(which(cpdb_db_KeggReactome$pathway %in% x))}))
    
    KeggReactome_ora_results$members_input_overlap <- str_replace_all(KeggReactome_ora_results$members_input_overlap,pattern = "/", replacement = ";")
  }
  # 
  # 
  # GO_BP <- rba_panther_enrich(
  #   genes = proteins_ENTREZ$UNIPROT,
  #   organism = 9606,
  #   annot_dataset = "GO:0008150",
  #   cutoff = 0.05,
  #   ref_genes = universe_phospho_ENTREZ$UNIPROT,
  #   ref_organism = 9606,test_type = "FISHER",correction = "FDR"
  # )
  # 
  # GO_BP <- data.frame(Pathway = GO_BP$result$term.label,
  #                     Source = rep(c("GO:BP"),length(GO_BP$result$term.label)),
  #                     padj_values = GO_BP$result$fdr)
  # 
  # 
  # GO_MF <- rba_panther_enrich(
  #   genes = proteins_ENTREZ$UNIPROT,
  #   organism = 9606,
  #   annot_dataset = "GO:0003674",
  #   cutoff = 0.05,
  #   ref_genes = universe_phospho_ENTREZ$UNIPROT,
  #   ref_organism = 9606,test_type = "FISHER",correction = "FDR"
  # )
  # 
  # GO_MF <- data.frame(Pathway = GO_MF$result$term.label,
  #                     Source = rep(c("GO:MF"),length(GO_MF$result$term.label)),
  #                     padj_values = GO_MF$result$fdr)
  # 
  
  
  complete_pathways <- KeggReactome_ora_results   
  
  complete_pathways <- complete_pathways[order(complete_pathways$padj_values),]
  
  # reactome_official <- rbind(GO_BP,GO_MF)
  
  save_name = paste(getwd(),"/complete_pathways.xlsx",sep="")
  if (k==1&file.exists(save_name)){file.remove(save_name)}
  if(nrow(complete_pathways)>0){
    write.xlsx(complete_pathways,save_name,sheetName=xls_sheet[k],append=TRUE,row.names = FALSE)
  }
  
  # save_name = paste(getwd(),"/reactome_GO.xlsx",sep="")
  # if (k==1&file.exists(save_name)){file.remove(save_name)}
  # if(nrow(reactome_official)>0){
  #   write.xlsx(reactome_official,save_name,sheetName=xls_sheet[k],append=TRUE)
  # }
  
}

return(1)

}