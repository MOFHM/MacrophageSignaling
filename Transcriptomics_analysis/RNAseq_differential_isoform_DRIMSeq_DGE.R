RNAseq_DEA_DTU <- function(working_dir,phenotype_names,experiments_folders,sample_info_fileName,Transcripts_fileName){

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# BiocManager::install("DESeq2")
# BiocManager::install("DRIMSeq")
# BiocManager::install("tximport")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("GenomicFeatures")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("magrittr")
# install.packages("readxl")
# install.packages("writexl")
# install.packages("readr")
# install.packages("XLConnect")


library(readxl)
#library(openxlsx)
#library(xlsx)
library(writexl)
library(tximport)
library(DRIMSeq)
library(GenomicFeatures)
library(readr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(DESeq2)
library(XLConnect)
library(AnnotationDbi)
library(magrittr)
library(org.Hs.eg.db)
library(limma)
library(edgeR)
library(stageR)


set.seed(21)
options(timeout = max(1000, getOption("timeout")))

#The same database that was used for Salmon alignment
gtf <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.chr_patch_hapl_scaff.annotation.gff3.gz"
txdb.filename <- "gencode.v33.annotation.sqlite"
txdb <- makeTxDbFromGFF(gtf)

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

setwd(working_dir)

for (i in 1:length(experiments_folders)){
  
  setwd(experiments_folders[i])
  #Read the experiment metadata
  samples <- read.table(file.path(getwd(), sample_info_fileName), header = TRUE,sep=',')
  
  #Read the files' names
  files <- file.path(getwd(), samples$sample_id)
  for (j in 1:length(files)){
    aux = list.files(files[j])
    ind <- grep(x = aux,pattern = Transcripts_fileName)
    files[j] <- paste(files[j],'/',aux[ind],sep="")
  }
  
  names(files) <- samples$sample_id
  
  
  #Import the Salmon data using the scaledTPM option that will change both the abundances (TMP) and the counts depending on the gene length
  #TRY dtuScaledTMP as well!!
  data <- tximport(files, type="salmon", txOut=TRUE,
                   countsFromAbundance="scaledTPM",tx2gene = tx2gene,countsCol = 5,lengthCol = 3,abundanceCol = 4)
  #Take only the measured data
  abundance <- data$counts
  abundance <- abundance[rowSums(abundance) > 0,]
  
  
  #DRIMSeq
  
  #Select only the transcripts that are known
  txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
  tab <- table(txdf$GENEID)
  txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
  
  txdf <- txdf[match(rownames(abundance),txdf$TXNAME),]
  
  #all(rownames(abundance) == txdf$TXNAME)
  
  #txdf <- tx2gene[match(rownames(abundance),tx2gene$TXNAME),]
  
  abundance_complete <- data.frame(gene_id=txdf$GENEID,
                                   feature_id=txdf$TXNAME,
                                   abundance)
  
  drimseq_results <- dmDSdata(counts=abundance_complete, samples=samples)
  #Filtering parameters
  n <- length(drimseq_results@samples$condition)
  n.small <- round(length(drimseq_results@samples$condition)/length(levels(drimseq_results@samples$condition))-1)
  filter_threshold = quantile(abundance)[3]
  #Filter the data
  drimseq_results <- dmFilter(drimseq_results,
                              min_samps_feature_expr=n.small, min_feature_expr=filter_threshold,
                              min_samps_feature_prop=n.small, min_feature_prop=0.05,
                              min_samps_gene_expr=n, min_gene_expr=filter_threshold)
  
  #Create the design matrix
  design_full <- model.matrix(~ 0 + condition, data=DRIMSeq::samples(drimseq_results))
  #colnames(design_full)
  
  #Run Drimseq analysis
  drimseq_results <- dmPrecision(drimseq_results, design=design_full, add_uniform=TRUE)
  drimseq_results <- dmFit(drimseq_results, design=design_full, add_uniform=TRUE)

  #Write the results
  dir.create(paste0(getwd(),"/DTU_DRIMSeq"))
  xls_name = paste(getwd(),'/DTU_DRIMSeq/FoldChange_DRIMSeq_Transcripts_Uniform.xlsx',sep="")
  file.remove(xls_name)
  xls_name_gene = paste(getwd(),'/DTU_DRIMSeq/FoldChange_DRIMSeq_Gene_Uniform.xlsx',sep="")
  file.remove(xls_name_gene)
  xls_name_stageR = paste(getwd(),'/DTU_DRIMSeq/FoldChange_DRIMSeq_StageR_Uniform.xlsx',sep="")
  file.remove(xls_name_stageR)
  xls_list <- list()
  xls_list_gene <- list()
  xls_list_stageR <- list()
  names_sheets <- c()
  analysis_names <- levels(DRIMSeq::samples(drimseq_results)$condition)
  
  #Extract the results for each comparison group on the transcript and on the gene level + do StageR analysis
  for (k in 1:length(colnames(design_full))){
    for (kk in (k+1):length(colnames(design_full))){
      if(kk<=length(colnames(design_full)) & kk!=k){
        con1 <- as.numeric(rep(0,length(colnames(design_full))))
        
        #First contrast - feature level
        con1[k] <- -1
        con1[kk] <- 1
        drimseq_results1 <- dmTest(drimseq_results, contrast = con1 )
        res1 <- DRIMSeq::results(drimseq_results1,level='feature')
        
        uniprot_names <- mapIds(org.Hs.eg.db, 
                                keys = gsub("\\..*","",res1$gene_id),
                                column = c("UNIPROT"),
                                keytype = "ENSEMBL",
                                multiVals = "first")
        
        symbol_names <- mapIds(org.Hs.eg.db, 
                               keys = gsub("\\..*","",res1$gene_id),
                               column = c("SYMBOL"),
                               keytype = "ENSEMBL",
                               multiVals = "first")
        
        #Instead of log2(FC) compute the log of proportions and log of cofficients
        a <- DRIMSeq::proportions(drimseq_results1)
        akk <- a[,grep( x = colnames(a), pattern = analysis_names[kk], ignore.case =TRUE)]
        akk <- rowMeans(akk)
        ak <- a[,grep( x = colnames(a), pattern = analysis_names[k], ignore.case =TRUE)]
        ak <- rowMeans(ak)
        log2_proprotions <- log2(akk/ak)
        a <- DRIMSeq::coefficients(drimseq_results1,level='feature')
        akk <- a[,grep( x = colnames(a), pattern = analysis_names[kk], ignore.case =TRUE)]
        ak <- a[,grep( x = colnames(a), pattern = analysis_names[k], ignore.case =TRUE)]
        log2_coefficents = log2(exp(akk-ak))
        
        
        res1 <- data.frame(Name_Gene = res1$gene_id, Name_Transcript = res1$feature_id,Name_Uniprot = uniprot_names, Name_Symbol = symbol_names, lr=res1$lr,
                           df = res1$df, pvalue = res1$pvalue, adj_pvalue = res1$adj_pvalue, log2_Proportions = log2_proprotions, log2_Coefficents = log2_coefficents)
        
        #First contrast - gene level
        res1_gene <- DRIMSeq::results(drimseq_results1,level='gene')
        
        uniprot_names <- mapIds(org.Hs.eg.db, 
                                keys = gsub("\\..*","",res1_gene$gene_id),
                                column = c("UNIPROT"),
                                keytype = "ENSEMBL",
                                multiVals = "first")
        
        symbol_names <- mapIds(org.Hs.eg.db, 
                               keys = gsub("\\..*","",res1_gene$gene_id),
                               column = c("SYMBOL"),
                               keytype = "ENSEMBL",
                               multiVals = "first")

        
        res1_gene <- data.frame(Name_Gene = res1_gene$gene_id,Name_Uniprot = uniprot_names, Name_Symbol = symbol_names, lr=res1_gene$lr,
                           df = res1_gene$df, pvalue = res1_gene$pvalue, adj_pvalue = res1_gene$adj_pvalue)
        
        #First part - StageR significant genes and transcripts extraction - according to vignette
        
        pScreen <- DRIMSeq::results(drimseq_results1)$pvalue
        strp <- function(x) substr(x,1,15)
        names(pScreen) <- strp(DRIMSeq::results(drimseq_results1)$gene_id)
        pConfirmation <- matrix(DRIMSeq::results(drimseq_results1, level = "feature")$pvalue, ncol = 1)
        rownames(pConfirmation) <- strp(DRIMSeq::results(drimseq_results1, level = "feature")$feature_id)
        tx2gene1 <- DRIMSeq::results(drimseq_results1, level = "feature")[, c("feature_id", "gene_id")]
        for (i in 1:2) tx2gene1[,i] <- strp(tx2gene1[,i])
        stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                              pScreenAdjusted = FALSE, tx2gene = tx2gene1)
        stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                         alpha = 0.05, allowNA = TRUE)
        padj_stageRObj <- getAdjustedPValues(stageRObj, order = TRUE,
                                   onlySignificantGenes = TRUE)
        
        uniprot_names <- mapIds(org.Hs.eg.db, 
                                keys = gsub("\\..*","",padj_stageRObj$geneID),
                                column = c("UNIPROT"),
                                keytype = "ENSEMBL",
                                multiVals = "first")
        
        symbol_names <- mapIds(org.Hs.eg.db, 
                               keys = gsub("\\..*","",padj_stageRObj$geneID),
                               column = c("SYMBOL"),
                               keytype = "ENSEMBL",
                               multiVals = "first")
        
        res1_stageR <- data.frame(Name_gene = padj_stageRObj$geneID, Name_Transcript = padj_stageRObj$txID, Name_Uniprot = uniprot_names, Name_Symbol = symbol_names,
                                     padj_gene = padj_stageRObj$gene, padj_transcript = padj_stageRObj$transcript)
        
        #Second contrast - feature level
        con2 <- as.numeric(rep(0,length(colnames(design_full))))
        con2[k] <- 1
        con2[kk] <- -1
        drimseq_results2 <- dmTest(drimseq_results, contrast = con2 )
        res2 <- DRIMSeq::results(drimseq_results2,level='feature')
        
        uniprot_names <- mapIds(org.Hs.eg.db, 
                                keys = gsub("\\..*","",res2$gene_id),
                                column = c("UNIPROT"),
                                keytype = "ENSEMBL",
                                multiVals = "first")
        
        symbol_names <- mapIds(org.Hs.eg.db, 
                               keys = gsub("\\..*","",res2$gene_id),
                               column = c("SYMBOL"),
                               keytype = "ENSEMBL",
                               multiVals = "first")
        
        a <- DRIMSeq::proportions(drimseq_results2)
        ak <- a[,grep( x = colnames(a), pattern = analysis_names[k], ignore.case =TRUE)]
        ak <- rowMeans(ak)
        akk <- a[,grep( x = colnames(a), pattern = analysis_names[kk], ignore.case =TRUE)]
        akk <- rowMeans(akk)
        log2_proprotions <- log2(ak/akk)
        a <- DRIMSeq::coefficients(drimseq_results2,level='feature')
        ak <- a[,grep( x = colnames(a), pattern = analysis_names[k], ignore.case =TRUE)]
        akk <- a[,grep( x = colnames(a), pattern = analysis_names[kk], ignore.case =TRUE)]
        log2_coefficents = log2(exp(ak-akk))
        
        res2 <- data.frame(Name_Gene = res2$gene_id, Name_Transcript = res2$feature_id,Name_Uniprot = uniprot_names, Name_Symbol = symbol_names, lr=res2$lr,
                           df = res2$df, pvalue = res2$pvalue, adj_pvalue = res2$adj_pvalue, log2_Proportions = log2_proprotions, log2_Coefficents = log2_coefficents)
        
        #Second contrast - gene level
        res2_gene <- DRIMSeq::results(drimseq_results2,level='gene')
        
        uniprot_names <- mapIds(org.Hs.eg.db, 
                                keys = gsub("\\..*","",res2_gene$gene_id),
                                column = c("UNIPROT"),
                                keytype = "ENSEMBL",
                                multiVals = "first")
        
        symbol_names <- mapIds(org.Hs.eg.db, 
                               keys = gsub("\\..*","",res2_gene$gene_id),
                               column = c("SYMBOL"),
                               keytype = "ENSEMBL",
                               multiVals = "first")
        
      

        res2_gene <- data.frame(Name_Gene = res2_gene$gene_id,Name_Uniprot = uniprot_names, Name_Symbol = symbol_names, lr=res2_gene$lr,
                                df = res2_gene$df, pvalue = res2_gene$pvalue, adj_pvalue = res2_gene$adj_pvalue)
        
        #Second part - StageR significant genes and transcripts extraction
        
        pScreen <- DRIMSeq::results(drimseq_results2)$pvalue
        strp <- function(x) substr(x,1,15)
        names(pScreen) <- strp(DRIMSeq::results(drimseq_results2)$gene_id)
        pConfirmation <- matrix(DRIMSeq::results(drimseq_results2, level = "feature")$pvalue, ncol = 1)
        rownames(pConfirmation) <- strp(DRIMSeq::results(drimseq_results2, level = "feature")$feature_id)
        tx2gene1 <- DRIMSeq::results(drimseq_results2, level = "feature")[, c("feature_id", "gene_id")]
        for (i in 1:2) tx2gene1[,i] <- strp(tx2gene1[,i])
        stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation,
                              pScreenAdjusted = FALSE, tx2gene = tx2gene1)
        stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu",
                                         alpha = 0.05, allowNA = TRUE)
        padj_stageRObj <- getAdjustedPValues(stageRObj, order = TRUE,
                                             onlySignificantGenes = TRUE)
        
        uniprot_names <- mapIds(org.Hs.eg.db, 
                                keys = gsub("\\..*","",padj_stageRObj$geneID),
                                column = c("UNIPROT"),
                                keytype = "ENSEMBL",
                                multiVals = "first")
        
        symbol_names <- mapIds(org.Hs.eg.db, 
                               keys = gsub("\\..*","",padj_stageRObj$geneID),
                               column = c("SYMBOL"),
                               keytype = "ENSEMBL",
                               multiVals = "first")
        
        res2_stageR <- data.frame(Name_gene = padj_stageRObj$geneID, Name_Transcript = padj_stageRObj$txID, Name_Uniprot = uniprot_names, Name_Symbol = symbol_names,
                                  padj_gene = padj_stageRObj$gene, padj_transcript = padj_stageRObj$transcript)
        
        
        
        names_sheets <- append(names_sheets,paste("Condition_",analysis_names[kk],"vs",analysis_names[k],sep=""))
        names_sheets <- append(names_sheets,paste("Condition_",analysis_names[k],"vs",analysis_names[kk],sep=""))
        xls_list <- append(xls_list,c(list(res1),list(res2)))
        xls_list_gene <- append(xls_list_gene,c(list(res1_gene),list(res2_gene)))
        xls_list_stageR <- append(xls_list_stageR,c(list(res1_stageR),list(res2_stageR)))
      }
    }
  }
  
  #names(xls_list) <- colnames(design_full)
  names(xls_list) <- names_sheets
  write_xlsx(xls_list,xls_name)#,overwrite=TRUE)
  
  names(xls_list_gene) <- names_sheets
  write_xlsx(xls_list_gene,xls_name_gene)#,overwrite=TRUE)
  
  names(xls_list_stageR) <- names_sheets
  write_xlsx(xls_list_stageR,xls_name_stageR)#,overwrite=TRUE)
  
  
  #DESeq2 analysis of the gene expression starting from transcripts
  
  data <- tximport(files, type="salmon", tx2gene=tx2gene)
  
  dds <- DESeqDataSetFromTximport(data,
                                     colData = samples,
                                     design = ~ condition)
  
  #Filter the data
  keep <- rowSums(counts(dds)) > 0
  dds <- dds[keep,]
  keep <- filterByExpr(dds, group = dds$condition, min.count = filter_threshold,min.total.count = filter_threshold,min.prop = 0.7,large.n = n.small)
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  analysis_names <- levels(factor(samples$condition))
  dir.create(paste0(getwd(),'/DGE_DESeq2'))
  xls_name <- paste(getwd(),'/DGE_DESeq2/FoldChange_DESeq2_Gene.xlsx',sep="")
  file.remove(xls_name)
  xls_list <- list()
  names_sheets <- c()
  
  #Write the results for each comparison group
  
  for (k in 1:length(analysis_names)){
    for (kk in (k+1):length(analysis_names)){
      if(kk<=length(analysis_names) & kk!=k){
        dds_res1 <- DESeq2::results(dds,contrast = c("condition",analysis_names[kk],analysis_names[k]))
        dds_res2 <- DESeq2::results(dds,contrast = c("condition",analysis_names[k],analysis_names[kk]))
        
        uniprot_names <- mapIds(org.Hs.eg.db, 
                                keys = as.character(gsub("\\..*","",dds_res2@rownames)),
                                column = c("UNIPROT"),
                                keytype = "ENSEMBL",
                                multiVals = "first")
        
        symbol_names <- mapIds(org.Hs.eg.db, 
                               keys = as.character(gsub("\\..*","",dds_res2@rownames)),
                               column = c("SYMBOL"),
                               keytype = "ENSEMBL",
                               multiVals = "first")
        
        dds_res1 = data.frame(Name_Gene = dds_res2@rownames, Name_Ensembl = dds_res1@rownames,Name_Uniprot = uniprot_names, Name_Symbol = symbol_names, baseMean=dds_res1@listData$baseMean,
                              log2FoldChange = dds_res1@listData$log2FoldChange,
                              lfcSE = dds_res1@listData$lfcSE, stat = dds_res1@listData$stat,
                              pvalue = dds_res1@listData$pvalue, padj = dds_res1@listData$padj)
        xls_list <- append(xls_list,list(dds_res1))
        

        
        uniprot_names <- mapIds(org.Hs.eg.db, 
                                keys = as.character(gsub("\\..*","",dds_res2@rownames)),
                                column = c("UNIPROT"),
                                keytype = "ENSEMBL",
                                multiVals = "first")
        
        symbol_names <- mapIds(org.Hs.eg.db, 
                               keys = as.character(gsub("\\..*","",dds_res2@rownames)),
                               column = c("SYMBOL"),
                               keytype = "ENSEMBL",
                               multiVals = "first")
        
        dds_res2 = data.frame(Name_Gene = dds_res2@rownames, Name_Ensembl = dds_res2@rownames,Name_Uniprot = uniprot_names, Name_Symbol = symbol_names, baseMean=dds_res2@listData$baseMean,
                              log2FoldChange = dds_res2@listData$log2FoldChange,
                              lfcSE = dds_res2@listData$lfcSE, stat = dds_res2@listData$stat,
                              pvalue = dds_res2@listData$pvalue, padj = dds_res2@listData$padj)
        xls_list <- append(xls_list,list(dds_res2))
        
        names_sheets <- append(names_sheets,paste("Condition_",analysis_names[kk],"vs",analysis_names[k],sep=""))
        names_sheets <- append(names_sheets,paste("Condition_",analysis_names[k],"vs",analysis_names[kk],sep=""))
      }
    }
  }
  
  names(xls_list) <- names_sheets
  write_xlsx(xls_list,xls_name)#,overwrite=TRUE)
  
  
}

}