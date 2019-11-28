#! /usr/bin/env R
options(scipen=999)

libs <- c("httr", "limma", "CAMERA", "DEFormats", "apeglm", "ggplot2", "tidyverse", "plyr", "dplyr", "biomaRt", "reshape2", "viper", "mixtools", "Biobase", "rgexf", "fgsea", "gtools", "tximport", "aracne.networks")
print("Loading libraries...")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = T)))})
print("Loaded")

runMsViper <- function(NETWORK, EXPRMAT, METADATA, TAG){

  regulonaracne <- suppressWarnings(aracne2regulon(afile = NETWORK, eset = EXPRMAT))

  ## create eset
  exprd <- read.table(EXPRMAT, row=1, header=T)
  pheno <- read_tsv(METADATA)
  if(sum(match(colnames(pheno),c("sampleID", "group")))!=3){
    print("Converting metadata colnames, ensure they are in sample, group order or msViper may fail!")
    colnames(pheno) <- c("sampleID", "group")
  }
  phenodf <- pheno %>%
            dplyr::mutate(sampleID = gsub("-",".",sampleID)) %>%
            dplyr::filter(sampleID %in% colnames(exprd)) %>%
            dplyr::rename(description=group) %>%
            as.data.frame()
  exprds <- exprd[,colnames(exprd) %in% phenodf$sampleID]
  phenos <- phenodf[phenodf$sampleID %in% colnames(exprds),]
  pData <- new("AnnotatedDataFrame",
                data=data.frame(row.names = colnames(exprds),
                                description = as.vector(phenos$description)))
  combined.eset <- ExpressionSet(assayData=as.matrix(exprds),
                                 phenoData=pData)

  ##generate signature for pairwise description levels
  LEVELS <- levels(pData$description)
  lapply(seq_along(LEVELS), function(f){

    level1 <- LEVELS[f]
    print(paste0("Working on: ", level1))

    ##require 10+ samples in Group
    if(table(pData$description %in% level1)[2] > 20){
      pairname <- paste0(gsub(" ","-",level1), "_vs_rest")
      signature <- viper::rowTtest(combined.eset, "description", level1)
      rownamesp <- rownames(signature$p.value)
      signature$p.value <- as.matrix(p.adjust(signature$p.value, method="BH"))
      rownames(signature$p.value) <- rownamesp
      signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
      nullmodel <- viper::ttestNull(combined.eset, "description", level1, LEVELS[! LEVELS %in% level1], per=5000, repos=TRUE, verbose=TRUE)

      ## msViper, leading edge analyses + bootstrapped
      mra <- msviper(signature, regulonaracne, nullmodel)
      mra <- viper::ledge(mra)
      mra <- viper::shadow(mra, regulators=0.01, shadow=0.01, verbose=FALSE, per=5000)
      mrac <- msviperCombinatorial(mra, regulators=0.01, verbose=FALSE)

      ##if statement to deny error from synergy analysis which otherwise exits from the script
      if(length(mrac$regulon)>length(mra$regulon)){
        mra <- msviperSynergy(mrac, verbose = FALSE)
        ##BH adjust p-values
        mra$es$q.value <- p.adjust(mra$es$p.value,method="BH")
      }
      if(length(mrac$regulon)==length(mra$regulon)){
        ##BH adjust p-values
        mra$es$q.value <- p.adjust(mra$es$p.value, method="BH")
      }

      ## outputs
      ##omit NaN as throws error
      mra$signature<-na.omit(mra$signature)
      pdf(paste0(TAG, ".", pairname, ".msViper.syn.results.pdf"))
        plot(mra, pval = mra$es$q.value, cex = .7)
      dev.off()

      mratab <- data.frame(NES=mra$es$nes,
                           Size=mra$es$size,
                           p.value=mra$es$p.value,
                           q.value=mra$es$q.value)

      write.table(mratab, paste0(TAG, ".", pairname, ".msViper.results.tsv"), sep="\t", quote=F, col.names=paste("gene",colnames(mratab), sep="\t"))

      retList <- list(mra, mratab)
      assignedName <- paste0(pairname)
      assign(assignedName, value=retList)
      saveFile <- paste0(TAG, ".", pairname, ".RData")
      save(list=assignedName, file=saveFile)
    }
  })
}

##filter out significant results (p-value)
filterSigRes <- function(RDATA, PVAL=NULL){

  if(is.null(PVAL)){
    PVAL <- 0.05
  }
  hmart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ens2ext <- as_tibble(getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), mart = hmart))

  load(RDATA)
  les <- ls()
  got <- get(les)
  signames <- names(got[[1]]$es$p.value[got[[1]]$es$p.value<PVAL])
  if(length(signames>0)){
    sigout <- ens2ext %>% dplyr::filter(ensembl_gene_id %in% signames | external_gene_name %in% signames)
    sigout$size <- got[[1]]$es$size[names(got[[1]]$es$size) %in% signames]
    sigout$pvalue <- got[[1]]$es$p.value[names(got[[1]]$es$p.value) %in% signames]
    sigout$qvalue <- got[[1]]$es$q.value[names(got[[1]]$es$q.value) %in% signames]
    outname <- gsub("RData", "sig.tsv", RDATA)
    write_tsv(sigout, path=outname)
  }}
