#! /usr/bin/env R
options(scipen=999)

libs <- c("sleuth", "httr", "limma", "CAMERA", "DEFormats", "apeglm", "ggplot2", "tidyverse", "plyr", "dplyr", "biomaRt", "reshape2", "viper", "mixtools", "Biobase", "rgexf", "fgsea", "gtools", "tximport", "aracne.networks")
print("Loading libraries...")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = T)))})
print("Loaded")

source("https://raw.githubusercontent.com/brucemoran/R/master/functions/expression/highest_exp_wide.func.R")

runMsViper <- function(NETWORK, EXPRMAT, METADATA, TAG){

  hmart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ens2ext <- as_tibble(getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), mart = hmart))
  save(ens2ext, file="ens2ext.RData")
  regulonaracne <- suppressWarnings(aracne2regulon(afile = NETWORK, eset = EXPRMAT))

  ## create eset
  exprd <- read_tsv(EXPRMAT)

  ##check geneID as exernal_gene_name, then convert to ensembl_gene_id
  exprdg <- left_join(ens2ext, exprd)
  exprdh <- highest_exp_wide(exprdg, "external_gene_name", "ensembl_gene_id", 3:dim(exprdg)[2])

  pheno <- read_tsv(METADATA)
  sumatch <- sum(match(colnames(pheno),c("sampleID", "group")))
  if(sumatch!=3 | is.na(sumatch)){
    print("Converting metadata colnames, ensure they are in sample, group order or msViper may fail!")
    colnames(pheno) <- c("sampleID", "group")
  }
  dotdash <- ifelse(length(strsplit(colnames(exprd)[3], "-")[[1]]),"-", ".")
  phenodf <- pheno %>%
            dplyr::mutate(sampleID = gsub("-", dotdash, sampleID)) %>%
            dplyr::filter(sampleID %in% colnames(exprd)) %>%
            dplyr::rename(description=group) %>%
            as.data.frame()
  exprds <- exprdh %>%
            dplyr::select(-external_gene_name) %>%
            as.data.frame() %>%
            column_to_rownames(., var = "ensembl_gene_id")
  phenos <- phenodf[phenodf$sampleID %in% colnames(exprds),]
  pData <- new("AnnotatedDataFrame",
                data=data.frame(row.names = colnames(exprds),
                                description = factor(phenos$description)))
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
      mra <- msviper(signature, regulonaracne, nullmodel, verbose=FALSE)
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
filterSigRes <- function(RDATA, ENS2EXT="ens2ext.RData", PVAL=NULL){

  load(ENS2EXT)

  if(is.null(PVAL)){
    PVAL <- 0.05
  }

  for(x in 1:length(RDATA)){
    print(paste0("Loading: ", RDATA[x]))
    load(RDATA[x])
  }

  les <- ls(pattern="vs")
  got <- get(les)
  names(got) <- les
  lapply(seq_along(got), function(ff){
  f <- got[[ff]]
  signames <- names(f$es$p.value[f$es$p.value<PVAL])
  if(length(signames>0)){
    sigout <- ens2ext %>% dplyr::filter(ensembl_gene_id %in% signames | external_gene_name %in% signames)
    sigout$size <- f$es$size[names(f$es$size) %in% signames]
    sigout$pvalue <- f$es$p.value[names(f$es$p.value) %in% signames]
    sigout$qvalue <- f$es$q.value[names(f$es$q.value) %in% signames]
    outname <- gsub("RData", "sig.tsv", grep("vs",RDATA,value=T)[ff])
    write_tsv(sigout, path=outname)
  }
  if(length(signames==0)){
    print(paste0("No significant genes found at p < ", PVAL))
  }
  })
}
