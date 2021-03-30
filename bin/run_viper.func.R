#! /usr/bin/env R
options(scipen=999)

# libs <- c("sleuth", "httr", "limma", "CAMERA", "DEFormats", "apeglm", "ggplot2", "tidyverse", "plyr", "dplyr", "biomaRt", "reshape2", "viper", "mixtools", "Biobase", "rgexf", "fgsea", "gtools", "tximport", "aracne.networks")
# print("Loading libraries...")
# libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = T)))})
# print("Loaded")

# source("https://raw.githubusercontent.com/brucemoran/R/master/functions/expression/highest_exp_wide.func.R")

#
# @param NETWORK ARACNe_AP network file
# @param EXPRMAT expression matrix in tab-separated format:
#                header line in format:
#                geneID sample_1 ... sampleID_n
#                all other lines in format
#                geneID_1 value_1 ... value_n etc.
# @param METADATA metadata to discriminate between cohorts for msViper analysis;
#                 tab-separated format; header line format:
#                 sample group
#                 all other lines:
#                 sampleID_1 group_0 etc.
# @param ENS_DATASET biomaRt::useEnsembl(dataset = <this>) default: human
# @param TAG string identifier of the run (default: "aracne")
# @return
#
parse_inputs <- function(NETWORK, EXPRMAT, METADATA, TAG, ENS_DATASET = "hsapiens_gene_ensembl"){

  ##process annotation data
  if(!file.exists("ens2ext.RData")){
    print("Annotating with biomaRt...")
    hmart <- biomaRt::useEnsembl(biomart = "ensembl",
                                 dataset = ENS_DATASET)
    ens2ext <- as_tibble(biomaRt::getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                                        mart = hmart))

    save(ens2ext, file="ens2ext.RData")
  } else {
    print("Loading annotations from biomaRt...")
    load("ens2ext.RData")
  }

  print("ARACNe-AP to regulon...")
  regulonaracne <- suppressWarnings(viper::aracne2regulon(afile = NETWORK,
                                                          eset = EXPRMAT))

  ## create eset
  print("Parsing expression matrix...")
  exprd <- readr::read_tsv(EXPRMAT)

  ##check geneID as exernal_gene_name, then convert to ensembl_gene_id
  exprdg <- dplyr::left_join(ens2ext, exprd)
  exprdh <- highest_exp_wides(wide_object = exprdg,
                             name_col = "external_gene_name",
                             id_col = "ensembl_gene_id",
                             value_cols = 3:dim(exprdg)[2])

  print("Processing metadata...")
  pheno <- readr::read_tsv(METADATA)
  sumatch <- sum(match(colnames(pheno), c("sampleID", "group")))

  if(sumatch != 3 | is.na(sumatch)){
    print("Converting metadata colnames, ensure they are in sample, group order or msViper may fail!")
    colnames(pheno) <- c("sampleID", "group")
  }

  dot_dash <- ifelse(length(strsplit(colnames(exprd)[3], "-")[[1]]),"-", ".")
  phenodf <- pheno %>%
             dplyr::mutate(sampleID = gsub("-", dot_dash, sampleID)) %>%
             dplyr::filter(sampleID %in% colnames(exprd)) %>%
             dplyr::rename(description = group) %>%
             as.data.frame()
  exprds <- exprdh %>%
            dplyr::select(-external_gene_name) %>%
            as.data.frame() %>%
            column_to_rownames(., var = "ensembl_gene_id")
  phenos <- phenodf[phenodf$sampleID %in% colnames(exprds),]
  print("Creating ExpressionSet...")
  pData <- new("AnnotatedDataFrame",
                data=data.frame(row.names = colnames(exprds),
                                description = factor(phenos$description)))
  combined.eset <- Biobase::ExpressionSet(assayData = as.matrix(exprds),
                                          phenoData = pData)
  print("Saving...")
  save(combined.eset, phenos, pData, exprds, exprdh, regulonaracne, file = paste0(TAG,".parse_inputs.RData"))
}
