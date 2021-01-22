#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  ------------------------------------------------------------------------------
                        NEXTFLOW ARACNe-AP AND msVIPER
  ------------------------------------------------------------------------------
  Usage:

  nextflow run brucemoran/aracne-ap

  Mandatory arguments:

    -profile        [str]       Configuration profile (required: singularity)
    --exprmat       [file]      expression matrix in tab-separated format:
                                header line in format:
                                geneID sample_1 ... sampleID_n
                                all other lines in format
                                geneID_1 values_1 ... values_n etc.
    --regulators    [file]      single geneID per line for genes of interest
                                acceptable geneIDs are ensembl_gene_id or
                                external_gene_name from Ensembl;
                                header line indicates name of geneIDs
    --metaData      [str]       metadata to discriminate between cohorts for
                                msViper analysis; tab-separated format;
                                header line format: sample group
                                all other lines: sampleID_1 group_0 etc.


  Optional arguments:

    --bootstraps     [int]       specify number of bootstraps (default: 100)
    --runID          [str]       label for the run, used to ID files (default: aracne)
    --outDir         [str]       output directory (default: ARACNe)
  """.stripIndet()
}

//Java task memory allocation via task.memory
javaTaskmem = { it.replace(" GB", "g") }

// 1.0: calculate threshold
Channel.fromPath("$params.exprmat", type: 'file', checkIfExists: true).set { exprmat_ch }
Channel.fromPath("$params.regulators", type: 'file', checkIfExists: true).set { regulators_ch }

process calc_thresh {

  publishDir "${params.outDir}/setup", mode: "copy", pattern: "*"

  input:
  file(exprmat) from exprmat_ch
  file(regulators) from regulators_ch

  output:
  tuple file(exprmat), file(regulators), file("${params.runID}.miThreshold_p1E-8.txt") into bootstrap

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  java ${taskmem} -jar /opt/miniconda/envs/jupyter_rnaseq/bin/aracne.jar \
    -e ${exprmat}\
    -o ${params.runID}.miThreshold_p1E-8.txt \
    --tfs ${regulators} \
    -p 1E-8 \
    --seed 843892 \
    --calculateThreshold
  """
}

// 2.0: bootstraps

process boot_strap {

  input:
  tuple file(exprmat), file(regulators), file(threshold) from bootstrap
  each s from 1..params.bootstraps

  output:
  file("${s}.bootstrap.txt") into consolidate

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  java ${taskmem} -jar /opt/miniconda/envs/jupyter_rnaseq/bin/aracne.jar \
    -e ${exprmat} \
    -o ${s}.bootstrap.txt \
    --tfs ${regulators} \
    -p 1E-8 \
    --seed ${s}
  """
}

/* 3.0: consolidate bootstraps
*/
process cons_boots {

  publishDir "${params.outDir}/setup", mode: "copy"

  input:
  file(bootstraps) from consolidate.collect()

  output:
  file("${params.tag}.network.txt") into aracne_network

  script:
  def taskmem = task.memory == null ? "" : "-Xmx" + javaTaskmem("${task.memory}")
  """
  java ${taskmem} -jar /opt/miniconda/envs/jupyter_rnaseq/bin/aracne.jar \
    -o ./ \
    --consolidate
  mv network.txt ${params.runID}".network.txt"
  """
}

// 4.0: Viper

Channel.fromPath("$params.metaData", type: 'file', checkIfExists: true).set { metafile_ch }

process viper {

  publishDir "${params.outDir}/msViper", mode: "copy"

  input:
  tuple file(exprmat), file(regulators) from viperdat
  file(network) from aracne_network
  file(metadata) from metafile_ch

  output:
  file('*') into msviper_complete

  script:
  """
  Rscript --vanilla \
    ${workflow.projectDir}/bin/run_viper.call.R \
    ${workflow.projectDir}/bin/run_viper.func.R \
    ${network} \
    ${exprmat} \
    ${metadata} \
    ${params.tag}
  """
}
