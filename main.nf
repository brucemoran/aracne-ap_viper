#!/usr/bin/env nextflow

/* Run ARACNe-AP and msViper
*/

params.help = ""

if (params.help) {
  log.info ''
  log.info '------------------------------'
  log.info 'NEXTFLOW ARACNe-AP AND msVIPER'
  log.info '------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run main.nf \\'
  log.info '  --inputCsv data/input.csv \\'
  log.info '  --metaData data/metadata.csv \\'
  log.info '  -profile standard,singularity'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '  --inputCsv     STRING      CSV format, requires header: \'expmat,regulators\'; these 3 entries are input files for ARACNe; exprmat is TSV format, row1 = geneID,sample_1,...,sample_n; rows2..-1 = geneID_1,values_1,...,values_n; regulators has header \'geneID\' as per exprmat row1, and one geneID per line for genes of interest (TFs, DEGs etc.)'
  log.info '  --metaData     STRING      CSV format with header: \'sample,group\'; group relates group information used to discriminate between samples for msViper analysis'
  log.info ''
  log.info 'Optional arguments:'
  log.info '  --tag     STRING      label for the run, used to ID files; baseDir name used if not specified here'
  exit 1
}

//tag
if(!params.tag){
  params.tag = file("${workflow.launchDir}").getBaseName()
}

/* 0.0: Print script, configs used to run to pipeline_info for posterity
*/
WRITEOUT = Channel.fromPath("$baseDir/{main.nf,nextflow.config,conf/*}")
process writeFiles {

  publishDir "pipeline_info", mode: "copy", pattern: "*"

  input:
  file(writefile) from WRITEOUT

  output:
  file('*') into written

}

/* 1.0: calculate threshold
*/
Channel.fromPath("$params.inputCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [file(row.exprmat), file(row.regulators)]
       .into { calcthresh; viperdat }

process calcThresh {

  publishDir "analysis/ARACNe", mode: "copy", pattern: "*"

  input:
  set file(exprmat), file(regulators) from calcthresh

  output:
  set file(exprmat), file(regulators), file("*.txt") into bootstrap

  script:
  """
  JVMEM=\$(echo ${task.memory} | cut -d " " -f 1)"G"
  java -jar /usr/local/ARACNe-AP/dist/aracne.jar \
            -e $exprmat\
            -o ./ \
            --tfs $regulators \
            -p 1E-8 \
            --seed 1 \
            --calculateThreshold
  """
}

/* 2.0: bootstraps
*/
process bootStrap {

  input:
  set file(exprmat), file(regulators), file(threshold) from bootstrap
  each s from 1..100

  output:
  set file("*.txt") into consolidate

  script:
  """
  JVMEM=\$(echo ${task.memory} | cut -d " " -f 1)"G"
  java -Xmx\$JVMEM -jar /usr/local/ARACNe-AP/dist/aracne.jar \
          -e $exprmat \
          -o ./ \
          --tfs $regulators \
          -p 1E-8 \
          --seed $s
  """
}

/* 3.0: consolidate bootstraps
*/
process consBoots {

  publishDir "analysis/ARACNe", mode: "copy", pattern: "*"

  input:
  set file(bootstraps) from consolidate.groupTuple()

  output:
  file("${params.tag}.network.txt") into aracne_network

  script:
  """
  JVMEM=\$(echo ${task.memory} | cut -d " " -f 1)"G"
  java -Xmx\$JVMEM -jar /usr/local/ARACNe-AP/dist/aracne.jar \
            -o ./ \
            --consolidate
  mv network.txt ${params.tag}".network.txt"
  """
}

/* 4.0: Viper
*/
METAFILE = Channel.fromPath("$params.metaData", type: 'file')

process viper {

  publishDir "analysis/msViper", mode: "copy", pattern: "*"

  input:
  set file(exprmat), file(regulators) from viperdat
  file(network) from aracne_network
  file(metadata) from METAFILE

  output:
  file('*') into msviper_complete

  script:
  """
  Rscript --vanilla \
    ${workflow.projectDir}/bin/run_viper.call.R \
    ${workflow.projectDir}/bin/run_viper.func.R \
    $network \
    $exprmat \
    $metadata \
    ${params.tag}
  """
}
