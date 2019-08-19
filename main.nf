#!/usr/bin/env nextflow

/* Run ARACNe-AP/$tag
*/

params.help = ""

if (params.help) {
  log.info ''
  log.info '------------------------------'
  log.info 'NEXTFLOW ARACNe-AP AND msVIPER'
  log.info '------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run ARACNe-AP/$tag.msviper.conda.nf \
                      --inputCsv data/input.csv \
                      --metaData data/metadata.csv \
                      -c "RNAseq.nextflow.config" \
                      -with-report "report.html" \
                      -with-timeline "timeline.html"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --inputCsv     STRING      CSV format, requires headers: expmat,regulators,seed; these 3 entries are files in data/; input requirements: tag is ID; exprmat is TSV format; row 1 = gene, sample_1, ..., sample_n; rows +2 geneID, values_1, ..., values_n; regulators has header \'gene\' equivalent to header in exprmat, one ID per line for genes of interest (TFs, DEGs etc.)'
  log.info '    --metaData     STRING      CSV format file: tag (matches tag in inputCsv file),metafile; metafile is CSV format with headers: sampleID,Group; Group relates group information used to discriminate between samples for Viper analysis'
  log.info 'Optional arguments:'
  log.info ''
  exit 1
}

/* 0.0: Print script, configs used to run to pipeline_info for posterity
*/
WRITEOUT = Channel.fromPath("$baseDir/{main.nf,nextflow.config,conf/*}")
process writeFiles {

  publishDir "$baseDir/pipeline_info/", mode: "copy", pattern: "*"

  input:
  file(writefile) from WRITEOUT

  script:
  """
  """
}

/* 1.0: calculate threshold
*/
Channel.fromPath("$params.inputCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [row.tag, file(row.exprmat), file(row.regulators)] }
       .into { calcthresh; viperdat }

process calcThresh {

  publishDir "$baseDir/analysis/$tag", mode: "copy", pattern: "*"

  input:
  set val(tag), file(exprmat), file(regulators) from calcthresh

  output:
  set val(tag), file(exprmat), file(regulators), file("*.txt") into bootstrap

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
  set val(tag), file(exprmat), file(regulators), file(threshold) from bootstrap
  each s from 1..100

  output:
  set val(tag), file("*.txt") into consolidate

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

  publishDir "$baseDir/analysis/$tag", mode: "copy", pattern: "*"

  input:
  set val(tag), file(bootstraps) from consolidate.groupTuple()

  output:
  set val(tag), file("${tag}.network.txt") into viperrun

  script:
  """
  JVMEM=\$(echo ${task.memory} | cut -d " " -f 1)"G"
  java -Xmx\$JVMEM -jar /usr/local/ARACNe-AP/dist/aracne.jar \
            -o ./ \
            --consolidate
  mv network.txt $tag".network.txt"
  """
}

/* 4.0: Viper
*/
Channel.fromPath("$params.metaData", type: 'file')
       .splitCsv( header: true )
       .map { row -> [file(row.metafile)] }
       .set { metafiles }

viperdat.join(viperrun).set{ viperin }

process viper {

  publishDir "$baseDir/analysis/$tag", mode: "copy", pattern: "*"
  publishDir "$baseDir/analysis", mode: "copy", pattern: "*sig.tsv"

  input:
  set val(tag), file(exprmat), file(regulators), file(network) from viperin
  each file(metafile) from metafiles


  output:
  file('*') into complete

  script:
  """
  Rscript --vanilla \
    ${params.binDir}/run_viper.call.R \
    ${params.binDir}/run_viper.func.R \
    $network \
    $exprmat \
    $metafile \
    $tag
  """
}
