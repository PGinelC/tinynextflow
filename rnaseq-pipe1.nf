#!/usr/bin/env nextflow

/*
 * The following pipeline parameters specify the reference genomes
 * and read pairs and can be provided as command line options
 * Run with this command: $ nextflow run rnaseq-pipe1.nf -with-docker
 */
params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )

    INDEX(params.transcriptome)
    FASTQC(read_pairs_ch)
    QUANT(INDEX.out, read_pairs_ch)
}

process INDEX {
    container 'combinelab/salmon:latest'

    tag "$transcriptome.simpleName"

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

process FASTQC {
    container 'staphb/fastqc:0.11.9'

    tag "FASTQC on $sample_id"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process QUANT {
    container 'combinelab/salmon:latest' 

    tag "$pair_id"
    publishDir params.outdir

    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path pair_id

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}