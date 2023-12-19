params.reads = "/Users/sbalgobin/Documents/Nexflow-RNAseq-Pipeline/nextflow-rnaseq-pipeline/data/gingiva_{1,2}.fastq"
params.transcriptome_file = "/Users/sbalgobin/Documents/Nexflow-RNAseq-Pipeline/nextflow-rnaseq-pipeline/data/feline_mouth.fasta"
params.multiqc = "/Users/sbalgobin/Documents/Nextflow-RNAseq-Pipeline/nextflow-rnaseq-pipeline/multiqc"
params.outdir = "../results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)

/*
 * Define the process to create a binary index
 * given the transcriptome file
 */
process CreateIndex {
    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t ${params.transcriptome_file} -i salmon_index
    """
}

/*
 * Define the RNA quantification process
 */
process RNAQuantification {
    tag "$pair_id"
    publishDir params.outdir, mode:'copy'

    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path "${pair_id}_quant"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o ${pair_id}_quant
    """
}

/*
 * Define the FastQC process
 */
process RunFastQC {
    tag "FASTQC on $pair_id"

    input:
    tuple val(pair_id), path(reads)

    output:
    path "fastqc_${pair_id}_logs"

    script:
    """
    mkdir fastqc_${pair_id}_logs
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${reads}
    """
}

/*
 * Define the MultiQC process
 */
process GenerateMultiQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    index_ch = CreateIndex()

    quant_ch = RNAQuantification(index_ch, read_pairs_ch)
    fastqc_ch = RunFastQC(read_pairs_ch)

    multiqc_ch_input = quant_ch.mix(fastqc_ch).collect()
    GenerateMultiQC(multiqc_ch_input)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
