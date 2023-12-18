params.reads = "/Users/sbalgobin/Documents/nextflow-rnaseq-pipeline/data/SSR24292984_{1,2}.fq"
params.transcriptome_file = "/Users/sbalgobin/Documents/nextflow-rnaseq-pipeline/data/SRR27213070.fasta"
params.multiqc = "/Users/sbalgobin/Documents/nextflow-rnaseq-pipeline/multiqc"
params.outdir = "results"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)

/*
 * define the INDEX process that creates a binary index
 * given the transcriptome file
 */
process index {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

workflow {
    index_ch = index(params.transcriptome_file)
}

index_ch.view()

Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .into { read_pairs_ch; read_pairs_ch2 }

read_pairs_ch.view()

process quantification {

    tag "$pair_id"
    publishDir params.outdir, mode:'copy'

    input:
    path index from index_ch
    tuple pair_id, path(reads) from read_pairs_ch

    output:
    path pair_id into quant_ch

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

process fastqc {
    tag "FASTQC on $pair_id"

    input:
    tuple pair_id, path(reads) from read_pairs_ch2

    output:
    path "fastqc_${pair_id}_logs" into fastqc_ch

    script:
    """
    mkdir fastqc_${pair_id}_logs
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${reads}
    """
}

multiqc_ch_input = quant_ch.mix(fastqc_ch).collect()

process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    path '*' from multiqc_ch_input

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}
