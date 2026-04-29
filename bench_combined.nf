nextflow.enable.dsl=2

include { RUN_METHYLSEQ } from './benchmarks/o2/run_methylseq_index'
include { POSTPROCESS_PAIRS } from './benchmarks/o2/combined_postprocess'

params.outdir = "results"
params.fasta  = null
params.fai    = null
params.input  = null
params.bam    = null
params.bwameth_index = null
params.inner_methylseq_extra = null
params.skip_inner_multiqc = false
params.methylseq_runner = null
params.methylseq_config = null

workflow {
    if (params.bam) {
        Channel.fromPath(params.bam).set { bam_ch }
    } else {
        Channel.fromPath(params.input).set { input_ch }
        Channel.fromPath(params.fasta).set { fasta_ch }
        Channel.fromPath(params.methylseq_config ?: "${projectDir}/methylseq.config").set { config_ch }

        def rows = file(params.input).readLines()
        def fastq_paths = rows.drop(1)
                             .findAll { it.trim() }
                             .collectMany { line ->
                                 def cols = line.split(',')
                                 [ cols[1].trim(), cols[2].trim() ]
                             }
                             .unique()

        Channel.value(fastq_paths.collect { file(it) }).set { fastq_ch }
        RUN_METHYLSEQ(input_ch, fasta_ch, config_ch, fastq_ch)
        bam_ch = RUN_METHYLSEQ.out.bam
    }

    POSTPROCESS_PAIRS(bam_ch, file(params.fasta), file(params.fai))
}
