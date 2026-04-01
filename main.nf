nextflow.enable.dsl=2

include { RUN_METHYLSEQ } from './modules/run_methylseq'
include { PARSE_PAIRS } from './modules/parse_pairs'
include { ANNOTATE_PAIRS } from './modules/annotate_pairs'
include { VALIDATE_PAIRS } from './modules/validate_pairs'

params.outdir = "results"
params.fasta  = null
params.fai    = null
params.input  = null
params.bam    = null

workflow {

    if (params.bam) {
        Channel.fromPath(params.bam).set { bam_ch }
    } else {
        Channel.fromPath(params.input).set { input_ch }
        Channel.fromPath(params.fasta).set { fasta_ch }
        Channel.fromPath("${projectDir}/methylseq.config").set { config_ch }

        def rows = file(params.input).readLines()
        def fastq_paths = rows.drop(1)
                             .findAll { it.trim() }
                             .collectMany { line ->
                                 def cols = line.split(',')
                                 [ cols[1].trim(), cols[2].trim() ]
                             }
                             .unique()

        Channel
            .value(fastq_paths.collect { file(it) })
            .set { fastq_ch }

        RUN_METHYLSEQ(input_ch, fasta_ch, config_ch, fastq_ch)
        bam_ch = RUN_METHYLSEQ.out.bam
    }
    
    PARSE_PAIRS(bam_ch, file(params.fai))
    ANNOTATE_PAIRS(PARSE_PAIRS.out.pairs, file(params.fasta), file(params.fai))
    VALIDATE_PAIRS(ANNOTATE_PAIRS.out.meth_pairs, file(params.fasta))
}
