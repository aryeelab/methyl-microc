nextflow.enable.dsl=2

include { SPLIT_FASTQ }        from './modules/split_fastq'
include { RUN_METHYLSEQ }      from './modules/run_methylseq'
include { PARSE_PAIRS }        from './modules/parse_pairs'
include { MERGE_DEDUP_PAIRS }  from './modules/merge_dedup_pairs'
include { ANNOTATE_PAIRS }     from './modules/annotate_pairs'
include { VALIDATE_PAIRS }     from './modules/validate_pairs'

params.outdir       = "results"
params.fasta        = null
params.fai          = null
params.input        = null
params.bam          = null
params.split_chunks = 4

workflow {

    if (params.bam) {
        Channel
            .fromPath(params.bam)
            .map { bam ->
                def sample_id = bam.name
                    .replaceFirst(/\.markdup\.sorted\.bam$/, '')
                    .replaceFirst(/\.bam$/, '')
                tuple(sample_id, 'merged', bam)
            }
            .set { bam_ch }

    } else {
        if (!params.input) {
            error "Please provide --input <samplesheet.csv>"
        }
        if (!params.fasta) {
            error "Please provide --fasta <reference.fa>"
        }

        Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    row['sample'].toString().trim(),
                    file(row['fastq_1'].toString().trim()),
                    file(row['fastq_2'].toString().trim())
                )
            }
            .set { samples_ch }

        SPLIT_FASTQ(samples_ch, params.split_chunks)

        chunk_fastq_ch = SPLIT_FASTQ.out.manifest
            .splitCsv(sep: '\t', header: false)
            .map { row ->
                tuple(
                    row[0].toString(),                  // sample_id
                    row[1].toString(),                  // chunk_id
                    file(row[2].toString()),            // chunk R1
                    file(row[3].toString())             // chunk R2
                )
            }

        RUN_METHYLSEQ(
            chunk_fastq_ch,
            file(params.fasta),
            file("${projectDir}/methylseq.config")
        )

        bam_ch = RUN_METHYLSEQ.out.bam
    }

    PARSE_PAIRS(bam_ch, file(params.fai))

    grouped_pairs_ch = PARSE_PAIRS.out.pairs
        .map { sample_id, chunk_id, pairs_file ->
            tuple(sample_id, pairs_file)
        }
        .groupTuple()

    MERGE_DEDUP_PAIRS(grouped_pairs_ch)

    ANNOTATE_PAIRS(
        MERGE_DEDUP_PAIRS.out.pairs,
        file(params.fasta),
        file(params.fai)
    )

    VALIDATE_PAIRS(
        ANNOTATE_PAIRS.out.meth_pairs,
        file(params.fasta)
    )
}
