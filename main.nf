nextflow.enable.dsl=2

include { READ_QC }            from './modules/read_qc'
include { SPLIT_FASTQ }        from './modules/split_fastq'
include { RUN_METHYLSEQ }      from './modules/run_methylseq'
include { MERGE_GLOBAL_BAM }   from './modules/merge_global_bam'
include { METHYLDACKEL_EXTRACT } from './modules/methyldackel_extract'
include { METHYLDACKEL_MBIAS }   from './modules/methyldackel_mbias'
include { PARSE_PAIRS }        from './modules/parse_pairs'
include { MERGE_DEDUP_PAIRS }  from './modules/merge_dedup_pairs'
include { BUILD_MATRICES }      from './modules/build_matrices'
include { PAIR_QC }            from './modules/pair_qc'
include { ANNOTATE_PAIRS }     from './modules/annotate_pairs'
include { METHYLATION_QC }     from './modules/methylation_qc'
include { VALIDATE_PAIRS }     from './modules/validate_pairs'

params.outdir          = "results"
params.fasta           = null
params.fai             = null
params.input           = null
params.bam             = null
params.reads_per_chunk       = 20_000_000
params.matrix_base_resolution = 1000
params.matrix_resolutions     = "nice"
params.methyl_quality_threshold = 20

workflow {

    if (params.bam) {

        /*
         * BAM-input mode:
         * treat the provided BAM as a single already-prepared BAM.
         * Do not create another merged/global BAM output.
         */
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

        /*
         * Read-level QC runs once per biological sample.
         */
        READ_QC(samples_ch)

        /*
         * FASTQs are split into artificial chunks for parallel processing.
         */
        SPLIT_FASTQ(samples_ch, params.reads_per_chunk)

        chunk_fastq_ch = SPLIT_FASTQ.out.manifest
            .splitCsv(sep: '\t', header: false)
            .map { row ->
                tuple(
                    row[0].toString(),
                    row[1].toString(),
                    file(row[2].toString()),
                    file(row[3].toString())
                )
            }

        /*
         * Each FASTQ chunk is processed independently by nf-core/methylseq.
         * Output:
         * tuple(sample_id, chunk_id, chunk.markdup.sorted.bam)
         */
        RUN_METHYLSEQ(
            chunk_fastq_ch,
            file(params.fasta),
            file("${projectDir}/methylseq.config")
        )

        bam_ch = RUN_METHYLSEQ.out.bam

        /*
         * Build one reusable BAM per biological sample.
         *
         * Remove chunk_id, then group only BAMs sharing the same sample_id.
         * This keeps multiple samples completely separate.
         */
        grouped_bam_ch = RUN_METHYLSEQ.out.bam
            .map { sample_id, chunk_id, bam ->
                tuple(sample_id, bam)
            }
            .groupTuple()

        /*
         * Per sample:
         *   merge chunks
         *   coordinate sort
         *   clear chunk-level duplicate marks
         *   normalize RG/library
         *   global Picard MarkDuplicates
         *   index
         *
         * Published to:
         *   results/bam/<sample>.markdup.sorted.bam
         *   results/bam/<sample>.markdup.sorted.bam.bai
         */
        MERGE_GLOBAL_BAM(grouped_bam_ch)

    /*
     * MethylDackel intermediate outputs.
     *
     * Side branch from the final globally merged and
     * globally duplicate-marked BAM.
     *
     * Existing pairs, matrices, and methylation branches
     * remain unchanged.
     */
    METHYLDACKEL_EXTRACT(
        MERGE_GLOBAL_BAM.out.bam,
        file(params.fasta),
        file(params.fai),
        params.methyl_quality_threshold
    )

    METHYLDACKEL_MBIAS(
        MERGE_GLOBAL_BAM.out.bam,
        file(params.fasta),
        file(params.fai)
    )
    }

    /*
     * Existing pair pipeline is unchanged.
     * Each BAM chunk is name-sorted inside PARSE_PAIRS before pairtools parse.
     */
    PARSE_PAIRS(
        bam_ch,
        file(params.fai)
    )

    /*
     * Group parsed pair chunks by biological sample.
     */
    grouped_pairs_ch = PARSE_PAIRS.out.pairs
        .map { sample_id, chunk_id, pairs_file ->
            tuple(sample_id, pairs_file)
        }
        .groupTuple()

    /*
     * Merge all pair chunks for each sample and perform global pair dedup.
     */
    MERGE_DEDUP_PAIRS(grouped_pairs_ch)

    /*
     * Build contact matrices from the final deduplicated pairs.
     *
     * Existing methylation processing is unchanged.
     *
     * Outputs per sample:
     *   results/matrices/<sample>/<sample>.cool
     *   results/matrices/<sample>/<sample>.mcool
     *   results/matrices/<sample>/<sample>.hic
     */
    BUILD_MATRICES(
        MERGE_DEDUP_PAIRS.out.pairs,
        file(params.fai),
        params.matrix_base_resolution,
        params.matrix_resolutions
    )

    /*
     * Pair-level QC.
     */
    PAIR_QC(
        MERGE_DEDUP_PAIRS.out.pairs
    )

    /*
     * Methylation annotation.
     */
    ANNOTATE_PAIRS(
        MERGE_DEDUP_PAIRS.out.pairs,
        file(params.fasta),
        file(params.fai)
    )

    /*
     * Methylation-level QC and validation.
     */
    METHYLATION_QC(
        ANNOTATE_PAIRS.out.meth_pairs
    )

    VALIDATE_PAIRS(
        ANNOTATE_PAIRS.out.meth_pairs,
        file(params.fasta)
    )
}
