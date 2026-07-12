process MERGE_DEDUP_PAIRS {

    tag "${sample_id}"

    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    tuple val(sample_id), path(pair_chunks)

    output:
    tuple val(sample_id), path("${sample_id}.pairs.gz"), emit: pairs
    tuple val(sample_id), path("${sample_id}.stats.txt"), emit: stats

    script:
    def chunk_list = pair_chunks.collect { it.name }.join(' ')

    """
    pairtools merge ${chunk_list} | \
        pairtools dedup \
            -o ${sample_id}.pairs.gz \
            --output-stats ${sample_id}.stats.txt
    """
}
