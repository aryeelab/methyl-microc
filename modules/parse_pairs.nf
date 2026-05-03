process PARSE_PAIRS {

    input:
    tuple val(sample_id), val(chunk_id), path(bam)
    path fai

    output:
    tuple val(sample_id), val(chunk_id), path("${sample_id}.${chunk_id}.sorted.pairs.gz"), emit: pairs

    script:
    def total_cpus    = task.cpus as int
    def parse_threads = Math.max(1, total_cpus.intdiv(2))
    def sort_threads  = Math.max(1, total_cpus.intdiv(2))

    """
    pairtools parse \\
        --min-mapq 30 \\
        --walks-policy 5unique \\
        --max-inter-align-gap 30 \\
        --drop-sam \\
        --add-columns pos5,pos3,cigar,seq \\
        --nproc-in ${parse_threads} \\
        --nproc-out ${parse_threads} \\
        --chroms-path $fai \\
        $bam | \\
        pairtools sort --nproc ${sort_threads} -o ${sample_id}.${chunk_id}.sorted.pairs.gz
    """
}
