process PARSE_PAIRS {

    input:
    tuple val(sample_id), val(chunk_id), path(bam)
    path fai

    output:
    tuple val(sample_id), val(chunk_id), path("${sample_id}.${chunk_id}.sorted.pairs.gz"), emit: pairs

    script:
    def total_cpus    = task.cpus as int
    def sort_bam_threads = Math.max(1, total_cpus.intdiv(2))
    def parse_threads    = Math.max(1, total_cpus.intdiv(4))
    def sort_pair_threads = Math.max(1, total_cpus.intdiv(4))

    """
    samtools sort -n \
        -@ ${sort_bam_threads} \
        -o ${sample_id}.${chunk_id}.namesorted.bam \
        $bam

    samtools view -H ${sample_id}.${chunk_id}.namesorted.bam | grep '^@HD'

    pairtools parse \
        --min-mapq 30 \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --drop-sam \
        --add-columns pos5,pos3,cigar,seq \
        --nproc-in ${parse_threads} \
        --nproc-out ${parse_threads} \
        --chroms-path $fai \
        ${sample_id}.${chunk_id}.namesorted.bam | \
        pairtools sort --nproc ${sort_pair_threads} -o ${sample_id}.${chunk_id}.sorted.pairs.gz
    """
}
