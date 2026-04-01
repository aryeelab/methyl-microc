process PARSE_PAIRS {

    publishDir "results/pairs", mode: 'copy'

    input:
    path bam
    path fai

    output:
    path "test_sample.pairs.gz", emit: pairs
    path "test_sample.stats.txt", emit: stats

    script:
    """
    pairtools parse \
        --min-mapq 30 \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --drop-sam \
        --add-columns pos5,pos3,cigar,seq \
        --nproc-in 8 \
        --nproc-out 8 \
        --chroms-path $fai \
        $bam | \
        pairtools sort --nproc 4 | \
        pairtools dedup -o test_sample.pairs.gz --output-stats test_sample.stats.txt
    """
}
