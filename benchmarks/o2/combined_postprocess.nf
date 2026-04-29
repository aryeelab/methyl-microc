process POSTPROCESS_PAIRS {

    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    path bam
    path fasta
    path fai

    output:
    path "*.pairs.gz", emit: pairs
    path "*.stats.txt", emit: stats
    path "*.meth.pairs.gz", emit: meth_pairs
    path "*.meth.pairs", emit: meth_pairs_plain

    script:
    def sample_id = bam.name
        .replaceFirst(/\.markdup\.sorted\.bam$/, '')
        .replaceFirst(/\.bam$/, '')

    def total_cpus = task.cpus as int
    def parse_threads = Math.max(1, total_cpus.intdiv(2))
    def sort_threads = Math.max(1, total_cpus.intdiv(4))

    """
    pairtools parse \
        --min-mapq 30 \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --drop-sam \
        --add-columns pos5,pos3,cigar,seq \
        --nproc-in ${parse_threads} \
        --nproc-out ${parse_threads} \
        --chroms-path $fai \
        $bam | \
        pairtools sort --nproc ${sort_threads} | \
        pairtools dedup -o ${sample_id}.pairs.gz --output-stats ${sample_id}.stats.txt

    python ${projectDir}/bin/annotate_pairs_methylation.py \
        --input ${sample_id}.pairs.gz \
        --fasta $fasta \
        --fai $fai \
        --output ${sample_id}.meth.pairs.gz

    gunzip -c ${sample_id}.meth.pairs.gz > ${sample_id}.meth.pairs

    python ${projectDir}/bin/validate_pairs_methylation.py \
        --pairs ${sample_id}.meth.pairs.gz \
        --fasta $fasta \
        --record 3 \
        >/dev/null 2>&1 \
    || {
        echo "VALIDATION FAILED !!!"
        echo "Please check the input pairs file and methylation annotation step."
        exit 0
    }
    """
}
