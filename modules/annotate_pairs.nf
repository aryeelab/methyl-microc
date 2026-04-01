process ANNOTATE_PAIRS {

    publishDir "${params.outdir}/pairs", mode: 'copy', pattern: "*.meth.pairs.gz"
    publishDir "${params.outdir}/pairs", mode: 'copy', pattern: "*.meth.pairs"

    input:
    path pairs
    path fasta
    path fai

    output:
    path "*.meth.pairs.gz", emit: meth_pairs
    path "*.meth.pairs",    emit: meth_pairs_plain

    script:
    def sample_id = pairs.name.replaceFirst(/\.pairs\.gz$/, '')

    """
    python ${projectDir}/bin/annotate_pairs_methylation.py \
        --input $pairs \
        --fasta $fasta \
        --fai $fai \
        --output ${sample_id}.meth.pairs.gz

    gunzip -c ${sample_id}.meth.pairs.gz > ${sample_id}.meth.pairs
    """
}
