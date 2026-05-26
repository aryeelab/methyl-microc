process METHYLATION_QC {

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path meth_pairs

    output:
    path "methylation_level_qc.html", emit: report

    script:
    """
    python ${projectDir}/bin/methylation_qc.py \
        --pairs ${meth_pairs} \
        --out methylation_level_qc.html \
        --frag-len 150 \
        --max-data-lines 10000000
    """
}
