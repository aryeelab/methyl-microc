process METHYLATION_QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(meth_pairs)

    output:
    tuple val(sample_id), path("methylation_level_qc.html"), emit: report

    script:
    """
    python ${projectDir}/bin/methylation_qc.py \
        --pairs ${meth_pairs} \
        --out methylation_level_qc.html \
        --frag-len 150 \
        --max-data-lines 10000000
    """
}
