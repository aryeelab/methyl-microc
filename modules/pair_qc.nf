process PAIR_QC {

    tag "${sample_id}"

    publishDir "${params.outdir}/qc/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(pairs)

    output:
    tuple val(sample_id), path("pair_level_qc.html"), emit: report

    script:
    """
    python ${projectDir}/bin/pair_level_qc.py \
        --pairs ${pairs} \
        --out pair_level_qc.html \
        --cis-distance 10000
    """
}
