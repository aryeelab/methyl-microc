process PAIR_QC {

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path pairs

    output:
    path "pair_level_qc.html", emit: report

    script:
    """
    python ${projectDir}/bin/pair_level_qc.py \
        --pairs ${pairs} \
        --out pair_level_qc.html \
        --cis-distance 10000
    """
}
