process READ_QC {

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "read_level_qc.html", emit: report

    script:
    """
    mkdir -p fastqc multiqc

    fastqc -t ${task.cpus} \
        -o fastqc \
        ${r1} ${r2}

    multiqc fastqc \
        -o multiqc \
        -n read_level_qc.html

    cp multiqc/read_level_qc.html read_level_qc.html
    """
}
