process METHYLDACKEL_EXTRACT {

    tag "${sample_id}"

    publishDir "${params.outdir}/methylDackelExtracts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai
    val methyl_quality_threshold

    output:
    tuple val(sample_id),
          path("${sample_id}_CHG.methylKit.gz"),
          path("${sample_id}_CHH.methylKit.gz"),
          path("${sample_id}_CpG.methylKit.gz"),
          emit: methylkits

    script:
    """
    set -euo pipefail

    echo "========================================"
    echo "MethylDackel extract"
    echo "Sample: ${sample_id}"
    echo "BAM: ${bam}"
    echo "Reference: ${fasta}"
    echo "Methyl quality threshold: ${methyl_quality_threshold}"
    echo "========================================"

    MethylDackel extract \
        --methylKit \
        -q ${methyl_quality_threshold} \
        --nOT 0,0,0,5 \
        --nOB 0,0,5,0 \
        -@ ${task.cpus} \
        --CHH \
        --CHG \
        -o ${sample_id} \
        ${fasta} \
        "${bam}"

    pigz -p ${task.cpus} \
        ${sample_id}_CpG.methylKit \
        ${sample_id}_CHG.methylKit \
        ${sample_id}_CHH.methylKit

    echo
    echo "Generated MethylDackel methylKit files:"
    ls -lh \
        ${sample_id}_CpG.methylKit.gz \
        ${sample_id}_CHG.methylKit.gz \
        ${sample_id}_CHH.methylKit.gz
    """
}
