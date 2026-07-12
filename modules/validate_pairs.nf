process VALIDATE_PAIRS {

    tag "${sample_id}"

    debug true

    input:
    tuple val(sample_id), path(pairs)
    path fasta

    output:
    tuple val(sample_id), stdout, emit: validation

    script:
    """
    python ${projectDir}/bin/validate_pairs_methylation.py \
        --pairs ${pairs} \
        --fasta ${fasta} \
        --record 3 \
        >/dev/null 2>&1 \
    || {
        echo "VALIDATION FAILED for sample: ${sample_id}"
        echo "Please check the input pairs file and methylation annotation step."
        exit 0
    }

    echo "VALIDATION PASSED for sample: ${sample_id}"
    """
}
