process VALIDATE_PAIRS {

    debug true

    input:
    path pairs
    path fasta

    output:
    stdout

    script:
    """
    python ${projectDir}/bin/validate_pairs_methylation.py \
        --pairs $pairs \
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
