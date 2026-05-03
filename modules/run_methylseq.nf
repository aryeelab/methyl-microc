process RUN_METHYLSEQ {

    input:
    tuple val(sample_id), val(chunk_id), path(r1), path(r2)
    path fasta
    path methyl_config

    output:
    tuple val(sample_id), val(chunk_id), path("results/**/deduplicated/*.markdup.sorted.bam"), emit: bam

    script:
    """
    cat > samplesheet.csv <<EOF
sample,fastq_1,fastq_2
${sample_id}_${chunk_id},${r1.name},${r2.name}
EOF

    ls -lh

    bash ${projectDir}/run_methyl_microc.sh \\
        --input samplesheet.csv \\
        --outdir results \\
        --fasta ${fasta.name}
    """
}
