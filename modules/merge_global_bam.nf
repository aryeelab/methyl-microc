process MERGE_GLOBAL_BAM {

    tag "${sample_id}"

    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_chunks)

    output:
    tuple val(sample_id),
          path("${sample_id}.markdup.sorted.bam"),
          path("${sample_id}.markdup.sorted.bam.bai"),
          emit: bam

    script:
    def bam_files = bam_chunks instanceof List ? bam_chunks : [bam_chunks]
    def bam_list  = bam_files.collect { it.name }.join(' ')
    def n_bams    = bam_files.size()

    """
    set -euo pipefail

    echo "Merging ${n_bams} BAM chunk(s) for sample ${sample_id}"

    # 1. Merge all artificial FASTQ-chunk BAMs.
    #    -u avoids unnecessary compression before the following sort.
    samtools merge \
        -@ ${task.cpus} \
        -u \
        merged.unsorted.bam \
        ${bam_list}

    # 2. Re-sort the merged BAM by genomic coordinate.
    samtools sort \
        -@ ${task.cpus} \
        -o merged.sorted.bam \
        merged.unsorted.bam

    # 3. Remove duplicate flags/tags inherited from chunk-level
    #    methylseq MarkDuplicates runs.
    #    Reads themselves are retained.
    samtools view \
        -@ ${task.cpus} \
        -b \
        --remove-flags DUP \
        -x DT \
        -o merged.sorted.cleardup.bam \
        merged.sorted.bam

    # 4. The chunks are artificial computational splits of one sample.
    #    Normalize read-group/library information so global duplicate
    #    marking treats them as one biological sample/library.
    picard AddOrReplaceReadGroups \
        INPUT=merged.sorted.cleardup.bam \
        OUTPUT=merged.sorted.cleardup.rg.bam \
        RGID="${sample_id}" \
        RGLB="${sample_id}" \
        RGPL="ILLUMINA" \
        RGPU="${sample_id}" \
        RGSM="${sample_id}" \
        SORT_ORDER=coordinate \
        CREATE_INDEX=false

    # 5. Recalculate duplicate status globally across the full sample.
    #    Duplicates are MARKED, not removed.
    picard MarkDuplicates \
        INPUT=merged.sorted.cleardup.rg.bam \
        OUTPUT=${sample_id}.markdup.sorted.bam \
        METRICS_FILE=${sample_id}.markdup.metrics.txt \
        REMOVE_DUPLICATES=false \
        ASSUME_SORT_ORDER=coordinate \
        CLEAR_DT=true \
        CREATE_INDEX=false

    # 6. Create BAM index.
    samtools index \
        -@ ${task.cpus} \
        ${sample_id}.markdup.sorted.bam

    # 7. Basic validation that BAM/index were written successfully.
    samtools quickcheck \
        ${sample_id}.markdup.sorted.bam

    echo "Final BAM:"
    ls -lh \
        ${sample_id}.markdup.sorted.bam \
        ${sample_id}.markdup.sorted.bam.bai
    """
}
