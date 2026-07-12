process ANNOTATE_PAIRS {

    tag "${sample_id}"

    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    tuple val(sample_id), path(pairs)
    path fasta
    path fai

    output:
    tuple val(sample_id), path("${sample_id}.meth.pairs.gz"), emit: meth_pairs

    script:
    def total_cpus      = task.cpus as int
    def annotate_jobs   = Math.max(1, total_cpus - 1)
    def pigz_threads    = Math.max(1, total_cpus)
    def lines_per_chunk = params.annotate_pairs_per_chunk ?: 5_000_000

    """
    set -euo pipefail

    mkdir -p annotate_chunks

    echo "[ANNOTATE_PAIRS] Sample: ${sample_id}" >&2
    echo "[ANNOTATE_PAIRS] Splitting ${pairs} into chunks of ${lines_per_chunk} pair lines" >&2

    pigz -dc ${pairs} | awk '
        /^#/ {
            print > "annotate_chunks/header.txt"
            next
        }
        {
            print
        }
    ' | split \
        -d \
        -a 6 \
        -l ${lines_per_chunk} \
        - \
        annotate_chunks/body_

    n_chunks=\$(find annotate_chunks -maxdepth 1 -type f -name 'body_[0-9]*' | wc -l)

    echo "[ANNOTATE_PAIRS] Number of chunks: \${n_chunks}" >&2

    if [[ "\${n_chunks}" -eq 0 ]]; then
        echo "ERROR: No annotation chunks were generated." >&2
        exit 1
    fi

    cat > annotate_one_chunk.sh <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

body="\$1"
base=\$(basename "\$body")

echo "[ANNOTATE_PAIRS] Start \${base}" >&2

cat annotate_chunks/header.txt "\$body" | \
python ${projectDir}/bin/annotate_pairs_methylation.py \
    --input - \
    --fasta ${fasta} \
    --fai ${fai} \
    --output - \
    --progress-every 0 \
    > "annotate_chunks/\${base}.meth.pairs"

echo "[ANNOTATE_PAIRS] Done \${base}" >&2
EOS

    chmod +x annotate_one_chunk.sh

    find annotate_chunks \
        -maxdepth 1 \
        -type f \
        -name 'body_[0-9]*' \
        -print0 | \
        sort -z | \
        xargs \
            -0 \
            -n 1 \
            -P ${annotate_jobs} \
            ./annotate_one_chunk.sh

    echo "[ANNOTATE_PAIRS] Concatenating annotated chunks" >&2

    first=1

    while IFS= read -r f; do
        if [[ "\${first}" -eq 1 ]]; then
            cat "\${f}"
            first=0
        else
            grep -v '^#' "\${f}"
        fi
    done < <(
        find annotate_chunks \
            -maxdepth 1 \
            -type f \
            -name 'body_*.meth.pairs' \
            | sort
    ) | pigz \
        -p ${pigz_threads} \
        > ${sample_id}.meth.pairs.gz

    pigz -t ${sample_id}.meth.pairs.gz

    echo "[ANNOTATE_PAIRS] Finished ${sample_id}.meth.pairs.gz" >&2
    """
}
