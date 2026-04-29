process RUN_METHYLSEQ {

    input:
    path samplesheet
    path fasta
    path methyl_config, stageAs: 'methylseq.config'
    path fastqs, stageAs: '*'

    output:
    path "results/**/deduplicated/*.markdup.sorted.bam", emit: bam

    script:
    def bwamethIndexArg = params.bwameth_index ? "--bwameth_index ${params.bwameth_index}" : ""
    def innerExtra = params.skip_inner_multiqc ? "--skip_multiqc" : (params.inner_methylseq_extra ?: "")
    def methylseqRunner = params.methylseq_runner ?: "${projectDir}/run_methylseq.sh"

    """
    python - <<'PY'
import csv
from pathlib import Path

infile = Path("${samplesheet.name}")
outfile = Path("samplesheet.staged.csv")

with infile.open(newline='') as fin, outfile.open("w", newline='') as fout:
    reader = csv.reader(fin)
    writer = csv.writer(fout, lineterminator='\\n')

    header = next(reader)
    writer.writerow(header)

    for row in reader:
        if not row or all(not x.strip() for x in row):
            continue
        row = [x.strip() for x in row]
        row[1] = Path(row[1]).name
        row[2] = Path(row[2]).name
        writer.writerow(row)
PY

    bash ${methylseqRunner} \
        --input samplesheet.staged.csv \
        --outdir results \
        --fasta ${fasta.name} \
        ${bwamethIndexArg} \
        ${innerExtra}
    """
}
