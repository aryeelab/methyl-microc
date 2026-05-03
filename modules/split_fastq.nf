process SPLIT_FASTQ {

    input:
    tuple val(sample_id), path(r1), path(r2)
    val n_chunks

    output:
    path "chunks/${sample_id}.chunks.tsv", emit: manifest
    path "chunks/*.fastq.gz", emit: chunk_fastqs

    script:
    """
    mkdir -p chunks

    python - <<'PY'
import gzip
from pathlib import Path

sample_id = "${sample_id}"
r1_path = Path("${r1}")
r2_path = Path("${r2}")
n_chunks = int("${n_chunks}")

outdir = Path("chunks")
outdir.mkdir(exist_ok=True)

writers = {}
manifest_rows = []

def get_writers(chunk_id):
    if chunk_id not in writers:
        r1_out = outdir / f"{sample_id}.{chunk_id}_R1.fastq.gz"
        r2_out = outdir / f"{sample_id}.{chunk_id}_R2.fastq.gz"
        writers[chunk_id] = (
            gzip.open(r1_out, "wt"),
            gzip.open(r2_out, "wt"),
            r1_out.resolve(),
            r2_out.resolve(),
        )
        manifest_rows.append((sample_id, chunk_id, str(r1_out.resolve()), str(r2_out.resolve())))
    return writers[chunk_id]

read_idx = 0

with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
    while True:
        rec1 = [f1.readline() for _ in range(4)]
        rec2 = [f2.readline() for _ in range(4)]

        if not rec1[0] and not rec2[0]:
            break

        if not all(rec1) or not all(rec2):
            raise RuntimeError("R1/R2 FASTQ files are truncated or out of sync")

        h1 = rec1[0].strip().split()[0]
        h2 = rec2[0].strip().split()[0]

        h1_core = h1[:-2] if h1.endswith('/1') else h1
        h2_core = h2[:-2] if h2.endswith('/2') else h2

        if h1_core != h2_core:
            raise RuntimeError(f"Read pair mismatch: {h1} vs {h2}")

        chunk_id = f"chunk{(read_idx % n_chunks) + 1:03d}"
        w1, w2, _, _ = get_writers(chunk_id)
        w1.write("".join(rec1))
        w2.write("".join(rec2))
        read_idx += 1

for w1, w2, _, _ in writers.values():
    w1.close()
    w2.close()

manifest = outdir / f"{sample_id}.chunks.tsv"
with manifest.open("w") as fout:
    for row in manifest_rows:
        fout.write("\\t".join(row) + "\\n")

print(f"Split {read_idx} read pairs from sample {sample_id} into {len(manifest_rows)} chunk(s).")
PY
    """
}
