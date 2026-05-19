process SPLIT_FASTQ {

    input:
    tuple val(sample_id), path(r1), path(r2)
    val reads_per_chunk

    output:
    path "chunks/${sample_id}.chunks.tsv", emit: manifest
    path "chunks/*.fastq.gz", emit: chunk_fastqs

    script:
    """
    mkdir -p chunks

    python - <<'PY'
import gzip
import subprocess
from pathlib import Path

sample_id = "${sample_id}"
r1_path = Path("${r1}")
r2_path = Path("${r2}")
reads_per_chunk = int("${reads_per_chunk}")
pigz_threads = max(1, int("${task.cpus}"))

if reads_per_chunk <= 0:
    raise ValueError("reads_per_chunk must be a positive integer")

outdir = Path("chunks")
outdir.mkdir(exist_ok=True)

manifest_rows = []

def open_pigz_writer(out_path):
    fout = open(out_path, "wb")
    proc = subprocess.Popen(
        ["pigz", "-p", str(pigz_threads), "-c"],
        stdin=subprocess.PIPE,
        stdout=fout
    )
    return proc, fout

def close_pigz_writer(proc, fout):
    if proc.stdin:
        proc.stdin.close()
    ret = proc.wait()
    fout.close()
    if ret != 0:
        raise RuntimeError(f"pigz failed with exit code {ret}")

read_idx = 0
current_chunk_id = None
r1_proc = r2_proc = None
r1_fout = r2_fout = None

with gzip.open(r1_path, "rb") as f1, gzip.open(r2_path, "rb") as f2:
    while True:
        rec1 = [f1.readline() for _ in range(4)]
        rec2 = [f2.readline() for _ in range(4)]

        if not rec1[0] and not rec2[0]:
            break

        if not all(rec1) or not all(rec2):
            raise RuntimeError("R1/R2 FASTQ files are truncated or out of sync")

        h1 = rec1[0].strip().split()[0]
        h2 = rec2[0].strip().split()[0]

        h1_core = h1[:-2] if h1.endswith(b'/1') else h1
        h2_core = h2[:-2] if h2.endswith(b'/2') else h2

        if h1_core != h2_core:
            raise RuntimeError(f"Read pair mismatch: {h1.decode(errors='ignore')} vs {h2.decode(errors='ignore')}")

        chunk_id = f"chunk{(read_idx // reads_per_chunk) + 1:03d}"

        if chunk_id != current_chunk_id:
            if current_chunk_id is not None:
                close_pigz_writer(r1_proc, r1_fout)
                close_pigz_writer(r2_proc, r2_fout)

            r1_out = outdir / f"{sample_id}.{chunk_id}_R1.fastq.gz"
            r2_out = outdir / f"{sample_id}.{chunk_id}_R2.fastq.gz"

            r1_proc, r1_fout = open_pigz_writer(r1_out)
            r2_proc, r2_fout = open_pigz_writer(r2_out)

            manifest_rows.append((sample_id, chunk_id, str(r1_out.resolve()), str(r2_out.resolve())))
            current_chunk_id = chunk_id

        r1_proc.stdin.write(b"".join(rec1))
        r2_proc.stdin.write(b"".join(rec2))

        read_idx += 1

if current_chunk_id is not None:
    close_pigz_writer(r1_proc, r1_fout)
    close_pigz_writer(r2_proc, r2_fout)

manifest = outdir / f"{sample_id}.chunks.tsv"
with manifest.open("w") as fout:
    for row in manifest_rows:
        fout.write("\\t".join(row) + "\\n")

print(f"Split {read_idx} read pairs from sample {sample_id} into {len(manifest_rows)} chunk(s), with up to {reads_per_chunk} read pairs per chunk.")
print(f"Used pigz with {pigz_threads} thread(s) per output stream.")
PY
    """
}
