#!/usr/bin/env bash
#SBATCH --partition=short
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=methyl_microc_io
#SBATCH --output=benchmarks/o2/logs/io_microbench_%j.out
#SBATCH --error=benchmarks/o2/logs/io_microbench_%j.err

set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
NFS_DIR="$ROOT/benchmarks/o2/io_nfs/${SLURM_JOB_ID}"
TMP_DIR="/tmp/${USER}/${SLURM_JOB_ID}/io_microbench"

mkdir -p "$NFS_DIR" "$TMP_DIR"

cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

run_test() {
    local label="$1"
    local dir="$2"
    local file="$dir/blob.bin"

    echo "## $label"
    df -h "$dir" | tail -1 || true
    stat -f -c 'fstype=%T' "$dir" || true

    echo "# sequential write 2GiB"
    /usr/bin/time -f "elapsed=%e user=%U sys=%S maxrss_kb=%M" \
        dd if=/dev/zero of="$file" bs=64M count=32 conv=fdatasync status=none

    echo "# sequential read 2GiB"
    /usr/bin/time -f "elapsed=%e user=%U sys=%S maxrss_kb=%M" \
        dd if="$file" of=/dev/null bs=64M status=none

    echo "# small file create 5000"
    mkdir -p "$dir/small"
    /usr/bin/time -f "elapsed=%e user=%U sys=%S maxrss_kb=%M" \
        bash -c 'for i in $(seq 1 5000); do printf "%s\n" "$i" > "$0/small/file_$i.txt"; done' "$dir"

    echo "# small file stat 5000"
    /usr/bin/time -f "elapsed=%e user=%U sys=%S maxrss_kb=%M" \
        find "$dir/small" -type f -exec stat {} \; >/dev/null
}

echo "host=$(hostname)"
echo "job=${SLURM_JOB_ID}"
echo "date=$(date -Is)"
echo "root=$ROOT"
echo "tmp_dir=$TMP_DIR"
echo "nfs_dir=$NFS_DIR"

run_test "node_local_tmp" "$TMP_DIR"
run_test "project_nfs" "$NFS_DIR"
