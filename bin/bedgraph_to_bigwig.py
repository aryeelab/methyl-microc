#!/usr/bin/env python3
"""
Convert bedGraph files to bigWig format using pyBigWig
Usage: python bedgraph_to_bigwig.py input.bedGraph output.bw chromosome_sizes.txt
"""

import sys
import pyBigWig
import os

def get_chromosome_sizes_from_fasta_index(fai_file):
    """Extract chromosome sizes from a FASTA index file"""
    from collections import OrderedDict
    chrom_sizes = OrderedDict()
    with open(fai_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom = parts[0]
                size = int(parts[1])
                chrom_sizes[chrom] = size
    return chrom_sizes

def convert_bedgraph_to_bigwig(bedgraph_file, bigwig_file, chrom_sizes):
    """Convert bedGraph to bigWig format"""

    # Create bigWig file
    bw = pyBigWig.open(bigwig_file, "w")

    # Add chromosome sizes (maintain order from the fai file)
    chrom_order = list(chrom_sizes.keys())
    bw.addHeader(list(chrom_sizes.items()))

    # Read bedGraph and group by chromosome
    chrom_entries = {}
    with open(bedgraph_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip track lines and empty lines
            if line.startswith('track') or line.startswith('#') or not line:
                continue

            parts = line.split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                value = float(parts[3])

                # Only add if chromosome exists in our sizes
                if chrom in chrom_sizes:
                    if chrom not in chrom_entries:
                        chrom_entries[chrom] = []
                    chrom_entries[chrom].append((start, end, value))

    # Process chromosomes in the order they appear in the fai file
    for chrom in chrom_order:
        if chrom in chrom_entries:
            # Sort entries for this chromosome by position
            chrom_entries[chrom].sort(key=lambda x: x[0])

            # Prepare data for pyBigWig
            chroms = [chrom] * len(chrom_entries[chrom])
            starts = [x[0] for x in chrom_entries[chrom]]
            ends = [x[1] for x in chrom_entries[chrom]]
            values = [x[2] for x in chrom_entries[chrom]]

            # Add entries for this chromosome
            bw.addEntries(chroms, starts, ends=ends, values=values)

    bw.close()
    print(f"Successfully converted {bedgraph_file} to {bigwig_file}")

def main():
    if len(sys.argv) != 4:
        print("Usage: python bedgraph_to_bigwig.py input.bedGraph output.bw reference.fa.fai")
        print("   or: python bedgraph_to_bigwig.py input.bedGraph output.bw chromosome_sizes.txt")
        sys.exit(1)
    
    bedgraph_file = sys.argv[1]
    bigwig_file = sys.argv[2]
    sizes_file = sys.argv[3]
    
    # Check if input files exist
    if not os.path.exists(bedgraph_file):
        print(f"Error: bedGraph file {bedgraph_file} not found")
        sys.exit(1)
    
    if not os.path.exists(sizes_file):
        print(f"Error: sizes file {sizes_file} not found")
        sys.exit(1)
    
    # Get chromosome sizes
    chrom_sizes = get_chromosome_sizes_from_fasta_index(sizes_file)
    
    if not chrom_sizes:
        print(f"Error: No chromosome sizes found in {sizes_file}")
        sys.exit(1)
    
    print(f"Found {len(chrom_sizes)} chromosomes")
    
    # Convert
    convert_bedgraph_to_bigwig(bedgraph_file, bigwig_file, chrom_sizes)

if __name__ == "__main__":
    main()
