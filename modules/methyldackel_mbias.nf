process METHYLDACKEL_MBIAS {

    tag "${sample_id}"

    publishDir "${params.outdir}/methylDackelExtracts/mbias", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta
    path fai

    output:
    path("*.svg"), emit: mbias_svg
    path("*.tsv"), emit: mbias_tsv

    script:
    """
    set -euo pipefail

    echo "========================================"
    echo "MethylDackel mbias"
    echo "Sample: ${sample_id}"
    echo "BAM: ${bam}"
    echo "Reference: ${fasta}"
    echo "========================================"

    echo -e "chr\\tcontext\\tstrand\\tRead\\tPosition\\tnMethylated\\tnUnmethylated\\tnMethylated(+dups)\\tnUnmethylated(+dups)" \
        > ${sample_id}.combined_mbias.tsv

    chrs=(\$(
        samtools view -H "${bam}" \
            | grep '@SQ' \
            | cut -f 2 \
            | sed 's/SN://' \
            | grep -v '_random' \
            | grep -v 'chrUn' \
            | sed 's/|/\\\\|/'
    ))

    for chr in "\${chrs[@]}"; do

        for context in CHH CHG CpG; do

            arg=''

            if [[ "\$context" == "CHH" ]]; then
                arg='--CHH --noCpG'
            elif [[ "\$context" == "CHG" ]]; then
                arg='--CHG --noCpG'
            fi

            join \
                -t \$'\\t' \
                -j1 \
                -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 \
                -a 1 \
                -e 0 \
                <(
                    MethylDackel mbias \
                        --noSVG \
                        \$arg \
                        -@ ${task.cpus} \
                        -r "\$chr" \
                        ${fasta} \
                        "${bam}" \
                    | tail -n +2 \
                    | awk '{print \$1"-"\$2"-"\$3"\\t"\$0}' \
                    | sort -k 1b,1
                ) \
                <(
                    MethylDackel mbias \
                        --noSVG \
                        --keepDupes \
                        -F 2816 \
                        \$arg \
                        -@ ${task.cpus} \
                        -r "\$chr" \
                        ${fasta} \
                        "${bam}" \
                    | tail -n +2 \
                    | awk '{print \$1"-"\$2"-"\$3"\\t"\$0}' \
                    | sort -k 1b,1
                ) \
                | sed "s/^/\${chr}\\t\${context}\\t/" \
                >> ${sample_id}.combined_mbias.tsv

        done
    done

    #
    # CHN SVG for positional/trimming bias
    #
    MethylDackel mbias \
        -@ ${task.cpus} \
        --noCpG \
        --CHH \
        --CHG \
        -r "\${chrs[0]}" \
        ${fasta} \
        "${bam}" \
        ${sample_id}_chn

    for f in *chn*.svg; do
        sed -i.bak \
            "s/Strand<\\/text>/Strand \$f \${chrs[0]} CHN <\\/text>/" \
            "\$f"
    done

    #
    # CpG SVG for positional/trimming bias
    #
    MethylDackel mbias \
        -@ ${task.cpus} \
        -r "\${chrs[0]}" \
        "${fasta}" \
        "${bam}" \
        ${sample_id}_cpg

    for f in *cpg*.svg; do
        sed -i.bak \
            "s/Strand<\\/text>/Strand \$f \${chrs[0]} CpG<\\/text>/" \
            "\$f"
    done

    echo
    echo "Generated MethylDackel mbias files:"
    ls -lh *.tsv *.svg
    """
}
