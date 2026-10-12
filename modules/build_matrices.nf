process BUILD_MATRICES {

    tag "${sample_id}"

    publishDir "${params.outdir}/matrices/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(pairs)
    path fai
    val matrix_base_resolution
    val matrix_resolutions

    output:
    tuple val(sample_id),
          path("${sample_id}.cool"),
          path("${sample_id}.mcool"),
          path("${sample_id}.hic"),
          emit: matrices

    script:
    def base_resolution = matrix_base_resolution as int
    def resolution_mode = matrix_resolutions.toString().trim()
    def max_resolution  = 10_000_000

    if (base_resolution <= 0) {
        throw new IllegalArgumentException(
            "matrix_base_resolution must be a positive integer"
        )
    }

    if (base_resolution > max_resolution) {
        throw new IllegalArgumentException(
            "matrix_base_resolution cannot be larger than 10 Mb"
        )
    }

    def requested_resolutions

    /*
     * DEFAULT MODE
     *
     * Nice steps relative to the base resolution,
     * stopping at 10 Mb.
     *
     * Default base = 1 kb:
     * 1 kb, 2 kb, 5 kb, 10 kb, 20 kb, 50 kb,
     * 100 kb, 200 kb, 500 kb, 1 Mb, 2 Mb, 5 Mb, 10 Mb
     */
    if (resolution_mode == 'nice') {

        requested_resolutions = []

        def scale = 1L

        while (true) {

            def candidates = [
                base_resolution * scale,
                base_resolution * 2L * scale,
                base_resolution * 5L * scale
            ]

            candidates.each { res ->
                if (res <= max_resolution) {
                    requested_resolutions << (res as int)
                }
            }

            if (base_resolution * 10L * scale > max_resolution) {
                break
            }

            scale *= 10L
        }

        requested_resolutions =
            requested_resolutions
                .unique()
                .sort()

    } else {

        /*
         * CUSTOM MODE
         */
        requested_resolutions = resolution_mode
            .split(',')
            .collect { it.trim() }
            .findAll { it }
            .collect { it as int }
            .unique()
            .sort()

        if (!requested_resolutions) {
            throw new IllegalArgumentException(
                "matrix_resolutions must be 'nice' or a comma-separated list of resolutions"
            )
        }

        if (!requested_resolutions.contains(base_resolution)) {
            throw new IllegalArgumentException(
                "Custom matrix_resolutions must include matrix_base_resolution (${base_resolution})"
            )
        }

        def invalid_resolutions = requested_resolutions.findAll {
            it < base_resolution ||
            it > max_resolution ||
            it % base_resolution != 0
        }

        if (invalid_resolutions) {
            throw new IllegalArgumentException(
                "All custom resolutions must be >= matrix_base_resolution, " +
                "<= 10 Mb, and integer multiples of matrix_base_resolution. " +
                "Invalid values: ${invalid_resolutions.join(',')}"
            )
        }
    }

    def resolution_args =
        "--resolutions " + requested_resolutions.join(' ')

    def expected_resolution_string =
        requested_resolutions.join(',')

    """
    set -euo pipefail

    mkdir -p tmp

    echo "========================================"
    echo "Building matrices for sample: ${sample_id}"
    echo "Original pairs: ${pairs}"
    echo "FAI: ${fai}"
    echo "Base resolution: ${base_resolution}"
    echo "Resolution mode: ${resolution_mode}"
    echo "Expected resolutions: ${expected_resolution_string}"
    echo "========================================"

    #
    # 1. Create matrix-ready pairs with FAI-ordered #chromsize header.
    #
    # IMPORTANT:
    # - original pairs file is NOT modified
    # - data rows are unchanged
    # - only #chromsize header lines are replaced/reordered
    #
    MATRIX_PAIRS="${sample_id}.matrix.pairs.gz"
    HEADER_OTHER="tmp/${sample_id}.header.other"
    HEADER_CHROMS="tmp/${sample_id}.header.chromsizes"
    HEADER_FINAL="tmp/${sample_id}.header.final"

    #
    # Keep all original header lines except #chromsize.
    #
    # Do NOT exit early from awk here because set -o pipefail would
    # interpret the resulting SIGPIPE from pigz as a process failure.
    #
    pigz -dc ${pairs} \
        | awk '
            /^#/ {
                if (\$0 !~ /^#chromsize:/) {
                    print
                }
            }
        ' > "\$HEADER_OTHER"

    #
    # Build #chromsize lines in authoritative FAI order.
    #
    awk 'BEGIN{OFS=" "} {print "#chromsize:", \$1, \$2}' ${fai} \
        > "\$HEADER_CHROMS"

    #
    # Reconstruct header:
    # insert FAI-ordered #chromsize lines immediately before #columns.
    #
    awk -v chromfile="\$HEADER_CHROMS" '
        /^#columns:/ {
            while ((getline line < chromfile) > 0) {
                print line
            }
            close(chromfile)
        }
        { print }
    ' "\$HEADER_OTHER" > "\$HEADER_FINAL"

    #
    # Write corrected header + unchanged data rows.
    #
    {
        cat "\$HEADER_FINAL"

        pigz -dc ${pairs} \
            | awk '!/^#/'
    } \
        | pigz -p ${task.cpus} \
        > "\$MATRIX_PAIRS"

    echo
    echo "Created temporary matrix-ready pairs:"
    ls -lh "\$MATRIX_PAIRS"

    #
    # 2. Sanity check:
    # decompressed data rows must be identical.
    #
    ORIGINAL_DATA_MD5=\$(
        pigz -dc ${pairs} \
            | awk '!/^#/' \
            | md5sum \
            | awk '{print \$1}'
    )

    MATRIX_DATA_MD5=\$(
        pigz -dc "\$MATRIX_PAIRS" \
            | awk '!/^#/' \
            | md5sum \
            | awk '{print \$1}'
    )

    echo "Original data MD5: \$ORIGINAL_DATA_MD5"
    echo "Matrix   data MD5: \$MATRIX_DATA_MD5"

    if [[ "\$ORIGINAL_DATA_MD5" != "\$MATRIX_DATA_MD5" ]]; then
        echo "ERROR: matrix-ready pairs data rows differ from original pairs."
        exit 1
    fi

    echo "Data rows are unchanged."

    #
    # 3. Build base-resolution .cool
    #
    hictk load \
        --format 4dn \
        --bin-size ${base_resolution} \
        --threads ${task.cpus} \
        --tmpdir tmp \
        --assume-sorted \
        --force \
        "\$MATRIX_PAIRS" \
        ${sample_id}.cool

    #
    # 4. Build multi-resolution .mcool
    #
    hictk zoomify \
        --copy-base-resolution \
        ${resolution_args} \
        --threads ${task.cpus} \
        --tmpdir tmp \
        --force \
        ${sample_id}.cool \
        ${sample_id}.mcool

    #
    # 5. Validate .mcool
    #
    hictk validate ${sample_id}.mcool \
        > ${sample_id}.mcool.validate.txt

    EXPECTED_RESOLUTIONS="${expected_resolution_string}"

    IFS=',' read -ra RESOLUTIONS <<< "\$EXPECTED_RESOLUTIONS"

    for RES in "\${RESOLUTIONS[@]}"; do

        if ! grep -Eq "\\"?\${RES}\\"?[[:space:]]*:" \
            ${sample_id}.mcool.validate.txt
        then
            echo "ERROR: requested resolution \${RES} is missing from ${sample_id}.mcool"
            echo
            echo "Expected resolutions:"
            echo "\$EXPECTED_RESOLUTIONS"
            echo
            echo "hictk validation output:"
            cat ${sample_id}.mcool.validate.txt
            exit 1
        fi

    done

    echo "All requested .mcool resolutions are present."

    #
    # 6. Build base-resolution .hic from the same matrix-ready pairs
    #
    hictk load \
        --format 4dn \
        --bin-size ${base_resolution} \
        --threads ${task.cpus} \
        --tmpdir tmp \
        --assume-sorted \
        --force \
        "\$MATRIX_PAIRS" \
        ${sample_id}.hic

    #
    # 7. Validate final .cool and .hic
    #
    hictk validate ${sample_id}.cool \
        > ${sample_id}.cool.validate.txt

    hictk validate ${sample_id}.hic \
        > ${sample_id}.hic.validate.txt

    echo
    echo "Generated matrices:"
    ls -lh \
        ${sample_id}.cool \
        ${sample_id}.mcool \
        ${sample_id}.hic

    echo
    echo "Matrix generation completed successfully."
    echo "Verified resolutions: ${expected_resolution_string}"
    """
}
