process ISESCAN {
    tag "$sample"
    label 'large'
    container 'quay.io/biocontainers/isescan:1.7.3--hab16a5f_1'
    cpus 2
    memory '16 GB'
    maxForks 1

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}.isescan.tsv"), path("${sample}.isescan.gff"), emit: records

    script:
    """
    isescan.py \\
        --seqfile "${assembly}" \\
        --output isescan_out \\
        --nthread ${task.cpus}

    if [[ -f "isescan_out/${assembly}.tsv" ]]; then
        cp "isescan_out/${assembly}.tsv" "${sample}.isescan.tsv"
    else
        touch "${sample}.isescan.tsv"
    fi
    if [[ -f "isescan_out/${assembly}.gff" ]]; then
        cp "isescan_out/${assembly}.gff" "${sample}.isescan.gff"
    else
        touch "${sample}.isescan.gff"
    fi
    """
}
