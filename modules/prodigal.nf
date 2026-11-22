process PRODIGAL {
    tag "$sample"
    label 'moderate'
    container 'quay.io/biocontainers/prodigal:2.6.3--h577a1d6_11'
    cpus 2
    memory '8 GB'
    maxForks 2

    input:
    tuple val(sample), path(assembly)

    output:
    tuple val(sample), path("${sample}.prodigal.gff"), path("${sample}.prodigal.faa"), path("${sample}.prodigal.fna"), emit: bundle

    script:
    """
    prodigal \\
        -i "${assembly}" \\
        -o "${sample}.prodigal.gff" \\
        -a "${sample}.prodigal.faa" \\
        -d "${sample}.prodigal.fna" \\
        -f gff \\
        -p meta
    """
}
