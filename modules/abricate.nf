process FORMAT_ABRICATE_DB {
    tag "$database"
    label 'light'
    container 'python:3.11-slim-bookworm'

    input:
    tuple val(database), path(reference)

    output:
    tuple val(database), path("${database}.abricate.fasta"), emit: formatted
    tuple val(database), path("${database}.reference_metadata.tsv"), emit: metadata

    script:
    """
    python3 "${projectDir}/bin/prepare_abricate_db.py" \\
        --database "${database}" \\
        --input "${reference}" \\
        --output "${database}.abricate.fasta" \\
        --metadata "${database}.reference_metadata.tsv"
    """
}

process PREPARE_ABRICATE_DB {
    tag "$database"
    label 'light'
    container 'quay.io/biocontainers/abricate:1.4.0--h05cac1d_0'
    memory '4 GB'
    maxForks 2

    input:
    tuple val(database), path(formatted_reference)

    output:
    tuple val(database), path("${database}.db.tar.gz"), emit: database

    script:
    """
    mkdir -p databases/${database}
    cp "${formatted_reference}" databases/${database}/sequences

    abricate --setupdb --datadir databases
    tar -czf "${database}.db.tar.gz" -C databases "${database}"
    """
}

process ABRICATE_SCAN {
    tag "${sample}:${database}"
    label 'light'
    container 'quay.io/biocontainers/abricate:1.4.0--h05cac1d_0'
    cpus 2
    memory '8 GB'
    maxForks 2

    input:
    tuple val(sample), path(assembly), val(database), path(database_archive)

    output:
    tuple val(sample), val(database), path("${sample}.${database}.abricate.tsv"), emit: hits

    script:
    """
    mkdir -p databases
    tar -xzf "${database_archive}" -C databases
    abricate \\
        --datadir databases \\
        --db "${database}" \\
        --minid ${params.min_identity} \\
        --mincov ${params.min_coverage} \\
        "${assembly}" > "${sample}.${database}.abricate.tsv"
    """
}
