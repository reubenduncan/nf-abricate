process PREPARE_ASSEMBLY {
    tag "$sample"
    label 'light'
    container 'python:3.11-bookworm'
    memory '4 GB'
    maxForks 4

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}.assembly.fasta"), emit: assembly

    script:
    """
    python3 "${projectDir}/bin/prepare_assembly.py" \\
        --input "${fasta}" \\
        --output "${sample}.assembly.fasta"
    """
}
