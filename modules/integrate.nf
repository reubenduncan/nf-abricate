process INTEGRATE_FEATURES {
    tag "$sample"
    label 'light'
    container 'python:3.11-slim-bookworm'
    memory '4 GB'
    maxForks 2

    input:
    tuple val(sample), path(assembly), path(prodigal_gff), path(prodigal_faa), path(prodigal_fna), path(resfinder_hits), path(vfdb_hits), path(isescan_tsv), path(isescan_gff)

    output:
    tuple val(sample), path("${sample}.contig_features.tsv"), emit: features
    tuple val(sample), path("${sample}.nearby_locus_pairs.tsv"), emit: pairs
    tuple val(sample), path("${sample}.hit_gene_context.tsv"), emit: context
    tuple val(sample), path("${sample}.summary.json"), emit: summary

    script:
    """
    python3 "${projectDir}/bin/integrate_features.py" \\
        --sample "${sample}" \\
        --assembly "${assembly}" \\
        --prodigal-gff "${prodigal_gff}" \\
        --resfinder "${resfinder_hits}" \\
        --vfdb "${vfdb_hits}" \\
        --isescan-tsv "${isescan_tsv}" \\
        --isescan-gff "${isescan_gff}" \\
        --max-distance ${params.max_pair_distance} \\
        --min-identity ${params.min_identity} \\
        --min-coverage ${params.min_coverage} \\
        --prefix "${sample}"
    """
}
