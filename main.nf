nextflow.enable.dsl=2

include { PRODIGAL } from './modules/prodigal'
include { FORMAT_ABRICATE_DB; PREPARE_ABRICATE_DB; ABRICATE_SCAN } from './modules/abricate'
include { ISESCAN } from './modules/isescan'
include { INTEGRATE_FEATURES } from './modules/integrate'
include { PREPARE_ASSEMBLY } from './modules/prepare_assembly'

workflow {
    main:
    if (!params.input) error 'Required: --input <CSV with sample,fasta columns>'
    if (!params.resfinder_fasta) error 'Required: --resfinder_fasta <pinned nucleotide FASTA>'
    if (!params.vfdb_fasta) error 'Required: --vfdb_fasta <pinned nucleotide FASTA>'
    if (!params.outdir) error 'Required: --outdir <output directory or s3:// URI>'

    assemblies = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fasta) error "Each input row needs sample and fasta: ${row}"
            tuple(row.sample.toString(), file(row.fasta.toString()))
        }

    references = Channel.of(
        tuple('resfinder', file(params.resfinder_fasta.toString())),
        tuple('vfdb', file(params.vfdb_fasta.toString()))
    )

    FORMAT_ABRICATE_DB(references)
    PREPARE_ABRICATE_DB(FORMAT_ABRICATE_DB.out.formatted)
    PREPARE_ASSEMBLY(assemblies)
    PRODIGAL(PREPARE_ASSEMBLY.out.assembly)
    ISESCAN(PREPARE_ASSEMBLY.out.assembly)

    assembly_files = PREPARE_ASSEMBLY.out.assembly.map { sample, assembly -> tuple(sample, assembly) }
    ABRICATE_SCAN(assembly_files.combine(PREPARE_ABRICATE_DB.out.database))

    resfinder_hits = ABRICATE_SCAN.out.hits
        .filter { sample, database, report -> database == 'resfinder' }
        .map { sample, database, report -> tuple(sample, report) }
    vfdb_hits = ABRICATE_SCAN.out.hits
        .filter { sample, database, report -> database == 'vfdb' }
        .map { sample, database, report -> tuple(sample, report) }

    hit_pair = resfinder_hits.join(vfdb_hits, by: 0)
        .map { sample, resfinder, vfdb -> tuple(sample, resfinder, vfdb) }

    annotation_bundle = PREPARE_ASSEMBLY.out.assembly
        .join(PRODIGAL.out.bundle, by: 0)
        .map { sample, assembly, gff, faa, fna -> tuple(sample, assembly, gff, faa, fna) }
    annotation_bundle
        .join(hit_pair, by: 0)
        .join(ISESCAN.out.records, by: 0)
        .map { sample, assembly, gff, faa, fna, resfinder, vfdb, is_tsv, is_gff ->
            tuple(sample, assembly, gff, faa, fna, resfinder, vfdb, is_tsv, is_gff)
        }
        .set { integration_inputs }

    INTEGRATE_FEATURES(integration_inputs)

    publish:
    prodigal_gff = PRODIGAL.out.bundle.map { sample, gff, faa, fna -> tuple(sample, gff) }
    prodigal_faa = PRODIGAL.out.bundle.map { sample, gff, faa, fna -> tuple(sample, faa) }
    prodigal_fna = PRODIGAL.out.bundle.map { sample, gff, faa, fna -> tuple(sample, fna) }
    abricate_hits = ABRICATE_SCAN.out.hits
    abricate_reference_metadata = FORMAT_ABRICATE_DB.out.metadata
    isescan_tsv = ISESCAN.out.records.map { sample, tsv, gff -> tuple(sample, tsv) }
    isescan_gff = ISESCAN.out.records.map { sample, tsv, gff -> tuple(sample, gff) }
    contig_features = INTEGRATE_FEATURES.out.features
    nearby_locus_pairs = INTEGRATE_FEATURES.out.pairs
    hit_gene_context = INTEGRATE_FEATURES.out.context
    sample_summaries = INTEGRATE_FEATURES.out.summary
}

output {
    prodigal_gff {
        path 'prodigal/gff'
    }
    prodigal_faa {
        path 'prodigal/proteins'
    }
    prodigal_fna {
        path 'prodigal/cds'
    }
    abricate_hits {
        path 'abricate/hits'
    }
    abricate_reference_metadata {
        path 'abricate/reference_metadata'
    }
    isescan_tsv {
        path 'isescan/tsv'
    }
    isescan_gff {
        path 'isescan/gff'
    }
    contig_features {
        path 'integrated/contig_features'
    }
    nearby_locus_pairs {
        path 'integrated/nearby_locus_pairs'
    }
    hit_gene_context {
        path 'integrated/hit_gene_context'
    }
    sample_summaries {
        path 'integrated/sample_summaries'
    }
}
