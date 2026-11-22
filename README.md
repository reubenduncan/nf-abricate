# nf-abricate

Contig-level screening and context integration for microbial assemblies. The workflow runs Prodigal, screens the nucleotide contigs against supplied ABRicate-compatible reference FASTAs, calls insertion sequences with ISEScan, and joins the features by contig coordinates.

The pipeline uses Nextflow DSL2 and the workflow-output syntax available from Nextflow 25.10. It is intended for metagenome assemblies and draft genomes, including contigs that were not retained in a final MAG catalogue.

## Workflow

1. Read a CSV with one assembly per row.
2. Run Prodigal in metagenome mode (`-p meta`) to predict CDS features and write GFF, protein FASTA, and CDS nucleotide FASTA.
3. Convert the supplied ResFinder and VFDB nucleotide FASTAs into ABRicate custom database format, index them, and retain a metadata table linking ABRicate labels to source headers.
4. Run ABRicate independently against each database. Default thresholds are 90% nucleotide identity and 80% reference coverage; both are configurable.
5. Run ISEScan on each assembly. ISEScan internally predicts proteins and searches its transposase profile models; its calls can be complete or partial.
6. Join Prodigal CDSs, AMR/VF matches, and ISEScan calls by contig and coordinates. Write nearby feature pairs and the Prodigal CDS context around each ABRicate hit.

## Inputs

Required parameters:

| Parameter | Description |
| --- | --- |
| `--input` | CSV file with `sample,fasta` columns. `fasta` can be a local path or an S3 URI supported by Nextflow. Each `sample` value must be unique. |
| `--resfinder_fasta` | Pinned ResFinder nucleotide FASTA. |
| `--vfdb_fasta` | Pinned VFDB nucleotide FASTA, such as Set A. |
| `--outdir` | Results directory or S3 URI. |

Optional parameters:

| Parameter | Default | Description |
| --- | ---: | --- |
| `--min_identity` | `90` | ABRicate minimum percent nucleotide identity. |
| `--min_coverage` | `80` | ABRicate minimum percent reference coverage. |
| `--max_pair_distance` | `10000` | Maximum number of intervening bases for a nearby-locus pair or hit-to-CDS context record. |

Example `samples.csv`:

```csv
sample,fasta
siteA_surface,assemblies/siteA_surface.fasta.gz
siteA_handpump,s3://example-bucket/assemblies/siteA_handpump.fasta.gz
```

Reference FASTAs are not bundled. Supply the exact pinned database files used for the analysis and keep their version, source URL, retrieval date, and checksums with the project records.

## Outputs

| Output | Contents |
| --- | --- |
| `prodigal/gff` | Prodigal CDS coordinates and gene-prediction attributes. |
| `prodigal/proteins` | Predicted proteins (`.faa`). |
| `prodigal/cds` | Predicted CDS nucleotide sequences (`.fna`). |
| `abricate/hits` | One raw tab-delimited ABRicate report per assembly and reference database, including identity, reference coverage, coordinates, strand, label, accession, and product fields. |
| `abricate/reference_metadata` | Mapping of each normalized custom database record to its original FASTA header and description. |
| `isescan/tsv` | ISEScan per-element records, including coordinates, family, strand, and complete/partial status. |
| `isescan/gff` | ISEScan GFF records. |
| `integrated/contig_features` | A unified table of Prodigal CDSs, ABRicate AMR/VF hits, and ISEScan elements. |
| `integrated/nearby_locus_pairs` | Pairs among AMR, VFDB, and IS features on the same contig within `max_pair_distance`, with gap, overlap, and orientation. |
| `integrated/hit_gene_context` | Prodigal CDSs within the distance threshold of each AMR/VFDB hit, including whether coordinates overlap. |
| `integrated/sample_summaries` | One JSON summary per input assembly. |

## Running

With Docker available locally:

```bash
nextflow run . -profile docker \
  --input samples.csv \
  --resfinder_fasta references/resfinder_all.fasta \
  --vfdb_fasta references/VFDB_setA_nt.fasta \
  --outdir results
```

On a managed compute environment that runs processes locally with Singularity, pass `conf/local-executor.config` as the launch configuration. On AWS Batch, use the platform compute environment's configured executor and container runtime.

## Software

- Prodigal 2.6.3 (`quay.io/biocontainers/prodigal:2.6.3--h577a1d6_11`)
- ABRicate 1.4.0 (`quay.io/biocontainers/abricate:1.4.0--h05cac1d_0`)
- ISEScan 1.7.3 (`quay.io/biocontainers/isescan:1.7.3--hab16a5f_1`)
- Python 3.11 (`python:3.11-slim-bookworm`) for reference formatting and feature integration

Tool processes use version-tagged public containers; the small integration steps use the official Python 3.11 slim image. The workflow retries tasks that exit with common termination signals twice, which helps recover from transient preemption; other task errors stop the run for diagnosis.

## Interpretation and limitations

- ABRicate is a nucleotide sequence match screen against the supplied reference set. It detects acquired-gene homologs represented in that set; it does not call resistance mutations or establish a clinical resistance phenotype. VFDB matches are homologous sequence evidence and do not alone establish pathogenicity or an intact virulence system.
- ISEScan identifies candidate insertion sequences through transposase profile matches and element-boundary features. Calls may be partial, and elements on contig edges can be truncated.
- A nearby-locus row establishes only that the reported features lie close together on the same assembled contig. It does not prove a complete mobile element, physical transfer, host identity, or sample-level presence. Short or fragmented contigs can split loci; coassembly can also complicate assignment to an individual library.
- Prodigal calls are gene predictions, not functional annotations. CDS proximity is reported as genomic context, not as a confirmed gene function.
- Results depend on the versions and contents of the supplied references. The workflow records the original FASTA headers but does not update or silently replace the databases.

## Citation

Please cite Prodigal, ABRicate, ISEScan, and the exact ResFinder/VFDB reference releases used. ISEScan publication: Xie and Tang, *Bioinformatics* 2017, 33(21):3340–3347, doi:10.1093/bioinformatics/btx433.
