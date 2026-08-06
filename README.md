# METAeDNA (beta version)

METAeDNA is an R package for environmental DNA metabarcoding workflows. It brings together paired-end read preprocessing, denoising, decontamination, reference-database construction, in-silico PCR, BLAST/VSEARCH wrappers, taxonomy querying, and diversity rarefaction/extrapolation in one reproducible analysis framework.

The package is designed for eDNA projects where users need to move from raw sequencing files to clean sequence tables, curated references, and taxonomic assignments while keeping intermediate results traceable.

## Features

- **Sequence preprocessing**: demultiplexing, paired-end merging, dereplication, denoising, replicate merging, and record counting.
- **Traceable FASTA metadata**: intermediate metadata can be stored in structured FASTA headers and converted back into tidy tables.
- **External tool wrappers**: convenient R wrappers for BLAST+ and VSEARCH with automatic executable detection.
- **Reference database preparation**: download and prepare NCBI BLAST databases; build local NCBI taxonomy SQLite databases.
- **Adaptive GenBank reference building**: construct marker-specific reference FASTA files online, locally, or with a hybrid Entrez + local BLAST workflow.
- **Taxonomic assignment support**: query local accession-to-taxid and taxonomy databases; infer LCA from BLAST hits.
- **In-silico PCR**: search primer binding sites and extract expected amplicons from reference sequences.
- **Decontamination and diversity checks**: identify potential contaminants and run rarefaction/extrapolation summaries.

## Installation

Install from GitHub:

```r
remotes::install_github("intactcircle/METAeDNA")
```

Or install from a local checkout:

```r
install.packages("remotes")
remotes::install_local("path/to/METAeDNA")
```

METAeDNA uses Bioconductor packages such as `Biostrings`. If needed, install Bioconductor first:

```r
install.packages("BiocManager")
BiocManager::install("Biostrings")
```

## System Requirements

- R >= 4.0.0
- C++17 compiler
- zlib
- Optional external tools:
  - BLAST+ for `blast_cmd()`, `blastn_cmd()`, and local reference building
  - VSEARCH for clustering and chimera-removal helpers

## Typical Workflow

```r
library(METAeDNA)

# 1. Merge paired-end reads
pairend(
  fastq_1 = "sample_R1.fastq.gz",
  fastq_2 = "sample_R2.fastq.gz",
  out_file = "merged.fastq",
  n_threads = 4
)

# 2. Dereplicate and denoise sequences
dereplicate(
  input_fasta = "merged.fastq",
  output_fasta = "dereplicated.fasta",
  n_threads = 4
)

denoise(
  input_fasta = "dereplicated.fasta",
  output_fasta = "denoised.fasta",
  n_threads = 8
)

# 3. Convert FASTA metadata into a sequence table
seq_table <- fasta_to_table("denoised.fasta", threads = 4)

# 4. Cluster similar sequences if needed
clustered <- vsearch_cluster(
  seq_table = seq_table,
  output_path = "cluster_result",
  identity = 0.97,
  threads = 8
)
```

Function arguments may differ depending on input format and analysis settings; see the package help pages for details.

## Reference Database Construction

METAeDNA can work with local NCBI taxonomy resources and BLAST databases.

### Build Local NCBI Taxonomy Database

```r
build_ncbi_taxonomy_database(
  taxonomy = TRUE,
  accession = TRUE,
  include_WGS = FALSE,
  taxdump_dir = "ncbi_taxdump",
  sql_path = "ncbi_taxonomy.sqlite",
  overwrite = FALSE
)
```

This creates a SQLite database that can be used by `get_taxid()` and `get_taxonomy()`.

### Download and Prepare BLAST Databases

```r
download_blast_db(
  db = "ITS_fungal",
  database_path = "blast_download"
)

prepare_blast_database(
  file_path = "blast_download",
  database_path = "blast_db",
  threads = 4
)
```

### Build Marker-Specific GenBank References

Online Entrez mode:

```r
ref <- build_gb_reference(
  mode = "online",
  taxon = "Saxifragales",
  gene = "ITS",
  email = "your.name@example.org",
  max_records = 500,
  output_dir = "Saxifragales_ITS_ref",
  threads = 2
)
```

Local BLAST expansion mode:

```r
ref <- build_gb_reference(
  mode = "local",
  seed_fasta = "seed.fasta",
  database = "blast_db/nt",
  output_dir = "local_reference",
  max_rounds = 3,
  round_seed_limit = 200,
  threads = 8
)
```

Hybrid mode first retrieves online seed sequences and then expands them through a local BLAST database:

```r
ref <- build_gb_reference(
  mode = "hybrid",
  taxon = "Saxifragales",
  gene = "ITS",
  database = "blast_db/nt",
  output_dir = "hybrid_reference",
  email = "your.name@example.org",
  threads = 8
)
```

## Taxonomy Utilities

Query accession numbers:

```r
taxid <- get_taxid(
  sql_path = "ncbi_taxonomy.sqlite",
  accessions = c("AB123456.1", "NM_001200025.1")
)
```

Query full lineages:

```r
lineage <- get_taxonomy(
  sql_path = "ncbi_taxonomy.sqlite",
  tax_ids = taxid$taxid
)
```

Infer LCA taxonomy from BLAST results:

```r
lca <- find_lca(
  df = blast_results,
  seqs_ID = "qseqid",
  taxid_colum = "staxids",
  accession_colum = "sseqid"
)
```

Convert GenBank FASTA headers into a DADA2-style reference format:

```r
format_to_dada(
  input_fasta_file = "genbank_raw.fasta",
  output_fasta_file = "dada_reference.fasta",
  sql_path = "ncbi_taxonomy.sqlite"
)
```

## In-Silico PCR

Search primer binding sites:

```r
sites <- find_primer_binding_sites(
  target_sequence = "reference.fasta",
  primer = "AGCTTAGGCTAC",
  max_mismatch = 2
)
```

Run in-silico PCR against reference sequences:

```r
amplicons <- in_silico(
  fasta_file = "reference.fasta",
  forward = "ACGTACGT",
  reverse = "TGCAAGCT",
  mis_match = 2,
  trim = TRUE,
  chunk_size = 5000
)
```

For multiple primer pairs, pass equal-length vectors to `forward` and `reverse`
and set `workers` in `in_silico()`.

## BLAST and VSEARCH Wrappers

Run BLASTN:

```r
blastn_cmd(
  input = "query.fasta",
  output = "blast_results.tsv",
  database = "blast_db/nt",
  max_target_seqs = 10,
  perc_identity = 90,
  qcov_hsp_perc = 90,
  threads = 8
)
```

Run VSEARCH clustering:

```r
vsearch_cluster(
  input_fasta = "seqs.fasta",
  output_path = "vsearch_cluster",
  identity = 0.97,
  mode = "cluster_fast",
  threads = 8
)
```

Remove chimeras:

```r
vsearch_remove_chimeras(
  input_fasta = "seqs.fasta",
  output_fasta = "seqs_nonchimeras.fasta",
  mode = "uchime3_denovo",
  output_mode = "nonchimeras",
  threads = 8
)
```

## Decontamination and Diversity

Potential contaminants can be screened with `decontamination()`, and sequencing-depth effects can be evaluated with rarefaction/extrapolation:

```r
rare <- rarefaction_extrapolation(
  temp_data = abundance_matrix,
  threads = 4,
  show_observed = TRUE,
  show_extrapolation = TRUE
)

rare$fig
```

## Function Overview

| Task | Functions |
| --- | --- |
| Paired-end and sequence preprocessing | `pairend()`, `demultiplex()`, `dereplicate()`, `denoise()`, `merge_replicate()`, `count_records()` |
| FASTA and table conversion | `fasta_to_table()`, `format_to_dada()` |
| BLAST wrappers and databases | `blast_cmd()`, `blastn_cmd()`, `download_blast_db()`, `prepare_blast_database()` |
| VSEARCH wrappers | `vsearch_cmd()`, `vsearch_cluster()`, `vsearch_remove_chimeras()` |
| Taxonomy databases | `build_ncbi_taxonomy_database()`, `build_col_taxonomy_database()`, `get_taxid()`, `get_taxonomy()`, `find_lca()` |
| Reference building | `build_gb_reference()`, `build_gb_reference_online()`, `build_gb_reference_local()` |
| In-silico PCR | `find_primer_binding_sites()`, `in_silico()` |
| Decontamination and diversity | `decontamination()`, `rarefaction_extrapolation()` |

## Notes

- Use absolute paths for large reference databases and external command-line tools when possible.
- NCBI downloads and accession databases can require substantial storage space.
- Online NCBI queries should include a contact email and conservative parallel settings.
- Local BLAST and VSEARCH functions require the corresponding executables to be installed and available in `PATH`, or supplied through the relevant executable-path argument.

## License

MIT. See [LICENSE](LICENSE) and [LICENSE.md](LICENSE.md).
