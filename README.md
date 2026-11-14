# METAeDNA

METAeDNA is an R package for environmental DNA (eDNA) metabarcoding analysis,
from paired-end sequencing reads to taxon abundance matrices and diversity
estimates. It provides tools for read merging, primer and tag assignment,
denoising, blank correction, reference database construction, and taxonomic
classification.

PCR reactions are identified by primer and nucleotide-tag combinations,
including technical replicates and blanks. The example library uses paired
8-bp tags before amplicons are pooled for library preparation
(Taberlet et al., 2018); other tag lengths can be specified in the reaction table.

## Contents

- [Workflow overview](#workflow-overview)
- [1. Sequence processing and reaction counts](#sequence-processing)
- [2. Reference databases and taxonomic assignment](#reference-databases-and-taxonomy)
- [3. Sequencing depth and community diversity](#sequencing-depth-and-diversity)

<a id="workflow-overview"></a>

## Workflow overview


The workflow is divided into three tutorials. Run them from the same
working directory; intermediate RDS files connect the steps.

| Tutorial                                                          | Main inputs                                                                    | Saved outputs                                                         |
|-------------------------------------------------------------------|--------------------------------------------------------------------------------|-----------------------------------------------------------------------|
| [1. Sequence processing](#sequence-processing)             | Paired FASTQ, reaction definitions, blank and replicate design                 | `data/workflow_tables.rds`: cleaned reaction counts and control table |
| [2. Reference databases and taxonomy](#reference-databases-and-taxonomy) | Species checklist, local nt or downloaded reference records, cleaned sequences | `data/taxonomy_results.rds`: classifications linked by `SEQ`          |
| [3. Sequencing depth and diversity](#sequencing-depth-and-diversity)   | Reaction counts, control table, classification results                         | Sample RRA, community matrix, depth summaries, diversity and plots    |

Part 1 ends at corrected reaction counts. Part 2 builds references and
assigns taxonomy. Part 3 selects classified sequences through `SEQ`,
calculates RRA within each reaction, and merges technical replicates.
Each reaction’s RRA denominator is its total retained classified reads.
Sequencing-depth summaries use reaction counts, not RRA. Part 2 also
accepts independent BLAST tables.

![METAeDNA workflow](vignettes/Workflow.svg)

### About METAeDNA

METAeDNA processes environmental DNA metabarcoding reads, assigns
taxonomy, and calculates taxon abundances and diversity.

Reaction labels are nucleotide tags used for tag-based multiplexing
during library preparation, not target-gene markers. The example uses
paired 8-bp tags to distinguish PCR reactions before pooling. Primer and
tag combinations retain the link to samples, technical replicates, and
blanks. Other tag lengths can be specified in the reaction table.

> **Installation**
>
> METAeDNA requires R 4.0 or later, a C++ compiler, and zlib. OpenMP is
> used for multithreading when supported by the compiler. Install the R
> package from the supplied source archive:

``` r
# install.packages("remotes") # Run once if needed.
# install from a local package file
remotes::install_local("METAeDNA_1.0.0.tar.gz", upgrade = "never")
# install from GitHub
library(METAeDNA)
packageVersion("METAeDNA")
```

### Configure external software

Software releases are available for each operating system:

-   [VSEARCH releases](https://github.com/torognes/vsearch/releases):
    chimera removal and sequence clustering.
-   [NCBI BLAST+
    downloads](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/):
    reference database construction and similarity searches.

Extract or install the software first. Give the folders containing the
executables: `vsearch` for VSEARCH, and `blastn`, `blastdbcmd`, and
`makeblastdb` for BLAST+. On Windows these files end in `.exe`. For
BLAST+, this is normally its `bin` folder.

The paths below identify the executable folders. A single path is
sufficient when both programs share a folder.

``` r
tool_dirs <- c("software/blast/bin", "software/vsearch/bin")
# A shared executable folder instead:
# tool_dirs <- "software/bin"

Sys.setenv(PATH = paste(c(tool_dirs, Sys.getenv("PATH")), collapse = .Platform$path.sep))
Sys.which(c("vsearch", "blastn", "blastdbcmd", "makeblastdb"))
```

This sets the search path for the current R session. An empty result
from `Sys.which()` means that executable was not found in the supplied
folders. If the commands are already on `PATH`, only the check is
needed.

<a id="sequence-processing"></a>

## 1. Sequence processing and reaction counts


This tutorial takes paired FASTQ reads through merging, demultiplexing,
sequence filtering and blank correction. Technical replicates remain
separate until classification and RRA calculation in part 3. Run the
steps in order from the same working directory. Package example inputs
are read with `system.file()`; output files are written under `data/`.

``` r
library(METAeDNA)
```

### Step 1. Prepare inputs and outputs

> **Preparation**
>
> 1.  **forward and reverse FASTQ files**
> 2.  **a reaction-definition table containing reaction tags and primer
>     information**
> 3.  **a table containing technical replicates, reaction IDs, and blank
>     controls (control table)**
> 4.  **setting outputs**

***Note: gzip-compressed FASTQ files can be read directly.***

> **Setting input files:**

``` r
forward_fastq <- system.file( "extdata", "trnL_40k_1.fq.gz", package = "METAeDNA") # Forward reads
reverse_fastq <- system.file( "extdata", "trnL_40k_2.fq.gz", package = "METAeDNA") # Reverse reads
sample_table <- system.file( "extdata", "trnL_table.txt", package = "METAeDNA")  # sample table
control_table <- system.file( "extdata", "control_table.csv", package = "METAeDNA") # control table
```

> **Reaction table**
>
> The reaction-definition table links each PCR reaction to its
> experiment, tags, and primers. `demultiplex()` uses these fields to
> assign reads to reactions. The supplied library uses tag-based
> multiplexing: pairs of 8-bp oligonucleotide tags at the 5’ ends of
> primers label individual PCR reactions before pooling for library
> preparation (Taberlet et al., 2018). These tags distinguish reactions,
> including technical replicates and blanks; they are not different
> target-gene markers. The software accepts tag lengths defined in the
> table, not only 8 bp. The format is:

| experiment | reaction_id  | reaction_tag      | forward_primer    | reverse_primer         |
|------------|--------------|-------------------|-------------------|------------------------|
| trnL_run1  | lake_s1_r1   | ATCGATCG:TGGCTGGC | GGGCAATCCTGAGCCAA | CCATTGAGTCTCTGCACCTATC |
| trnL_run1  | lake_s1_r2   | GCTAGCTA:ACCGACCG | GGGCAATCCTGAGCCAA | CCATTGAGTCTCTGCACCTATC |
| trnL_run1  | pcr_blank_r1 | TTAATTAA:GGCCGGCC | GGGCAATCCTGAGCCAA | CCATTGAGTCTCTGCACCTATC |

Column meanings:

-   `experiment` identifies a library, sequencing run, or experimental
    group;
-   `reaction_id` identifies one PCR reaction within an experiment;
-   `reaction_tag` contains forward\_<tag:reverse_tag>;
-   `forward_primer` and `reverse_primer` contain nucleotide or IUPAC
    codes.

Inspect the table before processing:

``` r
sample_info <- read.delim(sample_table, sep = '\t', check.names = FALSE)
names(sample_info)
head(sample_info)
```

> **Control table and technical replicates**
>
> The control table links each biological sample to its associated
> controls, including but not limited to filtration, extraction, and PCR
> blanks, with the control type specified in a dedicated column.

For each sequence and blank type, decontamination finds the blank count
linked to the target sample.

Corrected counts are limited to zero, so subtraction cannot produce
negative abundance. remove_status controls whether per-PCR denoising
status columns are removed from the returned table.

| sample_id | reaction_id | group | category | EB  | FB  | PB  |
|-----------|-------------|-------|----------|-----|-----|-----|
| lake_1    | PCR1        | run1  | sample   | EB1 | FB1 | PB1 |
| lake_2    | PCR2        | run1  | sample   | EB1 | FB1 | PB1 |
| EB1       | PCR_EB      | run1  | EB       | NA  | NA  | NA  |
| FB1       | PCR_FB      | run1  | FB       | NA  | NA  | NA  |
| PB1       | PCR_PB      | run1  | PB       | NA  | NA  | NA  |

Each row represents one PCR reaction.

Column meanings:

-   `sample_id` is the biological sample or control name
-   `reaction_id` identifies its PCR
-   `group` identifies the library
-   `category` states whether the row is a sample or a control.
-   Additional columns, such as `EB`, `FB`, and `PB`, specify the
    controls associated with each sample and can be added or modified as
    needed for different experimental designs.

``` r
control_table <- read.csv(control_table, check.names = FALSE)
head(control_table)
```

------------------------------------------------------------------------

### Step 2. Merge paired reads

> `pairend()` evaluates candidate overlaps of at least `min_overlap`
> bases and merges paired reads. A quality-aware consensus is produced
> for the selected overlap. At a mismatch, the base with the higher
> Phred-derived correctness probability is retained; a tie is resolved
> using the estimated background base frequencies.

``` r
# Create the output directory and save merged paired-end reads as FASTQ.
dir.create("data", recursive = TRUE, showWarnings = FALSE)
paired_file <- "data/01_paired.fastq" # output pair-end file

METAeDNA::pairend(
  fastq_1 = forward_fastq,
  fastq_2 = reverse_fastq,
  out_file = paired_file,
  min_overlap = 10L,
  max_mismatches = 5L,
  # compress_output =  TRUE, # Compression is optional. For sequences, compression can save storage but increases processing time.
  # compress_level = 3,
  alpha = 1,
  beta = -1,
  n_threads = 8L,
  batch_size = 100000L
)

METAeDNA::count_records(paired_file)
```

> Paired-end reads are merged using a PEAR-based overlap scoring method
> [(Zhang et al., 2014)](https://doi.org/10.1093/bioinformatics/btt593).
> Candidate overlaps are evaluated using base quality scores and
> empirical nucleotide frequencies. The overlap with the highest
> assembly score is retained and accepted only when its PEAR overlap
> significance is $p < 0.01$.

*Note:*

-   Background frequencies of A, C, G, and T are estimated from the
    current processing batch. Forward reads and reverse-complemented
    reverse reads are combined so that frequencies correspond to the
    merged-read orientation. The number of read pairs used for each
    estimate is controlled by `batch_size`. Equal frequencies of 0.25
    are used when frequencies cannot be estimated.

-   `min_overlap` specifies the minimum overlap length. `max_mismatches`
    specifies the maximum absolute number of mismatched bases permitted
    in the selected overlap; the default is 3. A read pair must satisfy
    both the PEAR significance test and the mismatch limit.

-   `n_threads` sets the number of processing threads.

-   Increasing `batch_size` can reduce input overhead but also increases
    memory use.

-   `reorder_fastq()` restores mate order when the paired files contain
    the same records in different orders:

``` r
METAeDNA::reorder_fastq(fastq_1 = forward_fastq, fastq_2 = reverse_fastq, output_fastq = "data/reverse_reads_reordered.fastq.gz")
```

------------------------------------------------------------------------

### Step 3. Demultiplex reads

> Demultiplexing searches for the configured `primer pair` and optional
> `reaction_tag`. Primer mismatches are controlled by `primer_mismatch`.
> `reaction_tag` must match exactly. Identified and unidentified reads
> are both retained.

The reaction-definition file may be comma-separated or tab-delimited.
For every read, demultiplex first searches for a compatible forward and
reverse primer pair. Primer recognition uses a bitwise Shift-Or search,
which supports IUPAC primer symbols and permits up to `primer_mismatch`
substitutions. A read is not assigned when a required primer is absent
or exceeds that mismatch limit.

When `with_tag = TRUE`, the sequence immediately associated with each
primer is also compared with the forward and reverse tags from
`reaction_tag`. Tag matching is exact: `primer_mismatch` does not allow
errors in tags. A successful record is trimmed to the sequence between
the matched primer boundaries and written as FASTA. Its header receives
a demultiplex JSON block containing the `experiment`, `reaction_id`,
observed `reaction_tag`, `primer` sequences, and forward/reverse
`mismatch` counts. Unassigned reads are written unchanged to
`unidentified_output` so that assignment failure can be diagnosed
instead of silently discarded.

``` r
# Reads assigned to a reaction, trimmed between the primers.
identified_file <- "data/02_identified.fasta"
# Reads without a matching primer/tag combination, retained for inspection.
unidentified_file <- "data/02_unidentified.fasta"

METAeDNA::demultiplex(
  fastq_file = paired_file,
  reaction_table = sample_table,
  identified_output = identified_file,
  unidentified_output = unidentified_file,
  n_threads = 8L,
  with_tag = TRUE,
  primer_mismatch = 2L
)

METAeDNA::count_records(identified_file)
METAeDNA::count_records(unidentified_file)
```

> *Note:* If almost all reads are unidentified, check:

-   whether the primers are in the expected orientation;
-   whether reaction tags are correct;
-   whether paired tags contain a colon;
-   whether with_tag matches the library design;
-   whether primer_mismatch is too strict.

------------------------------------------------------------------------

### Step 4. Collapse identical sequences

> `collapse_sequences()` stores one representative of each nucleotide
> sequence and aggregates its counts by experiment and sample.

The function reads the metadata created by pairend and demultiplex.
Records with **identical** nucleotide sequences are collapsed even when
they came from different samples. Their total abundance and the
abundance within every experiment–sample combination are summed. This
reduces file size and avoids performing the later similarity search
repeatedly on identical reads.

``` r
# One record per distinct sequence, with counts combined across reactions.
unique_file <- "data/03_unique.fasta"

METAeDNA::collapse_sequences(
  fasta_file = identified_file,
  output_fasta = unique_file,
  n_threads = 8L,
  batch_size = 100000L
)

METAeDNA::count_records(unique_file)
```

------------------------------------------------------------------------

### Step 5. Denoise sequences

> Denoising first applies sequence-length and total-count filters. It
> then builds a local neighbor graph from substitutions and indels.
> Lower-abundance sequences connected to more abundant neighbors are
> classified separately in each PCR. This PCR-aware abundance and
> neighborhood strategy is based on the algorithm distributed with
> `OBITools`.

`sequence_min_L` and `sequence_max_L` set the expected amplicon-length
range.

``` r
# Sequences retained after abundance filtering and denoising.
denoised_file <- "data/04_denoised.fasta"

METAeDNA::denoise(
  fasta_file = unique_file,
  output_fasta = denoised_file,
  sequence_min_L = 10L,
  sequence_max_L = 1000L,
  min_count = 10L,
  threshold = 0.5,
  max_diff = 1L,
  max_indel = 0L,
  keep_internal = FALSE,
  n_threads = 8L,
  enable_avx2 = TRUE
)

METAeDNA::count_records(denoised_file)
```

Two sequences are neighbors when they co-occur in a PCR and satisfy the
`max_diff` substitution limit or the max_indel indel limit. For a
neighboring pair with abundances $C_{\mathrm{high}}$ and
$C_{\mathrm{low}}$, an abundance relationship is created when:

$$
\frac{C_{\mathrm{low}}}{C_{\mathrm{high}}}<\mathrm{threshold}.
$$

The graph is evaluated separately for every experiment–sample
combination. Sequences are classified as `head`, `internal`, or
`singleton` according to their local relationships [(Boyer et al.,
2016)](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12428). The
`keep_internal` parameter controls whether internal variants remain in
the output. Sequence comparisons use bit-encoded operations;
`enable_avx2` allows AVX2 acceleration on compatible processors and
falls back when it is unavailable.

------------------------------------------------------------------------

### Step 6. Remove chimeras (optional)

> PCR chimeras arise when an incompletely extended molecule anneals to a
> different template and is extended again. Such artificial recombinants
> can inflate apparent diversity. The `vsearch_remove_chimeras` function
> calls the UCHIME implementation in VSEARCH and writes non-chimeric
> sequences to a new FASTA file.

Install VSEARCH and make sure its executable is on PATH, or provide its
directory with vsearch_exe_path. Keep the pre-chimera file because
parameter or software changes can alter the result.

``` r
# Denoised sequences retained after chimera removal.
nonchimera_file <- "data/05_nonchimera.fasta"

METAeDNA::vsearch_remove_chimeras(
  fasta_file = denoised_file,
  output_fasta =  nonchimera_file,
  vsearch_exe_path = NULL # Use the executable found on PATH.
)

METAeDNA::count_records(nonchimera_file)
```

------------------------------------------------------------------------

### Step 7. Convert FASTA to a sequence table

> `fasta_to_table()` converts FASTA sequences and header metadata into a
> table with one row per sequence and experiment, and count columns for
> PCR reactions.

``` r
sequence_table <- METAeDNA::fasta_to_table(
  fasta_file = nonchimera_file,
  n_threads = 8L
)

dim(sequence_table)
head(sequence_table)
```

The table contains `sequence` `identifiers`, `nucleotide sequences`,
`experiments`, `samples`, `counts`, and `denoising status` extracted
from FASTA JSON metadata.

| Column        | Content                                                        |
|---------------|----------------------------------------------------------------|
| `experiment`  | Experiment or library identifier                               |
| `ID`          | Sequence identifier                                            |
| `SEQ`         | Nucleotide sequence; the key used to join classification later |
| `count_TAG`   | Integer read counts for one reaction                           |
| `total_count` | Total read count reported for the sequence                     |

After replicate merging, `count_TAG` columns are replaced by sample RRA
columns. For multiple experiments, both tables must include
`experiment`.

------------------------------------------------------------------------

### Step 8. Decontamination

`decontamination()` uses blank reactions to correct every sample
reaction. `control_names` selects the control columns to use, for
example `"FB"` for filtration blanks or `c("EB", "FB", "PB")` for
extraction, filtration, and PCR blanks. Each selected column contains
the `sample_id` of the associated blank. With `method = "maximum"`, the
largest associated blank count is subtracted and negative results are
set to zero.

The example contains ten sequences and ten reactions:

``` r
decontamination_counts <- matrix(
  (seq_len(100) * 7) %% 31,
  nrow = 10,
  dimnames = list(NULL, paste0("count_R", 1:10))
)
mini_sequences <- data.frame(
  SEQ = paste0("sequence_", 1:10),
  decontamination_counts,
  check.names = FALSE
)
mini_controls <- data.frame(
  group = "G1",
  sample_id = c(paste0("lake_", 1:5), paste0("blank_", 1:5)),
  reaction_id = paste0("R", 1:10),
  category = c(rep("sample", 5), rep("FB", 5)),
  FB = c(paste0("blank_", 1:5), rep(NA, 5))
)
mini_clean <- METAeDNA::decontamination(
  sequence_table = mini_sequences,
  control_table = mini_controls,
  control_names = "FB",
  method = "maximum"
)
mini_clean
```

This creates a 10-sequence × 10-reaction table. Reactions R1–R5 are
samples, and R6–R10 are their corresponding filtration blanks. Blank
reaction columns are removed after correction. For the complete dataset:

``` r
decontamination_table <- METAeDNA::decontamination(sequence_table = sequence_table,
  control_table = control_table,
  control_names = c("EB", "FB", "PB"))
```

------------------------------------------------------------------------

### Step 9. Sequence clustering (Optional)

> `vsearch_cluster()` calls `VSEARCH` to cluster sequences at the
> requested identity. With `identity` = 0.97, sequences meeting the 97%
> identity criterion are represented by a cluster. Clustering is
> optional: retain denoised sequence variants when ASV-level resolution
> is required.

``` r
cluster_input <- decontamination_table
cluster_input$SEQ_ID <- paste0("seq_", seq_len(nrow(cluster_input)))
clustered_table <- METAeDNA::vsearch_cluster(
  sequence_table = cluster_input,
  output_path =  "data/clustered_sequences",
  identity = 0.98
)
```

VSEARCH also implements UNOISE3 through mode = “cluster_unoise”. Unlike
the identity modes, UNOISE uses both sequence distance and the abundance
ratio between a low-abundance sequence and a potential centroid.
Therefore identity is not used in this mode. For `sequence_table` input,
`vsearch_cluster()` writes the required size annotations from
total_count automatically.

------------------------------------------------------------------------

### Technical-replicate merging: calculation example

> Technical PCR replicates reduce stochastic amplification effects but
> must be combined before sample-level analysis. merge_replicate
> supports three strategies: mean_rra averages relative read abundance,
> max_rra keeps its maximum, and both_occur requires occurrence in at
> least repeat_times replicates.

Each `reaction_id` identifies one PCR reaction. Reactions with the same
`experiment`, `group`, and `sample_id` are treated as technical
replicates.

The following artificial counts illustrate replicate merging. In the
analysis workflow, merging takes place after classification in [part
3](#sequencing-depth-and-diversity), not at this stage.

``` r
replicate_counts <- matrix(
  (seq_len(100) * 11) %% 37,
  nrow = 10,
  dimnames = list(NULL, paste0("count_R", 1:10))
)
mini_replicates <- data.frame(
  SEQ = paste0("sequence_", 1:10),
  replicate_counts,
  check.names = FALSE
)
replicate_design <- data.frame(
  group = "G1",
  sample_id = rep(paste0("lake_", 1:5), each = 2),
  reaction_id = paste0("R", 1:10),
  category = "sample"
)
mini_merged <- METAeDNA::merge_replicate(
  sequence_table = mini_replicates,
  control_table = replicate_design,
  method = "mean_rra",
  remove_all_zero = FALSE
)
mini_merged
```

This creates 10 reaction columns for five samples. Each pair of
reactions is normalized separately and then averaged into one sample
column.

The RRA table has one column per biological sample. Empty reactions
contribute zero. Use `replicate_number = 3` to divide by three even when
fewer reactions are listed; the default uses the number of reactions
assigned to each sample.

### Step 10. Save reaction counts

Reaction counts are retained for classification, sequencing-depth
calculations, and subsequent RRA calculation. The optional clustering
example above does not replace the main ASV table. To analyse clusters
instead, use the clustered count table consistently for merging,
classification, and depth calculations.

``` r
saveRDS(
  list(sequence_table = sequence_table,
       decontamination_table = decontamination_table,
       control_table = control_table),
  "data/workflow_tables.rds"
)
```

Continue with [reference databases and
taxonomy](#reference-databases-and-taxonomy).

### References

Taberlet P, Bonin A, Zinger L, Coissac E. (2018). *Environmental DNA:
For Biodiversity Research and Monitoring*. Oxford University Press,
Chapter 6, “DNA amplification and multiplexing”.
<https://doi.org/10.1093/oso/9780198767220.001.0001>.

Boyer F, Mercier C, Bonin A, Le Bras Y, Taberlet P, Coissac E. (2016).
obitools: a unix-inspired software package for DNA metabarcoding.
*Molecular Ecology Resources*, 16, 176–182.

Zhang J, Kobert K, Flouri T, Stamatakis A. (2014). PEAR: a fast and
accurate Illumina Paired-End reAd mergeR. *Bioinformatics*, 30, 614–620.

<a id="reference-databases-and-taxonomy"></a>

## 2. Reference databases and taxonomic assignment


Build taxonomy resources and reference sequences, then classify the
sequences saved in [part 1](#sequence-processing). All paths are
relative to the same working directory. Database construction can also
be run on its own with a species checklist, independently of the
sequencing workflow.

### Step 1. Load sequence tables

``` r
library(METAeDNA)
library(Biostrings)
library(dplyr)

workflow_tables <- readRDS("data/workflow_tables.rds")
decontamination_table <- workflow_tables$decontamination_table
```

Paths are introduced in the steps that use them and are relative to the
current working directory. Keep the same working directory throughout
the tutorial.

### Step 2. Build the NCBI taxonomy database

The NCBI taxonomy database is separate from the BLAST nt database. The
nt database stores nucleotide sequences for similarity searches; the
taxonomy database stores names, taxonomic nodes, obsolete-taxid
redirects, and optionally accession-to-taxid mappings.

`build_ncbi_taxonomy_database()` downloads NCBI taxonomy files and
imports names and lineages into SQLite. The database supplies TaxIDs and
lineage information for reference retrieval and BLAST annotation.

| Argument           | NCBI file                   | SQLite content                            | Required for                                 |
|--------------------|-----------------------------|-------------------------------------------|----------------------------------------------|
| taxonomy = TRUE    | taxdump.tar.gz              | nodes, names, name aliases, merged taxids | name and lineage queries, COL–NCBI crosswalk |
| accession = TRUE   | nucl_gb.accession2taxid.gz  | GenBank nucleotide accession to taxid     | converting BLAST accessions to taxids        |
| include_WGS = TRUE | nucl_wgs.accession2taxid.gz | WGS accession to taxid                    | BLAST results containing WGS accessions      |

The source files come from the [NCBI taxonomy FTP
directory](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) and its
[accession2taxid
directory](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/).
They are retained in taxdump_dir and reused on later runs.

``` r
ncbi_dir <- "data/reference/ncbi_taxonomy" # Downloaded NCBI taxonomy source files.
taxonomy_db <- "data/reference/ncbi_taxonomy/taxonomy.sqlite" # SQLite database reused for lineage queries.
dir.create(ncbi_dir, recursive = TRUE, showWarnings = FALSE)

METAeDNA::build_ncbi_taxonomy_database(
  taxonomy = TRUE,
  accession = FALSE,
  include_WGS = FALSE,
  taxdump_dir = ncbi_dir,
  sql_path = taxonomy_db,
  download = TRUE,
  resume = TRUE,
  overwrite = TRUE,
  overwrite_scope = "database" # overwrite the entire SQL database.
)
```

***Query lineages:***

``` r
# Retrieve queries by taxid; taxid can be specified as a integer vector.
METAeDNA::get_taxonomy(
  sql_path = taxonomy_db,
  taxids = c(9606, 1779274)
)
# Retrieve queries by name; query names can be specified as a character vector.
METAeDNA::get_taxonomy(
  sql_path = taxonomy_db,
  taxon_names = c("Homo sapiens", "Acer acutum")
)
```

### Step 3. Add accession-to-TaxID mappings

Add a `accession-taxid` table to current SQL database. It can require
more than 50 GB of free space for the downloaded file, SQLite content,
and indexes. Set `include_WGS = TRUE` only when WGS accessions are
required because it downloads and imports an additional large file.

With `download = FALSE`, the required NCBI files must already exist in
`taxdump_dir`.

``` r
METAeDNA::build_ncbi_taxonomy_database(
  taxonomy = FALSE,
  accession = TRUE,
  include_WGS = FALSE,
  taxdump_dir = ncbi_dir,
  sql_path = taxonomy_db,
  download = TRUE,
  resume = TRUE,
  overwrite = TRUE, # allows an existing target table to be replaced if a table with the same name is already present in the SQL database.
  overwrite_scope = "selected", # limits the overwrite operation to the selected table(s). Other existing tables in the SQL database are preserved and remain unchanged.
  chunk_size = 1000000L,
  accession_index_columns = c("accession", "accession.version")
)
```

***Accession-taxid mapping:***

``` r
# convert accession to taxid
accession_taxids <- METAeDNA::get_taxid(
  sql_path = taxonomy_db,
  accessions = c("NC_045512.2", "NC_000001.11")
)
# get lineage
METAeDNA::get_taxonomy(
  sql_path = taxonomy_db,
  taxids = accession_taxids$taxid
)
```

### Step 4. Build the COL name database

Catalogue of Life (COL) is used here to standardize names before
querying NCBI. Download a COL Data Package (ColDP) release and extract
NameUsage.tsv. Use a dated annual or monthly release for reproducible
work and record its release date. The [COL download
page](https://www.catalogueoflife.org/data/download) provides full
releases.

For a quick installation, the following downloads the current COL Base
Release in ColDP format. Replace the `latest` archive with a dated
release for a published analysis, because the contents of the latest
release change over time. The Extended Release can be obtained by using
`xr_latest_coldp.zip` instead.

``` r
col_zip <- "data/reference/COL/latest_coldp.zip" # Downloaded COL archive.
col_db <- "data/reference/COL/catalogue_of_life.sqlite" # SQLite database reused for name matching.
dir.create("data/reference/COL", recursive = TRUE, showWarnings = FALSE)

download.file(
  "https://download.checklistbank.org/col/latest_coldp.zip",
  destfile = col_zip,
  mode = "wb"
)
unzip(col_zip, exdir = "data/reference/COL") # Folder containing the extracted COL tables.
```

Build the local SQLite database:

``` r
METAeDNA::build_col_taxonomy_database(
  col_data = "data/reference/COL/NameUsage.tsv", # Name-usage table extracted from the COL archive.
  db_path = col_db,
  overwrite = TRUE
)
```

The database retains accepted names, synonyms, authorship, rank,
identifiers, and lineage fields. Exact matching is preferred; optional
`fuzzy` matching should be inspected because a short edit distance does
not prove taxonomic identity.

***Run a query:***

``` r
col_matches <- METAeDNA::col_match_species(
  c("Quercus robur", "Felis leo"),
  db_path = col_db,
  fuzzy = TRUE,
  max_dist = 1
)
# Query the built-in Chinese aquatic plant checklist
aquatic_names <- read.csv(
  system.file("extdata", "aquatic_species.csv", package = "METAeDNA") # Built-in aquatic species checklist.
)

col_matches <- METAeDNA::col_match_species(
  aquatic_names$Latin,
  db_path = col_db,
  fuzzy = TRUE,
  max_dist = 1
)
col_matches
```

------------------------------------------------------------------------

### Step 5. Reference database construction

> The built-in aquatic plant checklist supplies species names for the
> trnL reference database. Other checklists can supply names in the same
> format.

``` r
aquatic_names_with_taxid <- METAeDNA::get_taxonomy(
  sql_path = taxonomy_db,
  taxon_names = aquatic_names$Latin
)
```

`col_ncbi_crosswalk()` links COL accepted names and synonyms to NCBI
TaxIDs.

``` r
name_crosswalk <- METAeDNA::col_ncbi_crosswalk(
  taxon_names = aquatic_names$Latin,
  col_db = col_db,
  ncbi_db = taxonomy_db
)
```

Two routes are available. A complete local nt database requires
substantial disk space and preparation time, but can be reused for
different taxa and markers without downloading the sequences again. A
marker-specific database is smaller and quicker to prepare, but its
contents are limited by the species list and marker query used during
retrieval.

------------------------------------------------------------------------

#### Route 1: use the local NCBI nt database

The compressed nt database is very large. Confirm available disk
capacity before starting. Download and extract it once:

``` r
blast_download_dir <- "data/reference/blast_downloads" # Archives downloaded below and then extracted.
nt_database <- "data/reference/blast_database/nt/nt" # Local BLAST database prefix, not a directory.

METAeDNA::download_blast_db(
  database_path = blast_download_dir,
  md5_threads = 4L
)
METAeDNA::prepare_blast_database(
  file_path = blast_download_dir,
  database_path = "data/reference/blast_database/nt", # Destination for extracted BLAST database files.
  n_threads = 8L
)
```

Use the `taxid` to select the required taxa, filter record titles by
marker annotation, and format the resulting FASTA as a smaller BLAST
database:

``` r
# Select records by TaxID and title from local nt, then format the FASTA for BLAST.
local_reference <- METAeDNA::build_gb_reference(
  mode = "local",
  database = nt_database,
  taxids = name_crosswalk$ncbi_taxid, #
  regex_expression = METAeDNA::marker_regex$trnL,
  output_dir = "data/reference/local_trnL", # Marker reference extracted from local nt.
  min_length = 10,
  max_length = 500,
  makeblastdb = TRUE,
  database_name = "trnL_reference"
)
local_reference$fasta
local_reference$blast_database
```

`taxids` at the order, family, and genus levels are also accepted. The
function automatically expands these higher-level `taxids` to include
the TaxIDs of their descendant taxa.

`regex_expression` is matched against nt record titles. Include the
common gene symbol and annotation wording used for the marker. The
complete nt database must include the NCBI taxonomy support files
required by `blastdbcmd -taxidlist`.

#### Route 2: search and download a marker-specific database

Use the same `taxids` (or `taxon_names`) to search NCBI Nucleotide. The
online function writes a FASTA file; `makeblastdb` then formats it as a
BLAST database:

``` r
# Narrow candidates on NCBI, check titles, download FASTA, then format for BLAST.
online_reference <- METAeDNA::build_gb_reference(
  mode = "online",
  taxids = name_crosswalk$ncbi_taxid %>% na.omit() %>% unique(),
  ncbi_query = METAeDNA::marker_entrez_query$trnL,
  # regex_expression = METAeDNA::marker_regex$trnL,
  output_dir = "data/reference/downloaded_trnL", # Reference sequences downloaded from NCBI.
  max_records = Inf,
  batch_size = 500,
  min_length = 10,
  max_length = 500,
  # makeblastdb = TRUE, # if TRUE, use the makeblastdb command from BLAST+ to build a database.
  database_name = "trnL_reference"
)
online_reference$fasta
online_reference$blast_database
online_reference$taxid_map
```

`ncbi_query` uses NCBI Entrez syntax to find candidates; the optional
`regex_expression` then checks their titles using the same rule as the
local route. Only titles passing both steps are sent to EFetch. The
Entrez query should be broad enough to include all records of interest:
the regex cannot recover records excluded earlier. For another marker,
change both the Entrez query and the title regex. To query names,
replace `taxids` with
`taxon_names = c("Quercus robur", "Quercus ilex")`.

`batch_size` is the number of accessions fetched per FASTA request (up
to 500). Titles are retrieved separately in batches of up to 10,000.
`max_records = Inf` checks all records returned by the Entrez query.
Supply an email address for NCBI E-utilities and, when available, an API
key.

`-parse_seqids` indexes the accession in each FASTA header; `-taxid_map`
assigns NCBI TaxIDs to those accessions. These are separate
`makeblastdb` options.

### Step 6. Test primers by in-silico PCR

> in_silico predicts which reference sequences can be amplified by a
> given primer pair. Predicted amplicons support primer screening and
> reference database construction. Predicted primer coverage does not
> replace experimental validation.

``` r
reference_fasta <- online_reference$fasta

predicted_amplicons <- METAeDNA::in_silico(
  fasta_file = reference_fasta,
  forward = "GGGCAATCCTGAGCCAA", # trnL g
  reverse = "CCATTGAGTCTCTGCACCTATC", # trnL h
  mis_match = 2,
  trim = TRUE,
  try_reverse = TRUE,
  mark_reverse = TRUE,
  chunk_size = 10000,
  show_progress = TRUE
)

length(predicted_amplicons)
Biostrings::width(predicted_amplicons)
```

With `trim = TRUE`, the returned DNAStringSet contains the predicted
amplicon between the primer boundaries. Set `trim = FALSE` to retain the
matched primers. `try_reverse = TRUE` also tests references stored in
the opposite orientation. With `mark_reverse = TRUE`, these amplicon
names end in \[REVERSE\]; set it to FALSE to retain the original
reference names. For several primer pairs, pass equal-length vectors to
forward and reverse. The result is a named list with one DNAStringSet
per primer pair:

> The function supports the simultaneous comparison of multiple primer
> pairs.

``` r
primer_comparison <- METAeDNA::in_silico(
  fasta_file = reference_fasta,
  forward = c("ATGTCACCACAAACAGAGACTAAAGC", "GGGCAATCCTGAGCCAA"),
  reverse = c("GTAAAATCAAGTCCACCRCG", "CCATTGAGTCTCTGCACCTATC"),
  mis_match = 2,
  trim = TRUE,
  n_threads = 4
)
```

***Note:*** Target sequence(s), find the primer-binding location:

``` r
METAeDNA::find_primer_binding_sites(
  target_sequence = reference_fasta,
  primer = "CCATTGAGTCTCTGCACCTATC",
  max_mismatch = 2
)
```

#### Build the trimmed BLAST database

Retain one amplicon per accession and remove the reverse-orientation
label before indexing sequence IDs. The online reference supplies the
accession-to-TaxID map.

``` r
trimmed_fasta <- "data/reference/downloaded_trnL/trimmed_trnL_reference.fas"
trimmed_db <- "data/reference/downloaded_trnL/trimmed_trnL_reference"

trimmed_seqs <- predicted_amplicons
names(trimmed_seqs) <- sub("\\[REVERSE\\]$", "", names(trimmed_seqs))
accessions <- sub(" .*$", "", names(trimmed_seqs))
trimmed_seqs <- trimmed_seqs[!duplicated(accessions)]

Biostrings::writeXStringSet(trimmed_seqs, trimmed_fasta)

METAeDNA::blast_cmd(
  base_cmd = "makeblastdb",
  args = c(
    "-in", trimmed_fasta,
    "-dbtype", "nucl",
    "-parse_seqids",
    "-taxid_map", online_reference$taxid_map,
    "-out", trimmed_db
  )
)
```

### Step 7. Search reference databases

Use one query per distinct cleaned sequence. `query_key` connects the
short FASTA IDs to `SEQ`, including sequences occurring in several
experiments.

``` r
query_sequences <- unique(decontamination_table$SEQ)
query_set <- DNAStringSet(query_sequences)
names(query_set) <- paste0("query_", seq_along(query_set))
query_fasta <- "data/reference/query_sequences.fasta"
writeXStringSet(query_set, query_fasta)
query_key <- data.frame(qseqid = names(query_set), SEQ = query_sequences)

nt_hits <- "data/reference/blast_nt.tsv"
marker_hits <- "data/reference/blast_marker.tsv"
trimmed_hits <- "data/reference/blast_trimmed.tsv"

METAeDNA::blastn_cmd(
  fasta_file = query_fasta, database = nt_database,
  output = nt_hits, max_target_seqs = 20, n_threads = 4
)
METAeDNA::blastn_cmd(
  fasta_file = query_fasta, database = online_reference$blast_database,
  output = marker_hits, max_target_seqs = 20, n_threads = 4
)
METAeDNA::blastn_cmd(
  fasta_file = query_fasta, database = trimmed_db,
  output = trimmed_hits, max_target_seqs = 20, n_threads = 4
)
```

These searches use full nt, the downloaded marker reference, and its
trimmed amplicons, respectively. Full nt can take much longer. Use
`local_reference$blast_database` instead of the online reference to
classify against the marker sequences extracted from local nt.

#### Run a BLAST command independently

`blast_cmd()` accepts an executable name and its command-line arguments.
It returns the exit status; the search results are written to `-out`.
The database argument is the database prefix, not its directory.

``` r
status <- METAeDNA::blast_cmd(
  base_cmd = "blastn",
  args = c(
    "-query", shQuote(query_fasta),
    "-db", shQuote(nt_database),
    "-out", shQuote(nt_hits),
    "-outfmt", shQuote("6 std qcovs staxids"),
    "-max_target_seqs", "20",
    "-num_threads", "4"
  )
)
```

For downstream `get_taxonomy()`, the BLAST table must have these 14
columns in this order, without a header when supplied as a file:

| Columns                            | Meaning                                      |
|------------------------------------|----------------------------------------------|
| `qseqid`, `sseqid`                 | Query ID and reference sequence ID           |
| `pident`, `length`                 | Percentage identity and alignment length     |
| `mismatch`, `gapopen`              | Mismatch and gap-opening counts              |
| `qstart`, `qend`, `sstart`, `send` | Query and reference alignment coordinates    |
| `evalue`, `bitscore`               | Alignment significance and score             |
| `qcovs`, `staxids`                 | Query coverage (%) and reference NCBI TaxIDs |

### Step 8. Retrieve taxonomy and assign LCA

`get_taxonomy()` reads a BLAST file or accepts a data frame with the
columns above. It adds NCBI lineage columns. Optional `col_path`
reconciles species names with COL accepted names and synonyms. Missing
TaxIDs are looked up by accession, which requires accession-to-TaxID
tables in `taxonomy_db`.

`find_lca()` accepts the resulting annotated **data frame**. The
thresholds below are illustrative; choose marker-appropriate values.
Omit `identity_threshold` to apply no additional rank-specific identity
cutoff.

``` r
identity_cutoffs <- c(species = 100, genus = 98, family = 95)

nt_lineages <- METAeDNA::get_taxonomy(
  sql_path = taxonomy_db, blast_result = nt_hits, col_path = col_db
)
marker_lineages <- METAeDNA::get_taxonomy(
  sql_path = taxonomy_db, blast_result = marker_hits, col_path = col_db
)
trimmed_lineages <- METAeDNA::get_taxonomy(
  sql_path = taxonomy_db, blast_result = trimmed_hits, col_path = col_db
)

nt_lca <- METAeDNA::find_lca(nt_lineages, identity_threshold = identity_cutoffs)
marker_lca <- METAeDNA::find_lca(marker_lineages, identity_threshold = identity_cutoffs)
trimmed_lca <- METAeDNA::find_lca(trimmed_lineages, identity_threshold = identity_cutoffs)

nt_taxonomy <- left_join(query_key, nt_lca, by = "qseqid")
marker_taxonomy <- left_join(query_key, marker_lca, by = "qseqid")
trimmed_taxonomy <- left_join(query_key, trimmed_lca, by = "qseqid")

nt_annotated <- left_join(decontamination_table, nt_taxonomy, by = "SEQ")
marker_annotated <- left_join(decontamination_table, marker_taxonomy, by = "SEQ")
trimmed_annotated <- left_join(decontamination_table, trimmed_taxonomy, by = "SEQ")
```

Queries with no hit retain missing taxonomy after the join.
Above-species LCA assignments receive labels such as `Quercus sp_1`;
these distinguish sequence groups but do not establish named species.

#### Input from a BLAST data frame

``` r
blast_columns <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs", "staxids"
) # an input data.frame should contain above columns.
blast_table <- read.delim(
  nt_hits, header = FALSE, col.names = blast_columns,
  quote = "", comment.char = "", check.names = FALSE
)
annotated_hits <- METAeDNA::get_taxonomy(
  sql_path = taxonomy_db, blast_result = blast_table
)
```

`get_taxonomy()` also accepts `taxids = c(9606, 3702)` or
`taxon_names = c("Homo sapiens", "Arabidopsis thaliana")`. Supply
exactly one of `taxids`, `taxon_names`, and `blast_result`.

#### Input from an annotated table to find_lca()

Each row represents one hit. Include `qseqid`, `sseqid`, `staxids`,
`pident`, `qcovs`, `evalue`, and all seven lineage columns. Unknown
ranks should be `NA`, rather than invented taxon names.

``` r
custom_hits <- data.frame(
  qseqid = rep(paste0("query_", 1:6), each = 2),
  sseqid = sprintf("example_%02d", 1:12),
  staxids = c(9606, 9606, 9598, 9597, 9606, 9598,
              9685, 9685, 9689, 9694, 9685, 9689),
  pident = rep(c(100, 99, 97, 100, 99, 97), each = 2),
  qcovs = 100, evalue = 1e-30,
  kingdom = "Metazoa", phylum = "Chordata", class = "Mammalia",
  order = rep(c("Primates", "Carnivora"), each = 6),
  family = rep(c("Hominidae", "Felidae"), each = 6),
  genus = c("Homo", "Homo", "Pan", "Pan", "Homo", "Pan",
            "Felis", "Felis", "Panthera", "Panthera", "Felis", "Panthera"),
  species = c(
    "Homo sapiens", "Homo sapiens",
    "Pan troglodytes", "Pan paniscus",
    "Homo sapiens", "Pan troglodytes",
    "Felis catus", "Felis catus",
    "Panthera leo", "Panthera tigris",
    "Felis catus", "Panthera leo"
  )
)
custom_lca <- METAeDNA::find_lca(custom_hits)
custom_lca
```

The sequence IDs in this example are illustrative. Each query has two
equally significant hits, so both are retained when calculating the LCA.

| Query     | Agreement between hits          | Expected assignment     |
|-----------|---------------------------------|-------------------------|
| `query_1` | Same species                    | Species: *Homo sapiens* |
| `query_2` | Different species of *Pan*      | Genus: *Pan*            |
| `query_3` | *Homo* and *Pan*                | Family: Hominidae       |
| `query_4` | Same species                    | Species: *Felis catus*  |
| `query_5` | Different species of *Panthera* | Genus: *Panthera*       |
| `query_6` | *Felis* and *Panthera*          | Family: Felidae         |

Genus and family assignments receive a `sp_` label followed by a number.
These examples use agreement between hits, without additional identity
thresholds.

With different input column names, specify `seqs_ID`, `taxid_column`,
`accession_column`, `identity_column`, and `coverage_column`. For a
different rank hierarchy, provide `tax_levels` from highest to lowest.
`identity_threshold` must be a named numeric vector, such as
`c(species = 100, genus = 98, family = 95)`.

### Step 9. Classify with DADA2

`format_to_dada()` converts the accession-labelled trimmed FASTA into a
DADA2 reference with semicolon-separated taxonomy in its headers. The
SQLite database must include accession-to-TaxID mappings. DADA2 is an
optional package.

``` r
dada_reference <- "data/reference/trnL_in_silico_dada2.fasta"
METAeDNA::format_to_dada(
  fasta_file = trimmed_fasta,
  output_fasta = dada_reference,
  sql_path = taxonomy_db
)
dada_assignment <- dada2::assignTaxonomy(
  seqs = query_sequences, refFasta = dada_reference,
  multithread = 4, tryRC = TRUE
)
dada_taxonomy <- as.data.frame(dada_assignment, stringsAsFactors = FALSE)
dada_taxonomy$SEQ <- rownames(dada_taxonomy)
dada_taxonomy$scientific_name <- gsub("_", " ", dada_taxonomy$Species)
dada_annotated <- left_join(decontamination_table, dada_taxonomy, by = "SEQ")
```

DADA2 joins directly by sequence. In the reference produced here,
`Species` contains the full scientific name. Missing species assignments
remain unassigned for the species-level matrix.

[sequencing depth and community diversity](#sequencing-depth-and-diversity).

<a id="sequencing-depth-and-diversity"></a>

## 3. Sequencing depth and community diversity


Use reaction counts for sequencing depth and rarefaction, and sample RRA
values for the community matrix and diversity. This tutorial reads the
outputs saved in [part 1](#sequence-processing) and [part
2](#reference-databases-and-taxonomy).

### Step 1. Load reaction counts and classifications

``` r
library(METAeDNA)
library(ggplot2)
library(tidyr)

workflow_tables <- readRDS("data/workflow_tables.rds")
taxonomy_results <- readRDS("data/taxonomy_results.rds")
decontamination_table <- workflow_tables$decontamination_table
control_table <- workflow_tables$control_table
trimmed_taxonomy <- taxonomy_results$trimmed_taxonomy
```

The classification table must contain `SEQ` and `scientific_name`, with
one classification per distinct sequence. To use another route, set
`trimmed_taxonomy <- taxonomy_results$nt_taxonomy`,
`taxonomy_results$marker_taxonomy`, or `taxonomy_results$dada_taxonomy`.

### Step 2. Calculate RRA and construct the community matrix

The classification table selects retained sequences from the corrected
reaction counts. `merge_replicate()` divides each retained count by the
total retained reads in that reaction, then averages the reaction RRA
values for each sample. Unassigned reads are excluded from both
numerator and denominator. `merge_to_matrix()` sums sequences sharing a
scientific name and returns samples in rows and taxa in columns.

`control_table` identifies sample-abundance columns from `group`,
`sample_id`, and, for multiple experiments, `experiment`; other numeric
columns such as `total_count` are excluded. `taxonomy_table` links each
`SEQ` to one `scientific_name`. The function sums abundances for
sequences assigned to the same name without recalculating RRA. Both
tables are required.

``` r
taxonomy_table <- trimmed_taxonomy[
  !is.na(trimmed_taxonomy$scientific_name) & nzchar(trimmed_taxonomy$scientific_name), ]
classified_table <- decontamination_table[
  decontamination_table$SEQ %in% taxonomy_table$SEQ, ]
merged_table <- METAeDNA::merge_replicate(
  sequence_table = classified_table,
  control_table = control_table,
  method = "mean_rra"
)
community_matrix <- METAeDNA::merge_to_matrix(
  merged_table = merged_table,
  control_table = control_table,
  taxonomy_table = taxonomy_table
)
dim(community_matrix)
community_matrix[seq_len(min(5, nrow(community_matrix))),
                 seq_len(min(5, ncol(community_matrix))), drop = FALSE]
```

The matrix contains relative abundances when the earlier
`merge_replicate()` call used `method = "mean_rra"`. Do not interpret
its row sums as read counts. Filtering classifications before merging
makes the denominator specific to the retained taxa. To use all
corrected reads as the denominator instead, replicates must be merged
before classification filtering. For DADA2, load
`taxonomy_results$dada_taxonomy` as the selected taxonomy table.

### Step 3. Inspect sequencing depth and rarefaction

Rarefaction uses the **reaction-level read counts** in
`decontamination_table`, before replicate merging. Taxonomy is matched
by `SEQ`. The plot shows rarefaction as solid lines, extrapolation as
dashed lines, and observed depths as points. An `experiment` column
automatically produces separate panels. Reused reaction IDs remain
separate across experiments. `facet_experiment = FALSE` pools
experiments instead.

``` r
depth_summary <- data.frame(
  reaction_id = sub("^count_", "", grep("^count_", names(decontamination_table), value = TRUE)),
  reads = colSums(decontamination_table[grep("^count_", names(decontamination_table))])
)
ggplot(depth_summary, aes(x = reorder(reaction_id, reads), y = reads)) +
  geom_col(fill = "#3B4992") +
  coord_flip() +
  labs(x = "Reaction", y = "Reads after decontamination") +
  theme_classic(base_size = 10)

rarefaction_result <- METAeDNA::rarefaction_extrapolation(
  sequence_table = decontamination_table,
  taxonomy_table = trimmed_taxonomy,
  n_threads = 4,
  show_observed = TRUE,
  show_extrapolation = TRUE
)
rarefaction_result$fig
```

### Step 4. Calculate sample diversity

`hill_number()` calculates richness (`q = 0`), effective Shannon
diversity (`q = 1`), and inverse Simpson diversity (`q = 2`). These are
all expressed as effective numbers of taxa, not raw Shannon entropy or
Simpson concentration.

``` r
diversity <- data.frame(
  sample_id = rownames(community_matrix),
  richness = METAeDNA::hill_number(community_matrix, q = 0),
  shannon_effective = METAeDNA::hill_number(community_matrix, q = 1),
  simpson_effective = METAeDNA::hill_number(community_matrix, q = 2)
)
diversity_long <- pivot_longer(diversity, -sample_id,
                                      names_to = "metric", values_to = "value")
ggplot(diversity_long, aes(x = reorder(sample_id, value), y = value)) +
  geom_col(fill = "#3B4992") +
  coord_flip() +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = "Sample", y = "Effective number of taxa") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank())
```

The rarefaction curves describe classified taxa; sequences without a
`scientific_name` are excluded from those curves. The depth bar plot
counts all retained reads. LCA labels above species level are sequence
groups, so richness computed from them should be reported as assigned
taxa rather than confirmed species richness.

### Step 5. Assess sampling-site coverage

`species_accumulation()` calculates expected taxon richness as the
number of sampled sites increases. Each row of `community_matrix`
represents one sample site after replicate merging. Repeated visits to
the same site must first be combined if the sampling unit is the site
rather than the visit.

The curve uses presence/absence and averages all subsets of the observed
sites analytically. It ends at the observed number of sites; it does not
predict unsampled taxa. A plateau indicates diminishing gains within
these sites, not proof of complete inventory.

``` r
site_design <- control_table[control_table$category == "sample", ]
if (!"experiment" %in% names(site_design)) {
  site_design$experiment <- unique(decontamination_table$experiment)
}
site_design$matrix_id <- paste(site_design$group, site_design$sample_id, sep = "_")
if (length(unique(site_design$experiment)) > 1) {
  site_design$matrix_id <- paste(site_design$experiment, site_design$matrix_id, sep = "_")
}
site_design <- unique(site_design[c("matrix_id", "experiment")])
site_experiment <- site_design$experiment[match(rownames(community_matrix), site_design$matrix_id)]
accumulation_result <- METAeDNA::species_accumulation(
  community_matrix, experiment = site_experiment
)
accumulation_result$fig
tail(accumulation_result$accumulation_data)
```

`additional_richness` records the expected gain from each added site. An
empty site counts toward sampling effort. Experiment groups are
calculated and displayed separately; their sites are not pooled.

### Step 6. Save analysis outputs

``` r
saveRDS(merged_table, "data/merged_rra.rds")
write.csv(accumulation_result$accumulation_data, "data/site_accumulation.csv", row.names = FALSE)
ggsave("data/site_accumulation.pdf", accumulation_result$fig, width = 158, height = 100, units = "mm")
write.csv(community_matrix, "data/community_matrix.csv")
write.csv(depth_summary, "data/sequencing_depth.csv", row.names = FALSE)
write.csv(diversity, "data/sample_diversity.csv", row.names = FALSE)
ggsave("data/rarefaction.pdf", rarefaction_result$fig, width = 158, height = 100, units = "mm")
```

## Reference

Taberlet P, Bonin A, Zinger L, Coissac E. (2018). *Environmental DNA: For
Biodiversity Research and Monitoring*. Oxford University Press, Chapter 6,
"DNA amplification and multiplexing".
[DOI](https://doi.org/10.1093/oso/9780198767220.001.0001).

## License

[MIT](LICENSE.md).

