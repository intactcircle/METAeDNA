# METAeDNA

METAeDNA is an R package for environmental DNA (eDNA) metabarcoding analysis,
from paired-end sequencing reads to taxon abundance matrices and diversity
estimates. It provides tools for read merging, primer and tag assignment,
denoising, blank correction, reference database construction, and taxonomic
classification.

Reaction assignment uses tag-based multiplexing, not multiple target-gene
markers. In the supplied library design, PCR reactions carry pairs of 8-bp
oligonucleotide tags at the 5' ends of the primers before amplicons are pooled
for library preparation (Taberlet et al., 2018). Primer and tag combinations
identify each reaction, including technical replicates and blanks. Tag length
is specified by the reaction table and is not restricted to 8 bp.


## Workflow overview

The workflow is divided into three step-by-step tutorials:

1. **Sequence processing:** merge forward and reverse reads, assign reactions,
   filter sequences, and correct blank counts while retaining reaction-level reads.
2. **Reference databases and taxonomy:** obtain reference sequences from local
   nt or NCBI online, optionally trim them by in-silico PCR, and classify queries
   with BLAST or DADA2.
3. **Sequencing depth and diversity:** select classified sequences, calculate
   reaction RRA and merge technical replicates, build sample-by-taxon matrices,
   plot rarefaction curves, and calculate Hill diversity.

![Workflow overview](vignettes/Workflow.svg)

[Overview](vignettes/00_Workflow_overview.Rmd) ·
[1. Sequence processing](vignettes/01_Sequence_processing.Rmd) ·
[2. Databases and taxonomy](vignettes/02_Reference_taxonomy.Rmd) ·
[3. Depth and diversity](vignettes/03_Community_analysis.Rmd)

## System Requirements

- R >= 4.0.0
- C++17 compiler
- zlib
- Optional external tools:
  - BLAST+ for `blast_cmd()`, `blastn_cmd()`, and local reference building.
  - VSEARCH for clustering and chimera-removal helpers.

## Installation

Install the supplied source archive locally; required dependencies are handled automatically:

```r
# install.packages("remotes") # Run once if needed.
remotes::install_local("METAeDNA_1.0.0.tar.gz", upgrade = "never")
library(METAeDNA)
```
Step-by-step tutorials are supplied in the source archive's `vignettes` folder.

## External software

[VSEARCH](https://github.com/torognes/vsearch/releases) is used for chimera
removal and optional sequence clustering.
[NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
is used for reference extraction, database construction, and similarity
searches. The release must match the operating system. The paths below identify
the folders containing the extracted executables:

```r
tool_paths <- c("software/blast/bin", "software/vsearch/bin")
Sys.setenv(PATH = paste(c(tool_paths, Sys.getenv("PATH")), collapse = .Platform$path.sep))
Sys.which(c("vsearch", "blastn", "blastdbcmd", "makeblastdb"))
```

For BLAST+, use its `bin` folder. This sets `PATH` for the current R session;
if the programs are already on `PATH`, only the availability check is needed.

### Notes

- Use absolute paths for large reference databases and external command-line tools when possible.
- NCBI downloads and accession databases can require substantial storage space.
- Online NCBI queries should include a contact email and use conservative parallel settings.
- Local BLAST and VSEARCH functions require the corresponding executables to be installed and available in `PATH`, or supplied through the relevant executable-path argument.

## Run the example: paired-end reads to corrected counts

The package includes forward and reverse reads, a reaction-definition table
containing tags and primers, and a control table linking reactions to samples,
technical replicates, and blanks. The following example uses these inputs and
writes intermediate files under `data/` in the current working directory.
VSEARCH must be available for the chimera-removal step.

```r
forward_fastq <- system.file("extdata", "trnL_40k_1.fq.gz", package = "METAeDNA")
reverse_fastq <- system.file("extdata", "trnL_40k_2.fq.gz", package = "METAeDNA")
reaction_table <- system.file("extdata", "trnL_table.txt", package = "METAeDNA")
control_table <- read.csv(system.file("extdata", "control_table.csv", package = "METAeDNA"))
dir.create("data", showWarnings = FALSE)

# Merge forward and reverse reads.
pairend(fastq_1 = forward_fastq, fastq_2 = reverse_fastq,
        out_file = "data/01_paired.fastq", min_overlap = 10L, n_threads = 4L)

# Assign reactions using primers and tags.
demultiplex(fastq_file = "data/01_paired.fastq", reaction_table = reaction_table,
            identified_output = "data/02_identified.fasta",
            unidentified_output = "data/02_unidentified.fasta",
            primer_mismatch = 2L, n_threads = 4L)

# Collapse identical sequences and denoise.
collapse_sequences(fasta_file = "data/02_identified.fasta",
                   output_fasta = "data/03_unique.fasta", n_threads = 4L)
denoise(fasta_file = "data/03_unique.fasta", output_fasta = "data/04_denoised.fasta",
        min_count = 10L, keep_internal = FALSE, n_threads = 4L)

# Remove chimeras and read reaction counts.
vsearch_remove_chimeras(fasta_file = "data/04_denoised.fasta",
                       output_fasta = "data/05_nonchimera.fasta")
sequence_table <- fasta_to_table(fasta_file = "data/05_nonchimera.fasta", n_threads = 4L)

# Correct blanks while retaining reaction-level counts.
decontamination_table <- decontamination(
  sequence_table, control_table, control_names = c("EB", "FB", "PB")
)

saveRDS(
  list(sequence_table = sequence_table, decontamination_table = decontamination_table,
       control_table = control_table),
  "data/workflow_tables.rds"
)
```

`decontamination_table` retains corrected integer counts for individual
reactions. After classification, retained sequences are selected by `SEQ`.
`merge_replicate()` then calculates each reaction's RRA using its total
retained classified reads as the denominator and averages technical replicates.
Reaction counts are used for sequencing depth and rarefaction; sample RRA is
used for the community matrix and diversity.

Continue with [reference databases and classification](vignettes/02_Reference_taxonomy.Rmd),
then [sequencing depth and diversity](vignettes/03_Community_analysis.Rmd).
Processing thresholds depend on the marker and library design.


## Reference

Taberlet P, Bonin A, Zinger L, Coissac E. (2018). *Environmental DNA: For
Biodiversity Research and Monitoring*. Oxford University Press, Chapter 6,
"DNA amplification and multiplexing".
[DOI](https://doi.org/10.1093/oso/9780198767220.001.0001).

## License

[MIT](LICENSE.md).

