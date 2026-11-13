2. Reference databases and taxonomic assignment
================

Build taxonomy resources and reference sequences, then classify the
sequences saved in [part 1](01_Sequence_processing.md). All paths are
relative to the same working directory. Database construction can also
be run on its own with a species checklist, independently of the
sequencing workflow.

## Step 1. Load sequence tables

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

## Step 2. Build the NCBI taxonomy database

The NCBI taxonomy database is separate from the BLAST nt database. The
nt database stores nucleotide sequences for similarity searches; the
taxonomy database stores names, taxonomic nodes, obsolete-taxid
redirects, and optionally accession-to-taxid mappings.

`build_ncbi_taxonomy_database()` downloads NCBI taxonomy files and
imports names and lineages into SQLite. The database supplies TaxIDs and
lineage information for reference retrieval and BLAST annotation.

| Argument | NCBI file | SQLite content | Required for |
|----|----|----|----|
| taxonomy = TRUE | taxdump.tar.gz | nodes, names, name aliases, merged taxids | name and lineage queries, COL–NCBI crosswalk |
| accession = TRUE | nucl_gb.accession2taxid.gz | GenBank nucleotide accession to taxid | converting BLAST accessions to taxids |
| include_WGS = TRUE | nucl_wgs.accession2taxid.gz | WGS accession to taxid | BLAST results containing WGS accessions |

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

## Step 3. Add accession-to-TaxID mappings

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

## Step 4. Build the COL name database

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

## Step 5. Reference database construction

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

### Route 1: use the local NCBI nt database

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

### Route 2: search and download a marker-specific database

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

## Step 6. Test primers by in-silico PCR

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

### Build the trimmed BLAST database

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

## Step 7. Search reference databases

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

### Run a BLAST command independently

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

| Columns | Meaning |
|----|----|
| `qseqid`, `sseqid` | Query ID and reference sequence ID |
| `pident`, `length` | Percentage identity and alignment length |
| `mismatch`, `gapopen` | Mismatch and gap-opening counts |
| `qstart`, `qend`, `sstart`, `send` | Query and reference alignment coordinates |
| `evalue`, `bitscore` | Alignment significance and score |
| `qcovs`, `staxids` | Query coverage (%) and reference NCBI TaxIDs |

## Step 8. Retrieve taxonomy and assign LCA

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

### Input from a BLAST data frame

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

### Input from an annotated table to find_lca()

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

## Step 9. Classify with DADA2

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

[sequencing depth and community diversity](03_Community_analysis.md).
