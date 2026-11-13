3. Sequencing depth and community diversity
================

Use reaction counts for sequencing depth and rarefaction, and sample RRA
values for the community matrix and diversity. This tutorial reads the
outputs saved in [part 1](01_Sequence_processing.md) and [part
2](02_Reference_taxonomy.md).

## Step 1. Load reaction counts and classifications

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

## Step 2. Calculate RRA and construct the community matrix

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

## Step 3. Inspect sequencing depth and rarefaction

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

## Step 4. Calculate sample diversity

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

## Step 5. Assess sampling-site coverage

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

## Step 6. Save analysis outputs

``` r
saveRDS(merged_table, "data/merged_rra.rds")
write.csv(accumulation_result$accumulation_data, "data/site_accumulation.csv", row.names = FALSE)
ggsave("data/site_accumulation.pdf", accumulation_result$fig, width = 158, height = 100, units = "mm")
write.csv(community_matrix, "data/community_matrix.csv")
write.csv(depth_summary, "data/sequencing_depth.csv", row.names = FALSE)
write.csv(diversity, "data/sample_diversity.csv", row.names = FALSE)
ggsave("data/rarefaction.pdf", rarefaction_result$fig, width = 158, height = 100, units = "mm")
```
