#' Decontaminate sequence table using blank samples
#'
#' This function removes contaminant sequences by comparing sample counts
#' against corresponding blank samples defined in a control table.
#' The procedure supports several removal strategies based on blank abundance.
#'
#' @param seq_table A data frame generated from FASTA, containing columns
#'   `experiment`, `ID`, `SEQ`, and multiple `count_*` columns.
#'
#' @param control_table A data frame describing sample–blank relationships.
#'   Must include columns `group`, `sample`, `tag`, and blank columns named
#'   in `control_names`.
#'
#' @param method The decontamination method. One of:
#'   \itemize{
#'     \item \code{"maximum"} — subtract the maximum blank count
#'     \item \code{"threshold"} — remove sequences when blank max < threshold
#'     \item \code{"mean"} — subtract the mean blank count
#'     \item \code{"both_threshold"} — remove only when all blank counts < threshold
#'   }
#'
#' @param control_names A character vector giving blank types (e.g. `c("FB","EB","PB")`).
#'
#' @param threshold Numeric threshold used in methods "threshold" and "both_threshold".
#'
#' @return The input \code{seq_table} with blank-corrected abundance columns.
#'
#' @details
#' For each blank type, the function identifies all associated blank samples
#' and maps them to target sample tags. Blank abundance is computed per sequence
#' and used to adjust sample counts. All adjustments are performed only on
#' sample columns corresponding to the affected blank group.
#'
#' @export
#'
#' @examples 1
decontamination <- function(seq_table, control_table, method = c("maximum", "threshold","mean","both_threshold"), control_names = NULL, threshold = 10, remove_status = TRUE){
  method <- match.arg(method)
  # control_names <- c("FB","EB","PB")
  control_table
  blank_map <- control_table %>%
    pivot_longer(cols = all_of(control_names),
                 names_to = "blank_type",
                 values_to = "blank_sample") %>%
    select(group,sample, tag, blank_type, blank_sample,category)

  for(bl in control_names){
    blank_tag_table <- blank_map %>% filter(sample %in% .$blank_sample, category  == bl) %>% distinct(group, sample, tag) %>% dplyr::rename(rm_tag = tag)
    de_blank_table_map <- blank_map %>% filter(blank_type  == bl, category == "sample")
    for(rm in blank_tag_table$sample %>% unique){
      temp_blank_tag_table <- blank_tag_table %>% filter(sample == rm)
      temp_de_table <- seq_table %>% select(experiment, ID, SEQ, any_of(paste0("count_", temp_blank_tag_table$rm_tag)))
      if(temp_de_table %>% select(contains("count")) %>% ncol){
        sample_tag_to_treat <- de_blank_table_map %>% filter(blank_sample == rm)
        if(method == "maximum"){
          remove_count <- temp_de_table %>% mutate(remove_c = temp_de_table %>% select(contains("count")) %>% apply(1,max))
          seq_table <- left_join(seq_table, remove_count %>% select(experiment , ID, SEQ, remove_c), by = join_by(experiment , ID, SEQ))
          seq_table <- seq_table %>% mutate(across(any_of(paste0("count_",sample_tag_to_treat$tag)), ~ pmax(.x - remove_c, 0)))

        }else if (method == "threshold"){

          remove_count <- temp_de_table %>% mutate(remove_c = as.numeric((temp_de_table %>% select(contains("count")) %>% apply(1,max) ) < threshold))
          seq_table <- left_join(seq_table, remove_count %>% select(experiment , ID, SEQ, remove_c), by = join_by(experiment , ID, SEQ))
          seq_table <- seq_table %>% mutate(across(any_of(paste0("count_",sample_tag_to_treat$tag)), ~ .x *remove_c))
        }else if (method == "both_threshold"){

          remove_count <- temp_de_table %>% mutate(remove_c = (temp_de_table %>% select(contains("count")) %>% apply(1,function(x, threshold){all(x<threshold) %>% as.numeric}, threshold = threshold) ))
          seq_table <- left_join(seq_table, remove_count %>% select(experiment , ID, SEQ, remove_c), by = join_by(experiment , ID, SEQ))
          # seq_table[6, paste0("count_", sample_tag_to_treat$tag[1])]
          seq_table <- seq_table %>% mutate(across(any_of(paste0("count_",sample_tag_to_treat$tag)), ~ .x *remove_c))

        }else if (method == "mean"){

          remove_count <- temp_de_table %>% mutate(remove_c = temp_de_table %>% select(contains("count")) %>% apply(1,function(x,threshold){mean(x) %>% round()}, threshold = threshold))
          seq_table <- left_join(seq_table, remove_count %>% select(experiment , ID, SEQ, remove_c), by = join_by(experiment , ID, SEQ))
          seq_table <- seq_table %>% mutate(across(any_of(paste0("count_",sample_tag_to_treat$tag)), ~ pmax(.x - remove_c, 0)))

        }
        seq_table <- seq_table %>% select(-remove_c, -any_of(paste0("count_", temp_blank_tag_table$rm_tag)))
      }
    }
  }
  if(remove_status){
    seq_table <- seq_table %>% select(-any_of(matches("status_")))
  }

  seq_table %>% return()
}
# test <- decontamination(seq_table, control_table, control_names = control_names)


#
#' Merge technical replicates using Relative Read Abundance (RRA)
#'
#' @description
#' This function merges technical replicates of PCR belonging to the same
#' biological sample, based on replicate definitions provided in `control_table`.
#' For each group–sample combination, replicate count columns (e.g., `count_A1`,
#' `count_A2`, `count_A3`) are first normalized to Relative Read Abundance (RRA),
#' and then merged using one of three strategies:
#'
#' \itemize{
#'   \item \code{"mean_rra"}: average RRA across replicates;
#'   \item \code{"max_rra"}: maximum RRA across replicates;
#'   \item \code{"both_occur"}: sequence must appear in at least
#'         \code{repeat_times} replicates; otherwise the merged abundance is set
#'         to zero.
#' }
#'
#' After merging, the original replicate count columns are removed and replaced
#' by a new column named \code{"group_sample"} (e.g., \code{"G1_S1"}).
#'
#'
#' @param seq_table A data frame containing per-sequence counts.
#'   Must include columns of the form \code{count_TAG}, where TAG corresponds
#'   to replicate identifiers in \code{control_table}.
#'
#' @param control_table A data frame defining the relationship between
#'   biological samples and technical replicates. Must contain at least:
#'   \itemize{
#'     \item \code{group}: biological group identifier;
#'     \item \code{sample}: biological sample name;
#'     \item \code{tag}: replicate tag (matching \code{count_TAG} in seq_table);
#'     \item \code{category}: sample type, where rows with
#'           \code{category == "sample"} define valid sequencing samples.
#'   }
#'
#' @param method Character string specifying the replicate-merging strategy.
#'   One of:
#'   \itemize{
#'     \item \code{"mean_rra"}
#'     \item \code{"max_rra"}
#'     \item \code{"both_occur"}
#'   }
#'
#' @param repeat_times Integer. For \code{"both_occur"} mode, a sequence must
#'   appear (abundance > 0) in at least this number of replicates to be retained.
#'
#' @param present_absent Logical. If \code{TRUE}, the merged replicate columns
#'   are converted into presence/absence (1/0).
#'
#' @param remove_all_zero Logical. If \code{TRUE} (default), sequences whose
#'   merged abundance is zero in all samples are removed.
#'
#'
#' @details
#' For each biological sample, the function:
#' \enumerate{
#'   \item Identifies replicate tags from \code{control_table};
#'   \item Extracts corresponding count columns from \code{seq_table};
#'   \item Normalizes each replicate column to RRA by dividing by its column sum;
#'   \item Merges replicates using the method specified by \code{method};
#'   \item Drops the original \code{count_} columns;
#'   \item Adds a merged abundance column named \code{group_sample}.
#' }
#'
#' This function operates at the sequence level: each row corresponds to an
#' individual sequence (OTU/ASV) and replicates are merged for that row.
#'
#'
#' @return
#' A modified version of \code{seq_table} containing one merged abundance column
#' per biological sample (e.g., \code{"G1_S1"}), with optional presence/absence
#' and zero-row removal performed.
#'
#'
#' @examples
#' \dontrun{
#' merged_seq <- merge_replicate(
#'     seq_table       = seq_tab,
#'     control_table   = ctrl_tab,
#'     method          = "both_occur",
#'     repeat_times    = 2,
#'     present_absent  = FALSE,
#'     remove_all_zero = TRUE
#' )
#' }
#'
#' @export
merge_replicate <- function(seq_table, control_table, method = c("mean_rra", "max_rra", "both_occur"), repeat_times = 1, present_absent = FALSE, remove_all_zero = T){
  method <- match.arg(method)
  rm_col <- (seq_table %>% select(all_of(contains("count_"))) %>% colSums()) == 0
  seq_table <- seq_table %>% select(-any_of(names(rm_col)[rm_col]))
  sample_table <- control_table %>% filter(category == "sample")
  sample_name <- c(sample_table$sample) %>% unique
  group_name <- sample_table$group %>% unique
  for(gn in group_name){
    for(sn in sample_name){
      merge_tags <- sample_table %>% filter(group == gn, sample == sn )
      if(nrow(merge_tags) == 0 ){
        seq_table <- seq_table %>% mutate(!!paste(gn,sn,sep = '_') := 0)
        next
        }
      temp_res <- seq_table %>% select(any_of(paste0("count_", merge_tags$tag)))
      if( temp_res %>% ncol){
        rra_data <-  temp_res %>% mutate(across(any_of(paste0("count_", merge_tags$tag)), ~ .x/sum(.x))) %>%
          {
            if(method == "mean_rra"){
               apply(.,1, mean)
            }else if (method == "max_rra"){
               apply(.,1, max)
            }else if (method == "both_occur"){
              temp_rra <- apply(.,1, mean)
              temp_tf_loc <- apply(.,1, function(x, repeat_times){as.numeric(length (which(x > 0)) >= repeat_times)} ,repeat_times = repeat_times) %>% as.logical()
              temp_rra[!temp_tf_loc] <- 0
              temp_rra
            }
          }

        seq_table <- seq_table %>% mutate(!!paste(gn,sn,sep = '_') :=  rra_data) %>%
          select(-any_of(paste0("count_", merge_tags$tag)))

      }else{
        seq_table <- seq_table %>% mutate(!!paste(gn,sn,sep = '_') := 0)
        next
      }

    }
  }
  if(remove_all_zero){
    seq_table <- seq_table[((seq_table %>% select(any_of(paste0(sample_table$group, "_", sample_table$sample))) %>% apply(1, sum)) >0),]
  }
  if(present_absent){
    seq_table <- seq_table %>% mutate(across(any_of(paste0(sample_table$group, "_", sample_table$sample)), ~ ifelse(.x >0, 1,0)))
  }
  return(seq_table)
}


#' Parallel computation and visualization of iNEXT rarefaction/extrapolation curves
#'
#' This function performs parallelized iNEXT computation for multiple samples
#' (columns) of an abundance matrix. For each sample, standardized diversity
#' estimates (qD) are computed based on rarefaction, extrapolation, and
#' coverage-based methods. The function returns both the combined iNEXT
#' estimates and a ggplot object for visualization.
#'
#' @param temp_data A matrix or data frame of species-by-sample abundance
#'   data. Each column represents a sample; columns with zero total abundance
#'   are removed automatically.
#' @param threads Integer. Number of parallel workers (future backend).
#' @param q Diversity order (Hill number). Default is 0.
#' @param endpoint Maximum sample size used for extrapolation in iNEXT.
#' @param nboot Number of bootstrap replicates used in iNEXT estimation.
#' @param x_limit Optional numeric vector of x-axis limits.
#' @param y_limit Optional numeric vector of y-axis limits.
#' @param show_legend Either \code{NULL} (no legend) or a string specifying
#'   legend position (e.g., \code{"right"}, \code{"bottom"}).
#' @param show_observed Logical; whether to plot observed richness points.
#' @param show_extrapolation Logical; whether to add extrapolation curves.
#' @param x_title X-axis label.
#' @param y_title Y-axis label.
#' @param based One of \code{"size"} or \code{"coverage"} specifying which
#'   estimator to visualize.
#'
#' @details
#' The function uses \pkg{future.apply} for parallel computation, where each
#' sample is processed independently via \code{iNEXT()}. Results are combined,
#' reshaped, and plotted using \pkg{ggplot2}.
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{inext_data}: A tidy data frame of iNEXT estimates.
#'   \item \code{fig}: A ggplot object of rarefaction/extrapolation curves.
#' }
#'
#' @importFrom iNEXT iNEXT
#'
#' @importFrom dplyr %>% select filter mutate across bind_rows bind_cols everything
#'
#' @importFrom ggplot2 ggplot geom_line geom_point scale_x_continuous
#'   scale_y_continuous labs theme element_text element_blank element_rect
#'
#' @importFrom grid unit
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_apply
#'
#' @importFrom purrr map_dfr
#'
#' @examples
#' \dontrun{
#' data <- matrix(rpois(1000, lambda = 3), ncol = 5)
#' # each colum represents a PCR sample
#' res <- iNEXT_parallel(data, threads = 4, show_observed = TRUE)
#' res$fig
#' }
#'
#' @export
rarefaction_extrapolation <- function(temp_data, threads = 1, q = 0, endpoint = 50000,nboot = 20,x_limit = NULL, y_limit = NULL, show_legend = NULL,show_observed = F, show_extrapolation =F, x_title = "Reads",y_title = "Species richness ",
                      based = c("size", "coverage")){

  plan(multisession, workers = threads)
  on.exit(plan(sequential))
  temp_data <- temp_data[(temp_data %>% colSums())>0]
  inext_data <- future_apply(temp_data, MARGIN = 2, function(x) {

    x <- x[x > 0]
    inext_data <- iNEXT(
      x,
      datatype = "abundance",
      q = q,
      endpoint = 150000,
      nboot = 20
    )
    bind_rows(
      bind_cols(inext_data$iNextEst$size_based, based = "size"),
      bind_cols(inext_data$iNextEst$coverage_based, based = "coverage")
    ) %>%
      return()
  }, future.seed = TRUE)
  plan(sequential)

  final_inext_data <- inext_data %>% map_dfr(~ ., .id = "tag")

  plt <- ggplot(data = final_inext_data %>% filter(based == based, Method == "Rarefaction"), aes(y = qD, x = m, colour = tag))+
    geom_line(alpha = 0.6, linewidth = 0.8)+
    {
      if(show_observed){
        geom_point(data = final_inext_data %>% filter(based == based, Method == "Observed"),alpha = 0.6)
      }
    }+
    {
      if(show_extrapolation){
        geom_line(data = final_inext_data %>% filter(based == based, Method == "Extrapolation"), linetype = "dashed",alpha = 0.6, linewidth = 0.8)
      }
    }+
    {
      if(is.null(x_limit)){
        scale_x_continuous(limits = x_limit)
      }
    }+
    {
    if(is.null(y_limit)){
      scale_y_continuous(limits = y_limit)
    }
    }+
    labs(x = x_title, y = y_title)+

    #
    theme(
      plot.title = element_text(hjust = 0.5, size = 16,face = "bold"),
      text = element_text(size = 16, face = "bold", family = "Source Han Serif SC"),
      legend.position = ifelse(is.null(show_legend), "none", show_legend),

      axis.text.x = element_text(angle =45,hjust = 1,vjust = 1),
      strip.text = element_blank(),
      strip.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      panel.spacing.y = unit(1, "lines")
    )
    return(list(inext_data = final_inext_data, fig = plt))
}

