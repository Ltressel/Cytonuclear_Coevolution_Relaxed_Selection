#!/usr/bin/env Rscript

########### Coordinated Evolutionary Rates in Oxidative Phosphorylation Complexes of Papilionoid Legumes
########### Main ERC analysis script
########### This script:
########### 1. loads trees and extracts root-to-tip + terminal branch lengths
########### 2. normalizes rates to the random nuclear tree
########### 3. runs ERC analyses from branch-length Excel files
########### 4. saves summary tables, bootstrap results, plot-ready data, and an .RData workspace

suppressPackageStartupMessages({
  library(ape)
  library(adephylo)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(boot)
  library(writexl)
})

options(stringsAsFactors = FALSE)

set.seed(123)

# ----------------------------- #
# User settings
# ----------------------------- #

tree_dir <- "."
analysis_dir <- "."
out_dir <- file.path(analysis_dir, "erc_output")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----------------------------- #
# Helper functions
# ----------------------------- #

get_terminal_brlen <- function(tr) {
  stopifnot(!is.null(tr$edge.length))
  term_idx <- which(tr$edge[, 2] <= length(tr$tip.label))
  tip_ids  <- tr$edge[term_idx, 2]
  vals     <- tr$edge.length[term_idx]
  names(vals) <- tr$tip.label[tip_ids]
  vals
}

named_vector_to_long <- function(x, code) {
  tibble(
    taxon = names(x),
    distance = as.numeric(x),
    code = code
  )
}

normalize_against_random <- function(df, random_code) {
  rand_df <- df %>%
    filter(code == random_code) %>%
    select(taxon, random_distance = distance)

  df %>%
    left_join(rand_df, by = "taxon") %>%
    mutate(normalized_distance = distance / random_distance)
}

safe_boot_ci <- function(boot_obj) {
  out <- tryCatch(
    boot.ci(boot_obj, type = "bca"),
    error = function(e) NULL
  )

  if (is.null(out) || is.null(out$bca)) {
    return(c(ci_low = NA_real_, ci_high = NA_real_))
  }

  c(ci_low = out$bca[4], ci_high = out$bca[5])
}

run_pairwise_erc <- function(df,
                             subset_col = NULL,
                             subset_value = NULL,
                             type_col,
                             x_group,
                             y_group,
                             value_col,
                             group_col,
                             analysis_name,
                             n_boot = 10000) {

  dat <- df

  if (!is.null(subset_col) && !is.null(subset_value)) {
    dat <- dat %>% filter(.data[[subset_col]] == subset_value)
  }

  xdat <- dat %>%
    filter(.data[[type_col]] == x_group) %>%
    arrange(Species)

  ydat <- dat %>%
    filter(.data[[type_col]] == y_group) %>%
    arrange(Species)

  common_species <- intersect(xdat$Species, ydat$Species)

  xdat <- xdat %>% filter(Species %in% common_species) %>% arrange(Species)
  ydat <- ydat %>% filter(Species %in% common_species) %>% arrange(Species)

  if (nrow(xdat) == 0 || nrow(ydat) == 0) {
    stop(paste("No overlapping species found for", analysis_name))
  }

  plot_df <- tibble(
    Species = xdat$Species,
    x_value = xdat[[value_col]],
    y_value = ydat[[value_col]],
    Group = xdat[[group_col]]
  )

  cor_test <- cor.test(plot_df$x_value, plot_df$y_value, method = "spearman", exact = FALSE)

  boot_fun <- function(data, indices) {
    d <- data[indices, , drop = FALSE]
    cor(d$x_value, d$y_value, method = "spearman")
  }

  boot_res <- boot(data = plot_df[, c("x_value", "y_value")], statistic = boot_fun, R = n_boot)
  ci_vals <- safe_boot_ci(boot_res)

  summary_tbl <- tibble(
    analysis_name = analysis_name,
    n_species = nrow(plot_df),
    spearman_rs = unname(cor_test$estimate),
    p_value = cor_test$p.value,
    mean_bootstrap_rs = mean(boot_res$t, na.rm = TRUE),
    ci_low = unname(ci_vals["ci_low"]),
    ci_high = unname(ci_vals["ci_high"])
  )

  list(
    summary = summary_tbl,
    plot_data = plot_df,
    bootstrap_values = tibble(
      analysis_name = analysis_name,
      bootstrap_rs = as.numeric(boot_res$t)
    )
  )
}

run_analysis_set <- function(file,
                             subset_col = NULL,
                             subset_values = NULL,
                             type_col,
                             pairings,
                             value_col,
                             group_col,
                             n_boot = 10000) {

  dat <- read_excel(file)
  summaries <- list()
  plot_data <- list()
  boot_data <- list()

  for (i in seq_len(nrow(pairings))) {
    subset_value <- if (!is.null(subset_values)) subset_values[i] else NULL

    res <- run_pairwise_erc(
      df = dat,
      subset_col = subset_col,
      subset_value = subset_value,
      type_col = type_col,
      x_group = pairings$x_group[i],
      y_group = pairings$y_group[i],
      value_col = value_col,
      group_col = group_col,
      analysis_name = pairings$analysis_name[i],
      n_boot = n_boot
    )

    summaries[[length(summaries) + 1]] <- res$summary
    plot_data[[length(plot_data) + 1]] <- res$plot_data %>%
      mutate(analysis_name = pairings$analysis_name[i])
    boot_data[[length(boot_data) + 1]] <- res$bootstrap_values
  }

  list(
    summary = bind_rows(summaries),
    plot_data = bind_rows(plot_data),
    bootstrap = bind_rows(boot_data)
  )
}

# ----------------------------- #
# Section 1: tree-based distances
# ----------------------------- #

tree_files <- c(
  all_mt_tree = "All_mt_genes_tree",
  all_n_mt_tree = "All_n_mt_tree",
  cc_tree = "cc_tree",
  gly_tree = "gly_tree",
  cr_tree = "cr_tree",
  cprp_tree = "cprp_tree",
  ci_mt_tree = "ci_mt_tree",
  ci_n_mt_tree = "ci_n_mt_tree",
  cii_mt_tree = "cii_mt_tree",
  cii_n_mt_tree = "cii_n_mt_tree",
  ciii_mt_tree = "ciii_mt_tree",
  ciii_n_mt_tree = "ciii_n_mt_tree",
  civ_mt_tree = "civ_mt_tree",
  civ_n_mt_tree = "civ_n_mt_tree",
  cv_mt_tree = "cv_mt_tree",
  cv_n_mt_tree = "cv_n_mt_tree",
  normalizing_tree = "nrand"
)

trees <- lapply(file.path(tree_dir, tree_files), read.tree)

root_to_tip_codes <- c(
  gly_tree = "gly.dist",
  all_mt_tree = "mt.dist",
  all_n_mt_tree = "nuc.dist",
  cc_tree = "cc.dist",
  cr_tree = "cr.dist",
  cprp_tree = "cprp.dist",
  ci_mt_tree = "ci_mt.dist",
  ci_n_mt_tree = "ci_n_mt.dist",
  cii_mt_tree = "cii_mt.dist",
  cii_n_mt_tree = "cii_n_mt.dist",
  ciii_mt_tree = "ciii_mt.dist",
  ciii_n_mt_tree = "ciii_n_mt.dist",
  civ_mt_tree = "civ_mt.dist",
  civ_n_mt_tree = "civ_n_mt.dist",
  cv_mt_tree = "cv_mt.dist",
  cv_n_mt_tree = "cv_n_mt.dist",
  normalizing_tree = "nrand.dist"
)

terminal_codes <- c(
  gly_tree = "gly.term",
  all_mt_tree = "mt.term",
  all_n_mt_tree = "nmt.term",
  cc_tree = "cc.term",
  cr_tree = "cr.term",
  cprp_tree = "cprp.term",
  ci_mt_tree = "ci_mt.term",
  ci_n_mt_tree = "ci_n_mt.term",
  cii_mt_tree = "cii_mt.term",
  cii_n_mt_tree = "cii_n_mt.term",
  ciii_mt_tree = "ciii_mt.term",
  ciii_n_mt_tree = "ciii_n_mt.term",
  civ_mt_tree = "civ_mt.term",
  civ_n_mt_tree = "civ_n_mt.term",
  cv_mt_tree = "cv_mt.term",
  cv_n_mt_tree = "cv_n_mt.term",
  normalizing_tree = "norm.term"
)

dist.dat <- bind_rows(lapply(names(root_to_tip_codes), function(nm) {
  named_vector_to_long(distRoot(trees[[nm]]), root_to_tip_codes[[nm]])
}))

term.dat <- bind_rows(lapply(names(terminal_codes), function(nm) {
  named_vector_to_long(get_terminal_brlen(trees[[nm]]), terminal_codes[[nm]])
}))

normalized_distances <- normalize_against_random(dist.dat, "nrand.dist")
normalized_terminal  <- normalize_against_random(term.dat, "norm.term")

write_xlsx(
  normalized_distances %>%
    select(taxon, code, distance, random_distance, normalized_distance),
  file.path(out_dir, "normalized_root_to_tip_distances.xlsx")
)

write_xlsx(
  normalized_terminal %>%
    select(taxon, code, distance, random_distance, normalized_distance),
  file.path(out_dir, "normalized_terminal_branch_lengths.xlsx")
)

# ----------------------------- #
# Section 2: ERC analyses
# ----------------------------- #

group_col <- "Other Legumes (OL) or 50-kb (50kb)"

# OXPHOS complexes: root-to-tip
complex_pairings <- tibble(
  analysis_name = c("Complex I", "Complex II", "Complex III", "Complex IV", "Complex V"),
  x_group = c("N-mt", "N-mt", "N-mt", "N-mt", "N-mt"),
  y_group = c("mt", "mt", "mt", "mt", "mt")
)

complexes_root <- run_analysis_set(
  file = file.path(analysis_dir, "mt_N_mt_OXPHOS_complexes_branch_lengths_FINAL.xlsx"),
  subset_col = "Complex",
  subset_values = c("CI", "CII", "CIII", "CIV", "CV"),
  type_col = "Type of gene group (mt or N-mt)",
  pairings = complex_pairings,
  value_col = "Normalized Root to Tip Branch Length",
  group_col = group_col
)

complexes_term <- run_analysis_set(
  file = file.path(analysis_dir, "mt_N_mt_OXPHOS_complex_branch_lengths_terminal_FINAL.xlsx"),
  subset_col = "Complex",
  subset_values = c("CI", "CII", "CIII", "CIV", "CV"),
  type_col = "Type of gene group (mt or N-mt)",
  pairings = complex_pairings %>% mutate(analysis_name = paste0(analysis_name, " (Terminal)")),
  value_col = "Normalized Terminal Branch Length",
  group_col = group_col
)

# Gene groups: root-to-tip
gene_pairings_root <- tribble(
  ~analysis_name,       ~file,                                       ~type_col,                          ~x_group, ~y_group, ~value_col,
  "N-mt OXPHOS",        "mt_N_mt_OXPHOS_branch_lengths.xlsx",        "Type of gene group (mt or N-mt)", "N-mt",   "mt",     "Normalized Root to Tip Branch Length",
  "Cell Cycle",         "mt_Cell_Cycle_branch_lengths.xlsx",         "Type of gene group (mt or CC)",   "CC",     "mt",     "Normalized Root to Tip Branch Length",
  "Glycolysis",         "mt_Glycolysis_branch_lengths.xlsx",         "Type of gene group (mt or GLY)",  "GLY",    "mt",     "Normalized Root to Tip Branch Length",
  "Cyto-ribo",          "mt_Cytosolic_ribosomal_branch_lengths.xlsx","Type of gene group (mt or CR)",   "CR",     "mt",     "Normalized Root to Tip Branch Length",
  "Plastid",            "mt_plastid_branch_lengths.xlsx",            "Type of gene group (mt or CpRP)", "CpRP",   "mt",     "Normalized Root to Tip Branch Length"
)

gene_pairings_term <- tribble(
  ~analysis_name,                 ~file,                                                ~type_col,                          ~x_group, ~y_group, ~value_col,
  "N-mt OXPHOS (Terminal)",       "mt_N_mt_OXPHOS_terminal_branch_lengths.xlsx",         "Type of gene group (mt or N-mt)", "N-mt",   "mt",     "Normalized Terminal Branch Length",
  "Cell Cycle (Terminal)",        "mt_Cell_Cycle_terminal_branch_lengths.xlsx",          "Type of gene group (mt or CC)",   "CC",     "mt",     "Normalized Terminal Branch Length",
  "Glycolysis (Terminal)",        "mt_Glycolysis_terminal_branch_lengths.xlsx",          "Type of gene group (mt or GLY)",  "GLY",    "mt",     "Normalized Terminal Branch Length",
  "Cyto-ribo (Terminal)",         "mt_Cytosolic_ribosomal_terminal_branch_lengths.xlsx", "Type of gene group (mt or CR)",   "CR",     "mt",     "Normalized Terminal Branch Length",
  "Plastid (Terminal)",           "mt_plastid_terminal_branch_lengths.xlsx",             "Type of gene group (mt or CpRP)", "CpRP",   "mt",     "Normalized Terminal Branch Length"
)

run_gene_row <- function(row_df) {
  run_pairwise_erc(
    df = read_excel(file.path(analysis_dir, row_df$file)),
    type_col = row_df$type_col,
    x_group = row_df$x_group,
    y_group = row_df$y_group,
    value_col = row_df$value_col,
    group_col = group_col,
    analysis_name = row_df$analysis_name,
    n_boot = 10000
  )
}

gene_root_results <- lapply(seq_len(nrow(gene_pairings_root)), function(i) run_gene_row(gene_pairings_root[i, ]))
gene_term_results <- lapply(seq_len(nrow(gene_pairings_term)), function(i) run_gene_row(gene_pairings_term[i, ]))

gene_root <- list(
  summary = bind_rows(lapply(gene_root_results, `[[`, "summary")),
  plot_data = bind_rows(lapply(gene_root_results, `[[`, "plot_data")),
  bootstrap = bind_rows(lapply(gene_root_results, `[[`, "bootstrap_values"))
)

gene_term <- list(
  summary = bind_rows(lapply(gene_term_results, `[[`, "summary")),
  plot_data = bind_rows(lapply(gene_term_results, `[[`, "plot_data")),
  bootstrap = bind_rows(lapply(gene_term_results, `[[`, "bootstrap_values"))
)

# ----------------------------- #
# Section 3: save outputs
# ----------------------------- #

summary_tables <- list(
  complexes_root_summary = complexes_root$summary,
  complexes_terminal_summary = complexes_term$summary,
  gene_groups_root_summary = gene_root$summary,
  gene_groups_terminal_summary = gene_term$summary
)

write_xlsx(summary_tables, file.path(out_dir, "ERC_summary_tables.xlsx"))

write.csv(complexes_root$plot_data, file.path(out_dir, "complexes_root_plot_data.csv"), row.names = FALSE)
write.csv(complexes_term$plot_data, file.path(out_dir, "complexes_terminal_plot_data.csv"), row.names = FALSE)
write.csv(gene_root$plot_data, file.path(out_dir, "gene_groups_root_plot_data.csv"), row.names = FALSE)
write.csv(gene_term$plot_data, file.path(out_dir, "gene_groups_terminal_plot_data.csv"), row.names = FALSE)

write.csv(complexes_root$bootstrap, file.path(out_dir, "complexes_root_bootstrap.csv"), row.names = FALSE)
write.csv(complexes_term$bootstrap, file.path(out_dir, "complexes_terminal_bootstrap.csv"), row.names = FALSE)
write.csv(gene_root$bootstrap, file.path(out_dir, "gene_groups_root_bootstrap.csv"), row.names = FALSE)
write.csv(gene_term$bootstrap, file.path(out_dir, "gene_groups_terminal_bootstrap.csv"), row.names = FALSE)

save(
  dist.dat,
  term.dat,
  normalized_distances,
  normalized_terminal,
  complexes_root,
  complexes_term,
  gene_root,
  gene_term,
  file = file.path(out_dir, "erc_results.RData")
)

message("ERC analysis complete.")
message("Outputs written to: ", normalizePath(out_dir))
