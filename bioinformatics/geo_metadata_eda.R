#!/usr/bin/env Rscript
#
# geo_metadata_eda.R
#
# Description:
#   This script downloads a GEO Series dataset by GSE ID, extracts and tidies
#   the sample metadata (phenotype data), and performs quick exploratory
#   data analysis (EDA). It automatically parses "characteristics_ch1"
#   fields into structured columns, repairs duplicate column names, and
#   generates summary text reports and simple plots (age histograms, sex barplots).
#
#
#
# Outputs:
#   - CSV file of parsed metadata per platform
#   - Summary text files describing sample composition
#   - Quick plots (age histograms, sex barplots) when applicable
#   - Top-level overview of platforms and sample counts
#
# Requirements:
#     GEOquery, BiocManager, dplyr, tidyr, tibble, stringr,
#     ggplot2, readr
#
# Author:
#   Yaro Markov (2025)


# Define this helper function
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
  library(tibble); library(ggplot2)
})

# Load packages
library(GEOquery)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readr)
library(stringr)

ensure_pkg <- function(pkgs) {
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install)) {
    if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
  }
}

ensure_pkg(c("GEOquery"))

# In your R notebook
gse_id  <- "GSE25837"         # replace with the GEO Series you want
outdir  <- "~/Desktop/GSE25837"          # output folder

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)



message("Downloading GEO series: ", gse_id)
# Try matrix first (fastest when available)
gse_obj <- tryCatch(
  GEOquery::getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE),
  error = function(e) { message("GSEMatrix failed, retrying with GSEMatrix=FALSE"); NULL }
)

if (is.null(gse_obj)) {
  gse_obj <- GEOquery::getGEO(gse_id, GSEMatrix = FALSE)
  # Build ExpressionSets per platform if needed
  # (Often not necessary—most Series have prebuilt matrices.)
}

# Standardize list of ExpressionSets
esets <-
  if (is.list(gse_obj) && all(sapply(gse_obj, inherits, "ExpressionSet"))) gse_obj else
    if (inherits(gse_obj, "ExpressionSet")) list(gse_obj) else
      stop("Could not obtain ExpressionSet objects. Check the GSE ID or GEO record structure.")

# ---- Robust helpers ----
make_unique_names <- function(x) make.unique(x, sep = "__")

parse_characteristics <- function(pheno_df) {
  # 1) Ensure unique names immediately
  ph <- tibble::as_tibble(pheno_df, .name_repair = "unique")
  
  # 2) Find characteristics columns
  ch_cols <- grep("^characteristics", names(ph), value = TRUE, ignore.case = TRUE)
  if (!length(ch_cols)) {
    names(ph) <- make_unique_names(names(ph))
    return(ph)
  }
  
  # 3) Long format, parse "key: value" pairs (robust to missing/blank)
  long <- ph |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(ch_cols),
      names_to = "char_col",
      values_to = "char_val"
    ) |>
    dplyr::mutate(char_val = as.character(char_val)) |>
    dplyr::filter(!is.na(char_val), char_val != "") |>
    dplyr::mutate(
      has_colon = stringr::str_detect(char_val, ":"),
      key = dplyr::if_else(
        has_colon,
        stringr::str_trim(stringr::str_remove(char_val, ":.*$")),
        char_col
      ),
      val = dplyr::if_else(
        has_colon,
        stringr::str_trim(stringr::str_replace(char_val, "^[^:]+:\\s*", "")),
        char_val
      )
    ) |>
    dplyr::group_by(.row, key) |>
    dplyr::mutate(inst = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(key = dplyr::if_else(inst > 1, paste0(key, "__", inst), key)) |>
    dplyr::select(.row, key, val)
  
  # 4) Wide format; if no parsed rows, return base
  if (nrow(long) == 0) {
    names(ph) <- make_unique_names(names(ph))
    return(ph)
  }
  
  wide <- tidyr::pivot_wider(
    long,
    id_cols = .row,
    names_from = key,
    values_from = val,
    values_fn = list(val = function(x) paste(unique(x), collapse = "; "))
  )
  
  out <- ph |>
    dplyr::mutate(.row = dplyr::row_number()) |>
    dplyr::left_join(wide, by = ".row") |>
    dplyr::select(-.row)
  
  # 5) Final guarantee: unique column names
  names(out) <- make_unique_names(names(out))
  out
}

# Quick EDA function with safe plotting
quick_eda <- function(meta, tag, outdir) {
  # Save a text summary
  summ_path <- file.path(outdir, paste0(tag, "_summary.txt"))
  sink(summ_path)
  cat("== Summary for", tag, "==\n")
  cat("N samples:", nrow(meta), "\n")
  cat("Columns:", paste(names(meta), collapse = ", "), "\n\n")
  
  # Show top informative categorical columns
  cat("-- Categorical columns (<= 20 unique, excluding IDs):\n")
  for (cn in names(meta)) {
    if (is.character(meta[[cn]]) || is.factor(meta[[cn]])) {
      u <- unique(meta[[cn]])
      if (length(u) > 1 && length(u) <= 20 && !grepl("geo_accession|sample|title|source_name", cn, ignore.case=TRUE)) {
        cat("\n", cn, "(", length(u), "levels )\n", sep="")
        print(table(meta[[cn]], useNA="ifany"))
      }
    }
  }
  
  # Basic missingness
  miss <- sapply(meta, function(x) sum(is.na(x) | x==""))
  cat("\n-- Missing values per column --\n")
  print(sort(miss, decreasing = TRUE))
  sink()
  
  # Try a couple of standard plots if columns exist
  # Age-like variable
  age_cols <- grep("age", names(meta), value = TRUE, ignore.case = TRUE)
  if (length(age_cols)) {
    ac <- age_cols[1]
    suppressWarnings({
      aa <- suppressWarnings(as.numeric(meta[[ac]]))
      p <- ggplot(data.frame(age=aa), aes(x=age)) +
        geom_histogram(bins = 30) + theme_minimal() +
        labs(title = paste(tag, "-", ac, "histogram"), x = ac, y = "Count")
      ggsave(file.path(outdir, paste0(tag, "_", ac, "_hist.png")), p, width=6, height=4, dpi=120)
    })
  }
  # Sex-like variable
  sex_cols <- grep("sex|gender", names(meta), value = TRUE, ignore.case = TRUE)
  if (length(sex_cols)) {
    sc <- sex_cols[1]
    p <- ggplot(meta, aes(x = .data[[sc]])) +
      geom_bar() + theme_minimal() +
      labs(title = paste(tag, "-", sc, "counts"), x = sc, y = "N")
    ggsave(file.path(outdir, paste0(tag, "_", sc, "_bar.png")), p, width=6, height=4, dpi=120)
  }
}




# Process each platform-specific ExpressionSet
for (i in seq_along(esets)) {
  eset <- esets[[i]]
  gpl  <- annotation(eset)
  tag  <- paste0(gse_id, "_", gpl %||% paste0("set", i))
  
  # Extract and tidy phenotype data (pData)
  pd <- Biobase::pData(eset) |>
    tibble::as_tibble(rownames = "geo_accession", .name_repair = "unique")
  pd <- parse_characteristics(pd)
  
  # Add simple derived columns
  pd <- pd |>
    mutate(across(everything(), ~ifelse(. %in% c("NA","NaN","null","None",""), NA, .))) |>
    mutate(across(where(is.character), ~na_if(., "")))
  
  # Save metadata
  out_csv <- file.path(outdir, paste0(tag, "_metadata.csv"))
  readr::write_csv(pd, out_csv)
  message("Saved metadata: ", out_csv)
  
  # Quick EDA
  quick_eda(pd, tag, outdir)
}

# Also save a top-level note about platforms & sizes
top_path <- file.path(outdir, paste0(gse_id, "_platforms_overview.txt"))
sink(top_path)
cat("== Platforms overview for", gse_id, "==\n")
for (i in seq_along(esets)) {
  eset <- esets[[i]]
  gpl  <- annotation(eset)
  cat(sprintf("[%d] GPL: %s | samples: %d | features: %d\n",
              i, gpl, ncol(eset), nrow(eset)))
}
sink()

message("Done. See outputs in: ", normalizePath(outdir))

