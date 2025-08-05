#!/usr/bin/env Rscript
library("tidyverse")
library("Seurat")
library("Matrix")
library("jsonlite")
library("transport")

args = commandArgs(trailingOnly=TRUE)
case_pos <- args[1]
ctrl_pos <- args[2]
outdir <- args[3]


# setwd("/home/projects/dtu_00062/people/sorsan/ob_anonymization_dataloss")
# case_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_case.rds"
# ctrl_pos <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/D1_ctrl.rds"
# outdir <- "out/data/D1/default/init_align/P1/param_0/second_align/M1/param_0/metrics/m1/default/D1.somefile.txt"

# startdir <- "/Users/jkd448/Documents/OmniBenchmark/Development/R_sandbox"
# setwd(startdir)
# case_pos <- paste0(startdir,"/D1_case.rds")
# ctrl_pos <- paste0(startdir,"/D1_ctrl.rds")
# outdir <- paste0(startdir,"/outdir/D1.somefile.json")

case_obj <- readRDS(case_pos)
ctrl_obj <- readRDS(ctrl_pos)

# Do something - write earth movers distance

# Normalize
case_obj <- NormalizeData(case_obj)
ctrl_obj <- NormalizeData(ctrl_obj)

# Convert Seurat objects to expression matrices
expr_orig <- case_obj[["RNA"]]$data
expr_anon <- ctrl_obj[["RNA"]]$data


# Find cells in common
cells_in_orig <- colnames(expr_orig)
cells_in_anon <- colnames(expr_anon)
cells_in_common <- intersect(colnames(expr_orig),colnames(expr_anon))

# Find genes with non-zero expression
nonz_genes_orig <- rownames(expr_orig)[rowSums(expr_orig) != 0]
nonz_genes_anon <- rownames(expr_anon)[rowSums(expr_anon) != 0]
nonz_genes_common <- intersect(nonz_genes_orig,nonz_genes_anon)

# Remove unique cells, and genes with zero-expression
expr_orig <- expr_orig[nonz_genes_common,cells_in_common]
expr_anon <- expr_anon[nonz_genes_common,cells_in_common]

# Calculate correlation genewise
expr_orig_centered <- expr_orig - rowMeans(expr_orig)
expr_anon_centered <- expr_anon - rowMeans(expr_anon)
numerator <- rowSums(expr_orig_centered * expr_anon_centered)
denom_orig <- sqrt(rowSums(expr_orig_centered^2))
denom_anon <- sqrt(rowSums(expr_anon_centered^2))
gw_corr <- numerator / (denom_orig * denom_anon)

# Remove NA values
gw_corr <- gw_corr[!is.na(gw_corr)]

rowSums(case_obj[["RNA"]]$counts[names(gw_corr)[gw_corr < 0.5],])

# Format output
json_obj <- list(
  module = "R4_correlation",
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
  metrics = list(
    correlation = list(
      gene_wise = list(
        mean = mean(gw_corr),
        median = median(gw_corr),
        std = sd(gw_corr),
        min = min(gw_corr),
        max = max(gw_corr),
        q25 = quantile(gw_corr, 0.25),
        q75 = quantile(gw_corr, 0.75)
      ),
      poorly_matched_genes = list(
        count = sum(gw_corr < 0.5),
        threshold = 0.5,
        gene_ids = names(gw_corr)[gw_corr < 0.5],
        emd_values = unname(gw_corr[gw_corr < 0.5])
      )
    )
  )
)

# Save output
write_json(json_obj, outdir, pretty = TRUE, auto_unbox = TRUE)
