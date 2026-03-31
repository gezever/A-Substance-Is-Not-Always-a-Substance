# ==============================================================================
# 05_overlap_lists.R
# Analysis 5: Overlap between regulatory lists (UpSet plot)
#
# PURPOSE
# -------
# Multiple European regulatory instruments address chemical substances, but
# they may list the same substance independently.  This analysis quantifies
# how much overlap exists across lists at the structural level (InChIKey),
# revealing which substances are of concern to multiple regulatory instruments
# simultaneously and which are source-specific.
#
# The UpSet plot is preferred over a Venn diagram because it scales to more
# than four sets while preserving set-intersection quantification.
#
# DATA PROVENANCE
# ---------------
# Input: `data/processed/all_substances.rds`
# Only entries with a known InChIKey are included; structure-undefined entries
# cannot be reliably matched across lists.  Sources with fewer than 10 unique
# InChIKeys are excluded to avoid visually dominant near-empty sets.
#
# OUTPUTS
# -------
# output/figures/Analysis_5_Overlap_between_lists_UpSet.pdf
# ==============================================================================

library(dplyr)
library(tidyr)
library(UpSetR)
library(here)
library(grid)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 5: UpSet plot of cross-list substance overlap
# ==============================================================================
# INTENT
# Identify which combinations of regulatory lists share the most substances.
# A high intersection across concern-level lists (SVHC, restriction, POPs)
# would indicate that multiple instruments converge on the same substances.
# Sources with very few InChIKeys are excluded to keep the plot readable;
# the threshold of ≥ 10 unique InChIKeys is a display heuristic.
# Own addition: minimum-size filter of n ≥ 10; not mandated by any framework.
# ==============================================================================

inchi_per_source <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(source, inchikey)

# Own addition: ≥ 10 InChIKey threshold to avoid near-empty bars in the plot
sources_sufficient <- inchi_per_source |>
  count(source) |>
  filter(n >= 10) |>
  pull(source)

upset_input <- inchi_per_source |>
  filter(source %in% sources_sufficient) |>
  mutate(present = 1L) |>
  pivot_wider(names_from  = source,
              values_from = present,
              values_fill = 0L) |>
  select(-inchikey) |>
  as.data.frame()

pdf(
  here("output", "figures", "Analysis_5_Overlap_between_lists_UpSet.pdf"),
  width = 14, height = 7
)
upset(
  upset_input,
  nsets           = ncol(upset_input),
  nintersects     = 30,
  order.by        = "freq",
  mb.ratio        = c(0.55, 0.45),
  main.bar.color  = "#4a90d9",
  sets.bar.color  = "#e07b3a",
  text.scale      = c(1.2, 1.1, 1, 1, 1.1, 1),
  mainbar.y.label = "Overlap (number of shared InChIKeys)",
  sets.x.label    = "InChIKeys per list"
)
grid.text(
  "Analysis 5: Overlap of structure-defined substances across lists",
  x = 0.65, y = 0.97,
  gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

message("05_overlap_lists.R: analysis 5 completed")
