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
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **UpSet plot** (Lex et al. 2014, *UpSet: Visualization of Intersecting Sets*, IEEE Trans. Vis. Comput. Graph. 20(12):1983–1992. DOI: 10.1109/TVCG.2014.2346248) | Multi-set intersection visualisation; preferred over Venn diagrams for more than four sets because it encodes intersection sizes as bar heights rather than area, and scales to any number of sets. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | Minimum ≥ 10 InChIKeys per source | Sources with fewer than 10 unique InChIKeys produce near-empty bars that are visually indistinguishable from zero and consume space without adding information.  The threshold is a display heuristic, not a substantive filter. |
# | `nintersects = 30` (top-30 intersections shown) | Showing all 2ⁿ − 1 possible intersections would be unreadable for n > 5; limiting to the 30 most frequent intersections captures the dominant overlap patterns while keeping the figure legible. |
#
# INTERPRETATION
# --------------
# The UpSet plot has two components:
#
# 1. **Set size bars (left side):** the total number of unique InChIKeys in
#    each regulatory source.  Large bars indicate broad-scope instruments;
#    small bars are narrow or structurally sparse lists.
#
# 2. **Intersection bars (top):** the number of InChIKeys shared by exactly
#    the combination of sources indicated by the filled dots below.  A single
#    filled dot = substance appears in only that source; multiple filled dots =
#    substance appears in all indicated sources simultaneously.
#
# **Reading patterns:**
# - **Tall single-source bars:** many substances regulated exclusively by one
#   instrument.  These are source-specific concerns not yet picked up by other
#   regulatory lists.
# - **Tall multi-source bars:** a substantial chemical core is co-regulated by
#   multiple instruments simultaneously.  These substances accumulate high
#   `list_score` in Analysis 14.
# - **Small intersections between large sets:** two large lists share few
#   substances despite their size — indicating distinct regulatory scope, not
#   structural coincidence.
#
# For pairwise similarity normalised by set size, see the Jaccard heatmap
# (Analysis 12a).  For the obligation-list subset only, see Analysis 12b.
#
# OUTPUTS
# -------
# output/figures/Analysis_5_Overlap_between_lists_UpSet.pdf
# output/tables/Analysis_5_set_sizes.csv
# output/tables/Analysis_5_membership.csv
# ==============================================================================

library(dplyr)
library(tidyr)
library(UpSetR)
library(here)
library(grid)
library(readr)

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
  "Overlap of structure-defined substances across lists",
  x = 0.65, y = 0.97,
  gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

# Columns (Analysis_5_set_sizes.csv):
#   source      — regulatory source identifier (only sources with >= 10 unique
#                 InChIKeys are included; others were excluded before plotting)
#   n_inchikeys — number of unique InChIKeys in this source; corresponds to
#                 the bar lengths on the left side of the UpSet plot
write_csv(
  inchi_per_source |>
    filter(source %in% sources_sufficient) |>
    count(source, name = "n_inchikeys"),
  here("output", "tables", "Analysis_5_set_sizes.csv")
)

# Columns (Analysis_5_membership.csv):
#   inchikey         — InChIKey of the substance; primary linking key
#   <source_1 … N>  — binary indicator (1 = present in source, 0 = absent)
#                     for each source that passed the >= 10 InChIKey threshold;
#                     the column names are the source identifiers.
#                     This is the exact matrix fed to UpSetR::upset().
write_csv(
  inchi_per_source |>
    filter(source %in% sources_sufficient) |>
    mutate(present = 1L) |>
    pivot_wider(names_from = source, values_from = present, values_fill = 0L),
  here("output", "tables", "Analysis_5_membership.csv")
)

message(sprintf(
  "Analysis 5: %d sources (>= 10 InChIKeys); top intersection shown up to 30",
  length(sources_sufficient)
))

message("05_overlap_lists.R: analysis 5 completed")
