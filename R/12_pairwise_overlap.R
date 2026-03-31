# ==============================================================================
# 12_pairwise_overlap.R
# Analysis 12: Pairwise source overlap — Jaccard heatmap and UpSet plot
#
# PURPOSE
# -------
# Analysis 5 visualised multi-way intersections across all sufficiently large
# sources.  This analysis focuses on *pairwise* overlap using the Jaccard
# index: for each pair of regulatory sources, it measures the fraction of
# shared InChIKeys relative to the combined pool.  A high Jaccard index
# between two sources means that substances regulated by one instrument are
# largely also regulated by the other — indicative of regulatory convergence
# or cross-referencing.
#
# A complementary UpSet plot (12b) restricts to the "obligation lists" (the
# core REACH/CLP/POPs instruments with binding legal consequences) to show
# how multi-list membership is distributed among the most consequential
# regulatory instruments.
#
# DATA PROVENANCE
# ---------------
# Input: `data/processed/all_substances.rds`
# Only structure-defined entries (InChIKey present) are included; the Jaccard
# index is computed on unique InChIKeys per source.  Sources without any
# InChIKey are excluded automatically (empty set → Jaccard = 0 against all).
#
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **Jaccard similarity index** (Jaccard 1912, *Bull. Soc. Vaud. Sci. Nat.* 37:547–579) | Pairwise overlap metric: J(A,B) = \|A∩B\| / \|A∪B\|; 1 = identical sets, 0 = no overlap.  Chosen because it normalises for set size, making sources of different lengths comparable. |
# | **Hierarchical agglomerative clustering** (Ward 1963; Murtagh & Legendre 2014, *J. Classif.* 31(3):274–295) | Orders heatmap rows/columns by Jaccard dissimilarity (1−J) so that similar sources appear adjacent, making cluster structure visible as diagonal blocks. |
# | **REACH Regulation (EC) No 1907/2006, Arts. 57–73; CLP Regulation (EC) No 1272/2008; EU POPs Regulation (EU) 2019/1021** | Legal basis for the `obligation_sources` selection in Analysis 12b: only instruments that impose binding legal consequences are included. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | Average linkage (`hclust(method = "average")`) | No framework mandates a specific linkage method; "average" (UPGMA) balances the extremes of complete linkage (sensitive to outliers) and single linkage (chaining tendency), making it a robust default for exploratory clustering. |
# | Jaccard computed on unique InChIKeys (not records) | Duplicate records for the same structure within a source would inflate the intersection count; deduplication ensures the index measures structural overlap, not administrative repetition. |
# | Text annotations suppressed for Jaccard < 0.02 | Values below 0.02 are visually indistinguishable from zero and clutter the heatmap; the threshold is aesthetic, not analytical. |
#
# INTERPRETATION
# --------------
# **Analysis 12a — Jaccard heatmap**
# The heatmap is symmetric; the diagonal is always 1.0 (a source is identical
# to itself).  Read the off-diagonal values:
#
# - **High Jaccard (> 0.5, dark blue):** two sources regulate largely the same
#   chemicals.  This may reflect legal cross-referencing (one instrument
#   explicitly references another), shared scope, or both.  Substances in this
#   cell appear on both lists and will accumulate high `list_score` in
#   Analysis 14.
# - **Moderate Jaccard (0.1–0.5):** partial overlap; both instruments share a
#   common core but each also contains source-specific substances.
# - **Near-zero Jaccard (< 0.05, white):** the instruments regulate largely
#   distinct chemical pools.  Cross-list harmonisation is irrelevant for these
#   pairs; each must be handled independently.
#
# Cluster blocks along the diagonal indicate groups of instruments with
# similar regulatory scope.  Sources that fall outside any cluster are
# structurally peripheral.
#
# **Analysis 12b — UpSet plot (obligation lists)**
# The intersection bars show how many InChIKeys appear in exactly the
# combination of obligation lists indicated below.  The tallest bar is the
# most common membership pattern.  Substances appearing in only one obligation
# list are in source-specific single-set bars; substances in multiple bars
# represent the regulatory core that multiple binding instruments converge on.
#
# An empty intersection (no bar) for a combination means no InChIKey is
# simultaneously on all those lists — either a regulatory gap or a structural
# coverage issue.
#
# OUTPUTS
# -------
# output/figures/Analysis_12a_Jaccard_heatmap.pdf
# output/figures/Analysis_12b_UpSet_obligation_lists.pdf
# output/tables/Analysis_12a_jaccard_matrix.csv
# output/tables/Analysis_12b_set_sizes.csv
# output/tables/Analysis_12b_membership.csv
# ==============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(UpSetR)
library(readr)
library(here)
library(grid)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 12a: Pairwise Jaccard heatmap
# ==============================================================================
# INTENT
# Identify pairs of regulatory sources that substantially overlap on
# structure-defined substances.  High Jaccard values may indicate that two
# instruments are addressing the same chemicals for similar regulatory
# purposes; near-zero values indicate distinct regulatory scope.
# The hierarchical clustering orders sources so that the cluster structure is
# immediately visible as diagonal blocks.
# ==============================================================================

inchi_src   <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(source, inchikey)

sources_vec <- sort(unique(inchi_src$source))

jaccard_mat <- outer(
  sources_vec, sources_vec,
  FUN = Vectorize(function(a, b) {
    set_a <- inchi_src$inchikey[inchi_src$source == a]
    set_b <- inchi_src$inchikey[inchi_src$source == b]
    n_int <- length(intersect(set_a, set_b))
    n_uni <- length(union(set_a, set_b))
    if (n_uni == 0L) 0 else n_int / n_uni
  })
)
rownames(jaccard_mat) <- sources_vec
colnames(jaccard_mat) <- sources_vec

# Hierarchical clustering on Jaccard dissimilarity (1 − J)
hc        <- hclust(as.dist(1 - jaccard_mat), method = "average")
src_order <- sources_vec[hc$order]

jaccard_long <- as.data.frame(jaccard_mat) |>
  rownames_to_column("source_a") |>
  pivot_longer(-source_a, names_to = "source_b", values_to = "jaccard") |>
  mutate(
    source_a = factor(source_a, levels = src_order),
    source_b = factor(source_b, levels = src_order)
  )

p12a <- ggplot(jaccard_long,
               aes(x = source_a, y = source_b, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(jaccard > 0.02, sprintf("%.2f", jaccard), "")),
    size = 2.5, colour = "grey20"
  ) +
  scale_fill_gradient(low = "white", high = "#084594",
                      limits = c(0, 1), name = "Jaccard") +
  coord_fixed() +
  labs(
    title    = "Analysis 12a: Pairwise Jaccard index between sources (on InChIKey)",
    subtitle = "Hierarchically clustered; 1 = complete overlap, 0 = no overlap",
    x        = NULL,
    y        = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    plot.subtitle   = element_text(colour = "grey40"),
    legend.position = "right"
  )

print(p12a)
ggsave(p12a,
       filename = here("output", "figures", "Analysis_12a_Jaccard_heatmap.pdf"),
       device = "pdf",
       height = 10, width = 11, units = "in")

# Columns:
#   source_a — first regulatory source in the pair
#   source_b — second regulatory source in the pair
#   jaccard  — Jaccard similarity index (0–1); 1 = identical InChIKey sets,
#              0 = no shared InChIKeys; diagonal entries (source_a == source_b)
#              are always 1.0; the matrix is symmetric (jaccard(a,b) == jaccard(b,a))
write_csv(
  jaccard_long |> mutate(source_a = as.character(source_a),
                         source_b = as.character(source_b)),
  here("output", "tables", "Analysis_12a_jaccard_matrix.csv")
)

# ==============================================================================
# Analysis 12b: UpSet plot restricted to obligation lists
# ==============================================================================
# INTENT
# The full UpSet plot (analysis 5) includes administrative and positive lists
# that have limited regulatory consequence.  Restricting to the seven
# obligation lists — those that impose a binding obligation or are a step
# towards one under REACH, CLP, POPs, or pesticide regulation — shows how
# multi-list membership distributes among the highest-consequence instruments.
# Sources selected on the basis of their legal standing under REACH Regulation
# (EC) No 1907/2006 (Arts. 57–73), CLP Regulation (EC) No 1272/2008 (Annex VI),
# and EU POPs Regulation (EU) 2019/1021.
# ==============================================================================

obligation_sources <- c(
  "restriction_list",      # REACH Art. 67-73 — binding use restrictions
  "candidate_list",        # REACH Art. 59 — SVHC identification
  "authorisation_list",    # REACH Art. 62 — prior authorisation required
  "pops_list",             # EU POPs Reg. (EU) 2019/1021 — elimination
  "eu_positive_list",      # Positive list — administrative oversight signal
  "harmonised_list",       # CLP Reg. (EC) No 1272/2008 Annex VI — harmonised classification
  "svhc_identification"    # REACH Art. 57-59 — SVHC formal identification step
)

upset_obl <- inchi_src |>
  filter(source %in% obligation_sources) |>
  mutate(present = 1L) |>
  tidyr::pivot_wider(names_from  = source,
                     values_from = present,
                     values_fill = 0L) |>
  select(-inchikey) |>
  as.data.frame()

pdf(here("output", "figures", "Analysis_12b_UpSet_obligation_lists.pdf"),
    width = 14, height = 7)
UpSetR::upset(
  upset_obl,
  nsets           = ncol(upset_obl),
  nintersects     = 25,
  order.by        = "freq",
  mb.ratio        = c(0.55, 0.45),
  main.bar.color  = "#4a90d9",
  sets.bar.color  = "#e07b3a",
  text.scale      = c(1.2, 1.1, 1, 1, 1.1, 1),
  mainbar.y.label = "Overlap (shared InChIKeys)",
  sets.x.label    = "InChIKeys per list"
)
grid.text(
  "Analysis 12b: Overlap of obligation lists (InChIKey)",
  x = 0.65, y = 0.97,
  gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

# Columns (Analysis_12b_set_sizes.csv):
#   source           — regulatory source identifier
#   n_inchikeys      — number of unique InChIKeys in that obligation list;
#                      corresponds to the bar lengths on the left side of the UpSet plot
#   n_non_structured — number of records in that obligation list without an InChIKey
#                      (no resolvable chemical structure); these entries are excluded
#                      from the UpSet plot but represent real regulatory scope
#   n_total          — n_inchikeys + n_non_structured; total records in the list
write_csv(
  all_substances |>
    filter(source %in% obligation_sources) |>
    group_by(source) |>
    summarise(
      n_inchikeys      = n_distinct(inchikey[!is.na(inchikey)]),
      n_non_structured = sum(is.na(inchikey)),
      n_total          = n(),
      .groups = "drop"
    ) |>
    arrange(desc(n_inchikeys)),
  here("output", "tables", "Analysis_12b_set_sizes.csv")
)

# Columns (Analysis_12b_membership.csv):
#   inchikey         — InChIKey of the substance (primary linking key)
#   <obligation_source_1 … N> — binary indicator (1 = present, 0 = absent)
#                 for each of the seven obligation sources; missing combinations
#                 are filled with 0.  This is the exact matrix fed to UpSetR.
write_csv(
  inchi_src |>
    filter(source %in% obligation_sources) |>
    mutate(present = 1L) |>
    pivot_wider(names_from  = source,
                values_from = present,
                values_fill = 0L),
  here("output", "tables", "Analysis_12b_membership.csv")
)

message("12_pairwise_overlap.R: analysis 12 completed")
