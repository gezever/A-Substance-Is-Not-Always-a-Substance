# ==============================================================================
# 07_network_visualisation.R
# Analysis 7: Source-source network — shared structure-defined substances
#
# PURPOSE
# -------
# Multiple European regulatory instruments independently list chemical
# substances.  This analysis visualises the *pairwise overlap structure* among
# regulatory sources as a network: each source is a node, and an edge connects
# two sources when they share at least `min_shared` structure-defined
# substances (InChIKey-linked).  Edge width encodes the number of shared
# substances; node size encodes the total number of structure-defined
# substances in that source.
#
# The network reveals which instruments form regulatory clusters (groups of
# lists that regulate largely the same chemicals), which sources are
# peripheral (few shared substances with others), and which sources are hubs
# (highly connected to many other instruments).
#
# Compared to the Jaccard heatmap (Analysis 12a), which shows exact pairwise
# similarity values, the network emphasises topological relationships and is
# easier to scan for overall structure.
#
# DATA PROVENANCE
# ---------------
# Input: `data/processed/all_substances.rds`
# Only entries with a known InChIKey are included; non-structure entries
# cannot be reliably matched across lists and are therefore excluded from the
# overlap count.
#
# METHODOLOGY
# -----------
# ## Established frameworks used as anchors
#
# | Framework | Role in methodology |
# |---|---|
# | **Fruchterman–Reingold / stress-majorisation layout** (Fruchterman & Reingold 1991, *Graph Drawing by Force-Directed Placement*, Softw. Pract. Exp. 21(11):1129–1164. DOI: 10.1002/spe.4380211102) | Force-directed layout; places highly connected nodes close together, revealing regulatory clusters without imposing a fixed geometry. |
#
# ## Own methodological additions
#
# | Choice | Justification |
# |---|---|
# | Minimum shared InChIKeys threshold (`min_shared`) | Pairs with very few shared substances (< `min_shared`) produce visually indistinguishable thin edges that add clutter without analytical value.  The threshold is documented in the plot subtitle so the reader knows what is excluded. |
# | Edge width scaled to log₁₀(n_shared + 1) | Raw counts span several orders of magnitude; log scaling prevents a few dominant edges from compressing all others to near-zero width. |
# | Node size scaled to log₁₀(n_inchikeys + 1) | Same rationale as edge width: list sizes vary widely across sources. |
#
# INTERPRETATION
# --------------
# The network has one node per regulatory source and one edge per pair that
# shares >= `min_shared` InChIKeys.  Four reading patterns matter:
#
# 1. **Isolated nodes (no edges)**
#    The source shares fewer than `min_shared` InChIKeys with every other
#    source.  This may mean the list is genuinely distinct in scope (e.g.,
#    a sector-specific instrument), or that its structural coverage is too low
#    to form matches (see Analysis 6 for coverage per source).
#
# 2. **Tight cluster (dense subgraph)**
#    Sources with many mutual edges and similar node sizes regulate largely
#    the same pool of chemicals.  This indicates regulatory convergence: the
#    same substance concerns multiple instruments simultaneously.  Substances
#    in these clusters will score high on `list_score` in Analysis 14.
#
# 3. **Hub node (many edges, large size)**
#    A source connected to most others is a broad-scope instrument.  Its
#    substances are widely regulated; removing it from scope would leave large
#    gaps in multi-list coverage.
#
# 4. **Thick vs. thin edges**
#    Edge width is proportional to log₁₀(shared InChIKeys + 1).  A thick edge
#    between two sources means a large absolute overlap; pairs of equal Jaccard
#    similarity but different sizes will show different edge widths.  For
#    normalised pairwise similarity, see the Jaccard heatmap (Analysis 12a).
#
# **Actionable reading:** sources in a tight cluster are candidates for
# harmonised regulatory treatment; isolated sources identify gaps where a
# substance could be on one list but overlooked by others.
#
# OUTPUTS
# -------
# output/figures/Analysis_7_Network_source_overlap.pdf
# output/tables/Analysis_7_edge_list.csv
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(igraph)
library(ggraph)
library(ggrepel)
library(readr)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 7: Source-source projected network
# ==============================================================================
# INTENT
# Project the bipartite substance–list structure onto the list dimension:
# two sources are connected if they share at least `min_shared` InChIKeys.
# This makes the graph legible (one node per source, ~10-15 nodes total) and
# directly answers "which regulatory instruments regulate the same chemicals?"
# ==============================================================================

inchi_src   <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(source, inchikey)

sources_vec <- sort(unique(inchi_src$source))

# Own addition: min_shared = 5 suppresses near-zero overlaps that add visual
# noise; the threshold is visible in the plot subtitle.
min_shared <- 5L

# Compute pairwise shared InChIKey counts for all source pairs
pairs    <- combn(sources_vec, 2, simplify = FALSE)
edge_rows <- lapply(pairs, function(pair) {
  n <- length(intersect(
    inchi_src$inchikey[inchi_src$source == pair[1]],
    inchi_src$inchikey[inchi_src$source == pair[2]]
  ))
  if (n >= min_shared)
    data.frame(from = pair[1], to = pair[2], n_shared = n)
})
edge_list <- do.call(rbind, Filter(Negate(is.null), edge_rows))

# Node attributes: total unique InChIKeys per source
node_df <- inchi_src |>
  count(source, name = "n_inchikeys") |>
  rename(name = source)

g <- graph_from_data_frame(edge_list, directed = FALSE, vertices = node_df)

set.seed(42)
p7 <- ggraph(g, layout = "stress") +
  geom_edge_link(
    aes(width = log10(n_shared + 1), alpha = log10(n_shared + 1)),
    colour = "#4a90d9"
  ) +
  geom_node_point(
    aes(size = log10(n_inchikeys + 1)),
    colour = "#e05c5c"
  ) +
  geom_node_text(
    aes(label = paste0(name, " (", n_inchikeys, ")")),
    repel = TRUE, size = 3, colour = "#222222", family = ""
  ) +
  scale_edge_width(range = c(0.4, 5), guide = "none") +
  scale_edge_alpha(range = c(0.25, 0.85), guide = "none") +
  scale_size(range = c(3, 12), guide = "none") +
  labs(
    title    = "Regulatory source network — shared structure-defined substances",
    subtitle = paste0(
      "Edge width proportional to log10(shared InChIKeys); ",
      "node size proportional to log10(list size); ",
      "edges shown for >= ", min_shared, " shared InChIKeys"
    )
  ) +
  theme_graph(base_size = 12, base_family = "") +
  theme(plot.subtitle = element_text(colour = "grey40", size = 9))

print(p7)
ggsave(p7,
       filename = here("output", "figures",
                       "Analysis_7_Network_source_overlap.pdf"),
       device = "pdf",
       height = 8, width = 11, units = "in")

# Columns (Analysis_7_edge_list.csv):
#   from       — first regulatory source in the pair
#   to         — second regulatory source in the pair
#   n_shared   — number of unique InChIKeys present in both sources;
#                only pairs with n_shared >= min_shared are included
write_csv(edge_list,
          here("output", "tables", "Analysis_7_edge_list.csv"))

message(sprintf(
  "Analysis 7: %d sources, %d edges (>= %d shared InChIKeys)",
  vcount(g), ecount(g), min_shared
))

message("07_network_visualisation.R: analysis 7 completed")
