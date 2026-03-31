# ==============================================================================
# 07_network_visualisation.R
# Analysis 7: Bipartite network visualisation (substance ↔ regulatory list)
#
# PURPOSE
# -------
# A substance that appears on multiple regulatory lists can be thought of as
# connected to those lists in a bipartite graph.  Visualising this network
# reveals structural patterns that tabular analyses obscure: which substances
# are regulatory hubs (many list memberships), which lists cluster together
# (many shared substances), and how fragmented or cohesive the regulatory
# landscape is overall.
#
# DATA PROVENANCE
# ---------------
# Input: `data/processed/all_substances.rds`
# Only substances present in ≥ 2 lists are included; single-list substances
# would form isolated nodes that add visual noise without analytical value.
# The analysis is restricted to structure-defined entries (InChIKey present)
# to ensure that "same substance on two lists" reflects genuine chemical
# identity and not a coincidental name match.
#
# METHODOLOGY
# -----------
# The graph is built with igraph (Csárdi & Nepusz 2006) and laid out with the
# Fruchterman–Reingold force-directed algorithm (seed fixed for
# reproducibility).  List nodes are labelled; substance nodes are unlabelled
# to avoid overplotting.
# Own addition: ≥ 2 list threshold and restriction to InChIKey-bearing entries
# are display heuristics, not framework requirements.
#
# OUTPUTS
# -------
# output/figures/Analysis_7_Network_visualisation—bipartite_graph.pdf
# ==============================================================================

library(dplyr)
library(ggplot2)
library(here)
library(igraph)
library(ggraph)
library(ggrepel)

# ------------------------------------------------------------------------------
# Load clean data
# ------------------------------------------------------------------------------

all_substances <- readRDS(here("data", "processed", "all_substances.rds"))

# ==============================================================================
# Analysis 7: Bipartite substance ↔ list network
# ==============================================================================
# INTENT
# Force-directed layout reveals list clusters (groups of lists that share many
# substances) and substance hubs (individual molecules that trigger obligations
# across many instruments simultaneously).  Hub substances in the network
# correspond to the highest-priority candidates from a regulatory perspective.
# Own addition: Fruchterman–Reingold layout with set.seed(42) for
# reproducibility; the specific seed value is arbitrary.
# ==============================================================================

# Restrict to substances present in ≥ 2 lists for interpretable graph density
inchi_multi <- all_substances |>
  filter(!is.na(inchikey)) |>
  distinct(inchikey, source) |>
  group_by(inchikey) |>
  filter(n() >= 2) |>
  ungroup()

edges <- inchi_multi |>
  select(from = inchikey, to = source)

node_substances <- data.frame(
  name = unique(edges$from),
  type = "substance",
  stringsAsFactors = FALSE
)
node_sources <- data.frame(
  name = unique(edges$to),
  type = "list",
  stringsAsFactors = FALSE
)
nodes <- rbind(node_substances, node_sources)

g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

set.seed(42)
p7 <- ggraph(g, layout = "fr") +
  geom_edge_link(alpha = 0.15, colour = "grey60") +
  geom_node_point(aes(colour = type, size = type)) +
  geom_node_text(
    aes(label = if_else(type == "list", name, NA_character_)),
    repel = TRUE, size = 3.5, colour = "#333333"
  ) +
  scale_colour_manual(values = c("substance" = "#4a90d9", "list" = "#e05c5c")) +
  scale_size_manual(values  = c("substance" = 1.5,        "list" = 5)) +
  labs(
    title    = "Bipartite network substance \u2194 regulatory list",
    subtitle = paste0(
      "Only substances present in \u2265 2 lists shown (n = ",
      nrow(node_substances), " substances, ", nrow(node_sources), " lists)"
    ),
    colour = NULL,
    size   = NULL
  ) +
  theme_graph(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.subtitle   = element_text(colour = "grey40")
  )

print(p7)
ggsave(p7,
       filename = here("output", "figures",
                       "Analysis_7_Network_visualisation\u2014bipartite_graph.pdf"),
       device = "pdf",
       height = 5, width = 10, units = "in")

message("07_network_visualisation.R: analysis 7 completed")
